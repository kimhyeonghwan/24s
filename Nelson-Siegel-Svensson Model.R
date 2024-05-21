library(tidyverse)

# https://rpubs.com/beniamino98/NelsonSiegelSvensonn

# (nelson_siegel) return NSS model rates for given parameters
# (optim_nelson_siegel) return optimal tibble for given term-structure, using finite difference method
# maturity : vector of year maturity
# term_sturucture : vector of real spot rate, order to maturity, omit %
# params : vector of NSS medel 6 parameters, beta0 to tau2

nelson_siegel <- function(maturity, params){
  
  # parameters 
  beta0 = params[1]
  beta1 = params[2]
  beta2 = params[3]
  beta3 = params[4]
  tau1  = params[5]
  tau2  = params[6]
  
  # function 
  beta0 + 
    beta1 * (1 - exp(-maturity/tau1)) / (maturity/tau1) + 
    beta2 * ((1 - exp(-maturity/tau1)) / (maturity/tau1) - exp(-maturity/tau1)) + 
    beta3 * ((1 - exp(-maturity/tau2)) / (maturity/tau2) - exp(-maturity/tau2)) 
  
}
optim_nelson_siegel <- function(term_structure = NULL, maturity = NULL, params = rep(0.1, 6),
                                contraint_threshold = 0.1, tau_threshold = 10, beta_threshold = 30){
  
  loss_function <- function(params){
    
    # parameters 
    beta0 = params[1]
    beta1 = params[2]
    beta2 = params[3]
    beta3 = params[4]
    tau1  = params[5]
    tau2  = params[6]
    
    # contraints
    contraint = abs(beta0 + beta1 - term_structure[1])
    
    if(contraint > contraint_threshold ){
      return(NA)
    }
    
    # if(tau1 > tau_threshold | tau2 > tau_threshold){
    #   return(NA)
    # }
    # 
    # if(beta2 > beta_threshold | beta3 > beta_threshold ){
    #   return(NA)
    # }
    
    sum(
      (term_structure - beta0 - 
         beta1 *  (1 - exp(-maturity/tau1)) / (maturity/tau1) -
         beta2 * ((1 - exp(-maturity/tau1)) / (maturity/tau1) - exp(-maturity/tau1)) - 
         beta3 * ((1 - exp(-maturity/tau2)) / (maturity/tau2) - exp(-maturity/tau2)) 
      )^2)
    
  }
  
  # set initial values for to beta0 and beta1 that respect the contraint 
  # beta0 is long-term rate & beta1 is short-term rate minus beta0
  params = c(term_structure[length(term_structure)],
             term_structure[1]-term_structure[length(term_structure)],
             params[3:6] )
  
  # beta0 & beta1 is equal to (short-term rate)/2
  # params = c(term_structure[1]/2,term_structure[1]/2,params[3:6] )
  
  # minimization of the loss function 
  optimization = optim(params, loss_function)
  
  # optimal parameters 
  optim_params = c(beta0 = optimization$par[1], 
                   beta1 = optimization$par[2], 
                   beta2 = optimization$par[3],
                   beta3 = optimization$par[4],
                   tau1  = optimization$par[5], 
                   tau2  = optimization$par[6])
  
  # fitted value from Nelson-Siegel base function 
  fitted_values = nelson_siegel(maturity, optim_params)
  
  # compute the mean square error
  mse_fit = sd(fitted_values - term_structure, na.rm = TRUE)
  
  # output 
  dplyr::tibble( 
    maturity = list(maturity),
    term_structure = list(term_structure),
    start_params = list(params),
    optim_params = list(optim_params),
    fit = list(fitted_values),
    mse = mse_fit
  )
} 

# (random_params) return random number follows an uniform distribution
# n.params: number of parameter to generate 
# params.min: for random generation, minimum parameter
# params.max: for random generation, maximum parameter
# seed: to control randomness 

random_params <- function(n.params = 1, params.min = 0, params.max = 5, seed = 1){
  
  set.seed(seed)
  
  runif(n.params, params.min, params.max)
  
}

# (calibrate_nelson_siegel) return n-times "optim_nelson_siegel" outputs
#                           using Monte-carlo simulation for random 3 parameters(beta2 to tau2)
# n: number of simulations 
# params.min: for random generation, minimum parameter
# params.max: for random generation, maximum parameter
# verbose: disply progress in the importation.

calibrate_nelson_siegel <- function(object, n = 100, params.min = -5, params.max = 5, 
                                    contraint_threshold = 0.01, verbose = TRUE ){
  
  # initialize the parameters 
  term_structure = object$term_structure[[1]]
  maturity = object$maturity[[1]]
  
  # list containing all the simulations 
  simulations = list()
  
  # safe version to avoid errors if we made many simulations
  safe_optim = purrr::safely(optim_nelson_siegel)
  
  for(i in 1:n){
    
    # generate a random seed
    random_seed = mean(random_params(10, 0, 100000, seed = i))
    
    # generate random parameters 
    random_params = random_params(6, params.min, params.max, seed = random_seed)
    
    simulations[[i]] = safe_optim(term_structure  = term_structure, maturity = maturity, params = random_params, contraint_threshold = contraint_threshold)$result
    
    if( verbose &(i%%50 == 0)){message("Simulations: ", i, "/", n)}
    
  }
  
  # unique dataset for all the simulation
  simulations_df = dplyr::bind_rows(simulations) 
  
  # add the initial object
  simulations_df = dplyr::bind_rows(object, simulations) 
  
  # index for the simulations 
  simulations_df = dplyr::mutate(simulations_df, n = 1:nrow(simulations_df)) 
  
  return(simulations_df)
  
}

# (optimal_params_ns) summarize of all above, return optimal result minimize MSE from n-times simulation.
# term_structure = NULL
# maturity = NULL
# n = 1000
# params.init = rep(0.1, 6)
# params.min = -5
# params.max = 5
# contraint_threshold = 0.1
# label = NULL (label for the plot)
# verbose = TRUE

optimal_params_ns <- function(term_structure = NULL, maturity = NULL, 
                              n = 1000, params.init = rep(0.1, 6), params.min = -5, params.max = 5,
                              contraint_threshold = 0.1, label = NULL, verbose = TRUE){
  
  # first fit 
  first_fit_ns = optim_nelson_siegel(term_structure = term_structure, maturity = maturity, params = params.init)
  
  # simulations 
  sim_fit_ns = calibrate_nelson_siegel(first_fit_ns, n = n, params.min = params.min, params.max = params.max, contraint_threshold = contraint_threshold)
  
  # best parameters 
  df_optim_params = sim_fit_ns[which(sim_fit_ns$mse == min(sim_fit_ns$mse, na.rm = TRUE)),]
  
  
  # setting the title of the plot 
  if(!is.null(label) & is.character(label)){
    
    plot_title = paste0("Fitted Nelson-Siegel-Svensonn vs Real Value ", "(", label, ")" )
    
  } else {
    
    plot_title = "Fitted Nelson-Siegel-Svensonn vs Real Value"
  }
  
  
  # plot of fitted vs real values 
  plot_df =  dplyr::inner_join(
    dplyr::tibble(
      t = maturity,
      pred = df_optim_params$fit[[1]]
    ),
    dplyr::tibble(
      t = maturity,
      real = term_structure
    ),
    by = "t"
  ) 
  
  # Plot Real value vs Fitted Values 
  plot_ns = plot_df %>%
    mutate(label = paste0("T = ", round(t, 3)))%>%
    ggplot()+
    geom_point(aes(t, pred), color = "red", size = 2, alpha = 0.8) +
    geom_point(aes(t, real), color = "black", alpha = 0.5)+
    geom_line(aes(t, pred), color = "red", size = 1) +
    geom_line(aes(t, real), color = "black", linetype = "dashed") +
    geom_label(aes(t+1, real-0.001, label = label), size = 1.5)+
    theme(axis.text.x = element_text(angle = 25, face = "bold", size = 7), 
          axis.text.y = element_text(face = "bold"), 
          axis.title  = element_text(face = "bold"),
          plot.title  = element_text(face = "bold"),
          plot.subtitle = element_text(face = "italic"),
          plot.caption = element_text(face = "italic"),
          panel.grid.major.x = element_line(colour="grey60", linetype="dotted"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(colour="grey60", linetype="dotted"),
          legend.text = element_text(face = "italic", size = 5),
          legend.title = element_text(face = "bold"),
          legend.position = "top" ) +
    scale_x_continuous(breaks=c(0, 0.5, 1,2,3,4,5,7,10,15, 20, 25, 30))+
    ggtitle(plot_title, subtitle = "Fitted Value in Red and Real Values in Black" )+
    xlab("Maturities")+
    ylab("") + 
    labs(caption = paste0("Mean Squared Error of the Fit: ", round(df_optim_params$mse, 6)))+
    theme_bw()
  
  
  
  # output 
  structure(
    list(
      optim_params = tibble( beta0=df_optim_params$optim_params[[1]][1],
                             beta1=df_optim_params$optim_params[[1]][2],
                             beta2=df_optim_params$optim_params[[1]][3],
                             beta3=df_optim_params$optim_params[[1]][4],
                             tau1=df_optim_params$optim_params[[1]][5],
                             tau2=df_optim_params$optim_params[[1]][6]),
      plot_ns = plot_ns,
      simulations = sim_fit_ns,
      df_optim = df_optim_params
      
    )
  )
  
}

# example

setwd("/Users/hwan/Desktop/Homepage/study_24spring")

uspar <- tibble()
uspar <- read_csv("investment_hw/usparyield.csv") %>% 
  arrange(Date) %>% 
  mutate(year=substr(Date,7,10),
         month=substr(Date,1,2)) %>% 
  group_by(year,month) %>% 
  slice(1) %>% 
  pivot_longer(cols = contains(" "),values_to = "par", names_to = "time") %>% 
  mutate(maturity=c(1,2,3,4,6,12,24,36,60,84,120,240,360)/12) %>% 
  ungroup()

ggplot(uspar,aes(x=maturity,y=par,colour=month))+
  geom_line()+
  theme_bw()

uspar_jan <- uspar %>% filter(month=="01")
uspar_feb <- uspar %>% filter(month=="02")
uspar_mar <- uspar %>% filter(month=="03")
uspar_apr <- uspar %>% filter(month=="04")
uspar_may <- uspar %>% filter(month=="05")

Jan <- optimal_params_ns(uspar_jan$par,uspar_jan$maturity,n=1000,contraint_threshold = 0.25)
Feb <- optimal_params_ns(uspar_feb$par,uspar_feb$maturity,n=1000,contraint_threshold = 0.25)
Mar <- optimal_params_ns(uspar_mar$par,uspar_mar$maturity,n=1000,contraint_threshold = 0.25)
Apr <- optimal_params_ns(uspar_apr$par,uspar_apr$maturity,n=1000,contraint_threshold = 0.25)
May <- optimal_params_ns(uspar_may$par,uspar_may$maturity,n=1000,contraint_threshold = 0.25)

Jan$plot_ns
Feb$plot_ns
Mar$plot_ns
Apr$plot_ns
May$plot_ns

optim_parameters <- Jan$optim_params %>% 
  union_all(Feb$optim_params) %>% 
  union_all(Mar$optim_params) %>% 
  union_all(Apr$optim_params) %>% 
  union_all(May$optim_params) %>% 
  mutate(month=c("Jan","Feb","Mar","Apr","May"))
optim_parameters

est_jun <- tibble(maturity=uspar_may$maturity,
                  par=May$df_optim$fit[[1]])

ggplot(est_jun,aes(x=maturity,y=par))+
  geom_line(colour="blue", size=1)+
  scale_x_continuous(breaks=c(0, 0.5, 1,2,3,4,5,7,10,15, 20, 25, 30))+
  ggtitle("Expected US par-yield curve in June 1, based on NSS in May 1.")+
  xlab("Maturities")+
  ylab("") + 
  theme_bw()
