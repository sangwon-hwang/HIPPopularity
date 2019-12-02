
    # make sure required packages are installed and loaded
    source("./code/requiredPackages.R")
    source("./code/functions-fitting-data.R")

    # functions
    # Evaluations
    forecast_MSE <- function(df01, df02) {
      mean((df01$magnitude[91:120] - df02$magnitude[91:120])^2)
    }  
    
    
    forecast_MAE <- function(df01, df02) {
      mean(abs(df01$magnitude[91:120] - df02))
    }
    
    
    forecast_accuracy <- function(df01, df02) {
      
      error_rate <- abs(df01$magnitude[91:120] - df02)/(df01$magnitude[91:120])
      temp_accuracy <- rep(1, times = 30) - error_rate
      temp_list <- which(0>temp_accuracy)
      n <- length(temp_list)
      
      for(i in 1:n){
        temp_index <- temp_list[i]
        temp_accuracy[temp_index] <- 0
      }
      
      temp_accuracy
    }
    
    
    forecast_score <- function(df01, df02, a1, a2) {
      distances <- df02 - df01$magnitude # estimated values - true values
      
      a1 <- a1
      a2 <- a2
      
      a1_list <- which(0>distances)
      a2_list <- which(0<=distances)
      
      a1_n <- length(a1_list)
      a2_n <- length(a2_list)
      
      exp_list_a1 <- c()
      exp_list_a2 <- c()
      
      for(i in 1:a1_n){
        temp_distance <- a1_list[i]
        exp_list_a1[i] <- exp(-temp_distance/a1) - 1
      }
      
      for(i in 1:a2_n){
        temp_distance <- a2_list[i]
        exp_list_a2[i] <- exp(temp_distance/a2) - 1
      }
      
      return(list(sum(exp_list_a1), sum(exp_list_a2)))
    }
    
    
    # data load
    # endogenous(views) config
    options(digits=14)
      # 감소구간
      observed_cl_dec <- read.csv(file='./data/1min_120_cl_times_3rd.csv')
      # 증가구간
      observed_cl_inc <- read.csv(file='./data/1min_120_cl_times_2nd.csv')
    # exogoneous(shares) config
      # 감소구간
      observed_cv_dec <- read.csv(file='./data/1min_120_cv_times_3rd.csv')
      # 증가구간
      observed_cv_inc <- read.csv(file='./data/1min_120_cv_times_2nd.csv')

    # list
      Lst <- list(1,2,3,4,5,6,7,8)
      HIP <- list(1,2,3,4,5,6,7,8)
    
    # share setup       
    # exogoneous event: case: decreasing
      # 0.20045197016606
      shares <- ( observed_cv_dec$magnitude + observed_cl_inc$magnitude * 1/17.1 ) 
      cor(observed_cl_dec$magnitude[1:120], shares[1:120], use='complete.obs', method='pearson')
      Lst[[1]] <- shares
      
      # 0.40864585590557
      shares <- ( observed_cv_dec$magnitude + observed_cl_inc$magnitude * 1/100 ) 
      cor(observed_cl_dec$magnitude[1:120], shares[1:120], use='complete.obs', method='pearson')
      Lst[[2]] <- shares
      
      # 0.60923578698778
      shares <- ( observed_cv_dec$magnitude + observed_cl_dec$magnitude * 1/20 ) 
      cor(observed_cl_dec$magnitude[1:120], shares[1:120], use='complete.obs', method='pearson')
      Lst[[3]] <- shares
      
      # 0.80069927612292
      shares <- ( observed_cv_dec$magnitude + observed_cl_dec$magnitude * 1/6.5 ) 
      cor(observed_cl_dec$magnitude[1:120], shares[1:120], use='complete.obs', method='pearson')
      Lst[[4]] <- shares
      
    # exogoneous event: case: increasing
      # 0.20299594548747
      shares <- ( observed_cv_inc$magnitude + observed_cl_dec$magnitude * 1/19 )
      cor(observed_cl_inc$magnitude[1:120], shares[1:120], use='complete.obs', method='pearson')
      Lst[[5]] <- shares
      
      # 0.40070615500896
      shares <- ( observed_cv_inc$magnitude + observed_cl_dec$magnitude * 1/230 )
      cor(observed_cl_inc$magnitude[1:120], shares[1:120], use='complete.obs', method='pearson')
      Lst[[6]] <- shares
      
      # 0.60048848461137    
      shares <- ( observed_cv_inc$magnitude + observed_cl_inc$magnitude * 1/28.9 ) 
      cor(observed_cl_inc$magnitude[1:120], shares[1:120], use='complete.obs', method='pearson')
      Lst[[7]] <- shares
      
      # 0.79396502231167
      shares <- ( observed_cv_inc$magnitude + observed_cl_inc$magnitude * 1/10 ) 
      cor(observed_cl_inc$magnitude[1:120], shares[1:120], use='complete.obs', method='pearson')
      Lst[[8]] <- shares
    
      
      for(i in 1:length(Lst)) {
      
      shares <- Lst[[i]]
      
      if( i/2 <= 2 ){
        views <- observed_cl_dec$magnitude[1:120]
      }
      else{
        views <- observed_cl_inc$magnitude[1:120]
      }
        
      # fitting  
      fitted_params <- fit_series(data_series = views[1:90], ext_infl = list(shares = shares[1:90]),
                                  lowerBound = c(gamma = 0, eta = 0, K = 0, beta = 0.1, c = -Inf, theta = 0, mu1 = 0),
                                  upperBound = c(gamma = Inf, eta = Inf, K = Inf, beta = 0.1, c = Inf, theta = Inf, mu1 = Inf) )
      
      # CIF 
      # based on the fitted parameters, use HIP to generate a fitted series
      fitted_counts <- generate_simulated_data(params = fitted_params$model$par, time = 89, ext_infl = list(shares = shares[1:120]) )$Count
      
      # forecast 
      forecasted_counts <- generate_simulated_data(params = fitted_params$model$par, time = 119, ext_infl = list(shares = shares[1:120]), prefix = views[1:90] )$Count[91:120]
      
      HIP[[i]] <- c(fitted_counts, forecasted_counts)
      
      rm(fitted_params)
      
      }

    hip_dec_0.2 <- data.frame(time=1:120, magnitude=HIP[[1]]); 	colnames(hip_dec_0.2) <- c("time","magnitude")
    hip_dec_0.4 <- data.frame(time=1:120, magnitude=HIP[[2]]); 	colnames(hip_dec_0.4) <- c("time","magnitude")
    hip_dec_0.6 <- data.frame(time=1:120, magnitude=HIP[[3]]); 	colnames(hip_dec_0.6) <- c("time","magnitude")
    hip_dec_0.8 <- data.frame(time=1:120, magnitude=HIP[[4]]); 	colnames(hip_dec_0.8) <- c("time","magnitude")

    hip_inc_0.2 <- data.frame(time=1:120, magnitude=HIP[[5]]); 	colnames(hip_inc_0.2) <- c("time","magnitude")
    hip_inc_0.4 <- data.frame(time=1:120, magnitude=HIP[[6]]); 	colnames(hip_inc_0.4) <- c("time","magnitude")
    hip_inc_0.6 <- data.frame(time=1:120, magnitude=HIP[[7]]); 	colnames(hip_inc_0.6) <- c("time","magnitude")
    hip_inc_0.8 <- data.frame(time=1:120, magnitude=HIP[[8]]); 	colnames(hip_inc_0.8) <- c("time","magnitude")
    
    # write files
    write.csv(hip_dec_0.2, './data/hip_dec_0.2.csv', row.names=TRUE)
    write.csv(hip_dec_0.4, './data/hip_dec_0.4.csv', row.names=TRUE)
    write.csv(hip_dec_0.6, './data/hip_dec_0.6.csv', row.names=TRUE)
    write.csv(hip_dec_0.8, './data/hip_dec_0.8.csv', row.names=TRUE)

    write.csv(hip_inc_0.2, './data/hip_inc_0.2.csv', row.names=TRUE)
    write.csv(hip_inc_0.4, './data/hip_inc_0.4.csv', row.names=TRUE)
    write.csv(hip_inc_0.6, './data/hip_inc_0.6.csv', row.names=TRUE)
    write.csv(hip_inc_0.8, './data/hip_inc_0.8.csv', row.names=TRUE)

