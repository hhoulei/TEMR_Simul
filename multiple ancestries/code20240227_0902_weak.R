library(MASS)
library(glmnet)
library(doMC)
library(doParallel)
library(foreach)
library(TwoSampleMR)

g=100

####data generation####
DataGeneration <- function(g,b,rhob,siga1,siga2,race) {
  
  betaXG <- data.frame(
    betaX1G = c(rnorm(5,0.2,0.05),abs(rnorm(95,0,0.01))),
    betaX2G = c(rnorm(5,0.2,0.05),abs(rnorm(95,0,0.01))),
    betaX3G = c(rnorm(5,0.2,0.05),abs(rnorm(95,0,0.01))),
    betaX4G = c(abs(rnorm(10,0,0.01)),rnorm(90,0.2,0.03))
  )
  
  seseby <- data.frame(
    seseby1 <- runif(g,0.01,0.13),
    seseby2 <- runif(g,0.02,0.12),
    seseby3 <- runif(g,0.02,0.25),
    seseby4 <- runif(g,0.0002,0.01)
  )
  
  sigma_b <- data.frame(
    sigma_b1 = seseby[,1]/betaXG[,1],
    sigma_b2 = seseby[,2]/betaXG[,2],
    sigma_b3 = seseby[,3]/betaXG[,3],
    sigma_b4 = seseby[,4]/betaXG[,4]
  )
  
  gamma <- data.frame(
    gamma1 = runif(g,siga1,siga2),
    gamma2 = runif(g,siga1,siga2),
    gamma3 = runif(g,siga1,siga2),
    gamma4 = runif(g,siga1,siga2))
  
  
  betaY <- matrix(0, nrow = g, ncol = race)
  sigma_matrices <- list() # 用于存储每个j循环生成的Sigma矩阵
  
  # 构建每对种族间的rho索引矩阵
  # 输入的rhob里面应该有6个元素
  rho_index <- 0
  rho_matrix <- matrix(0, nrow = race, ncol = race)
  for(r1 in 1:(race-1)){
    for(r2 in (r1+1):race){
      rho_index <- rho_index + 1
      rho_matrix[r1, r2] <- rhob[rho_index]
      rho_matrix[r2, r1] <- rhob[rho_index] # rho值是对称的
    }
  }
  
  for(j in 1:g){
    Sigma <- matrix(0, nrow = race, ncol = race)
    for(r1 in 1:race){
      for(r2 in 1:race){
        if(r1 == r2){
          Sigma[r1, r2] <- (seseby[j,r1])^2
        } else {
          Sigma[r1, r2] <- rho_matrix[r1, r2] * seseby[j,r1] * seseby[j,r2]
        }
      }
    }
    
    #cat('j=',j,'\n')
    betaY_j <- mvrnorm(1, mu = as.numeric(b*betaXG[j,]+gamma[j,]), Sigma = Sigma)
    betaY[j, ] <- betaY_j
    
    Sigma_mrtp <- matrix(0, nrow = race, ncol = race)
    for(r1 in 1:race){
      for(r2 in 1:race){
        if(r1 == r2){
          Sigma_mrtp[r1, r2] <- (sigma_b[j,r1])^2
        } else {
          Sigma_mrtp[r1, r2] <- rho_matrix[r1, r2] * sigma_b[j,r1] * sigma_b[j,r2]
        }
      }
    }
    
    sigma_matrices[[j]] <- Sigma_mrtp # 存储Sigma矩阵
  }
  
  beta <- (betaY/betaXG)
  betat <- (betaY - gamma)/betaXG
  
  
  # 动态构建data.frame
  data_list <- list()
  for(i in 1:ncol(beta)) {
    data_list[[paste0("betat", i)]] <- betat[, i]
    data_list[[paste0("sebetamrtp", i)]] <- sigma_b[,i]
    data_list[[paste0("sebeta", i)]] <- seseby[,i]
    data_list[[paste0("gamma", i)]] <- gamma[,i]
    data_list[[paste0("betaXG", i)]] <- betaXG[, i]
    data_list[[paste0("betaY", i)]] <- betaY[, i]
  }
  
  data_sum <- as.data.frame(data_list)
  
  return(list(data_sum = data_sum, sigma_matrices = sigma_matrices))
}


####bootstrap####
Bootstrap <- function(n_bootstrap,betat, betat_know, known_b, sigma_b_all, rhob, race, gb,initial_mle,b){
  
  # 初始化自举结果存储
  bootstrap_results <- matrix(0, nrow = n_bootstrap, ncol = 2)  # 假设我们存储参数估计值和p值
  
  # 进行自举
  for (i in 1:n_bootstrap) {
    success <- FALSE
    attempt <- 0
    
    while (!success && attempt < (maxRetries/2)) {
      attempt <- attempt + 1
      
      # 生成重采样数据集
      bootstrap_indices <- sample(length(sigma_b_all), size = 100, replace = TRUE)
      
      betat_n <- betat[bootstrap_indices]
      betat_know_n <- betat_know[bootstrap_indices,]
      sigma_b_all_n <- lapply(bootstrap_indices, function(idx) sigma_b_all[[idx]]) 
      
      # 使用MLikelihood函数计算参数估计值
      res_mle <- MLikelihood(betat_n, betat_know_n, known_b, sigma_b_all_n, rhob, race, gb,b)
      
      if (!is.null(res_mle)) {
        success <- TRUE
        # 存储结果
        bootstrap_results[i, ] <- res_mle}else {
          bootstrap_results[i, ] <- c(NA,NA)}
    }
  }
  
  
  # 计算置信区间
  bootstrap_results <- na.omit(bootstrap_results)
  confidence_interval <- quantile(bootstrap_results[,1], probs = c(0.025, 0.975))
  if(initial_mle >= confidence_interval[1] && initial_mle <= confidence_interval[2]){
    coverage_rate <- 1
  }else{
    coverage_rate <- 0
  }
  
  return(coverage_rate)
}

coverCI <- function(beta,se,initial_b){
  u_CI <- beta + se*1.96
  l_CI <- beta - se*1.96
  if(initial_b >= l_CI && initial_b <= u_CI){
    coverage_rate <- 1
  }else{
    coverage_rate <- 0
  }
  
  return(coverage_rate)
}

############### MRTP #############
MLikelihood <- function(beta, beta_know,known_b, sigma_b_all, rhob, race, gb, b){
  
  gamma_know <- matrix(rep(0,300),ncol = 100)
  gamma <- rep(0,100)
  bbe <- b[2:4]
  
  # 定义 nll 函数
  nll <- function(bbe1) {
    
    sigma_bb <- list()
    pp <- NULL
    det_sigma <- NULL
    for (i in 1:gb) {
      sigma_a0 <- as.numeric(sigma_b_all[[i]][1, 1])
      sigma_a1 <- sigma_b_all[[i]][2:race, 1]
      sigma_aa <- sigma_b_all[[i]][2:race, 2:race]
      inv_sigma <- ginv(sigma_aa)
      sigma_bb <- t(sigma_a1) %*% inv_sigma %*% sigma_a1
      det_sigma[i] <- as.numeric(sigma_a0-sigma_bb)
      
      diff <- beta[i] - bbe1 - gamma[i] - sigma_a1 %*% inv_sigma %*% t(beta_know[i,] - known_b - gamma_know[,i])
      pp[i] <- as.numeric(diff * diff / det_sigma[i] )
    }
    
    sum(gb * log(2 * pi) + 0.5 * log(abs(det_sigma)) + 0.5 * pp)
  }
  
  # 初始参数
  #bba <- rep(0, (race - 1))
  
  tryFirstMethod <- function() {
    
    fit <- stats4::mle(minuslog=nll, 
                       start=list(bbe1=0))
    res <- fit@coef
    
    fit1 <- stats4::mle(minuslog=nll, 
                        start=list(bbe1=0),
                        fixed = list(bbe1=0))
    
    stat_beta1=  2 * (fit@min - fit1@min)
    pvalue_beta1 = pchisq(-stat_beta1,1,lower.tail=F)
    
    return(c(res, pvalue_beta1))
  }
  
  result <- tryCatch(
    tryFirstMethod(),
    error = function(e) NULL  # 如果第二种方法也失败，则跳过
  )
  
  if (is.null(result)) {
    cat("Both methods failed for the current iteration.\n")
    return(NULL)  # 可能需要根据主函数的期望返回值进行调整
  }
  
  return(result)
}


######COMP###########
Comp <- function(x, gb,b,rhob,siga1,siga2,race, maxRetries) {
  cat('x=', x, '\n')
  
  success <- FALSE
  attempt <- 0
  
  while (!success && attempt < maxRetries) {
    attempt <- attempt + 1
    
    fdata <- DataGeneration(gb,b,rhob,siga1,siga2,race)
    
    betat <- t(fdata$data_sum[grep("^betat\\d+$", colnames(fdata$data_sum), value = TRUE)[1]])
    betat_know <- fdata$data_sum[grep("^betat\\d+$", colnames(fdata$data_sum), value = TRUE)[2:race]]
    sigma_b_all <- fdata$sigma_matrices
    known_b <- b[2:4]
    
    res_mle <- MLikelihood(betat, betat_know, known_b, sigma_b_all, rhob, race, gb,b)
    coverage_TEMR <- Bootstrap(n_bootstrap,betat, betat_know, known_b,
                               sigma_b_all, rhob, race, gb,b[1],b)
    
    if (!is.null(res_mle)) {
      success <- TRUE
      
      ####Other MR
      betaX1G <- t(fdata$data_sum[grep("^betaXG\\d+$", colnames(fdata$data_sum), value = TRUE)[-race]])
      betaX1G <- betaX1G[1,]
      betaY1G <- t(fdata$data_sum[grep("^betaY\\d+$", colnames(fdata$data_sum), value = TRUE)[-race]])
      betaY1G <- betaY1G[1,]
      sigma_b <- t(fdata$data_sum[grep("^sebeta\\d+$", colnames(fdata$data_sum), value = TRUE)[-race]])
      sigma_b1 <- sigma_b[1,]
      sigma_b2 <- sigma_b[2,]
      
      obj <- MendelianRandomization::mr_input(bx=betaX1G,bxse=sigma_b2,
                                              by=betaY1G,byse = sigma_b1)
      mr_ivw_re <- MendelianRandomization::mr_ivw(obj)
      mr_simple_median_re <- mr_simple_median(betaX1G,betaY1G,sigma_b2,sigma_b1)
      mr_weighted_median_re <- mr_weighted_median(betaX1G,betaY1G,sigma_b2,sigma_b1)
      mr_penalised_weighted_median_re <- mr_penalised_weighted_median(betaX1G,betaY1G,sigma_b2,sigma_b1)
      
      coverage_ivw <- coverCI(mr_ivw_re$Estimate, mr_ivw_re@`StdError`,b[1])
      coverage_simple_median <- coverCI(mr_simple_median_re$b, mr_simple_median_re$se,b[1])
      coverage_weighted_median <- coverCI(mr_weighted_median_re$b, mr_weighted_median_re$se,b[1])
      coverage_penalised_weighted_median <- coverCI(mr_penalised_weighted_median_re$b, 
                                                    mr_penalised_weighted_median_re$se,b[1])
      
      
      MR_dat <- data.frame(
        beta.exposure=betaX1G,
        beta.outcome=betaY1G, 
        se.exposure=sigma_b2,
        se.outcome=sigma_b1
      )
      mr_mode_re <- mr_mode(MR_dat, mode_method = "all")
      simple_mode_coverage <- coverCI(mr_mode_re[1,5],mr_mode_re[1,6],b[1])
      weighted_mode_coverage <- coverCI(mr_mode_re[2,5],mr_mode_re[2,6],b[1])
      penalised_mode_coverage <- coverCI(mr_mode_re[3,5], mr_mode_re[3,6],b[1])
      
      res_all <- c(mr_ivw_re@Estimate, mr_ivw_re@StdError, mr_ivw_re@Pvalue, coverage_ivw,res_mle, coverage_TEMR,
                   mr_simple_median_re$b, mr_weighted_median_re$b, mr_penalised_weighted_median_re$b, 
                   mr_simple_median_re$se, mr_weighted_median_re$se, mr_penalised_weighted_median_re$se,
                   mr_simple_median_re$pval, mr_weighted_median_re$pval, mr_penalised_weighted_median_re$pval,
                   coverage_simple_median,coverage_simple_median,coverage_penalised_weighted_median,
                   mr_mode_re[1,5], mr_mode_re[2,5], mr_mode_re[3,5],
                   mr_mode_re[1,6], mr_mode_re[2,6], mr_mode_re[3,6],
                   mr_mode_re[1,9], mr_mode_re[2,9], mr_mode_re[3,9],
                   simple_mode_coverage,weighted_mode_coverage,penalised_mode_coverage,
                   gb,b, rhob, race, siga1, siga2)
      
      res_all <- matrix(res_all ,nrow=1)
      
    } else {
      cat("Both methods failed for the current iteration, retrying...\n")
    }
  }
  
  if (!success) {
    return(NULL)  # Or handle the failure appropriately
  }
  
  return(res_all)
}


#############MAIN############
main <- function(NN,x,gb,b,rhob,siga1,siga2,race, maxRetries, mc){
  
  registerDoMC(mc)
  
  tt1 <- foreach(x=1:NN,
                 .combine=rbind,
                 .errorhandling = "remove",
                 .packages = c("MASS","glmnet","TwoSampleMR")) %dopar% {
                   Comp(x,gb,b,rhob,siga1,siga2,race, maxRetries)
                 }
  return(tt1)
}


##########SIMULATION##########
# 示例参数
# race = 4
# Nb <- c(100, 100, 100, 100) # N1-N5
# g=gb <- 10
# pb <- c(0.1, 0.2, 0.3, 0.4) # p1-p5
# betab <- c(1, 2, 3, 4) # beta1-beta5
# sigma_xb <- c(1, 2, 3, 4) # sigma_x1-sigma_x5
# bb <- rep(0,4)
# rhob <- seq(0.1,0.6,0.1)
# siga1=0
# siga2=0.2

#参数设置
NN=1000
maxRetries=10
n_bootstrap=100
mc=20
race = 4
gb <- 100

bb_all <- list(rep(0,race),
               rep(0.05,race)) # beta1-beta4相等
rho_bb <- list(rep(0.1,sum(1:(race-1))),
               rep(0.9,sum(1:(race-1))))

siga1b <- c(0)
siga2b <- c(0)

######## 构建ccname向量##########
# 定义beta和rho的个数
n_beta <- race - 1
n_rho <- sum(1:(race - 1))

# 构建b和rho的字符串向量
b_strings <- paste0("b", 1:race)
rho_strings <- paste0("rho", 1:n_rho)

ccname <- c('IVW_beta1','IVW_var1','IVW_pval1','IVW_coverage',
            
            'TEMR_beta1','TEMR_pvalue_beta1','TEMR_coverage',
            
            'simple_median_beta1','weighted_median_beta1','penalised_median_beta1',
            'simple_median_se1','weighted_median_se1','penalised_median_se1',
            'simple_median_pval1','weighted_median_pval1','penalised_median_pval1',
            'simple_median_coverage','weighted_median_coverage','penalised_median_coverage',
            
            'simple_mode_beta1','weighted_mode_beta1','penalised_mode_beta1',
            'simple_mode_se1','weighted_mode_se1','penalised_mode_se1',
            'simple_mode_pval1','weighted_mode_pval1','penalised_mode_pval1',
            'simple_mode_coverage','weighted_mode_coverage','penalised_mode_coverage',
            'g',b_strings, rho_strings, 'race', "siga1", "siga2")

# result <- NULL
# for (x in 1:100) {
#   res <- Comp(x,gb,b,rhob,siga1,siga2,race, maxRetries)
#   result <- rbind(result,res)
# }


#########遍历##########
m=n=t=p=q=1

for(m in 1:length(bb_all)){
  for(t in 1:length(rho_bb)){
    for(p in 1:length(siga1b)){
      b=bb_all[[m]]
      rhob=rho_bb[[t]]
      siga1=siga1b[p]
      siga2=siga2b[p]
      result <- main(NN,x, gb,b,rhob,siga1,siga2,race, maxRetries, mc)
      #result <- Comp(x,gb,b,rhob,siga1,siga2,race, maxRetries)
      result <- as.data.frame(result)
      colnames(result) <- ccname
      write.csv(result,paste0("MC","_g_",gb,"_bb_",b[1],"_rho_b_",rhob[6],
                              "_siga1_",siga1,"_siga2_",siga2,"_race_",race,
                              "_0902_weak.csv"))
      
    }
  }
} 
