setwd("/home/houl/20240120/202401200828-SNP")
# setwd("/home/jihanbing/MRTP-20240120")
#setwd("C://Users//Lenovo//Desktop//异质人群MR//20231227")

library(readr)
library(MASS)
library(doParallel)
library(foreach)
library(parallel)
library(doMC)
library(foreach)
# library(AnchorRegression)
library(TwoSampleMR)
library(MendelianRandomization)

##########data generation##########
DataGeneration <- function(g,b1,b2,rho_b,siga1,siga2){
  
  
  # betaX1G <- c(rnorm(5,0.2,0.05),abs(rnorm(95,0,0.01)))
  # betaX2G <- c(abs(rnorm(10,0,0.01)),rnorm(90,0.2,0.03))
  
  betaX1G <- rnorm(g,0.2,0.05)
  betaX2G <- rnorm(g,0.2,0.03)
  
  seseby1 <- runif(g,0.01,0.13)
  seseby2 <- runif(g,0.0002,0.01)
  
  sigma_b1 <- seseby1/betaX1G
  sigma_b2 <- seseby2/betaX2G
  
  gamma1 <- runif(g,siga1,siga2)
  gamma2 <- runif(g,siga1,siga2)
  
  betaY1G <- NULL
  betaY2G <- NULL
  for(j in 1:g){
    
    betaY <- mvrnorm(1,mu = c(b1*betaX1G[j]+gamma1[j],
                              b2*betaX2G[j]+gamma2[j]),
                     Sigma = matrix(c(seseby1[j]^2,
                                      seseby1[j]*seseby2[j]*rho_b,
                                      seseby1[j]*seseby2[j]*rho_b,
                                      seseby2[j]^2),
                                    nrow=2))
    betaY1G <- c(betaY1G,betaY[1])
    betaY2G <- c(betaY2G,betaY[2])
  }
  
  beta1 <- betaY1G/betaX1G
  beta2 <- betaY2G/betaX2G
  betat1 <- (betaY1G - gamma1)/betaX1G
  betat2 <- (betaY2G - gamma2)/betaX2G
  
  
  data_sum <- data.frame(betaY1G=betaY1G,
                         betaY2G=betaY2G,
                         betaX1G=betaX1G,
                         betaX2G=betaX2G,
                         sebeta1=seseby1,
                         sebeta2=seseby2,
                         sebeta1mrtp=sigma_b1,
                         sebeta2mrtp=sigma_b2,
                         gamma1=gamma1,
                         gamma2=gamma2,
                         betat1=betat1,
                         betat2=betat2,
                         rho_b=rho_b)
  
  return(data_sum)
}

####bootstrap####
Bootstrap <- function(n_bootstrap,dataset,initial_mle){
  
  # 初始化自举结果存储
  bootstrap_results <- matrix(0, nrow = n_bootstrap, ncol = 2)  # 假设我们存储参数估计值和p值
  
  # 进行自举
  for (i in 1:n_bootstrap) {
    # 生成重采样数据集
    dataset_boot <- dataset[sample(nrow(dataset), size = 100, replace = TRUE), ]
    
    # 使用MLikelihood函数计算参数估计值
    res_mle <- MLikelihood(dataset_boot$betat1, dataset_boot$betat2, b2, dataset_boot$sebeta1, dataset_boot$sebeta2, rho_b, dataset_boot$gamma1, dataset_boot$gamma2)
    
    # 存储结果
    bootstrap_results[i, ] <- res_mle
  }
  
  # 计算置信区间
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
MLikelihood <- function(betat1,betat2,b2,sigma_b1,sigma_b2,rho_b,gamma1,gamma2){
  
  bb2=b2
  ggamma1=0
  ggamma2=0
  
  nll <- function(bb1){
    det_sigma = (1-rho_b^2)*(sigma_b1^2)
    pp = (betat1-bb1-ggamma1-rho_b*sigma_b1*(betat2-bb2-ggamma2)/sigma_b2)^2
    sum(log(2*pi)+0.5*log(det_sigma)+0.5*pp/det_sigma)
  }
  
  fit <- stats4::mle(minuslog=nll, 
                     start=list(bb1=0))
  res <- fit@coef
  
  fit1 <- stats4::mle(minuslog=nll, 
                      start=list(bb1=0),
                      fixed = list(bb1=0) )
  
  stat_beta1=  2 * (fit@min - fit1@min)
  pvalue_beta1 = pchisq(-stat_beta1,1,lower.tail=F)
  
  return(c(res,pvalue_beta1))
}


######COMP###########
Comp <- function(x,g,b1,b2,rho_b,siga1,siga2){
  cat('x=',x,'\n')
  dataset <- DataGeneration(g,b1,b2,rho_b,siga1,siga2)
  
  #####MRTP
  sigma_b1 <- dataset$sebeta1mrtp
  sigma_b2 <- dataset$sebeta2mrtp
  gamma1 <- dataset$gamma1
  gamma2 <- dataset$gamma2
  
  betat1 <- dataset$betat1
  betat2 <- dataset$betat2
  
  res_mle <- MLikelihood(betat1,betat2,b2,sigma_b1,sigma_b2,rho_b,gamma1,gamma2)
  res_mle
  
  coverage_TEMR <- Bootstrap(n_bootstrap,dataset,b1)
  
  ####Other MR
  obj <- MendelianRandomization::mr_input(bx=dataset$betaX1G,bxse=dataset$sebeta2,
                                          by=dataset$betaY1G,byse = dataset$sebeta1)
  mr_ivw_re <- MendelianRandomization::mr_ivw(obj)
  # mr_ivw_re <- mr_ivw(dataset$betaX1G,dataset$betaY1G,dataset$sebeta2,dataset$sebeta1)
  coverage_ivw <- coverCI(mr_ivw_re$Estimate, mr_ivw_re@`StdError`,b1)
  mr_simple_median_re <- mr_simple_median(dataset$betaX1G,dataset$betaY1G,dataset$sebeta2,dataset$sebeta1)
  coverage_simple_median <- coverCI(mr_simple_median_re$b, mr_simple_median_re$se,b1)
  mr_weighted_median_re <- mr_weighted_median(dataset$betaX1G,dataset$betaY1G,dataset$sebeta2,dataset$sebeta1)
  coverage_weighted_median <- coverCI(mr_weighted_median_re$b, mr_weighted_median_re$se,b1)
  mr_penalised_weighted_median_re <- mr_penalised_weighted_median(dataset$betaX1G,dataset$betaY1G,dataset$sebeta2,dataset$sebeta1)
  coverage_penalised_weighted_median <- coverCI(mr_penalised_weighted_median_re$b, mr_penalised_weighted_median_re$se,b1)
  
  MR_dat <- data.frame(
    beta.exposure=dataset$betaX1G,
    beta.outcome=dataset$betaY1G, 
    se.exposure=dataset$sebeta2,
    se.outcome=dataset$sebeta1
  )
  mr_mode_re <- mr_mode(MR_dat, mode_method = "all")
  simple_mode_coverage <- coverCI(mr_mode_re[1,5],mr_mode_re[1,6],b1)
  weighted_mode_coverage <- coverCI(mr_mode_re[2,5],mr_mode_re[2,6],b1)
  penalised_mode_coverage <- coverCI(mr_mode_re[3,5], mr_mode_re[3,6],b1)
  
  
  res_all <- c(mr_ivw_re@Estimate, mr_ivw_re@StdError, mr_ivw_re@Pvalue, coverage_ivw, res_mle, coverage_TEMR,
               mr_simple_median_re$b, mr_weighted_median_re$b, mr_penalised_weighted_median_re$b, 
               mr_simple_median_re$se, mr_weighted_median_re$se, mr_penalised_weighted_median_re$se,
               mr_simple_median_re$pval, mr_weighted_median_re$pval, mr_penalised_weighted_median_re$pval,
               coverage_simple_median,coverage_simple_median,coverage_penalised_weighted_median,
               mr_mode_re[1,5], mr_mode_re[2,5], mr_mode_re[3,5],
               mr_mode_re[1,6], mr_mode_re[2,6], mr_mode_re[3,6],
               mr_mode_re[1,9], mr_mode_re[2,9], mr_mode_re[3,9],
               simple_mode_coverage,weighted_mode_coverage,penalised_mode_coverage,
               g,b1,b2,rho_b,siga1,siga2)
  
  res_all <- matrix(res_all ,nrow=1)
  
  return(res_all)
}

#############MAIN############
main <- function(NN,x,g,b1,b2,rho_b,siga1,siga2,mc){
  
  registerDoMC(mc)
  
  tt1 <- foreach(x=1:NN,
                 .combine=rbind,
                 .errorhandling = "remove",
                 .packages = c("MASS")) %dopar% {
                   Comp(x,g,b1,b2,rho_b,siga1,siga2)
                 }
  
  # cl <- makeCluster(mc)
  # registerDoParallel(cl)
  # 
  # tt1 <- foreach(x=1:NN,
  #                .combine=rbind,
  #                .packages=c("MASS","AnchorRegression")) %dopar% {
  #                  Comp(x,betaX1G,betaX2G,g,b1,b2,rho_b,siga1,siga2)
  #                }
  # stopCluster(cl)
  
  
  return(tt1)
}


##########SIMULATION##########
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
            
            'g','b1','b2','rho_b',"siga1","siga2")

# result1 <- NULL
# for (x in 1:1000) {
# 
#   result <- Comp(x,betaX1G,betaX2G,g,b1,b2,rho_b,siga1,siga2)
#   result1 <- rbind(result1,result)
# }
# result2 <- data.frame(result1)
# colnames(result2) <- ccname

NN=1000
mc=20
N1b=N3b=3000
N2b=300000
N4b=300000

n_bootstrap <- 1000
bb1 <- c(0,0.05)
bb2 <- c(0,0.05)
# bb1 <- c(0,0.05,0.1,0.15,0.2)
# bb2 <- c(0,0.05,0.1,0.15,0.2)
rho_bb <- c(0.1,0.9)
siga1b <- c(0,0)
siga2b <- c(0,0.01)
gb=c(25,50,200)


#m=n=t=p=q=p=1
result_mean <- NULL
for(p in 1:length(siga1b)){
  for(n in 1:length(bb1)){
      for(t in 1:length(rho_bb)){
        for (q in 1:length(gb)) {
          b1=bb1[n]
          b2=bb2[n]
          rho_b=rho_bb[t]
          siga1=siga1b[p]
          siga2=siga2b[p]
          g=gb[q]
          result <- main(NN,x,g,b1,b2,rho_b,siga1,siga2,mc)
          #result <- Comp(x,betaX1G,betaX2G,g,b1,b2,rho_b,siga1,siga2)
          # result <- NULL
          # for (x in 1:1000) {
          #   
          #   result1 <- Comp(x,betaX1G,betaX2G,g,b1,b2,rho_b,siga1,siga2)
          #   result <- rbind(result1,result)
          # }
          result <- as.data.frame(result)
          colnames(result) <- ccname
          write.csv(result,paste0("MC","_N2_",N2b,"_N4_",N4b,"_g_",g,
                                  "_beta1_",b1,"_beta2_",b2,"_rho_b_",rho_b,
                                  "_siga1_",siga1,"_siga2_",siga2,"_0830.csv"))
      } 
    }
  }
}
