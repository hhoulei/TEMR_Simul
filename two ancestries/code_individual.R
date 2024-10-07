setwd("/home/houl/20240120/202401200925")

library(readr)
library(MASS)
library(doParallel)
library(foreach)
library(parallel)
library(doMC)
library(foreach)
library(TwoSampleMR)
library(MendelianRandomization)

Wald_Ratio <- function(betaXG,betaYG,sebetaXG,sebetaYG){
  WR <- betaYG/betaXG
  varWR_2 <- (betaYG^2)*(sebetaXG^2)/(betaXG^4)+
    (sebetaYG^2)/(betaXG^2)
  
  return(data.frame(WR=WR,
                    varWR_2=varWR_2))
}

####bootstrap####
Bootstrap <- function(n_bootstrap,beta2b,rho_b,dataset,initial_mle){
  
  # 初始化自举结果存储
  bootstrap_results <- matrix(0, nrow = n_bootstrap, ncol = 2)  # 假设我们存储参数估计值和p值
  
  # 进行自举
  for (i in 1:n_bootstrap) {
    # 生成重采样数据集
    dataset_boot <- dataset[sample(length(dataset$beta1), size = 100, replace = TRUE), ]
    
    # 使用MLikelihood函数计算参数估计值
    res_mle <- MLikelihood(dataset_boot$beta1, dataset_boot$beta2, beta2b, 
                           dataset_boot$se1, dataset_boot$se2, rho_b)
    
    # 存储结果
    bootstrap_results[i, ] <- res_mle
  }
  
  se <- sd(bootstrap_results[,1])
  
  # 计算置信区间
  confidence_interval <- quantile(bootstrap_results[,1], probs = c(0.025, 0.975))
  if(initial_mle >= confidence_interval[1] && initial_mle <= confidence_interval[2]){
    coverage_rate <- 1
  }else{
    coverage_rate <- 0
  }
  
  return(c(se,coverage_rate))
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
MLikelihood <- function(betat1,betat2,beta2b,sigma_b1,sigma_b2,rho){
  
  bb2=beta2b
  ggamma1=0
  ggamma2=0
  
  nll <- function(bb1){
    det_sigma = (1-rho^2)*(sigma_b1^2)
    pp = (betat1-bb1-ggamma1-rho*sigma_b1*(betat2-bb2-ggamma2)/sigma_b2)^2
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

DataGenerator_wald <- function(N1,N2,g,p1,p2,beta1b,beta2b,rho_b,sigma_x1b,sigma_x2b){
  
  G1 <- NULL
  for(i in 1:g){
    Gg <- rbinom(N1,2,p1)
    G1 <- cbind(G1,Gg)
  }
  
  G2 <- NULL
  for(i in 1:g){
    Gg <- rbinom(N2,2,p2)
    G2 <- cbind(G2,Gg)
  }
  
  betaX1G <- rnorm(g,0.2,0.05)
  betaX2G <- rnorm(g,0.2,0.03)
  
  U1 <- rnorm(N1,0,0.05)
  U2 <- rnorm(N2,0,0.05)
  
  X1 <- G1 %*% betaX1G + rnorm(N1,0,0.15) + U1
  X2 <- G2 %*% betaX2G + rnorm(N2,0,0.15) + U2
  
  sigma_beta <- matrix(c(sigma_x1b^2,rho_b*sigma_x1b*sigma_x2b,
                         rho_b*sigma_x1b*sigma_x2b,sigma_x2b^2),nrow=2)
  beta <- mvrnorm(n=g, c(beta1b, beta2b),  sigma_beta)
  beta1 <- beta[,1]
  beta2 <- beta[,2]
  
  Y1 <- G1 %*% (betaX1G * beta1) + rnorm(N1,0,0.01) + U1 
  Y2 <- G2 %*% (betaX2G * beta2) + rnorm(N1,0,0.01) + U1 
  
  re_beta <- mean(sd(X1)/sd(Y1)*beta1)
  datasimul <- list(G1=G1,
                    G2=G2,
                    X1=X1,
                    X2=X2,
                    Y1=Y1,
                    Y2=Y2,
                    re_beta=re_beta)
  return(datasimul)
}

Comp_wald <- function(N1b,N2b,gb,p1b,p2b,beta1b,beta2b,
                      siga1,siga2,rho_b,sigma_x1b,sigma_x2b){
  
  fdata1 <- DataGenerator_wald(N1b,N2b,gb,p1b,p2b,beta1b,beta2b,rho_b,sigma_x1b,sigma_x2b)
  
  betaX1G <- NULL
  betaX2G <- NULL
  betaY1G <- NULL
  betaY2G <- NULL
  seX1G <- NULL
  seX2G <- NULL
  seY1G <- NULL
  seY2G <- NULL
  for(j in 1:gb){
    FF1 <- lm(fdata1$X1~fdata1$G1[,j])
    betaX1G <- c(betaX1G,FF1$coef[2])
    seX1G <- c(seX1G,summary(FF1)$coef[2,2])
    FF2 <- lm(fdata1$X2~fdata1$G2[,j])
    betaX2G <- c(betaX2G,FF2$coef[2])
    seX2G <- c(seX2G,summary(FF2)$coef[2,2])
    FF3 <- lm(fdata1$Y1~fdata1$G1[,j])
    betaY1G <- c(betaY1G,FF3$coef[2])
    seY1G <- c(seY1G,summary(FF3)$coef[2,2])
    FF4 <- lm(fdata1$Y2~fdata1$G2[,j])
    betaY2G <- c(betaY2G,FF4$coef[2])
    seY2G <- c(seY2G,summary(FF4)$coef[2,2])
  }
  
  res_wald1 <- Wald_Ratio(betaX1G,
                          betaY1G,
                          seX1G,
                          seY1G)
  res_wald2 <- Wald_Ratio(betaX2G,
                          betaY2G,
                          seX2G,
                          seY2G)
  
  rho_w <- cor(res_wald1$WR,res_wald2$WR)
  rho_w
  re_beta <- sd(betaX1G)/sd(betaY1G)*beta1b
                #fdata1$re_beta
  
  gamma1 <- runif(gb,siga1,siga2)
  gamma2 <- runif(gb,siga1,siga2)
  
  betaY1G <- betaY1G + gamma1
  betaY2G <- betaY2G + gamma2
  
  return(data.frame(beta1=res_wald1$WR,
                    se1=sqrt(res_wald1$varWR_2),
                    beta2=res_wald2$WR,
                    se2=sqrt(res_wald2$varWR_2),
                    betaX1G=betaX1G,
                    betaY1G=betaY1G,
                    seX1G=seX1G,
                    seY1G=seY1G,
                    seX2G=seX2G,
                    seY2G=seY2G,
                    rho_w=rho_w,
                    re_beta=re_beta))
}


######COMP###########
Comp <- function(x,N1b,N2b,gb,p1b,p2b,beta1b,beta2b,rho_b,siga1,siga2,sigma_x1b,sigma_x2b){
  cat('x=',x,'\n')
  dataset <- Comp_wald(N1b,N2b,gb,p1b,p2b,beta1b,beta2b,
                       siga1,siga2,rho_b,sigma_x1b,sigma_x2b)
  re_beta <- mean(dataset$re_beta)
  
  sigma_b1 <- dataset$se1
  sigma_b2 <- dataset$se2
  
  #####MRTP
  
  betat1 <- dataset$beta1
  betat2 <- dataset$beta2
  
  # hist(sigma_b1[!sigma_b1>1])
  # hist(sigma_b2)
  # hist(betat1)
  # hist(betat2)
  # 
  # res_mle <- MLikelihood(betat1,betat2,beta2b,sigma_b1,
  #                        sigma_b2,mean(dataset$rho_w))
  res_mle <- MLikelihood(betat1,betat2,beta2b,rep(sd(betat1),gb),
                         rep(sd(betat2),gb),mean(dataset$rho_w))
  res_mle
  
  coverage_TEMR <- Bootstrap(n_bootstrap,beta2b,mean(dataset$rho_w),dataset,beta1b)
  coverage_TEMR
  
  ####Other MR
  obj <- MendelianRandomization::mr_input(bx=dataset$betaX1G,bxse=dataset$seX1G,
                                          by=dataset$betaY1G,byse=dataset$seY1G)
  mr_ivw_re <- MendelianRandomization::mr_ivw(obj)
  coverage_ivw <- coverCI(mr_ivw_re$Estimate, mr_ivw_re@`StdError`,beta1b)
  c(mr_ivw_re@Estimate, mr_ivw_re@StdError, mr_ivw_re@Pvalue)
  
  mr_simple_median_re <- mr_simple_median(dataset$betaX1G,dataset$betaY1G,
                                          dataset$seX1G,dataset$seY1G)
  coverage_simple_median <- coverCI(mr_simple_median_re$b, mr_simple_median_re$se,beta1b)
  mr_weighted_median_re <- mr_weighted_median(dataset$betaX1G,dataset$betaY1G,
                                              dataset$seX1G,dataset$seY1G)
  coverage_weighted_median <- coverCI(mr_weighted_median_re$b, mr_weighted_median_re$se,beta1b)
  mr_penalised_weighted_median_re <- mr_penalised_weighted_median(dataset$betaX1G,dataset$betaY1G,
                                                                  dataset$seX1G,dataset$seY1G)
  coverage_penalised_weighted_median <- coverCI(mr_penalised_weighted_median_re$b, mr_penalised_weighted_median_re$se,beta1b)
  
  # mr_simple_mode_re <-MendelianRandomization::mr_mbe(obj,weighting = "unweighted")
  # coverage_simple_mode <- coverCI(mr_simple_mode_re@Estimate, mr_simple_mode_re@StdError,beta1b)
  # mr_weighted_mode_re <-MendelianRandomization::mr_mbe(obj,weighting = "weighted")
  # coverage_weighted_mode <- coverCI(mr_weighted_mode_re@Estimate, mr_weighted_mode_re@StdError,beta1b)
  # mr_penalised_weighted_mode_re <- MendelianRandomization::mr_mbe(obj,weighting = "weighted",stderror = "delta")
  # coverage_penalised_weighted_mode <- coverCI(mr_penalised_weighted_mode_re@Estimate, mr_penalised_weighted_mode_re@StdError,beta1b)
  
  MR_dat <- data.frame(
    beta.exposure=dataset$betaX1G,
    beta.outcome=dataset$betaY1G, 
    se.exposure=dataset$seX1G,
    se.outcome=dataset$seY1G
  )
  mr_mode_re <- mr_mode(MR_dat, mode_method = "all")
  simple_mode_coverage <- coverCI(mr_mode_re[1,5],mr_mode_re[1,6],beta1b)
  weighted_mode_coverage <- coverCI(mr_mode_re[2,5],mr_mode_re[2,6],beta1b)
  penalised_mode_coverage <- coverCI(mr_mode_re[3,5], mr_mode_re[3,6],beta1b)
  
  
  res_all <- c(mr_ivw_re@Estimate, mr_ivw_re@StdError, mr_ivw_re@Pvalue, coverage_ivw, 
               res_mle, coverage_TEMR,
               mr_simple_median_re$b, mr_weighted_median_re$b, mr_penalised_weighted_median_re$b, 
               mr_simple_median_re$se, mr_weighted_median_re$se, mr_penalised_weighted_median_re$se,
               mr_simple_median_re$pval, mr_weighted_median_re$pval, mr_penalised_weighted_median_re$pval,
               coverage_simple_median,coverage_simple_median,coverage_penalised_weighted_median,
               # mr_simple_mode_re@Estimate, mr_weighted_mode_re$Estimate, mr_penalised_weighted_mode_re$Estimate, 
               # mr_simple_mode_re$StdError, mr_weighted_mode_re$StdError, mr_penalised_weighted_mode_re$StdError,
               # mr_simple_mode_re$Pvalue, mr_weighted_mode_re$Pvalue, mr_penalised_weighted_mode_re$Pvalue,
               mr_mode_re[1,5], mr_mode_re[2,5], mr_mode_re[3,5],
               mr_mode_re[1,6], mr_mode_re[2,6], mr_mode_re[3,6],
               mr_mode_re[1,9], mr_mode_re[2,9], mr_mode_re[3,9],
               simple_mode_coverage,weighted_mode_coverage,penalised_mode_coverage,
               gb,beta1b,beta2b,rho_b,mean(dataset$rho_w),siga1,siga2,
               sigma_x1b,sigma_x2b,re_beta)
  res_all
  res_all <- matrix(res_all ,nrow=1)
  res_all <- as.data.frame(res_all)
  colnames(res_all) <- ccname
  res_all
  
  write_csv(res_all,file =paste0("MC","_N1_",N1b,"_N2_",N2b,"_g_",gb,
                                 "_beta1_",beta1b,"_beta2_",beta2b,"_rho_b_",rho_b,
                                 "_siga1_",siga1,"_siga2_",siga2,"_sigma_x1b_",sigma_x1b,
                                 "_sigma_x2b_",sigma_x2b,"_0928-5-b.csv"), append = TRUE)
  
  # return(res_all)
}


#############MAIN############
main <- function(NN,x,N1b,N2b,gb,p1b,p2b,beta1b,beta2b,rho_b,siga1,siga2,sigma_x1b,sigma_x2b,mc){
  
  registerDoMC(mc)
  
  foreach(x=1:NN,
          .combine=rbind,
          .packages = c("MASS")) %dopar% {
            Comp(x,N1b,N2b,gb,p1b,p2b,beta1b,beta2b,rho_b,siga1,siga2,sigma_x1b,sigma_x2b)
          }
}

##########SIMULATION##########
ccname <- c('IVW_beta1','IVW_se1','IVW_pval1','IVW_coverage',
            
            'TEMR_beta1','TEMR_pvalue_beta1','TEMR_se1','TEMR_coverage',
            
            'simple_median_beta1','weighted_median_beta1','penalised_median_beta1',
            'simple_median_se1','weighted_median_se1','penalised_median_se1',
            'simple_median_pval1','weighted_median_pval1','penalised_median_pval1',
            'simple_median_coverage','weighted_median_coverage','penalised_median_coverage',
            
            'simple_mode_beta1','weighted_mode_beta1','penalised_mode_beta1',
            'simple_mode_se1','weighted_mode_se1','penalised_mode_se1',
            'simple_mode_pval1','weighted_mode_pval1','penalised_mode_pval1',
            'simple_mode_coverage','weighted_mode_coverage','penalised_mode_coverage',
            
            'g','b1','b2','rho_b',"rho_w","siga1","siga2","sigma_x1b","sigma_x2b","re_beta")



NN=1000
gb=100
mc=50
N2b=50000
p1b=0.2
p2b=0.3
n_bootstrap=1000
sigma_x1b=0.01
sigma_x2b=0.005

N1bb=c(1000,3000)
beta1bb <- c(0,0.002)
beta2bb <- c(0,0.002)
rho_bb <- c(0.99,0.5,0.4)
siga1b <- c(0,0)
siga2b <- c(0,0.002)

m=n=t=p=z=1
for(t in 1:length(rho_bb)){
  for(p in 1:length(siga1b)){
    for(n in 1:length(beta1bb)){
      for(m in 1:length(beta2bb)){
        for(z in 1:length(N1bb)){
          # if(!(n==1 && m==1)){
            N1b=N1bb[z]
            beta1b=beta1bb[n]
            beta2b=beta2bb[m]
            rho_b=rho_bb[t]
            siga1=siga1b[p]
            siga2=siga2b[p]
            main(NN,x,N1b,N2b,gb,p1b,p2b,beta1b,beta2b,rho_b,siga1,siga2,sigma_x1b,sigma_x2b,mc)
            
          # }
          # write_csv(result,file =paste0("MC","_N1_",N1b,"_N2_",N2b,"_g_",gb,
          #                         "_beta1_",beta1b,"_beta2_",beta2b,"_rho_b_",rho_b,
          #                         "_siga1_",siga1,"_siga2_",siga2,"_0926.csv"), append = TRUE)
          
        }
      }
    }
  }
}

