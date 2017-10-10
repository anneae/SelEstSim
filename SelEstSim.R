
##### SelEstSim.R ##############################################################
# By Anne Eaton
#
# The SelEstSim function performs the simulations for the paper and 
# saves the performance metrics that we report. 
#
# The function usage: SelEstSim(nsim, nmarkers, npatients, pwisecorr=0)
# where nsim is the number of simulations, nmarkers is the number of markers and
# must be >2 (the first two markers have a true effect and the other markers 
# have no true effect), npatients is the number of patients and must be 
# >nmarkers, and pwisecorr is the pairwise correlation between each pair of 
# marker levels. 
#
# The output is a nsim x 13 x (8+nmarkers) array, where the second dimesion 
# (columns) represents method (TrueVar, Full, Back, UFilt, LASSO_1SE, ENet_1SE, 
# Ridge_1SE, LASSO, ENet, Ridge, LASSO_ref, ENet_ref, Ridge_ref, and the third 
# dimension represents performance metric (NullModel, DidNotConv, CorrectModel,
# CalibSlope, Cindex, RMSE, P25PredProb, P75PredProb, Beta1, Beta2, Beta3, ...).
# NullModel, DidNotConv and Correct model are 0/1 variables where 1 represents 
# fitting a null model, trying to fit a model that did not converge, and 
# selecting the correct model (i.e. setting only beta1 and beta2 not equal to 
# zero).
#
# The methods are as described in the paper except the ENet method is fit using
# glmnet with alpha = 0.5 to be consistent with the other methods. (In the paper
# ENet referred to simultaneously optimizing over alpha and lambda with the 
# opt2D function.) 
# 
# Example - the first 3 scenarios in our paper, 100 simulations each: 
# fifty_patients<- SelEstSim(100, 12, 50)
# onehundred_patients<- SelEstSim(100, 12, 100)
# fivehundred_patients<- SelEstSim(100, 12, 500)
################################################################################

SelEstSim<-function(nsim, nmarkers, npatients, pwisecorr=0){
    require(glmnet)
    require(rms) 
    require(MASS) 
    require(pROC)   
    if (nmarkers<3) stop("nmarkers must be >2")
    if (npatients<nmarkers) stop("npatients must be <= nmarkers")
    if (pwisecorr<0 | pwisecorr>1) stop("pwisecorr must be between 0 and 1")

    method_names<-c('TrueVar','Full','Back','UFilt',
                    'LASSO_1SE','ENet_1SE','Ridge_1SE',
                    'LASSO','ENet','Ridge',
                    'LASSO_ref','ENet_ref','Ridge_ref')
    beta_names<-'Beta1'
    for (j in 2:nmarkers) beta_names<-c(beta_names, paste('Beta', j, sep=''))
    result_names<-c('NullModel', 'DidNotConv', 'CorrectModel',
                    'CalibSlope', 'Cindex', 'RMSE','P25PredProb','P75PredProb',
                    beta_names)
    retval <- array(data     = NA,
                    dim      = c(nsim, 13, 8+nmarkers),
                    dimnames = list(1:nsim, method_names,result_names))
    
    sigma<-matrix(pwisecorr, nrow = nmarkers, ncol=nmarkers)
    diag(sigma)<-1
    sigma<-sigma*(2.550817^2)
    
    test_data<-data.frame(mvrnorm(n=10000, mu=rep(4.9,nmarkers), Sigma = sigma))
    test_data$lp<- -3.5536+0.3774*test_data$X1+0.2504*test_data$X2
    test_data$outcome<-rbinom(10000, size = 1, 
                              prob = exp(test_data$lp)/(1+exp(test_data$lp)))
    test_data<-test_data[,-which(names(test_data) =='lp')]                   
    
    for (i in 1:nsim){
        dat<-data.frame(mvrnorm(n=npatients, mu=rep(4.9,nmarkers), Sigma=sigma))
        dat$lp<- -3.5536+0.3774*dat$X1+0.2504*dat$X2
        dat$outcome<-rbinom(npatients,size=1, prob= exp(dat$lp)/(1+exp(dat$lp)))
        dat<-dat[,-which(names(dat) =='lp')]
        
        # True variable model
        fit <- tryCatch(glm(outcome~ X1+X2, data = dat, family = binomial),
                        error=function(e) e, warning=function(w) w)
        if (is(fit, 'warning')) retval[i,1,1:2] <-c(F,T)
        else if (max(abs(fit$coefficients[2:3])<10^-14)) retval[i,1,1:2]<-c(T,F)
        else {
            pred<-predict(fit,as.data.frame(test_data[,1:2]),type = 'link')
            retval[i, 1, 1:3]<-c(F,F,T)
            retval[i, 1, 4:8]<-calc_perf(pred, test_data)
            retval[i, 1, 9:10]<-fit$coefficients[2:3]
            retval[i, 1, 11:(8+nmarkers)]<-0
        }
        # full model
        fit <- tryCatch(glm(outcome~., data = dat, family = binomial),
                        error=function(e) e, warning=function(w) w)
        if (is(fit, 'warning')) retval[i,2,1:2] <-c(F,T)
        else if (max(abs(fit$coefficients[-1])<10^-14)) retval[i,2,1:2]<-c(T,F)
        else {
            pred<-predict(fit,as.data.frame(test_data[,1:nmarkers]),type='link')
            retval[i, 2, 1:3]<-c(F,F,F)
            retval[i, 2, 4:8]<-calc_perf(pred, test_data)
            retval[i, 2, 9:(8+nmarkers)]<-fit$coefficients[2:(nmarkers+1)]
        }
        #Stepwise selection  
        if (is(fit, 'warning')) retval[i,3,1:2] <-c(F,T)
        else {
            fit<-step(fit, direction = "backward", trace = 0)
            if (length(fit$coefficients)==0 |
                max(abs(fit$coefficients[-1]))<10^-14) retval[i,3,1:2]<-c(T,F)
            else{ pred<-predict(fit,as.data.frame(test_data[,1:nmarkers]),
                              type='link')
                retval[i, 3, 1:3]<-c(F,F,F)
                retval[i, 3, 4:8]<-calc_perf(pred, test_data)
                retval[i, 3, 9:(8+nmarkers)]<-0
                for (j in 1:nmarkers) { 
                    pos<-which(names(fit$coefficients)==paste('X',j,sep=''))
                    if (length(pos)>0) retval[i,3,8+j]<-fit$coefficients[pos]
                }
                if (length(names(fit$coefficients)) == 3)
                    retval[i,3,3]<-min(names(fit$coefficients)[2:3]
                                       ==c('X1','X2'))
            }
        }
        # univariate filtering
        sig_on_uni<-rep(F, nmarkers)
        for (j in 1:nmarkers){ sig_on_uni[j]<-summary(glm(dat$outcome~ dat[,j], 
                                       family = binomial))$coefficients[2,4]<.05
        }
        if (sum(sig_on_uni)==0) retval[i, 4, 1:3]<-c(T,F, F)
        else { reduced_dat<-dat[,c(sig_on_uni, T)]
            fit <- tryCatch(glm(outcome~ ., data = reduced_dat, 
                family = binomial), error=function(e) e, warning=function(w) w)
            if (is(fit, 'warning')) retval[i,4,1:2] <-c(F,T)
            else {pred<-predict(fit,as.data.frame(test_data[,1:nmarkers]),
                                type='link')
                retval[i, 4, 1:3]<-c(F,F,F)
                retval[i, 4, 4:8]<-calc_perf(pred, test_data)
                retval[i, 4, 9:(8+nmarkers)]<-0
                for (j in 1:nmarkers) { 
                    pos<-which(names(fit$coefficients)==paste('X',j,sep=''))
                    if (length(pos)>0) retval[i,4,8+j]<-fit$coefficients[pos]
                }
                if (length(names(fit$coefficients)) == 3)
                    retval[i,4,3]<-min(names(fit$coefficients)[2:3]
                                   ==c('X1','X2'))
            }
        }
        # LASSO 
        res<-enet(dat, test_data, 1, '1se', refit= T)
        retval[i, 5, 1:(8+nmarkers)] <- res[[1]]
        retval[i, 11, 1:(8+nmarkers)] <- res[[2]]
        retval[i, 8, 1:(8+nmarkers)] <- enet(dat, test_data, 1, 'min')
        # ENET 
        res <- enet(dat, test_data, .5, '1se', refit= T)
        retval[i, 6, 1:(8+nmarkers)] <- res[[1]]
        retval[i, 12, 1:(8+nmarkers)] <- res[[2]]
        retval[i, 9, 1:(8+nmarkers)] <- enet(dat, test_data, .5, 'min')
        # Ridgelike
        res <- enet(dat, test_data, .2, '1se', refit= T)
        retval[i, 7, 1:(8+nmarkers)] <- res[[1]]
        retval[i, 13, 1:(8+nmarkers)] <- res[[2]]
        retval[i, 10, 1:(8+nmarkers)] <- enet(dat, test_data, .2, 'min')
    }
    return(retval)
}
 
# Does penalized regression (with or without refitting) and grabs the numbers 
enet<-function(dat, test_data, alpha, lambda.choice, refit = F){
    nmarkers<-dim(dat)[2]-1
    enet.return<-rep(NA, (8+nmarkers))
    cvfit <- cv.glmnet(as.matrix(dat[,1:nmarkers]),dat$outcome, 
                       family='binomial', alpha = alpha) 
    if (lambda.choice == '1se') sel.lam <- cvfit$lambda.1se
    if (lambda.choice == 'min') sel.lam <- cvfit$lambda.min
    fit<-glmnet(as.matrix(dat[,1:nmarkers]),dat$outcome, 
                family='binomial', alpha = alpha, lambda = sel.lam) 
    if (max(abs(fit$beta))<=10^-15) {enet.return[1:2]<-c(T,F)
        refit.return <- enet.return}
    else {enet.return[1:3]<-c(F, F, F)
        if (fit$beta[1]>0&fit$beta[2]>0&sum(fit$beta[3:nmarkers])==0) 
            enet.return[3]<-T
        pred<-predict(cvfit,  as.matrix(test_data[,1:nmarkers]),
                  s = sel.lam, type = "link")   
        enet.return[4:8]<-calc_perf(pred, test_data)
        enet.return[9:(8+nmarkers)]<-fit$beta[1:nmarkers]
        if (refit == T) {refit.return<-rep(NA, (8+nmarkers))
            refit.return[1:3]<-enet.return[1:3]
            dat_reduced<-dat[,c((fit$beta!=0)[1:nmarkers],T)]
            refi <- tryCatch(glm(outcome~ ., data = dat_reduced, 
                family = binomial), error=function(e) e, warning=function(w) w)
            if (is(refi, 'warning')) refit.return[2] <-T
            else {pred<-predict(refi,as.data.frame(test_data[,1:nmarkers]),
                                type='link')
                refit.return[4:8]<-calc_perf(pred, test_data)
                refit.return[9:(8+nmarkers)]<-0
                for (j in 1:nmarkers) {
                    pos<-which(names(refi$coefficients)==paste('X',j,sep=''))
                    if (length(pos)>0) refit.return[8+j]<-refi$coefficients[pos]
                }
            }
        }
    }
    if (refit == F) return(enet.return)
    if (refit == T) return(list(enet.return, refit.return))
}

# Takes predictions (one per test dataset patient, linear predictor scale)
# and calculates calib slope, c index, RMSE, and IQRs of pred probabilities
calc_perf<-function(pred, test_data){ 
    pred2<- exp(pred)/(1+exp(pred))
    return(c(lrm(test_data$outcome~pred)$coefficients[2],
        auc(roc(as.numeric(test_data$outcome)~as.numeric(pred))),
        sqrt(sum((as.numeric(test_data$outcome)-as.numeric(pred2))^2)/10000),
        quantile(pred2, probs = c(.25, .75))))
}