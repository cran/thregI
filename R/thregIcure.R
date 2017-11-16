## The following is the definition of the thregIcure function for interval-censored data with cure propolation
"thregIcure" <-
function (formula,data)
{
 ### Read in all arguments
 cl <- match.call()
 indx <- match(c("formula", "data"), names(cl), nomatch=0)
 if (indx[1] ==0) stop("A formula argument is required")
 mf<- cl[c(1, indx)]
 f <- Formula(formula)
 f1<-formula(f,lhs=1)
 f1<-Formula(f1)
 mf[[1]] <- as.name("model.frame")
 mf$formula <- if(missing(data)) terms(f1) else terms(f1, data=data)
 mf$formula <- f1
 mf <- eval(mf, parent.frame())
 if (nrow(mf) ==0) stop("No (non-missing) observations")
 Terms <- terms(mf)
 Y <- model.extract(mf, "response")
 if (!inherits(Y, "Surv")) stop("Response must be a survival object")
 type <- attr(Y, "type")
 if (type!='interval') stop(paste("thregI package only support \"", type, "\" survival data", sep=''))
 f2<-formula(f1, lhs = 0)
 if (length(f2[[2]])!=3) stop(paste("Predictors for both lny0 and mu should be specified"))
 x_lny0<-model.matrix(f1, data, rhs=1)
 x_mu<-model.matrix(f1, data, rhs=2)
 x_lamda<-model.matrix(f1, data, rhs=3)
 k<-dim(x_lamda)[2]
 left<-Y[,1]
 right<-Y[,2]
 delta<-Y[,3]  #delta=1:event, delta=3:left/interval delta=0:right
               #no delta=2 for left=NA since our dateset no NA
 ####################################
 #(0,R]---->delta=3, left=0, right=R
 #(a,b]---->delta=3, left=a, right=b
 #(L,Inf]-->delta=0, left=L, right=1
 #(a,a]---->delta=1, left=a, right=1
 ###################################
 delta3=matrix(0,length(left),1)
 # fix the exact times .............................................................
 for (i in 1 :length(left))
 {
   if (delta[i]==1) {right[i]=left[i]  #right=a instead of 1
                     delta3[i]=1}      #delta3=1 for exact time
 }
 #right_max=10*max(right[right!=Inf])
 right_max=10*max(right)
 delta1=matrix(0,length(left),1)
 delta2=matrix(0,length(left),1)
 for (i in 1 :length(left)){
   #if (right[i]==1) { right[i]=right_max }
   if (delta[i]==0) {right[i]=right_max}
   else if (left[i]==0) { delta1[i]=1 }   #delta1=1 for left censoring
   else if (delta[i]==3){ delta2[i]=1 }   #delta2=1 for interval censoring
 }                                        #delta1=delta2=delta3=0 for right censoring
####################################
#(0,R]---->delta1=1, left=0, right=R
#(a,b]---->delta2=1, left=a, right=b
#(L,Inf]-->delta1=delta2=delta3=0, left=L, right=Inf
#(a,a]---->delta3=1, left=a, right=a
###################################

 lny0<-function(para_lny0){x_lny0%*%para_lny0}
 mu<-function(para_mu){x_mu%*%para_mu}
 lamda<-function(para_lamda){x_lamda%*%para_lamda}
 d<-function(para)
 {
   para_lny0=para[1:length(dimnames(x_lny0)[[2]])]
   para_mu=para[(length(dimnames(x_lny0)[[2]])+1):(length(dimnames(x_lny0)[[2]])+length(dimnames(x_mu)[[2]]))]
   -mu(para_mu)/exp(lny0(para_lny0))
 }
 v<-function(para)
 {
   para_lny0=para[1:length(dimnames(x_lny0)[[2]])]
   exp(-2*lny0(para_lny0))
 }
 su<-function(para){
   pnorm((1-d(para)*left)/sqrt(v(para)*left))-exp(2*d(para)/v(para))*pnorm(-(1+d(para)*left)/sqrt(v(para)*left))
 }

 sv<-function(para){
   pnorm((1-d(para)*right)/sqrt(v(para)*right))-exp(2*d(para)/v(para))*pnorm(-(1+d(para)*right)/sqrt(v(para)*right))
 }

 cu<-function(para){
   para_lamda=para[(length(dimnames(x_lny0)[[2]])+length(dimnames(x_mu)[[2]])+1):(length(dimnames(x_lny0)[[2]])+length(dimnames(x_mu)[[2]])+length(dimnames(x_lamda)[[2]]))]
   exp(lamda(para_lamda))/(1+exp(lamda(para_lamda)))
 }

 logdf<-function(para){
      -.5*(log(2*pi*v(para)*(right^3))+(d(para)*right-1)^2/(v(para)*right))
 }

 logf0<-function(para) {
 -sum(delta1*log(1-sv(para)))-sum(delta2*log(su(para)-sv(para)), na.rm = TRUE)-sum((1-delta1-delta2-delta3)*log(su(para)))-sum(delta3*logdf(para))
 }

 p0<-rep(0,(length(dimnames(x_lny0)[[2]])+length(dimnames(x_mu)[[2]])))
 est0<-nlm(logf0, p0, hessian = TRUE)
 loglik0 = (-1)*est0$minimum

 logf<-function(para) {
   #-sum(delta1*log(cu(para)*(1-sv(para))), na.rm = TRUE)-sum(delta2*log(cu(para)*(su(para)-sv(para))), na.rm = TRUE)-sum((1-delta1-delta2-delta3)*log((1-cu(para))+cu(para)*su(para)), na.rm = TRUE)-sum(delta3*logdf(para), na.rm = TRUE)-sum(delta3*log(cu(para)), na.rm = TRUE)
   -sum(delta1*log(cu(para)*(1-sv(para))))-sum(delta2*log(cu(para)*(su(para)-sv(para))), na.rm = TRUE)-sum((1-delta1-delta2-delta3)*log((1-cu(para))+cu(para)*su(para)))-sum(delta3*logdf(para))-sum(delta3*log(cu(para)))
 }
 p<-rep(0,(length(dimnames(x_lny0)[[2]])+length(dimnames(x_mu)[[2]])+length(dimnames(x_lamda)[[2]])))
 est<-nlm(logf, p, hessian = TRUE) #nlm: minimize function logf
 names(est$estimate) <-c(paste("lny0:",dimnames(x_lny0)[[2]]),paste("mu:",dimnames(x_mu)[[2]]),paste("logit(p):",dimnames(x_lamda)[[2]]))
 loglik = (-1)*est$minimum
 #Pvalue = 1-pchisq(2*(loglik-loglik0),df=1)
 #print(paste('Goodness of fit test: p-value =',Pvalue,sep=" "))

 fit<-list(coefficients  = est$estimate,
           var    = solve(est$hessian),
           loglik = loglik,
           AIC    = (-2)*loglik+2*(length(dimnames(x_lny0)[[2]])+length(dimnames(x_mu)[[2]])+length(dimnames(x_lamda)[[2]])),
           Pvalue = 1-pchisq(2*(loglik-loglik0),df=k),
           iter   = est$iterations,
           call   = cl,
           mf     = mf,
           lny0   = dimnames(x_lny0)[[2]],
           mu     = dimnames(x_mu)[[2]],
           lamda  = dimnames(x_lamda)[[2]])
 class(fit) <- 'thregIcure'
 fit
}
