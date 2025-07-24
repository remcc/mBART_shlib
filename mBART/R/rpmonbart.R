rpmonbart=function(
x.train,y.train, x.test=matrix(0.0,0,0),
sigest=NA, sigdf=3, sigquant=.90,
k=2.0,
power=.8, base=.25, #regular bart is 2,.95
sigmaf=NA,
lambda=NA,
fmean = mean(y.train),
ntree=200L,
ndpost=1000L, nskip=100L,
mgsize=50L,
nkeeptrain=ndpost,nkeeptest=ndpost,
nkeeptestmean=ndpost,
nkeeptreedraws=ndpost,
printevery=50,
seed=99L,
mc.cores=2L
)
{
   cat("***** in rpmonbart\n")
##   require(parallel)
##   require(Rcpp)

   p=ncol(x.train)
   n=nrow(x.train)

   #--------------------------------------------------
   #prior
   nu=sigdf
   if(is.na(lambda)) {
      if(is.na(sigest)) {
         if(p < n) {
            df = data.frame(x.train,y.train)
            lmf = lm(y.train~.,df)
            sigest = summary(lmf)$sigma
         } else {
            sigest = sd(y.train)
         }
      }
      qchi = qchisq(1.0-sigquant,nu)
      lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
   }

   if(is.na(sigmaf)) {
      tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree));
   } else {
      tau = sigmaf/sqrt(ntree)
   }
   #--------------------------------------------------
   RNGkind("L'Ecuyer-CMRG")
   set.seed(seed)
   mc.reset.stream()

   mc.ndpost <- (ndpost %/% mc.cores)
   cat("ndpost: ",ndpost,"\n")
   cat("mc.ndpost: " ,mc.ndpost,"\n")

   mc.nkeeptrain = (nkeeptrain %/% mc.cores)
   mc.nkeeptest = (nkeeptest %/% mc.cores)
   mc.nkeeptestmean = (nkeeptestmean %/% mc.cores)
   mc.nkeeptreedraws = (nkeeptreedraws %/% mc.cores)

   cat("mc.cores: ",mc.cores,"\n")
   for(i in 1:mc.cores) {
      mcparallel(
         monbart(x.train,y.train, x.test,
               sigest=NA,
               sigdf=nu,
               sigquant=sigquant,
               k=2,
               power=power,
               base=base,
               sigmaf=tau*sqrt(ntree),
               lambda=lambda,
               fmean=fmean,
               ntree=ntree,
               ndpost=mc.ndpost,
               nskip=nskip,
               mgsize=mgsize,
               nkeeptrain = mc.nkeeptrain,
               nkeeptest = mc.nkeeptest,
               nkeeptestmean = mc.nkeeptestmean,
               nkeeptreedraws = 0,
               printevery=printevery
         ),
         silent=(i!=1)
      )
   }

   post.list <- mccollect()
   cat("length of post.list: ", length(post.list),"\n")
   np = nrow(x.test)

   if(nkeeptrain) yhat.train.ret=post.list[[1]]$yhat.train
   ##burnL = vector("list",mc.cores)
   ##burnL[[1]] = post.list[[1]]$sigma[1:nskip]
   burnL = cbind(post.list[[1]]$sigma)
   ##sigma=c(post.list[[1]]$sigma[1:mc.ndpost+nskip])
   sigma = post.list[[1]]$sigma[-(1:nskip)]
   if(np) yhat.test.ret=post.list[[1]]$yhat.test
   yhat.train.mean = post.list[[1]]$yhat.train.mean
   if(np) yhat.test.mean = post.list[[1]]$yhat.test.mean
   if(mc.cores>1) {
   for(i in 2:mc.cores) {
      if(nkeeptrain) yhat.train.ret = rbind(yhat.train.ret,post.list[[i]]$yhat.train)
      if(np) yhat.test.ret = rbind(yhat.test.ret,post.list[[i]]$yhat.test)
      ##burnL[[i]]=post.list[[i]]$sigma[1:nskip]
      burnL=cbind(burnL, post.list[[i]]$sigma)
      sigma =c(sigma, post.list[[i]]$sigma[-(1:nskip)])
      ##sigma =c(sigma,post.list[[i]]$sigma[1:mc.ndpost+nskip])
      yhat.train.mean = yhat.train.mean + post.list[[i]]$yhat.train.mean
      if(np) yhat.test.mean = yhat.test.mean + post.list[[i]]$yhat.test.mean
   }
   }

   ret = list()
   if(nkeeptrain) ret$yhat.train = yhat.train.ret
   ret$yhat.train.mean = yhat.train.mean/mc.cores
   if(np) {
      if(nkeeptest) ret$yhat.test = yhat.test.ret
      if(nkeeptestmean) ret$yhat.test.mean = yhat.test.mean/mc.cores
   }
   ret$sigma.=sigma
   ret$sigma = burnL

   return(ret)
}
