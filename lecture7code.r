
y=rnorm(100,2,4)#在均值为2，方差为4的正态分布中抽取100个样本
#############Metropolis-Hastings
exam8.5M=function(y,n,a,b){
  th=matrix(0,ncol=2,nrow=n)
  m=length(y)
  th[1,]=c(a,b)
  pth=matrix(0,ncol=2,nrow=n)
  accept=NULL 
  for(i in 1:(n-1)){
    Q<-runif(2,-0.2,0.2) # 提议密度U[-1,1] 
    pth[i+1,1]<-th[i,1]+Q[1] #均值建议取值
    pth[i+1,2]<-th[i,2]+Q[2] #方差建议取值
    U<-runif(1,0,1) #控制值
    aa=((th[i,2]/pth[i+1,2])^(m+2))*exp(-(sum((y-pth[i+1,1])^2))/(2*(pth[i+1,2])^2)+
                                          (sum((y-th[i,1])^2))/(2*(th[i,2])^2))
    alpha<-min(aa,1)
    if (U>=alpha){th[i+1,1]=th[i,1]
                  th[i+1,2]=th[i,2]
                  accept[i]=0} 
    if (U<alpha){th[i+1,1]= pth[i+1,1]
                 th[i+1,2]= pth[i+1,2]
                 accept[i]=1} 
  }
  list(mu=th[,1],sig=th[,2],accept=accept)
}

fitm=exam8.5M(y,5000,1,1)
par(mfrow=c(1,2))
plot(fitm$mu,xlab="迭代次数",ylab="均值")#画均值迭代图
plot(fitm$sig,xlab="迭代次数",ylab="标准差")
#比较贝叶斯和频率学派
mean(fitm$mu[1000:5000])
mean(fitm$sig[1000:5000])
mean(y)##频率
sd(y)##频率
sum(fitm$accept)


########################Gibbs抽样
#编写函数
exam8.5g=function(y,n,a,b){
  th=matrix(0,ncol=2,nrow=n)
  m=length(y)
  th[1,]=c(a,b)
  for(i in 2:n){
    th[i,1]=rnorm(1,mean(y),sqrt(th[i-1,2]/m))#公式（1）
    th[i,2]=1/rgamma(1,m/2,sum((y-th[i,1])^2)/2)##公式（2）
  }
  list(mu=th[,1],sig=(th[,2])^0.5)
}

#应用
fit=exam8.5g(y,5000,1,1)
par(mfrow=c(1,2))
plot(fit$mu,xlab="迭代次数",ylab="均值")#画均值迭代图
plot(fit$sig,xlab="迭代次数",ylab="标准差")
#比较贝叶斯和频率学派
mean(fit$mu[1000:5000])
mean(fit$sig[1000:5000])
mean(y)##频率
sd(y)##频率


#比较两者方法

par(mfrow=c(2,2))
plot(fitm$mu,xlab="迭代次数",ylab="均值",main="Metropolis-Hastings")
plot(fitm$sig,xlab="迭代次数",ylab="标准差",main="Metropolis-Hastings")
plot(fit$mu,xlab="迭代次数",ylab="均值",main="Gibbs")
plot(fit$sig,xlab="迭代次数",ylab="标准差",main="Gibbs")


######################################################################自编logit函数

#######################函数
logitM=function(x,y,n,a,b){
  beta=matrix(0,ncol=2,nrow=n)
  m=length(y)
  beta[1,]=c(a,b)
  pbeta=matrix(0,ncol=2,nrow=n)
  accept=NULL 
  for(i in 1:(n-1)){
    Q<-runif(2,-1,1) # 提议密度U[-1,1] 
    pbeta[i+1,1]<-beta[i,1]+Q[1] #均值建议取值
    pbeta[i+1,2]<-beta[i,2]+Q[2] #方差建议取值
    U<-runif(1,0,1) #控制值
    cur.eta<-beta[i,1]+beta[i,2]*x
    prop.eta<-pbeta[i+1,1]+pbeta[i+1,2]*x
    aa=(sum(y*prop.eta-log(1+exp(prop.eta)))-sum(y*cur.eta-log(1+exp(cur.eta)))
        +sum(dnorm(pbeta[i+1,1], 0,10,log=TRUE))
        +sum(dnorm(pbeta[i+1,2], 0,10,log=TRUE))
        -sum(dnorm(beta[i,1], 0,10,log=TRUE))
        -sum(dnorm(beta[i,2], 0,10,log=TRUE)))##取对数的
    aaa=exp(aa)    
    alpha<-min(aaa,1)
    if (U>=alpha){beta[i+1,1]=beta[i,1]
                  beta[i+1,2]=beta[i,2]
                  accept[i]=0} 
    if (U<alpha){beta[i+1,1]= pbeta[i+1,1]
                 beta[i+1,2]= pbeta[i+1,2]
                 accept[i]=1} 
  }
  list(beta0=beta[,1],beta1=beta[,2],accept=accept)
}

wais=read.table("D:\\wais.txt",header=TRUE)
x<-wais[,1]
y<-wais[,2]
fitl=logitM(x=wais[,1],y=wais[,2],8000,0,0)
par(mfrow=c(1,2))
plot(fitl$beta0,xlab="迭代次数",ylab="beta0")
plot(fitl$beta1,xlab="迭代次数",ylab="beta1")
#比较贝叶斯和频率学派
mean(fitl$beta0[2000:8000])#2.4040
mean(fitl$beta1[2000:8000])# -0.3235 
sum(fitl$accept)###接受率偏低不好


  
####与经典的比较一下
fitc=glm(y~x,family=binomial(link="logit"))
summary(fitc)

