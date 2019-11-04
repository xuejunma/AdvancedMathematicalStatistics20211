
y=rnorm(100,2,4)#�ھ�ֵΪ2������Ϊ4����̬�ֲ��г�ȡ100������
#############Metropolis-Hastings
exam8.5M=function(y,n,a,b){
  th=matrix(0,ncol=2,nrow=n)
  m=length(y)
  th[1,]=c(a,b)
  pth=matrix(0,ncol=2,nrow=n)
  accept=NULL 
  for(i in 1:(n-1)){
    Q<-runif(2,-0.2,0.2) # �����ܶ�U[-1,1] 
    pth[i+1,1]<-th[i,1]+Q[1] #��ֵ����ȡֵ
    pth[i+1,2]<-th[i,2]+Q[2] #�����ȡֵ
    U<-runif(1,0,1) #����ֵ
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
plot(fitm$mu,xlab="��������",ylab="��ֵ")#����ֵ����ͼ
plot(fitm$sig,xlab="��������",ylab="��׼��")
#�Ƚϱ�Ҷ˹��Ƶ��ѧ��
mean(fitm$mu[1000:5000])
mean(fitm$sig[1000:5000])
mean(y)##Ƶ��
sd(y)##Ƶ��
sum(fitm$accept)


########################Gibbs����
#��д����
exam8.5g=function(y,n,a,b){
  th=matrix(0,ncol=2,nrow=n)
  m=length(y)
  th[1,]=c(a,b)
  for(i in 2:n){
    th[i,1]=rnorm(1,mean(y),sqrt(th[i-1,2]/m))#��ʽ��1��
    th[i,2]=1/rgamma(1,m/2,sum((y-th[i,1])^2)/2)##��ʽ��2��
  }
  list(mu=th[,1],sig=(th[,2])^0.5)
}

#Ӧ��
fit=exam8.5g(y,5000,1,1)
par(mfrow=c(1,2))
plot(fit$mu,xlab="��������",ylab="��ֵ")#����ֵ����ͼ
plot(fit$sig,xlab="��������",ylab="��׼��")
#�Ƚϱ�Ҷ˹��Ƶ��ѧ��
mean(fit$mu[1000:5000])
mean(fit$sig[1000:5000])
mean(y)##Ƶ��
sd(y)##Ƶ��


#�Ƚ����߷���

par(mfrow=c(2,2))
plot(fitm$mu,xlab="��������",ylab="��ֵ",main="Metropolis-Hastings")
plot(fitm$sig,xlab="��������",ylab="��׼��",main="Metropolis-Hastings")
plot(fit$mu,xlab="��������",ylab="��ֵ",main="Gibbs")
plot(fit$sig,xlab="��������",ylab="��׼��",main="Gibbs")


######################################################################�Ա�logit����

#######################����
logitM=function(x,y,n,a,b){
  beta=matrix(0,ncol=2,nrow=n)
  m=length(y)
  beta[1,]=c(a,b)
  pbeta=matrix(0,ncol=2,nrow=n)
  accept=NULL 
  for(i in 1:(n-1)){
    Q<-runif(2,-1,1) # �����ܶ�U[-1,1] 
    pbeta[i+1,1]<-beta[i,1]+Q[1] #��ֵ����ȡֵ
    pbeta[i+1,2]<-beta[i,2]+Q[2] #�����ȡֵ
    U<-runif(1,0,1) #����ֵ
    cur.eta<-beta[i,1]+beta[i,2]*x
    prop.eta<-pbeta[i+1,1]+pbeta[i+1,2]*x
    aa=(sum(y*prop.eta-log(1+exp(prop.eta)))-sum(y*cur.eta-log(1+exp(cur.eta)))
        +sum(dnorm(pbeta[i+1,1], 0,10,log=TRUE))
        +sum(dnorm(pbeta[i+1,2], 0,10,log=TRUE))
        -sum(dnorm(beta[i,1], 0,10,log=TRUE))
        -sum(dnorm(beta[i,2], 0,10,log=TRUE)))##ȡ������
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
plot(fitl$beta0,xlab="��������",ylab="beta0")
plot(fitl$beta1,xlab="��������",ylab="beta1")
#�Ƚϱ�Ҷ˹��Ƶ��ѧ��
mean(fitl$beta0[2000:8000])#2.4040
mean(fitl$beta1[2000:8000])# -0.3235 
sum(fitl$accept)###������ƫ�Ͳ���


  
####�뾭��ıȽ�һ��
fitc=glm(y~x,family=binomial(link="logit"))
summary(fitc)

