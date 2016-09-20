###########################################################################
###########################################################################
###########################################################################
#
# 	  		FONCTIONS POUR ESTIMER UNE DENSITE  
#
###########################################################################
###########################################################################
###########################################################################






# Cette fonction trace le graphe de l'histogramme à m classes
# des observations x1,...,xN. L'argument x de la fonction est
# le vecteur qui contient les observations, tandis que m désigne
# le nombre de classes à utiliser. L'argument add indique si
# l'histogramme doit être superposé sur un graphe existant ou
# être tracé sur un nouveau repère orthogonal. 

histogram=function(x,m,add,...)
{
a=min(x)
b=max(x)
d=(b-a)
h=d/m
hatf=1:m
n=length(x)
for (j in 1:m){hatf[j]=sum(((j-1)*h<=x-a)*(x-a<j*h))/(n*h)}

xleft=a-h+(1:m)*h
xright=xleft+h
ybottom=(1:m)*0
ytop=hatf

if (add==F)
{
plot(c(a-h,xleft,b),c(0,hatf,0),type="n",xlab="Les classes", ylab="Estimateur de densité",...)
rect(xleft, ybottom, xright, ytop, col = "cyan", border = "darkblue", lwd = 1)
}
else {rect(xleft, ybottom, xright, ytop, col = "cyan", border = "darkblue", lwd = 1)}

}

###########################################################################
###########################################################################
###########################################################################
#
# 			VALIDATION CROISEE POUR LES HISTOGRAMMES 
#
###########################################################################
###########################################################################
###########################################################################

CV_hist=function(x)
{
n=length(x)
a=min(x)
b=max(x)

N=round(n/5);

b=b+(b-a)/N
m_CV=1

J0=2/(n-1)

J=1:N;

for (m in 2:(N+1))
	{
	h=(b-a)/m
	hatp=1:m
	A=(1:m)%*%t((1:n)*0+1)
	xx=((1:m)*0+1)%*%t((x-a)/h)
	hatp=rowSums(((A-1)<=xx)*(xx<A))/n
	J[m-1]=2-(n+1)*sum(hatp^2)
	remove(hatp)
	J[m-1]=J[m-1]/((n-1)*h)
	if (J[m-1]<J0) {m_CV=m; J0=J[m-1]}
	}
op=par(mfcol=c(1,2),pty="m",omi=c(0,0,0,0))

plot(2:(N+1),J,type='l',lwd=2,col='darkred',main='La courbe de la fonction de validation croisée',,xlab='nb de classes',ylab='CV')

h=(b-a)/m_CV
hatf=1:m_CV
n=length(x)
m=m_CV
for (j in 1:m_CV){hatf[j]=sum(((j-1)*h<=x-a)*(x-a<j*h))/(n*h)}
xleft=a-h+(1:m)*h
xright=xleft+h
ybottom=(1:m)*0
ytop=hatf
plot(c(a-h/n,xleft,b),c(0,hatf,0),type="n",xlab="Les classes", ylab="Estimateur de densité",main="Histogramme avec le nb de classes optimal")
rect(xleft, ybottom, xright, ytop, col = "cyan", border = "darkblue", lwd = 1)
par(op)

return(m_CV)
}


#x=galaxies
#CV_hist(x)

###########################################################################
###########################################################################
###########################################################################
#
#   				Estimateur à noyau  
#
###########################################################################
###########################################################################
###########################################################################

KernelEst=function(x,h,kern,...)
{
n=length(x)
a=min(x)
b=max(x)
a=a-(b-a)
b=b+(b-a)
tt=(a+(b-a)*(1:500)/500)%*%t((1:n)*0+1)
xx=((1:500)*0+1)%*%t(x)
z=(xx-tt)/h

if (kern=='EP')
	{
	hatf=(1/(n*h))*as.vector(rowSums(3*(1-z^2)*(abs(z)<=1)/4))
	}
if (kern=='Rect')
	{
	hatf=(1/(n*h))*as.vector(rowSums((abs(z)<=1)/2))
	}
if (kern=='Tri')
	{
	hatf=(1/(n*h))*as.vector(rowSums((1-abs(z))*(abs(z)<=1)))
	}

if (kern=='Gaus')
	{
	hatf=(1/(n*h))*as.vector(rowSums(dnorm(z,0,sd=1/2.2)))
	}
if (kern=='sinc')
	{
	hatf=(1/(pi*n*h))*as.vector(rowSums(sin(z)/z))
	}


plot(a+(b-a)*(1:500)/500,hatf,type='l',xlab='',ylab='',col='darkred',lwd=2,...)	
hatf
}

#x=rnorm(400)
#ff=KernelEst(x,0.8,"EP",ylim=c(-0.05,0.4))
#h=10*0.8
#curve(dnorm,col='blue',lwd=2,add=T)

###########################################################################
###########################################################################
###########################################################################
#
#   					VALIDATION CROISEE  
#
###########################################################################
###########################################################################
###########################################################################
CV_kern=function(x)
{
n=length(x)
a=min(x)
b=max(x)
b=b+(b-a)/n
h_CV=(b-a)/n

J0=20*dnorm(0)/(n*h_CV)
xx=((1:n)*0+1)%*%t((x))-t(((1:n)*0+1)%*%t((x)))

if (n<100)
	N=round(n/2)
else
	N=round(n/4)
end;

J=1:N;

for (m in 1:N)
	{
	h=m*(b-a)/N
	J[m]=2*dnorm(0)/(n*h)+(1/(n^2*h))*sum(dnorm(xx/h,0,sqrt(2))-2*dnorm(xx/h))
	if (J[m]<J0) {h_CV=h; J0=J[m]}
	}
op <- par(mfcol=c(1,2),pty="m",omi=c(0,0,0,0))

plot((1:N)*(b-a)/N,J,type='l',lwd=2,col='darkred',main='La courbe de la fonction de validation croisée',,xlab='fenêtre',ylab='CV')

n=length(x)
a=min(x)
b=max(x)
a=a-(b-a)
b=b+(b-a)
tt=(a+(b-a)*(1:500)/500)%*%t((1:n)*0+1)
xx=((1:500)*0+1)%*%t(x)
z=(xx-tt)/h_CV
hatf=(1/(n*h_CV))*as.vector(rowSums(dnorm(z)))
plot(a+(b-a)*(1:500)/500,hatf,type='l',col='darkred',lwd=2,main='Estimateur à noyau avec la fenêtre optimale',xlab="",ylab="")
par(op)

return(2.2*h_CV)
}


# CV_kern(rnorm(200))