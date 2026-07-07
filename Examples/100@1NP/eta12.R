## On Leonardo?:
# module load profile/bioinf
# module load R
## to install, in R:

#Lorentzian function
lorentzian <- function(A,lambda,omega)
{
 A*lambda/(lambda*lambda+omega*omega)
}

#Loss function
loss <- function(A,lambda,omega)
{
 A*omega/(lambda*lambda+omega*omega)
}

#x11()


df <- scan("tau1.out", list("a",0,0))
A1 <- as.numeric(unlist(df[2]))
lambda1 <- as.numeric(unlist(df[3]))
nf1=length(A1)-1
df <- scan("tau2.out", list("a",0,0))
A2 <- as.numeric(unlist(df[2]))
lambda2 <- as.numeric(unlist(df[3]))
nf2=length(A2)-1

Gamma=1000*lambda1[1]/2.
nstep=10000
step=Gamma/nstep
omega <- seq(step,Gamma,by=step)
#plot(omega,lorentzian(A[1],lambda[1],omega),xlim=c(0,Gamma))
#Sys.sleep(5)

omega0=step
sum=0.
sum1=0.
for(i in 1:nf1)
{
# sum = sum+loss(A[i],lambda[i],omega0)
 sum = sum+lorentzian(A1[i],lambda1[i],omega0)
 sum1 = sum1 + A1[i]/lambda1[i]
#  cat(x[j],y[j],"\n")
}
y0_1 <- sum
y1_1 <- sum1
sum=0.
sum1=0.
for(i in 1:nf2)
{
# sum = sum+loss(A[i],lambda[i],omega0)
 sum = sum+lorentzian(A2[i],lambda2[i],omega0)
 sum1 = sum1 + A2[i]/lambda2[i]
#  cat(x[j],y[j],"\n")
}
y0_2 <- sum
y1_2 <- sum1

#y1 = y1/y2
#print(y0_1)
#print(y1_1)
#print(y0_2)
#print(y1_2)

y2_1=sum(A1)
y2_2=sum(A2)
A1 = A1/y2_1
A2 = A2/y2_2
y0_1=sum(A1)
y0_2=sum(A2)
print(y0_1)
print(y0_2)
WW1=694329146.466667
WW2=1436360531.266000
alpha=WW1/WW2
w1=1.
w2=1.
A1 = w1*A1
A2 = w2*A2

x <- rep(0.,nstep)
y <- rep(0.,nstep)
for(j in seq_along(omega))
{
# x[j] = log(omega[j])
 x[j] = omega[j]
}
for(j in seq_along(omega))
{
 sum=0.
 for(i in 1:nf1) 
 {
# y[j] <- log(lorentzian(A[1],lambda[1],omega[j]))
#  sum = sum+lorentzian(A[i],lambda[i],omega[j])
  sum = sum+loss(A1[i],lambda1[i],omega[j])
  y[j] <- sum
#  cat(x[j],y[j],"\n")
 }
 for(i in 1:nf2) 
 {
# y[j] <- log(lorentzian(A[1],lambda[1],omega[j]))
#  sum = sum+lorentzian(A[i],lambda[i],omega[j])
  sum = sum+loss(A2[i],lambda2[i],omega[j])
  y[j] <- sum
#  cat(x[j],y[j],"\n")
 }
}

#plot.new()
#plot.window(xlim=c(0,1), ylim=c(5,10))
#abline(a=6, b=3)
#axis(1)
#axis(2)
#title(main="eta_R(omega)")
#title(xlab="omega")
#title(ylab="eta")
#box()

 out <- data.frame(x,y)
 write.table(out, file = "plot.dat", append = FALSE, quote = FALSE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)

# To exit from R:
#q()
