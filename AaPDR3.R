#Package: solve ode
library("deSolve")
#Define model without transporter 
model <- function (t, x, params) {
  P <- x[1]
  M <- x[2]
  
  r <- params["r"]
  C <- params["C"]
  d <- params["d"]
  k <- params["k"]
  
  
  
  dPdt <- (r-d*M)*P*((C-P)/C) 
  dMdt <- k*P
  dxdt <- c(dPdt, dMdt)
  list(dxdt)
}

params <- c(r=0.31271, C=100, d=0.001, k=0.0252)#h g
xstart<- c(P=2*10^-12, M=0) #2pg

times<-seq(from=0, to=120, by=1)
out<-ode(
  func = model,
  y = xstart,
  times = times,
  parms = params
)
out.df<-data.frame(out)
#Santalol 220.356 g/mol
###################################################################
#withtransporter
model2 <- function (t, x, params) {
  P <- x[1]
  M <- x[2]
  E <- x[3]
  
  r <- params["r"]
  C <- params["C"]
  d <- params["d"]
  k <- params["k"]
  Km <- params["Km"]
  Vmas <- params["Vmas"]
  
  
  #Michaelis Menten equation MATE
  #dE/dt=Vmax*[M]/(Km+[M])
  dPdt <- (r-d*M)*P*((C-P)/C) 
  dMdt <- k*P - (Vmas*M)/(Km+M)
  dEdt <- (Vmas*M)/(Km+M) #metabolite export by transporter
  
  dxdt <- c(dPdt, dMdt, dEdt)
  list(dxdt)
}
params <- c(r=0.31271, C=100, d=0.001, k=0.0252,Km=1.398*10^-8 ,Vmas=1.069*10^-6)
xstart<- c(P=2*10^-12, M=0, E=0)

times<-seq(from=0, to=120, by=1)
out<-ode(
  func = model2,
  y = xstart,
  times = times,
  parms = params
)
out.df1<-data.frame(out)
out.df1$Total<-out.df1$M+out.df1$E

biomassdf<-data.frame(out.df$time,out.df$P,out.df1$P)
colnames(biomassdf)<-c("Time","WithoutTransporter","WithTransporter")
library(ggplot2)
p<-ggplot(biomassdf,aes(x=Time))+
  geom_line(aes(y=WithoutTransporter,color="Without transporter"),linewidth=1)+
  geom_line(aes(y=WithTransporter,color="With transporter"),linewidth=1)+
  scale_colour_manual("", 
                      breaks = c("Without transporter", "With transporter"),
                      values = c("#D81B60", "#2224f0")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  theme_classic()+
  xlab("Time(h)")+
  ylab("Biomass(g)")+
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.text.y =  element_text(size = 12),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 12))
p
#####################################################
santaloldf<-data.frame(out.df$time,out.df$M,out.df1$Total)
colnames(santaloldf)<-c("Time","WithoutTransporter","WithTransporter")
s<-ggplot(santaloldf,aes(x=Time))+
  geom_line(aes(y=WithoutTransporter,color="Without transporter"),linewidth=1)+
  geom_line(aes(y=WithTransporter,color="With transporter"),linewidth=1)+
  scale_colour_manual("", 
                      breaks = c("Without transporter", "With transporter"),
                      values = c("#D81B60", "#2224f0")) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  theme_classic()+
  xlab("Time(h)")+
  ylab("Total Santalol(g)")+
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.text.y =  element_text(size = 12),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 12))
s


