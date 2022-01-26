psmc.result<-function(file,i.iteration=25,mu=1e-8,s=100,g=1)
{
  X<-scan(file=file,what="",sep="\n",quiet=TRUE)
  
  START<-grep("^RD",X)
  END<-grep("^//",X)
  
  X<-X[START[i.iteration+1]:END[i.iteration+1]]
  
  TR<-grep("^TR",X,value=TRUE)
  RS<-grep("^RS",X,value=TRUE)
  
  write(TR,"temp.psmc.result")
  theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
  N0<-theta0/4/mu/s
  
  write(RS,"temp.psmc.result")
  a<-read.table("temp.psmc.result")
  Generation<-as.numeric(2*N0*a[,3])
  Ne<-as.numeric(N0*a[,4])
  
  file.remove("temp.psmc.result")
  
  n.points<-length(Ne)
  YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
              Generation[n.points])*g
  Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
        Ne[n.points])
  
  data.frame(YearsAgo,Ne)
}

plotPsmc.allPops<-function(keywords, label, legend.names,
                           save.as="png", height=7, width=12,
                           #i.iteration=25, mu=2.5e-8, s=100, g=25,
                           #ylim=c(0,90000), xlim=c(7500,10000000),
                           i.iteration=25, mu=1.9e-10, s=100, g=0.0073,
                           ylim=c(0,12500000), xlim=c(275,750000),
                           col=rainbow(length(keywords)))
{
  n<-length(keywords)
  files<-grep("\\.psmc$",dir(),value=TRUE)
  #dev.new(height=height,width=width)
  par()
  plot(1,1,
       ylim=ylim,xlim=xlim,
       log="x",type="b",
       main=label[1],ylab="Ne",xlab="",
       bty="n", axes=FALSE)
  #axis(1, at=c(10000, 7500000), labels = c("recent", "ancient"))
  axis(1, at=c(1000, 10000, 100000), labels = c("1'000","10'000","100'000"))
  #axis(2, at=c(0, 75000), labels = c(0,"high"))
  axis(2, at=c(0, 10500000), labels = c(0,"high"))
  for(i in 1:n)
  {
    subfiles<-grep(paste("^",keywords[i],sep=""),files,value=TRUE)
    n.sub<-length(subfiles)
    for(i.sub in 1:n.sub)
    {
      lines(psmc.result(subfiles[i.sub],i.iteration,mu,s,g),
            type="l",col=col[i],lwd=2)
    }
  }
  legend("topright",legend=c(expression(italic("H. magnivora")), expression(italic("H. tvaerminnensis")* " - Northern"), expression(italic("H. tvaerminnensis")* " - Southern")),col=c("orangered2", "steelblue4", "steelblue1"),lty=1,lwd=2,bty = "n")
  #savePlot(filename=paste(label[1],save.as,sep="."),type=save.as)
  #dev.off()
}

setwd("~/pascal/Master/Master_publication/data/psmc")

png("Samples_PacBio.png", width=12, height=7, units="in", res=150)
plotPsmc.allPops(c("BE-OM-2.0.08.psmc",
                   "BE-T1-2.0.35.psmc",
                   "BE-WH2-11.0.34.psmc",
                   "CY-PA3-2.0.08.psmc",
                   "ES-DO1-1.0.18.psmc",
                   "ES-OM-5.0.14.psmc",
                   "FI-BR1-39-6.0.23.psmc",
                   "FI-FHS2-11-3.0.53.psmc",
                   "FI-OER-3-3.0.05.psmc",
                   "IL-G-3.0.07.psmc",
                   "IL-PS-6.0.07.psmc",
                   "IL-SH-4.0.08.psmc",
                   "RU-BAYA1-1.0.07.psmc",
                   "RU-TY1-3.0.13.psmc",
                   "RU-TY2-2.0.10.psmc",
                   "RU-ZB1-1.0.13.psmc",
                   "SE-G4-20.0.09.psmc",
                   "SE-H1-4.0.15.psmc"), label = "", legend.names = c("BE-OM-2",
                                                                     "BE-T1-2",
                                                                     "BE-WH2-11",
                                                                     "south",
                                                                     "south",
                                                                     "south",
                                                                     "north",
                                                                     "north",
                                                                     "north",
                                                                     "south",
                                                                     "south",
                                                                     "south",
                                                                     "north",
                                                                     "north - russia",
                                                                     "north - russia",
                                                                     "north - russia",
                                                                     "north",
                                                                     "north"), col = c("orangered2",
                                                                                       "orangered2",
                                                                                       "orangered2",
                                                                                       "steelblue1",
                                                                                       "steelblue1",
                                                                                       "steelblue1",
                                                                                       "steelblue4",
                                                                                       "steelblue4",
                                                                                       "steelblue4",
                                                                                       "steelblue1",
                                                                                       "steelblue1",
                                                                                       "steelblue1",
                                                                                       "steelblue4",
                                                                                       "steelblue4",
                                                                                       "steelblue4",
                                                                                       "steelblue4",
                                                                                       "steelblue4",
                                                                                       "steelblue4"))
dev.off()
