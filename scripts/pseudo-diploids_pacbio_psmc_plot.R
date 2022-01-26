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
                           #ylim=c(0,100000), xlim=c(7500,10000000),
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
       main=label[1],ylab="Ne",xlab="Time in the past [y]               ",
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
  legend("topright",legend=legend.names,col=col,lty=1,lwd=2,bty = "n", y.intersp=0.5)
  #savePlot(filename=paste(label[1],save.as,sep="."),type=save.as)
  #dev.off()
}

setwd("~/pascal/Master/Master_publication/data/psmc/pseudo-diploid")


plotPsmc.allPops(c("BE-OM-2_BE-T1-2.1a2a.psmc",
                   "FI-BR1-39-6_IL-SH-4.1a2a.psmc",
                   "FI-BR1-39-6_RU-BAYA1-1.1a2a.psmc",
                   "IL-SH-4_ES-DO1-1.1a2a.psmc",
                   "RU-BAYA1-1_BE-OM-2.1a2a.psmc",
                   "ES-DO1-1_BE-OM-2.1a2a.psmc"
), label = "", legend.names = c(
  expression(italic("H.m.")* " & "* italic("H.m.")),
  expression(italic("H.t.")* "-N & "* italic("H.t.")* "-S"),
  expression(italic("H.t.")* "-N & "* italic("H.t.")* "-N"),
  expression(italic("H.t.")* "-S & "* italic("H.t.")* "-S"),
  expression(italic("H.t.")* "-N & "* italic("H.m.")), 
  expression(italic("H.t.")* "-S & "* italic("H.m."))),col = c("orangered2",
                                      "blue",
                                      "steelblue4",
                                      "steelblue1",
                                      "green",
                                      "springgreen4"))
