
#Taken from moveNT 
sim_mov<-function(type=c("2states", "OU"), npatches=5, ratio=5, nswitch=150, ncore=200,spacecore=200, seq_visit=sample(1:npatches, nswitch, replace=T),
                  stepDist= "gamma", angleDist = "vm",  stepPar = c(0.5,3,1,5), anglePar = c(pi,0,0.5,2), s=diag(40,2), grph=F) {
  
  coordx<-sample(seq(0,20,2), npatches, replace=F)*spacecore
  coordy<-sample(seq(0,20,2), npatches, replace=F)*spacecore
  nmig=ncore/ratio
  out<-data.frame()
  for (i in 1:(nswitch-1)){
    
    if(type=="2states") {
      core<-moveHMM::simData(nbAnimals=1,nbStates=2,stepDist=stepDist,angleDist=angleDist,stepPar=stepPar, anglePar=anglePar,zeroInflation=F,obsPerAnimal=ncore)
      corex<-core$x+coordx[seq_visit[i]]
      corey<-core$y+coordy[seq_visit[i]]
      Corri1<-rep(2, ncore)
    }
    
    if(type=="OU") {
      core<-adehabitatLT::simm.mou(date=1:ncore, b=c(coordx[seq_visit[i]],coordy[seq_visit[i]]), s=s)
      corex<-ld(core)$x
      corey<-ld(core)$y
      Corri1<-rep(2, ncore)
    }
    
    if(seq_visit[i] != seq_visit[i+1]) {
      mig<-adehabitatLT::simm.bb(date=1:nmig, begin=c(tail(corex,1), tail(corey,1)), end=rnorm(2, c(coordx[seq_visit[i+1]],coordy[seq_visit[i+1]]), sd=25))
      Corri2<-rep(1, nmig)
      sub<-cbind(c(corex, ld(mig)$x), c(corey, ld(mig)$y), c(Corri1, Corri2))
      
    }
    if(seq_visit[i] == seq_visit[i+1]) {
      sub<-cbind(corex, corey, Corri1)
      colnames(sub)<-c("V1", "V2", "V3")
    }
    out<-rbind(out, sub)
  }
  names(out)<-c("x", "y", "Corri")
  out<-adehabitatLT::as.ltraj(out[,1:2], as.POSIXct(1:nrow(out), origin = "1960-01-01", tz="GMT"), id="id", infolocs=data.frame(out$Corri))
  if(grph==T) {plot(out)}
  return(out)
}


#First main function, calculate time in and out of each pixels
traj2timing<-function(mov, res=100, grid=NULL) {
   mov<-adehabitatLT::ld(mov)
  mov[,13]<-1:nrow(mov)
  tt<-sp::SpatialPoints(mov[,1:2])
  tt1<-apply(coordinates(tt), 2, min)
  tt2<-apply(coordinates(tt), 2, max)
  if(is.null(grid)){ras<-raster::raster(xmn=floor(tt1[1])-res, ymn=floor(tt1[2])-res,xmx=ceiling(tt2[1])+res, ymx=ceiling(tt2[2])+res, res=res)}
  if(!is.null(grid)){ras<-raster::crop(grid, tt)}
  values(ras)<-1:ncell(ras)
  mov$pix_start<-raster::extract(ras,tt)
  
  timing<-list()
  timing[[ncell(ras)]]<-NA
  mm<-max(table(mov$pix_start))
  timing<-lapply(timing, function(x) data.frame(time_in=as.POSIXct(rep(NA,mm)), time_out=as.POSIXct(rep(NA,mm))))
  
  timing[[mov$pix_start[1]]]$time_in[2]<-mov$date[1] 
  
 for (i in 2:nrow(mov)) {
  #for (i in 2:30) {
  if(mov$pix_start[(i-1)]!=mov$pix_start[i]) {
     
  a<-sum(!is.na(timing[[mov$pix_start[i-1]]]$time_out))+2
  b<-sum(!is.na(timing[[mov$pix_start[i]]]$time_in))+2
  
  #Maybe improve here - right now, linear interpolation - assuming from the center 
  if(mov$dist[i-1]<=res) {fract<-0.5} 
  if(mov$dist[i-1]>res) {fract<-res/mov$dist[i-1]/2}
  
   timing[[mov$pix_start[i-1]]]$time_out[a]<-mov$date[i-1]+mov$dt[i-1]*fract #### Adjust here this is giving the whole step interval dt even if the patch
   timing[[mov$pix_start[i]]]$time_in[b]<-mov$date[i]-mov$dt[i-1]*fract ### Substract something based on the step length
  }
  }
  timing<-lapply(timing, na.omit)
  return(list(timing, ras))
}
  
### Intermediary functions 
      #add the time columns
      time_diff<-function(x, unit="secs") {
        x$total_time<-as.numeric(difftime(x$time_out, x$time_in, units=unit))
        return(x)}
      
      
      #add an interval 4th columns 
      time_interval<-function(x, unit="secs") {
        time2<-c(x$time_in[-1], NA)
        x$time_interval<-as.numeric(difftime(time2,x$time_out, units=unit))
        return(x)
      }
      
      #frequency (nrow)
      freq_visit<-function(x) {return(nrow(x))}
      
      #total time (sum 3rd column)
      total_duration<-function(x) {return(ifelse(nrow(x)>0, sum(x$total_time, na.rm=T), 0))}
      
      #avg duration (mean 3rd column)
      mean_duration<-function(x) {return(ifelse(nrow(x)>0, mean(x$total_time, na.rm=T), 0))}
      
      #cv duration (cv 3rd column)
      cv_duration<-function(x) {return(ifelse(nrow(x)>0, cv(x$total_time, na.rm=T), 0))}
      sd_duration<-function(x) {return(ifelse(nrow(x)>0, sd(x$total_time, na.rm=T), 0))}
      
      
      #avg interval
      mean_interval<-function(x) {return(ifelse(nrow(x)>0, mean(x$time_interval, na.rm=T), 0))}
      #cv interval 
      cv_interval<-function(x) {return(ifelse(nrow(x)>0, cv(x$time_interval, na.rm=T), 0))}
      
      #sd interval 
      sd_interval<-function(x) {return(ifelse(nrow(x)>0, sd(x$time_interval, na.rm=T), 0))}

#Second big function, use the timing and convert to raster with the different statistics 
timing2stack<-function(timing_ls, unit_time="secs") {
  timing<-timing_ls[[1]]
  ras<-timing_ls[[2]]
  timing<-lapply(timing, function(x) time_diff(x, unit=unit_time))
  timing<-lapply(timing, function(x) time_interval(x, unit=unit_time))
  grid<-stack(ras,ras,ras,ras,ras,ras,ras,ras)
  values(grid[[1]])<-unlist(lapply(timing, freq_visit))
  values(grid[[2]])<-unlist(lapply(timing, total_duration))
  values(grid[[3]])<-unlist(lapply(timing, mean_duration))
  values(grid[[4]])<-unlist(lapply(timing, cv_duration))
  values(grid[[5]])<-unlist(lapply(timing, sd_duration))
  values(grid[[6]])<-unlist(lapply(timing, mean_interval))
  values(grid[[7]])<-unlist(lapply(timing, cv_interval))
  values(grid[[8]])<-unlist(lapply(timing, sd_interval))
    names(grid)<- c("Number visits", "Total duration", "Mean duration", "CV duration", "SD duration", "Mean interval", "CV interval", "SD interval")
  return(grid)
  }

#Useful function that test different grid size and see how it influences the time spent in a pixel (peak on the figure prob best resolution).  
res_test<-function(mov, res_seq=c(50,100,150), unit_time="secs") {
  ls<-pbapply::pblapply(res_seq, function (x) traj2timing(mov, x))
  timing<-lapply(ls, function(x) x[[1]])
  timing<-lapply(1:length(timing), function(y) lapply(timing[[y]], function(x) time_diff(x, unit=unit_time)))
  rt<-lapply(1:length(timing), function(y) unlist(lapply(timing[[y]], total_duration)))
  rt_var<-lapply(rt, function(x) log(cv(x[x>0])))
  #rt_var<-lapply(rt, function(x) log(var(x)))
  plot(res_seq, rt_var, xlab="Resolution", ylab="CV RT", type="l")
}
  
# Add a function that could exclude very short excursion (like when an animal might just have left for a short time)? Maybe not required  

#### Function for clustering 

clust_use<-function(stck,col=c(1,2,3,4,6,7),nb_clust=1:5, min_fix=3 ) {
data<-values(stck)  
data<-data[,col]
data<-scale(data[data[,1]>min_fix,]) ## Need to be visited at least twice to be considered or 3 times for CV/SD interval included 
clust<-Mclust(data, G=nb_clust)
print(clust$parameters$mean)
ras<-stck[[1]]
values(ras)[values(ras)>min_fix]<-clust$classification
plot(ras)
out<-list(clust, ras)
return(out)
}



#Examples: 
library(adehabitatLT)
library(raster)
library(mclust)

#Simulate movement (or use own data based on adehabitatLT package - one individual at a time!)
mov<-sim_mov(type="OU", npatches=3, grph=T)  

# Test different resolution of grid - higher values probably better pixel size  
res_test(mov, res_seq=seq(20, 400, 40))  
  
#Extract timing of in and out of each pixel 
timing_ls<-traj2timing(mov, res=50, grid=NULL)

#Produce a stack imagery  - we need to figure if 
stck<-timing2stack(timing_ls) 
plot(stck)


#Clustering
test<-clust_use(stck)

### Example with actual locs 
data(puechcirc)
puechcirc
res_test(na.omit(puechcirc[1]), res_seq=seq(100, 400, 40))  

#Extract timing of in and out of each pixel 
timing_ls<-traj2timing(na.omit(puechcirc[1]), res=120, grid=NULL)

#Produce a stack imagery  - we need to figure if 
stck<-timing2stack(timing_ls) 
plot(stck)

#Clustering
test<-clust_use(stck) # One cluster dectected 
