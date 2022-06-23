ncd2<-function(m=m3, plot=TRUE){
  m[order(as.numeric(POS))]-> m
  windows_dt <- data.table(Start = m$POS[1], End=m$POS[nrow(m)], Mid=m$Mid[1], W=m$W[1], MaxMaf=m$MaxMaf[1], MidMaf=m$MidMaf[1],Sim=m$Sim[1], Truth=m$Truth[1], Stat=m$Stat[1])
  res<-bind_cols(windows_dt, m %>% summarise(N_Sites=n(), N_FD=sum(FD), N_Pol=sum(SegSite), IS=N_FD+N_Pol, NCD2=sqrt(sum((MAF-0.5)^2)/IS), NCD2_0.4=sqrt(sum((MAF-0.4)^2)/IS), NCD2_0.3=sqrt(sum((MAF-0.3)^2)/IS)))
  setDT(res)
  res[,NCD2_0.4_std:=((NCD2_0.4)/0.4)*0.5] #normalized version
  res[,NCD2_0.3_std:=((NCD2_0.3)/0.3)*0.5] #normalized version
  res[,Start:=NULL][, End:=NULL][, N_Pol:=NULL]
  if(plot == TRUE){
    m[,MAF:=round(MAF,2)]
    p<-SFS_func(X=m[SegSite==T], type="MAF")
    return(p)
  }else{
    return(res)  
  }
}

#f<-function(x, tf=0.5, nfd=2){sum(c((x-tf)^2, rep(0, nfd)))/(length(x)+nfd)}
#f2<-function(f)integrate(Vectorize(f), lower=min(x), upper=max(x))

ncd1<-function(m=m3, plot=TRUE){
  m[order(as.numeric(POS))]-> m
  m<-m[SegSite==T] #only SNPs
  windows_dt <- data.table(Start = m$POS[1], End=m$POS[nrow(m)], Mid=m$Mid[1], W=m$W[1], MaxMaf=m$MaxMaf[1], MidMaf=m$MidMaf[1],Sim=m$Sim[1], Truth=m$Truth[1], Stat=m$Stat[1])
  res<-bind_cols(windows_dt, m %>% summarise(N_Sites=n(), N_FD=0, N_Pol=sum(SegSite), IS=N_FD+N_Pol, NCD1=sqrt(sum((MAF-0.5)^2)/IS), NCD1_0.4=sqrt(sum((MAF-0.4)^2)/IS), NCD1_0.3=sqrt(sum((MAF-0.3)^2)/IS)))
  setDT(res)
  res[,NCD1_0.4_std:=((NCD1_0.4)/0.4)*0.5] #normalized version
  res[,NCD1_0.3_std:=((NCD1_0.3)/0.3)*0.5] #normalized version
  res[,Start:=NULL][, End:=NULL][, N_Pol:=NULL]
  if(plot == TRUE){
    m[,MAF:=round(MAF,2)]
    p<-SFS_func(X=m[SegSite==T], type="MAF")
    return(p)
  }else{
    return(res)  
  }
}
NCD2 <- function(X,  W = 100000, S = 50000, fd=T) {
  
  X[order(as.numeric(POS))]-> X
  windows_dt <-
    data.table(POS = seq(from=1, to=X[nrow(X), POS], S))[
      , POS2 := POS + W][
        -length(POS)]; #create dt with POS, then add POS2, then remove last row.
  
  setkey(windows_dt, POS, POS2) #set two keys makes processing faster.
  X[, POS2 := POS]
  X<-X[,.(CHR,POS,POS2,MAF,ID,FD, Seg1)]
  setkey(X, POS,POS2)
  X_windows <-
    foverlaps(X, windows_dt, type = "within", nomatch = 0L)
  X_windows<-X_windows[,.(CHR,POS,POS2,MAF,Seg1,FD, ID)]
  #this is not ideal because counts end position but it's okay as long as i document this well.
  X_windows[ , Window := .GRP, by = .(POS, POS2)][
    order(Window,POS)]
  X_windows[,N_Sites:=.N,by = Window][,N_SNP:= sum(Seg1==T),by = Window]
  if(fd==T){
    X_windows[, N_FD:=sum(FD==T), by=Window]
    X_windows[, IS:=N_FD+N_SNP, by=Window]
    X_windows[,PtoD:= N_SNP/(N_FD+1), by=Window]
    unique(X_windows, by="Window")-> X2_windows
    
    X2_windows[,MAF:=NULL][,DAF:=NULL][,FD:=NULL][,Seg1:=NULL]
    setkey(X2_windows, Window)
    setkey(X_windows, Window)
    X_NCD<-X_windows[Seg1==T|FD==T][,NCD2_tf0.5 := sqrt(sum((MAF-0.5)^2)/IS), by=Window]
    X_NCD[,NCD2_tf0.45 := sqrt(sum((MAF-0.45)^2)/IS), by=Window]
    X_NCD[,NCD2_tf0.4 := sqrt(sum((MAF-0.4)^2)/IS), by=Window]
    X_NCD[,NCD2_tf0.35 := sqrt(sum((MAF-0.35)^2)/IS), by=Window]
    X_NCD[,NCD2_tf0.3 := sqrt(sum((MAF-0.3)^2)/IS), by=Window]
    X_NCD[,MaxMAf:=round(max(MAF),2), by=Window]
    
  }else{
    X_windows[, N_FD:=NA]
    X_windows[, IS:=N_SNP]
    X_windows[,PtoD:= NA]
    unique(X_windows, by="Window")-> X2_windows
    X2_windows[,MAF:=NULL][,DAF:=NULL][,FD:=NULL][,Seg1:=NULL]
    setkey(X2_windows, Window)
    setkey(X_windows, Window)
    
    X_NCD<-X_windows[Seg1==T][,NCD1_tf0.5 := sqrt(sum((MAF-0.5)^2)/IS), by=Window]
    X_NCD[,NCD1_tf0.45 := sqrt(sum((MAF-0.45)^2)/IS), by=Window]
    X_NCD[,NCD1_tf0.4 := sqrt(sum((MAF-0.4)^2)/IS), by=Window]
    X_NCD[,NCD1_tf0.35 := sqrt(sum((MAF-0.35)^2)/IS), by=Window]
    X_NCD[,NCD1_tf0.3 := sqrt(sum((MAF-0.3)^2)/IS), by=Window]
  }
  setkey(X_NCD, Window)
  unique(X_NCD, by="Window")-> X_NCD
  X_NCD[,W:=W][,S:=S][,MidPos:=POS+S]
  win<-sample(unique(X_windows[,Window]),1) #sample a random window
  Y<-X_windows[Window==win]
  p<-SFS_func(X=Y, type="DAF")
  return(list(NCD=X_NCD, plot=p))
}
