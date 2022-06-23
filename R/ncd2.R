#' @export
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
