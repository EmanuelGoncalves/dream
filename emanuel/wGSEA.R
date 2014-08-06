wGSEA<-function(ranking, norm_express, signature, p=1, display=FALSE, returnRS=FALSE, significance=FALSE, trial=1000){
  signature<-unique(intersect(signature,ranking))
  
  HITS<-is.element(ranking,signature)+0
  R<-norm_express*HITS
  
  hitCases<-cumsum(abs(R)^p)
  NR<-sum(abs(R)^p)
  
  
  missCases<-cumsum(1-HITS)
  
  
  N<-length(ranking)
  N_Nh<-length(ranking)-length(signature)
  
  Phit<-hitCases/NR
  Pmiss<-missCases/N_Nh
  
  m<-max(abs(Phit-Pmiss))
  t<-which(abs(Phit-Pmiss)==m)
  
  if (length(t)>1){t<-t[1]}
  
  peak<-t
  ES<-Phit[t]-Pmiss[t]
  RS<-Phit-Pmiss
  
  names(ES)<-NULL
  
  if (display){
    if (ES>=0){c<-"red"}else{c<-"green"}
    plot(0:N,c(0,Phit-Pmiss),col=c,type="l",xlim=c(0,N),ylim=c(-(abs(ES)+0.5*(abs(ES))),abs(ES)+0.5*(abs(ES))),xaxs="i",bty="l",axes=FALSE,
         xlab="Gene Rank Position",ylab="Running Sum")
    par(new=TRUE)
    plot(0:N,rep(0,N+1),col='gray',type="l",new=FALSE,xlab="",ylab="",ylim=c(-(abs(ES)+0.5*(abs(ES))),abs(ES)+0.5*(abs(ES))))
    axis(side=2)   
  }  
  
  
  P<-NA
  if(significance) {
    EMPES<-computeSimpleEMPES(ranking,norm_express,signature,trial)
    P<-(ES>=0)*(length(which(EMPES>=ES))/trial)+(ES<0)*(length(which(EMPES<=ES))/trial)
  }
  
  if (returnRS) {
    POSITIONS<-which(HITS==1)
    names(POSITIONS)<-ranking[which(HITS==1)]
    
    POSITIONS<-POSITIONS[order(names(POSITIONS))]
    names(POSITIONS)<-names(POSITIONS)[order(names(POSITIONS))]
    
    return(list(ES=ES,RS=RS,POSITIONS=POSITIONS,PEAK=t))
  } else {
    if (significance) {
      return(list(ES=ES,P=P))
    } else {
      return(ES)
    }
  }
}


computeSimpleEMPES<-function(ranking,exp_value_profile,signature,trials){
	
	ngenes<-length(ranking)
	siglen<-length(intersect(signature,ranking))
	
	ES<-rep(NA,trials)
	
	for (i in 1:trials){
		
		shuffled_signature<-ranking[sample(1:ngenes,siglen)]
		tmp<-wGSEA(ranking,exp_value_profile,shuffled_signature,display=FALSE,significance=FALSE)
		ES[i]<-tmp
    }
	
	return(ES)
	
}