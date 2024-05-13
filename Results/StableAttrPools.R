# Auxiliary analysis R codes for attractor cases of BioNet or Simulation (Fig. 4,5)

# StablePoint: loaded "*_stab.txt" data.
# Period：loaded "*_step.txt" data.
CounterAttrCase_Bionet<-function(StablePoint,Period){# Fig.4
  pools=cbind.data.frame(StablePoint,',',Period);
  ttt=apply(pools,1,paste,sep ="", collapse = "")# Combine strings
  ttt[Period[,1]==-1024]="ComplexPool,-1024";
  hzhz<-as.matrix(table(ttt));# Each string (state) counter.
  pers=unlist(strsplit(rownames(hzhz),','))
  pers=as.numeric(pers[2*c(1:nrow(hzhz))]);
  res=cbind.data.frame(hzhz,pers);
  colnames(res)=c("case","period");
  return(res);
}

# Period：loaded "*_step.txt" data.
CounterAttrCase_Simulation<-function(Period){# Fig.5
  fps=sum(Period[,1]==1);
  sma=sum((Period[,1]>1)&(Period[,1]<1024))
  lar=sum(Period[,1]==-1024);
  res=c(fps,sma,lar)/nrow(Period);
  names(res)=c("FP","small","large");
  return(res);
}