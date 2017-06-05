

#tr_c_f="PTMtopographer/bin/can_sites_properties_first_k.tsv"
#tr_s="PTMtopographer/bin/can_sites_states_first_k.tsv"
#test_c_f="PTMtopographer/bin/can_sites_properties_second_k.tsv"
#test_d_f="PTMtopographer/bin/decoy_sites_properties_second_k.tsv"
#out_c="candidate_second_k_rf.tsv"
#out_d="decoy_second_k_rf.tsv"


library(randomForest)
library(data.table)

randomforest=function(tr_c_f,tr_s,test_c_f, test_d_f,out_c,out_d)
{
  
  
  pspx=fread(tr_c_f)
  
  mypspx=as.data.frame(pspx)
  
  
    
  mypspx=mypspx[,-28]
  
  shrinkx=mypspx[,1:27]
  
  pspy=fread(tr_s)
  mypspy=as.data.frame(pspy)
  
  
  
  unmid=which(mypspy[[1]]==0)
  mypspy[[1]][unmid]="unmodified"
  mypspy[[1]][-unmid]="modified"
  
  
  shrinky=as.data.frame(mypspy)
  
  shrinky=factor(shrinky[[1]])
  
  y=relevel(shrinky, ref = "unmodified")
  
  
  ##################
  #data to predict
  
  expx=fread(test_c_f)
  myexpx=as.data.frame(expx)
  
  myexpx=myexpx[,-28]
  
  mynewexpx=as.matrix(myexpx[,1:27])
  
  #####
  
  decoyexpx=fread(test_d_f)
  decoymyexpx=as.data.frame(decoyexpx)
  
  decoymyexpx=decoymyexpx[,-28]
  
  decoymynewexpx=as.matrix(decoymyexpx[,1:27])
  
  
  set.seed(123)
  
  
  
  myrfmodel=randomForest(as.matrix(shrinkx),y) 
  
  
  mypredict=predict(myrfmodel,mynewexpx,type="prob")
  
  mydecoy=predict(myrfmodel,decoymynewexpx,type="prob")
  
  write.table(mypredict[,2],out_c,quote=F,row.names = F)
  
  write.table(mydecoy[,2],out_d,quote=F,row.names = F)
  
  
}

  

randomforest(tr_c_f,tr_s,test_c_f, test_d_f,out_c,out_d)
  
  
  
  
  
  
  
  
  
