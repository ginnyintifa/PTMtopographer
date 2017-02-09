

tr_c_f="/Users/ginny/A_anal_PTM/PTM_COM/bin/decoy_t/second/can_sites_properties.tsv"
tr_s="/Users/ginny/A_anal_PTM/PTM_COM/bin/decoy_t/second/can_sites_states.tsv"
test_c_f="/Users/ginny/PTMtopographer_NOV/bin/can_sites_properties_test.tsv"
test_d_f="/Users/ginny/PTMtopographer_NOV/bin/decoy_sites_properties_test.tsv"
out_c="candidate_left_t.tsv"
out_d="decoy_left_t.tsv"


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
  
  
  
  
  
  
  
  
  
