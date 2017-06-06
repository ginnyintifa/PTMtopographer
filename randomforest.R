

#tr_c_f="PTMtopographer/bin/can_sites_properties_first_k.tsv"
#tr_s="PTMtopographer/bin/can_sites_states_first_k.tsv"
#test_c_f="PTMtopographer/bin/can_sites_properties_second_k.tsv"
#test_d_f="PTMtopographer/bin/decoy_sites_properties_second_k.tsv"
#out_c="candidate_second_k_rf.tsv"
#out_d="decoy_second_k_rf.tsv"


#check if the user has these packages installed.

pkg=installed.packages()

if (!"randomForest" %in% pkg)
install.packages("randomForest")

if (!"data.table" %in% pkg)
install.packages("data.table")





library(randomForest)
library(data.table)

randomforest=function(tr_c_f,tr_s,test_c_f, test_d_f,out_c,out_d)
{
  
  
  pspx=fread(tr_c_f)
  
  mypspx=as.data.frame(pspx)
  
  
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
  
  mynewexpx=as.matrix(myexpx[,1:27])
  
  #####
  
  decoyexpx=fread(test_d_f)
  decoymyexpx=as.data.frame(decoyexpx)
  
  
  decoymynewexpx=as.matrix(decoymyexpx[,1:27])
  
  
  set.seed(123)
  
  cat("get input frames prepared")
  cat("\n")
  
  
  myrfmodel=randomForest(as.matrix(shrinkx),y) 
  
  
  cat("get model trained")
  cat("\n")
  
  
  mypredict=predict(myrfmodel,mynewexpx,type="prob")
  
  mydecoy=predict(myrfmodel,decoymynewexpx,type="prob")
  
  
  cat("get prediction")
  cat("\n")
  
  
  
  write.table(mypredict[,2],out_c,quote=F,row.names = F)
  
  write.table(mydecoy[,2],out_d,quote=F,row.names = F)
  
  
}

  

#randomforest(tr_c_f,tr_s,test_c_f, test_d_f,out_c,out_d)
  
  
  
  
  
  
  
  
  
