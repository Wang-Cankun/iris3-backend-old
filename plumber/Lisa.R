



# for lisa co:co expression gene module 
run_Lisa<-function(co,species="human"){
  system("rm -rf /fs/project/PAS1475/Xiaoying/Lisa1/ct*")
  CT<-unique(unlist(strsplit(names(co),split="_"))[seq(1,length(unlist(strsplit(names(co),split="_"))),2)])
  for (i in CT){
    system(paste0("mkdir /fs/project/PAS1475/Xiaoying/Lisa1/",i))
  }
  for (i in (1:length(co))){
    ct <- unlist(strsplit(names(co[i]),split="_"))[1]
    write.table(co[[i]],paste0("/fs/project/PAS1475/Xiaoying/Lisa1/",ct,"/",names(co[i])),quote=F,sep="\t",row.names = F,
                col.names = F)
  }
  if (species=="human"){
    system(paste0("sh /fs/project/PAS1475/Xiaoying/Lisa1/d.sh ",length(CT)))
  }else{
    system(paste0("sh /fs/project/PAS1475/Xiaoying/Lisa1/mouse.sh ",length(CT)))
  }
  
  return (CT)
}
CT<-run_Lisa(co,"human")

tf_pval_0.05 <- list()
for (ct in (1:length(CT))){
  all_res <- list.files(paste0("/fs/project/PAS1475/Xiaoying/Lisa1/ct",ct,"_results"),full.names = T)
  res_list <- list()
  top_list <- list()
  for (i in 1:length(all_res)) {
    this_name <- basename(all_res[i])
    this_name <- sub(".lisa.tsv","",this_name)
    res_list[[i]] <- read.table(paste0(all_res[i]), sep = "\t", header = T)
    tf_pval_0.05[[paste(ct,i,sep="_")]]<- unique(res_list[[i]][,3][res_list[[i]][,12]<0.05])
    #print(tf_pval_0.05)
    #print(length(tf_pval_0.05))
  }
}

str(tf_pval_0.05)
