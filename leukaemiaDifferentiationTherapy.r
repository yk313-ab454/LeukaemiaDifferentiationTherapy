library("stringr")
strGEOID="GSE48558"
chip="GPL6244"

leukaemiaDifferentiationTherapy<-function(){
  
  strGct_DB= paste0("./data/",strGEOID,".gct");
  strCls=paste0("./data/cls_",strGEOID,".cls");
  strPlatform=paste0("./data/",chip,".chip");
  
  for(strDirection in c("up","down")){
    strFolderName = paste0("output/" , strGEOID) ;
    strFolderName_dir = paste0( strFolderName , "/" , strDirection);
    strGMT=paste0("./data/cMapGenes_",strDirection,"_ngenes_20_preprocessed_names_cmapNames.gmt");
    strCMD = paste0("java -cp ./gsea/gsea2-2.2.4.jar -Xmx3130m xtools.gsea.Gsea -res " , strGct_DB ," -cls " , strCls , "#disease"  , "_versus_control" ,  " -gmx " , strGMT , " -collapse true -mode Max_probe -norm meandiv -nperm 0 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis -metric Signal2Noise -sort real -order descending -chip " , strPlatform , " -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 0 -rnd_seed timestamp -save_rnd_lists false -zip_report false -out " , strFolderName_dir , " -gui false -set_min 1")
    system(strCMD)
  }
  mergScores();
}
##
mergScores<-function(){
  
  dfUp <- file.info(list.files("./output/GSE48558/up/", full.names = T));
  upFolder=rownames(dfUp)[which.max(dfUp$mtime)][1]
  
  dfdown <- file.info(list.files("./output/GSE48558/down/", full.names = T));
  downFolder=rownames(dfdown)[which.max(dfdown$mtime)][1]
  
  file_up_h=  dir(upFolder,pattern = "gsea_report_for_control_.*xls",full.names = TRUE);
  file_up_d=dir(upFolder,pattern = "gsea_report_for_disease_.*xls",full.names = TRUE);
  file_down_h= dir(downFolder,pattern = "gsea_report_for_control_.*xls",full.names = TRUE);
  file_down_d=  dir(downFolder,pattern = "gsea_report_for_disease_.*xls",full.names = TRUE);
  
  tbl_up_h=read.table(file_up_h,header=TRUE,sep="\t")
  tbl_up_d=read.table(file_up_d,header=TRUE,sep="\t")
  tbl_up= rbind(tbl_up_h,tbl_up_d);
  tbl_up=tbl_up[order(tbl_up[,1]),];
  
  tbl_down_h=read.table(file_down_h,header=TRUE,sep="\t")
  tbl_down_d=read.table(file_down_d,header=TRUE,sep="\t")
  tbl_down=rbind(tbl_down_h,tbl_down_d);
  tbl_down=tbl_down[order(tbl_down[,1]),];
  
  score_Up=tbl_up[,"ES"];
  score_down=tbl_down[,"ES"];
  score=(as.numeric(score_Up) - as.numeric(score_down))/2;
  sp= str_split( tbl_up[,"NAME"],"_");
  instance_ID = unlist(lapply(sp, function(x)return(x[1])));
  cmap_name = unlist(lapply(sp, function(x)return(x[2])));
  
  tbl_results = cbind(instance_ID,cmap_name,score_Up,score_down,score);
  tbl_results = tbl_results[order(score),];
  rank=1:nrow(tbl_results);
  cID=rep(NA,nrow(tbl_results));
  tbl_results=cbind(rank,tbl_results,cID);
  colnames(tbl_results)= c("rank","instance_id","cmap name","score_up","scoreDown","score","cID");
  write.table(tbl_results,paste0("./Results/",strGEOID,".txt"),sep="\t",row.names = FALSE,quote = FALSE);

}
######
leukaemiaDifferentiationTherapy()