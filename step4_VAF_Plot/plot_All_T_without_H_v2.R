library(stringr)


setwd("/Users/ql4/Documents/work2022/work05/work_0524_mito_SNPs")
df_info<-read.csv("dog_ID.txt")
df<-read.csv("VAF_matrix_added_category_v2.1.txt",header=T,sep="\t",check.names = FALSE)

df_HT1<-read.csv("Tumor_no_host.txt",header=T)

pdf("VAF_T_without_H_v2.pdf", height = 8, width = 8)
par(mfrow = c(1, 1), mar = c(3,3,3,2), oma = c(0,0,4,0))

for (i in df_HT1$X0){
  
  dog_id<-strtoi(str_extract(i,regex("\\d+")))
  dog_info<-paste(i,df_info[which(df_info$numberID==dog_id),"info"],sep="  ")
  
  
  T_col_i_label<-paste(i,"label",sep="_")
  T_col_i_cat<-paste(i,"category",sep="_")
  T_col_i_germ<-paste(i,"germ",sep="_")
  
  
  df1<-df[which(df[[i]]>0),c("POS",i,T_col_i_label,T_col_i_cat,T_col_i_germ)]
  df_unsomatic_T<-df1[which(df1[[T_col_i_germ]]!="somatic"),]
  #df_unsomatic_T<-df_unsomatic[which(df_unsomatic[[i]]>0),]
  
  df_somatic_T<-df1[which(df1[[T_col_i_germ]]=="somatic"),]
  #df_somatic_T<-df_somatic[which(df_somatic[[i]]>0),]
  
  pch_i<-20
  lwd_i<-1
  
  plot(df_unsomatic_T$POS, df_unsomatic_T[[i]], col = alpha(df_unsomatic_T[[T_col_i_label]],1), pch = pch_i, cex = 1.4, ylim = c(0, 1),xlim=c(0,17000),xlab="Position",ylab="VAF")
  if(length(df_somatic_T)!=0){
    print("Tumour Somatic")
    points(df_somatic_T$POS, df_somatic_T[[i]], col = alpha(df_somatic_T[[T_col_i_label]],0.4), pch = pch_i, cex = 1.4)
    
  }
  abline(h=c(0.5),lty=2,lwd=lwd_i,col="grey")
  mtext(dog_info, outer = TRUE)

}

dev.off()

