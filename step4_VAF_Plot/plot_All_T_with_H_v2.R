library(stringr)


setwd("/Users/ql4/Documents/work2022/work05/work_0524_mito_SNPs")
df_info<-read.csv("dog_ID.txt")
df<-read.csv("VAF_matrix_added_category_v2.1.txt",header=T,sep="\t",check.names = FALSE)

df_HT<-read.csv("Tumor_has_host.txt",header=T)

pdf("VAF_T_with_H_v2.pdf", height = 11.75, width = 8.25)
par(mfrow = c(2, 1), mar = c(3,3,3,2), oma = c(0,0,4,0))

for (i in df_HT$X0){
  
  dog_id<-strtoi(str_extract(i,regex("\\d+")))
  dog_info<-paste(i,df_info[which(df_info$numberID==dog_id),"info"],sep="  ")
  
  
  T_col_i_label<-paste(i,"label",sep="_")
  T_col_i_cat<-paste(i,"category",sep="_")
  T_col_i_germ<-paste(i,"germ",sep="_")
  
  H_col_i<-df_HT[which(df_HT$X0==i),"X1"]
  H_col_i_label<-paste(H_col_i,"label",sep="_")
  
  
  df1<-df[which((df[[H_col_i]]>0)|(df[[i]]>0)),c("POS",H_col_i,i,T_col_i_label,T_col_i_cat,T_col_i_germ,H_col_i_label)]
  df_unsomatic<-df1[which(df1[[T_col_i_germ]]!="somatic"),]
  df_unsomatic_H<-df_unsomatic[which(df_unsomatic[[H_col_i]]>0),]
  df_unsomatic_T<-df_unsomatic[which(df_unsomatic[[i]]>0),]
  
  df_somatic<-df1[which(df1[[T_col_i_germ]]=="somatic"),]
  df_somatic_H<-df_somatic[which(df_somatic[[H_col_i]]>0),]
  df_somatic_T<-df_somatic[which(df_somatic[[i]]>0),]
  
  
  print(H_col_i)
  print(H_col_i_label)
  
  pch_i<-20
  lwd_i<-1
  
  plot(df_unsomatic_H$POS, df_unsomatic_H[[H_col_i]], col = df_unsomatic_H[[H_col_i_label]], pch = pch_i, cex = 1.4, ylim = c(0, 1),xlim=c(0,17000),xlab="",ylab="VAF")
  if(length(df_somatic_H)!=0){
    print("Host Somatic")
    points(df_somatic_H$POS, df_somatic_H[[H_col_i]], col = alpha(df_somatic_H[[H_col_i_label]],0.4), pch = pch_i, cex = 1.4)
    
  }
  abline(h=c(0.1),lty=2,lwd=lwd_i,col="grey")
  mtext(dog_info, outer = TRUE)
  
  
  plot(df_unsomatic_T$POS, df_unsomatic_T[[i]], col = alpha(df_unsomatic_T[[T_col_i_label]],1), pch = pch_i, cex = 1.4, ylim = c(0, 1),xlim=c(0,17000),xlab="Position",ylab="VAF")
  if(length(df_somatic_T)!=0){
    print("Tumour Somatic")
    points(df_somatic_T$POS, df_somatic_T[[i]], col = alpha(df_somatic_T[[T_col_i_label]],0.4), pch = pch_i, cex = 1.4)
    
  }
  abline(h=c(0.1,0.9),lty=2,lwd=lwd_i,col="grey")
}

dev.off()

