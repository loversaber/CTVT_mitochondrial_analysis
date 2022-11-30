import pandas as pd
import re

df_samples=pd.read_csv("numberID_samples_distribution.txt",sep="\t",header=0)
df_samples=df_samples.copy()
df_samples["numberID"]=df_samples["numberID"].astype(str)
df_all=pd.read_excel("2022_05_23_samples.xlsx",sheet_name="CTVT sample database COMPLETE")
df_all=df_all.copy()
#df_all["numberID"]=df_all["EPM_number"].str.split(r",|\s+",expand=True)[0]
l=[]
for index,row in df_all.iterrows():
	#print(row["EPM_number"])
	#if re.match(",",str(row["EPM_number"])):
	if "," in str(row["EPM_number"]):
		k=re.split(",",str(row["EPM_number"]))[0]
		print(row["EPM_number"])
	#elif re.match(r" \(",str(row["EPM_number"])):
	elif " (" in str(row["EPM_number"]):
		k=re.split(" \(",str(row["EPM_number"]))[0]
		print(row["EPM_number"])
	else:
		k=str(row["EPM_number"])
	l.append(k)
df_all["numberID"]=l

df_all[["numberID"]].to_csv("All_dog_ID.csv",sep=",",header=False,index=False)
df_all1=df_all[["numberID","Town","Region","Country","Month_collected","Year_collected"]]

print(df_samples.head(5))
print(df_all.head(5))

print(len(df_samples),len(df_all))
df_end=pd.merge(df_samples,df_all1,on=["numberID"],how="left")
print(len(df_end))
df_end.to_csv("meta_matrix.csv",sep=",",header=True,index=False)
