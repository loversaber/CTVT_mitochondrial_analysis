from sys import argv
import re
import pandas as pd
from Bio import SeqIO

def get_MT_ref_fa():
    fasta=SeqIO.parse("/nfs/users/nfs_q/ql4/lustre_ql4/data/reference_sequences/CanFam3.1/CanFam3.1.70_MT.fa","fasta")#16727
    for rec in fasta:
        if rec.id=="MT":
            print(len(rec.seq))
            print(rec.seq)
            seq=str(rec.seq)
    return seq

def generate_seqs_RaxML(df_VAF,samples_RAxML):#generate sequences for all 890 samples
    df_info=pd.read_csv("dog_ID3.txt")##use new modified one
    opt_file=open("RaxML_MSA_v1.fa","w")
    seq_id00=">REF"
    opt_file.write(seq_id00+"\n")
    opt_file.write(mito_ref_seq+"\n")
    cols5=list(df_VAF.columns[:5])
    for sample_i in samples_RAxML:
        sample_i_label=sample_i+"_label"
        df_tmp0=df_VAF[cols5+[sample_i_label]]
        df_tmp1=df_tmp0[df_tmp0[sample_i_label]!="grey"].reset_index(drop=True)
        seq0=mito_ref_seq
        dog_id=re.findall(r"\d+",sample_i)[0]
        sample_info=df_info.loc[df_info["numberID"]==int(dog_id),"info"].values[0]
        for index,row in df_tmp1.iterrows():
            pos=row["POS"]
            ref=row["REF"]
            alt=row["ALT"]
            pos_index=pos-1
            if seq0[pos_index]!=ref:
                print(f"#Not Match:{pos} {ref} {alt} {seq0[pos_index]}")
            seq0=seq0[:pos_index]+alt+seq0[pos_index+1:]
        check_num_subs=[i for n,i in enumerate(seq0) if i!=mito_ref_seq[n]]
        n_subs=len(df_tmp1)
        print(f"{sample_i}\t{n_subs}\t{len(check_num_subs)}")
        seq_id=">"+sample_i.split("_label")[0]+"|"+str(n_subs)+"|"+sample_info
        opt_file.write(seq_id+"\n")
        opt_file.write(seq0+"\n")
    opt_file.close()

def check_project(p,d):
    if p in d:
        if "H" in d[p]:
            Hsamples=d[p]["H"]
            nH=len(Hsamples)
            p_H_samples="|".join(Hsamples)
        else:
            Hsamples=[]
            nH=0
            p_H_samples=""
        if "T" in d[p]:
            Tsamples=d[p]["T"]
            nT=len(Tsamples)
            p_T_samples="|".join(Tsamples)
        else:
            nT=0
            p_T_samples=""
    else:
        Hsamples=[]
        nH=0
        p_H_samples=""
        nT=0
        p_T_samples=""
    if nH>=1 and nT>=1:
        syb="yes"
    elif nH==0 and nT>=1:
        syb="only_Tumor"
    elif nH>=1 and nT==0:
        syb="only_Host"
    else:
        syb="no"
    return p_H_samples,p_T_samples,nH,nT,Hsamples,syb

def get_N(dog_id,d_number_samples):#samples and number of samples in p1 p2
    d_project=d_number_samples[dog_id]
    p1_H_samples,p1_T_samples,p1_nH,p1_nT,p1_Hsamples_list,syb1=check_project("p1",d_project)
    p2_H_samples,p2_T_samples,p2_nH,p2_nT,p2_Hsamples_list,syb2=check_project("p2",d_project)
    N_H=p1_nH+p2_nH
    N_T=p1_nT+p2_nT
    H_samples_list=p1_Hsamples_list+p2_Hsamples_list
    opt_line=f"{dog_id}\t{p1_H_samples}\t{p1_T_samples}\t{p1_nH}\t{p1_nT}\t{p2_H_samples}\t{p2_T_samples}\t{p2_nH}\t{p2_nT}\t{N_H}\t{N_T}\n"
    return opt_line,H_samples_list,syb1,syb2

def find_host_tumor_pair(samples_name):
    opt_numberID_samples=open("numberID_samples_distribution.txt","w")
    opt_numberID_samples.write("\t".join(["numberID","Host_p1","Tumers_p1","nH_p1","nT_p1","Host_p2","Tumers_p2","nH_p2","nT_p2","nH","nT"])+"\n")
    l_already_check=[]
    d_number_samples={}
    num_id_project_check=[]
    for i in samples_name:
        #name=i.split("-Dog")[0]
        name=i
        print(name)
        sample_number_id=re.findall(r"\d+",name)[0]
        if sample_number_id not in l_already_check:#add the dog
            l_already_check.append(sample_number_id)
            d_number_samples[sample_number_id]={}
        if "p1" in name:
            project="p1"
        elif "p2" in name:
            project="p2"
        num_id_proj=f"{sample_number_id}_{project}"
        if num_id_proj not in num_id_project_check:#add project
            num_id_project_check.append(num_id_proj)
            d_number_samples[sample_number_id][project]={}
        host_or_tumor=re.findall(r"T|H|A|B",name)[0]
        if host_or_tumor=="A":
            host_or_tumor="H"
            print("unknown Host or Tumor",i)
        elif host_or_tumor=="B":
            host_or_tumor="T"
            print("unknown Host or Tumor",i)
        #if host_or_tumor=="A" or host_or_tumor=="B":
        print(sample_number_id,host_or_tumor,name)
        if host_or_tumor not in d_number_samples[sample_number_id][project]:#add the sample
            d_number_samples[sample_number_id][project][host_or_tumor]=[name]
        else:
            d_number_samples[sample_number_id][project][host_or_tumor].append(name)
    ###check what samples don't have host 
    l_H_samples=[]
    l_H_T_yes=[]
    l_T_only=[]
    l_H_T_yes_H=[]
    l_H_only=[]
    for sample_number in d_number_samples:
        opt_line,H_samples_list,p1_syb,p2_syb=get_N(sample_number,d_number_samples)#syb for checking H T all exist in P1 or P2
        opt_numberID_samples.write(opt_line)
        if p1_syb=="yes":
            n_T=len(d_number_samples[sample_number]["p1"]["T"])
            l_H_T_yes.extend(d_number_samples[sample_number]["p1"]["T"])
            l_H_T_yes_H.extend([d_number_samples[sample_number]["p1"]["H"][0]]*n_T)
        elif p1_syb=="only_Tumor":
            l_T_only.extend(d_number_samples[sample_number]["p1"]["T"])
        elif p1_syb=="only_Host":
            l_H_only.extend(d_number_samples[sample_number]["p1"]["H"])
        if p2_syb=="yes":
            n_T=len(d_number_samples[sample_number]["p2"]["T"])
            l_H_T_yes.extend(d_number_samples[sample_number]["p2"]["T"])
            l_H_T_yes_H.extend([d_number_samples[sample_number]["p2"]["H"][0]]*n_T)
        elif p2_syb=="only_Tumor":
            l_T_only.extend(d_number_samples[sample_number]["p2"]["T"])
        elif p2_syb=="only_Host":
            l_H_only.extend(d_number_samples[sample_number]["p2"]["H"])
        l_H_samples.extend(H_samples_list)
    print(d_number_samples)
    opt_numberID_samples.close()
    return d_number_samples,l_H_samples,l_H_T_yes,l_T_only,l_H_T_yes_H,l_H_only

def find_low_cov_host(l_H_samples,df_VAF):#p1 and p2 seperately
    df_VAF=df_VAF.copy()
    #df_VAF["ID"]=df_VAF["POS"].astype(str)+"_"+df_VAF["REF"]+"_"+df_VAF["ALT"]
    df_VAF.index=df_VAF["ID"]
    df_VAF_H=df_VAF[l_H_samples]
    df_VAF_H=df_VAF_H.copy()
    l_H_samples_NV=[]
    df_VAF_H=df_VAF_H.copy()
    for col in df_VAF_H.columns:
        df_VAF_H=df_VAF_H.copy()
        df_VAF_H[col+"_NV"]=df_VAF_H[col].str.split("|",expand=True)[1].astype(int)
        l_H_samples_NV.append(col+"_NV")
    df_H_NV=df_VAF_H[l_H_samples_NV]
    df_H_NV.to_csv("H_NV.txt",sep="\t",index=False,header=True)
    d={}
    low_cov_H=[]
    for i in l_H_samples:
        avg_cov=df_VAF[i].str.split("|",expand=True)[2].astype(int).mean()
        if avg_cov<20:
            low_cov_H.append(i)
    return low_cov_H,df_H_NV

def transfer(df,clade,germ_som):
    c=df["Mut"].str.split(">",expand=True)
    df["ID"]=df["POS"].astype(str)+"_"+c[0]+"_"+c[1]
    df["clade"]=clade
    df["germ_som"]=germ_som
    return df[["ID","clade","germ_som"]]

def generate_germline_somatic_clade_defining():
    c1_germ=transfer(pd.read_csv("clade1_germline",sep=r"\s+",names=["POS","Mut"]),"clade1","germline")
    c1_som =transfer(pd.read_csv("clade1_somatic",sep=r"\s+",names=["POS","Mut"]),"clade1","somatic")
    df_c1=pd.concat([c1_germ,c1_som]).reset_index(drop=True)
    #
    c2_germ=transfer(pd.read_csv("clade2_germline",sep=r"\s+",names=["POS","Mut"]),"clade2","germline")
    df_c2=c2_germ.reset_index(drop=True)
    #
    c3_germ=transfer(pd.read_csv("clade3_germline",sep=r"\s+",names=["POS","Mut"]),"clade3","germline")
    df_c3=c3_germ.reset_index(drop=True)
    #
    c4_germ=transfer(pd.read_csv("clade4_germline",sep=r"\s+",names=["POS","Mut"]),"clade4","germline")
    c4_som=transfer(pd.read_csv("clade4_somatic",sep=r"\s+",names=["POS","Mut"]),"clade4","somatic")
    df_c4=pd.concat([c4_germ,c4_som]).reset_index(drop=True)
    #
    c5_germ=transfer(pd.read_csv("clade5_germline",sep=r"\s+",names=["POS","Mut"]),"clade5","germline")
    c5_som=transfer(pd.read_csv("clade5_somatic",sep=r"\s+",names=["POS","Mut"]),"clade5","somatic")
    df_c5=pd.concat([c5_germ,c5_som]).reset_index(drop=True)
    #
    c6_germ=transfer(pd.read_csv("clade6_germline",sep=r"\s+",names=["POS","Mut"]),"clade6","germline")
    c6_som=transfer(pd.read_csv("clade6_somatic",sep=r"\s+",names=["POS","Mut"]),"clade6","somatic")
    df_c6=pd.concat([c6_germ,c6_som]).reset_index(drop=True)
    #
    c7_germ=pd.read_csv("clade7_germline",sep=r"\s+",names=["POS","Mut"])
    df_c7=transfer(c7_germ,"clade7","germline").reset_index(drop=True)
    #
    df_ALL=pd.concat([df_c1,df_c2,df_c3,df_c4,df_c5,df_c6,df_c7]).reset_index(drop=True)
    l_c1=df_c1["ID"].to_list()
    l_c2=df_c2["ID"].to_list()
    l_c3=df_c3["ID"].to_list()
    l_c4=df_c4["ID"].to_list()
    l_c5=df_c5["ID"].to_list()
    l_c6=df_c6["ID"].to_list()
    l_c7=df_c7["ID"].to_list()
    unique_c1=set(l_c1)-set(l_c2+l_c3+l_c4+l_c5+l_c6+l_c7)
    unique_c2=set(l_c2)-set(l_c1+l_c3+l_c4+l_c5+l_c6+l_c7)
    unique_c3=set(l_c3)-set(l_c2+l_c1+l_c4+l_c5+l_c6+l_c7)
    unique_c4=set(l_c4)-set(l_c2+l_c3+l_c1+l_c5+l_c6+l_c7)
    unique_c5=set(l_c5)-set(l_c2+l_c3+l_c4+l_c1+l_c6+l_c7)
    unique_c6=set(l_c6)-set(l_c2+l_c3+l_c4+l_c5+l_c1+l_c7)
    unique_c7=set(l_c7)-set(l_c2+l_c3+l_c4+l_c5+l_c6+l_c1)
    d={1:unique_c1,2:unique_c2,3:unique_c3,4:unique_c4,5:unique_c5,6:unique_c6,7:unique_c7}
    print("#Clade unique Sub:")
    print(d)
    d1={}
    for clade in d:
        for sub_id in d[clade]:
            d1[sub_id]=clade
    l_ALL=list(unique_c1)+list(unique_c2)+list(unique_c3)+list(unique_c4)+list(unique_c5)+list(unique_c6)+list(unique_c7)
    df_END=df_ALL[df_ALL["ID"].isin(l_ALL)].reset_index(drop=True)
    df_END.to_csv("Unique_CladeX_Substitution.txt",sep="\t",header=True,index=False)
    return d1,df_END

def check_cladeX(Sub_ID,d1,df_END):
    d_clade_color={1:"red",2:"blue",3:"green",4:"yellow",5:"purple",6:"brown",7:"cyan"}
    df_sub_clade=df_END[df_END["ID"]==Sub_ID]
    if not df_sub_clade.empty:
        clade= d1[Sub_ID]
        clade_color=d_clade_color[clade]
        germ_som=df_sub_clade["germ_som"].values[0]
    else:
        clade= 0 #not clade define
        clade_color="black"#not clade define substitution
        germ_som="NULL"
    return clade,clade_color,germ_som

def substitution_category(d_H_T,df_VAF,df_VAF1,l_H_samples,l_H_T_yes,l_T_only,l_H_T_yes_H,l_H_only):
    df_VAF=df_VAF.copy()
    df_VAF=df_VAF.reset_index(drop=True)
    cols5= df_VAF.columns[5:]

    df_VAF["ID"]=df_VAF["POS"].astype(str)+"_"+df_VAF["REF"]+"_"+df_VAF["ALT"]
    #l_low_cov_H,df_H_NV=find_low_cov_host(l_H_samples,df_VAF)
    d_cladeX_unique,df_cladeX_unique=generate_germline_somatic_clade_defining()
    d_T_label={}
    d_T_category={}
    d_T_germ={}
    samples_RAxML=l_H_T_yes+l_T_only+list(set(l_H_T_yes_H))+l_H_only
    for sample_T in samples_RAxML:
        T_label=sample_T+"_label"
        cat_label=sample_T+"_category"
        germ_soma_label=sample_T+"_germ"
        d_T_label[T_label]=[]
        d_T_category[cat_label]=[]
        d_T_germ[germ_soma_label]=[]
    for index,row in df_VAF.iterrows():
        print("#INDEX:",index)
        Sub_ID=row["ID"]#f"{row['POS']}_{row['REF']}_{row['ALT']}"
        mark_clade0,clade_color0,germ_soma0=check_cladeX(Sub_ID,d_cladeX_unique,df_cladeX_unique)
        for sample_H in list(set(l_H_T_yes_H))+l_H_only:####should be careful
            H_vaf_list1=row[sample_H].split("|")
            H_nAlt1=int(H_vaf_list1[1])
            H_vaf1=float(H_vaf_list1[0])
            if H_nAlt1>=3 and H_vaf1>=0.1:
                H_clade_i=mark_clade0
                H_clade_color=clade_color0
                H_germ=germ_soma0
            else:
                H_clade_i=-1
                H_clade_color="grey"
                H_germ="NULL"
            H_label_key=sample_H+"_label"
            cat_label_key=sample_H+"_category"
            germ_label=sample_H+"_germ"
            d_T_label[H_label_key].append(H_clade_color)#label)
            d_T_category[cat_label_key].append(H_clade_i)
            d_T_germ[germ_label].append(H_germ)
            if sample_H=="1665H-2-Dog_p2" and row["POS"]==2683:
                print("MATCH",row["POS"],H_vaf_list1,mark_clade0,clade_color0)
                print(sample_H,row[sample_H],H_clade_color,H_clade_i)

        for sample_T in l_H_T_yes:
            project_i=re.findall(r"p1|p2",sample_T)[0]
            number_id=re.findall(r"\d+",sample_T)[0]
            sample_H=d_H_T[number_id][project_i]["H"][0]
            H_vaf_list=row[sample_H].split("|")
            T_vaf_list=row[sample_T].split("|")
            H_vaf=float(H_vaf_list[0])
            T_vaf=float(T_vaf_list[0])
            H_nAlt=int(H_vaf_list[1])
            T_nAlt=int(T_vaf_list[1])
            H_nTotal=int(H_vaf_list[1])
            H_clade_color=d_T_label[sample_H+"_label"]
            if sample_T=="1016T-Dog_p1" and Sub_ID=="15931_A_G":
                print("Check Sub:")
                print(row[sample_T])
                print(T_nAlt,T_vaf)
                print(len(H_clade_color))
                print(sample_H,H_clade_color,T_nAlt,T_vaf)
                print(H_clade_color[index])
                print(row[sample_H])
            if H_clade_color[index]=="grey":#absence in Host
                if T_nAlt>=3 and T_vaf>=0.1: #presence in Tumor
                    label=mark_clade0#"pass"
                    clade_color=clade_color0
                    germ_mark=germ_soma0
                else:
                    label=-1#"deleted"#
                    clade_color="grey"
                    germ_mark="NULL"
            else:# presence in Host
                if T_nAlt>=3 and T_vaf>=0.9:
                    label=mark_clade0
                    clade_color=clade_color0
                    germ_mark=germ_soma0
                else:
                    label=-1#"deleted"#
                    clade_color="grey"
                    germ_mark="NULL"
            T_label=sample_T+"_label"
            cat_label=sample_T+"_category"
            germ_label1=sample_T+"_germ"
            d_T_label[T_label].append(clade_color)#label)
            d_T_category[cat_label].append(label)
            d_T_germ[germ_label1].append(germ_mark)

        for sample_T1 in l_T_only:
            #print(f"#pro only tumor: {sample_T1}")
            T_vaf_list1=row[sample_T1].split("|")
            T_vaf1=float(T_vaf_list1[0])
            T_nAlt1=int(T_vaf_list1[1])
            if T_nAlt1>=3 and T_vaf1>=0.5:
                label1=mark_clade0#"pass"
                clade_color1=clade_color0
                germ_mark1=germ_soma0
            else:
                label1=-1#"deleted"
                clade_color1="grey"
                germ_mark1="NULL"
            T_label1=sample_T1+"_label"
            cat_label1=sample_T1+"_category"
            germ_label2=sample_T1+"_germ"
            d_T_label[T_label1].append(clade_color1)
            d_T_category[cat_label1].append(label1)
            d_T_germ[germ_label2].append(germ_mark1)
    print("!!!!Add the color column and mark column to dataframe")
    #df_VAF=df_VAF1.copy()
    print(len(df_VAF))
    for T_label in d_T_label:
        df_VAF=df_VAF.copy()
        print(T_label,len(d_T_label[T_label]))
        df_VAF[T_label]=d_T_label[T_label]
    for cat_label in d_T_category:
        df_VAF=df_VAF.copy()
        df_VAF[cat_label]=d_T_category[cat_label]
    for germ_lable in d_T_germ:
        df_VAF=df_VAF.copy()
        df_VAF[germ_lable]=d_T_germ[germ_lable]
    
    #df_VAF.to_csv("VAF_matrix_added_category_RAxML.txt",sep="\t",header=True,index=False)
    df_VAF2=df_VAF1.copy()
    for T_label in d_T_label:
        df_VAF2=df_VAF2.copy()
        df_VAF2[T_label]=d_T_label[T_label]
    for cat_label in d_T_category:
        df_VAF2=df_VAF2.copy()
        df_VAF2[cat_label]=d_T_category[cat_label]
    for germ_lable in d_T_germ:
        df_VAF2=df_VAF2.copy()
        df_VAF2[germ_lable]=d_T_germ[germ_lable]
    df_VAF2.to_csv("VAF_matrix_added_category_v2.1.txt",sep="\t",header=True,index=False)

    df_VAF=df_VAF.copy()
    df_VAF["ID"]=df_VAF["POS"].astype(str)+"_"+df_VAF["REF"]+"_"+df_VAF["ALT"]

    for sample_i in l_T_only+l_H_only+list(set(l_H_T_yes_H)):
        color_i=sample_i+"_label"
        category_i=sample_i+"_category"
        df_tmp=df_VAF[["ID",sample_i,color_i,category_i]]
        df_tmp=df_tmp.copy()
        df_tmp[sample_i+"_vaf"]=df_tmp[sample_i].str.split("|",expand=True)[0].astype(float)
        df_tmp1=df_tmp[df_tmp[sample_i+"_vaf"]>0]
        df_tmp1=df_tmp1[["ID",sample_i,color_i,category_i]]
        df_tmp1.to_csv(sample_i+".tmp_vaf",sep="\t",header=True,index=False)
    for sample_T in l_H_T_yes:
        project_i=re.findall(r"p1|p2",sample_T)[0]
        number_id=re.findall(r"\d+",sample_T)[0]
        sample_H=d_H_T[number_id][project_i]["H"][0]
        color_i=sample_T+"_label"
        category_i=sample_T+"_category"
        color_i1=sample_H+"_label"
        category_i1=sample_H+"_category"
        df_tmp=df_VAF[["ID",sample_T,color_i,category_i,sample_H,color_i1,category_i1]]
        df_tmp=df_tmp.copy()
        df_tmp[sample_T+"_vaf"]=df_tmp[sample_T].str.split("|",expand=True)[0].astype(float)
        #df_tmp[sample_H+"_vaf"]=df_tmp[sample_H].str.split("|",expand=True)[0].astype(float)
        df_tmp1=df_tmp[df_tmp[sample_T+"_vaf"]>0]
        df_tmp1=df_tmp1[["ID",sample_T,color_i,category_i,sample_H,color_i1,category_i1]]
        df_tmp1.to_csv(sample_T+".tmp_vaf",sep="\t",header=True,index=False)
    return df_VAF2,samples_RAxML

def pro_vcf(vcf_file):
    opt_vaf_mtx=open("VAF_matrix.txt","w")
    opt_vaf_mtx1=open("VAF_matrix1.txt","w")
    #opt_cov_mtx=open("Cov_matrix.txt","w")
    with open(vcf_file,"r") as f:
        d_sample_index={}
        for line in f:
            if line.startswith("#CHROM"):
                line_list=line.strip().split("\t")
                cols=line_list[:5]
                samples_name=line_list[9:]
                d_number_samples_H_T,l_H_samples,l_H_T_yes,l_T_only,l_H_T_yes_H,l_H_only=find_host_tumor_pair(samples_name)
                print("#l_H_T_yes:")
                print(l_H_T_yes)
                df_H_T_all_have=pd.DataFrame(l_H_T_yes)
                print(len(l_H_T_yes),len(l_H_T_yes_H))
                df_H_T_all_have[1]=l_H_T_yes_H
                df_H_T_all_have[[0,1]].to_csv("Tumor_has_host.txt",header=True,index=False)
                print("#Only Tumors:")
                print(l_T_only)
                df_T_only=pd.DataFrame(l_T_only)
                df_T_only[0].to_csv("Tumor_no_host.txt",header=True,index=False)
                opt_vaf_mtx.write("\t".join(cols+samples_name)+"\n")
                opt_vaf_mtx1.write("\t".join(cols+samples_name)+"\n")
                for n,sample_i in enumerate(samples_name):
                    d_sample_index[sample_i]=n
            elif not line.startswith("#"):
                line_list=line.strip().split("\t")
                POS=int(line_list[1])
                if POS<16129 or POS>16430:
                    cols1=line_list[:5]
                    format_info=line_list[8]
                    samples_info=line_list[9:]
                    l_vaf=[]
                    l_vaf1=[]
                    for info_i in samples_info:
                        kk=dict(zip(format_info.split(":"),info_i.split(":")))
                        nTotal=int(kk["NR"])
                        nMut=int(kk["NV"])
                        if nTotal!=0:
                            vaf=nMut/nTotal
                        else:
                            vaf=0
                        item_i=f"{vaf}|{nMut}|{nTotal}"
                        item_i1=f"{vaf}"
                        l_vaf.append(item_i)
                        l_vaf1.append(item_i1)
                    opt_vaf_mtx.write("\t".join(cols1+l_vaf)+"\n")
                    opt_vaf_mtx1.write("\t".join(cols1+l_vaf1)+"\n")
                else:
                    print(f"#SNP in low complexity region: {POS}")
    opt_vaf_mtx.close()
    opt_vaf_mtx1.close()
    df_VAF=pd.read_csv("VAF_matrix.txt",sep="\t",header=0)
    df_VAF1=pd.read_csv("VAF_matrix1.txt",sep="\t",header=0)
    df_VAF,samples_RAxML=substitution_category(d_number_samples_H_T,df_VAF,df_VAF1,l_H_samples,l_H_T_yes,l_T_only,l_H_T_yes_H,l_H_only)
    ###generate the MSA fasta file 
    generate_seqs_RaxML(df_VAF,samples_RAxML)
    print("#l_H_T_yes_H:",len(l_H_T_yes_H))
    print("#list(set(l_H_T_yes_H)):",len(list(set(l_H_T_yes_H))))
    print(f"samples_RAxML:{len(samples_RAxML)}")

if __name__=="__main__":
    mito_ref_seq=get_MT_ref_fa()
    vcf_file=argv[1]
    pro_vcf(vcf_file)



