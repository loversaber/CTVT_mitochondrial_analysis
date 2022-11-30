from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import re
from sys import argv

opt_fa=open("RaxML_MSA_v2.3.fa","w")
seq3="ref_with_coyote2_MSA.fa"
fil="RaxML_MSA_v1.fa"
with open(seq3, "r") as handle, open(fil,"r") as fil900:
    for record in SeqIO.parse(handle, "fasta"):
        print(record.id,len(record.seq))
        seq0=record.seq
        gap_index=[n for n,i in enumerate(seq0) if i=="-"]
        print(len(gap_index))
        opt_fa.write(">"+record.id+"\n")
        opt_fa.write(str(record.seq)+"\n")
        if record.id=="REF":
            kk=gap_index
            ref_seq=str(record.seq)
    print(kk)
    print(len(ref_seq))
    recs_list=[]
    for line in fil900:
        if line.startswith(">"):
            opt_fa.write(line)
            #print(line)
            nSub0=int(line.split("|")[1])
        else:
            seq=line.strip()
            for index_i in kk:
                #seq=seq[:index_i]+"-"+seq[index_i+1:]#not replace just add
                seq=seq[:index_i]+"-"+seq[index_i:]
            #print(len(seq))
            nSub=len([i for n,i in enumerate(seq) if i!=ref_seq[n]])
            #print(nSub)
            if nSub0!=nSub:
                print("No")
            else:
                print("Yes")
            opt_fa.write(seq+"\n")
opt_fa.close()


