import re
import pandas as pd

files451="project_1630_451samples"
files447="project_1630_447samples"

l451=pd.read_csv(files451,header=None)[0].to_list()

l447=[]
with open(files447,"r") as f:
    for i in f:
        i1=i.strip().split("/")[-1].split(".v1.sample.dupmarked.bam")[0]
        i2=i.strip().split("/")[-2]
        if i1!=i2:
            print(i1,i2,i1==i2)
        l447.append(i1)

l4=set(l451)-set(l447)
print(l4)
