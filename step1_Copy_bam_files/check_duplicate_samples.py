import re
import pandas as pd

files451="project_1630_451samples"
files447="project_1531_453samples"

l451=pd.read_csv(files451,header=None)[0].to_list()
deleted_4samples=['1740T-Dog', '1607H-Dog', '1562B-Dog', '1663T-Dog']
l_1630=[i for i in l451 if i not in deleted_4samples]

l_1531=pd.read_csv(files447,header=None)[0].to_list()
l4=[i for i in l_1630 if i in l_1531]

print(len(l4),len(l_1531),len(l_1630))

import glob

all1630=glob.glob("/nfs/irods-cgp-*/intproj/1630/sample/*/*.bam")
all1630=[i.split("/")[-1] for i in all1630]
print(len(all1630))
print(all1630[:5])
all1531=glob.glob("/nfs/irods-cgp-*/intproj/1531/sample/*/*.bam")
all1531=[i.split("/")[-1] for i in all1531]
print(len(all1531))
print(all1531[:5])
a=[i for i in all1630 if i in all1531]
print(a[:5])
print(len(a))
