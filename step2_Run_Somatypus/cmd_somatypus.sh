export MODULEPATH=/software/CGP/modules/modulefiles:$MODULEPATH
module load /software/CGP/modules/modulefiles/vcftools/0.1.16
#VCFTOOLS_PATH=/software/CGP/modules/installs/vcftools/vcftools-0.1.16/bin
VCFTOOLS=$(which vcftools)
VCFTOOLS_PATH=$(dirname $VCFTOOLS)

# Define a specialised PATH for working with legacy platypus
PLATYPUS_PATH=/nfs/dog_n_devil/adrian/software/somatypus/src:/nfs/dog_n_devil/adrian/software/somatypus/utils:/nfs/dog_n_devil/adrian/software/platypus/bin:/usr/local/lsf/9.1/linux2.6-glibc2.3-x86_64/bin:/opt/renci/icommands//bin:/usr/local/lsf/9.1/linux2.6-glibc2.3-x86_64/etc:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/software/badger/bin:/software/oracle-ic-11.2:/software/bin:/software/CGP/bin:$VCFTOOLS_PATH
PLATYPUS_LDFLAGS=-L/nfs/dog_n_devil/adrian/software/platypus/htslib-1.2.1
PLATYPUS_LD_LIBRARY_PATH=/nfs/dog_n_devil/adrian/software/platypus/htslib-1.2.1

# Run somatypus
somatypus=/nfs/users/nfs_k/kg8/lustre_ms/projects/ctvt_horizontal_transfer/scripts/somatypus/src/somatypus

PATH=$PLATYPUS_PATH \
LDFLAGS=$PLATYPUS_LD_FLAGS \
LD_LIBRARY_PATH=$PLATYPUS_LD_LIBRARY_PATH \

python3=/nfs/users/nfs_q/ql4/anaconda3/bin/python3
script_path=/nfs/users/nfs_q/ql4/script/cmd_index_error_files
new_samtools_path=/software/CGP/modules/installs/samtools/samtools-1.14/bin
old_samtools_path=/nfs/users/nfs_q/ql4/lustre_ql4/tools/samtools-1.2/samtools

bam_files_path=/nfs/users/nfs_q/ql4/lustre_ql4/work/ctvt_mito/bam_files/bam_files_900

somatypus_opt_path=/nfs/users/nfs_q/ql4/lustre_ql4/work/ctvt_mito/run_somatypus/somatypus_result

ref_path=/nfs/users/nfs_q/ql4/lustre_ql4/data/reference_sequences/CanFam3.1/Canis_familiaris.CanFam3.1.70.dna.toplevel.fa
mito_region=/nfs/users/nfs_q/ql4/lustre_ql4/work/ctvt_mito/test_run_somatypus/mitochondrial_region.txt

$somatypus -i $bam_files_path -g $ref_path -o $somatypus_opt_path -r $mito_region -c 8 -p '--rmsmqThreshold=20 --qdThreshold=5 --maxReads=50000000 --bufferSize=2500'


