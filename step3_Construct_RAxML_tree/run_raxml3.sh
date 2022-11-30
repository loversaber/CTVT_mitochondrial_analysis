fa_dir=/nfs/users/nfs_q/ql4/lustre_ql4/work/ctvt_mito/run_somatypus/somatypus_result/check_vaf_result

/nfs/users/nfs_q/ql4/anaconda3/bin/raxmlHPC-PTHREADS-AVX2 -s $fa_dir/RaxML_MSA_v2.4.fa -f a -n raxml3 -m GTRGAMMA -c 4 -T 8 -x 12345 -p 12345 -# autoMRE
