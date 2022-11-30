samtools=/software/CGP/modules/installs/samtools/samtools-1.14/bin/samtools


for i in `ls ../$1/*.bam`
do
	a=${i%%.bam}
	b=${a##*/}
	d=`echo $b|cut -d"." -f1`
	c=$d"_"$2
	echo $c
	new_bam_file=$c".new.bam"
	cmd="$samtools addreplacerg -@ 20 -r "ID:$c" -r "SM:$c" -r "LB:$c" -r "PL:ILLUMINA" -o $new_bam_file $i"
	echo $cmd
	$cmd
	$samtools index -@ 20 -b $new_bam_file
	#ls -alth $i
done
