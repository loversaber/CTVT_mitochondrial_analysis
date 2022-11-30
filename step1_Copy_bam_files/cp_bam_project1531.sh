for i in `ls /nfs/irods-cgp*/intproj/1531/sample/*/*.bam`
do
	name=`echo $i | cut -d "/" -f 7`
	echo $name
	ls -alth $i
	cp $i ./
done
