#!/bin/bash
fqdir1=/media/data7/hxy/PM1_lncRNA/HK/fq
fqdir2=/media/data7/hxy/PM1_illumima_202101/ANNO_XS01KF2020120029_PM-XS01KF2020120029-02_AHF3JTCCX2_2021-01-08/Rawdata/
outdir=/media/data7/hxy/PM1_lncRNA/subSam
thread=40

# conda
source ~/miniconda3/etc/profile.d/conda.sh
source ~/miniconda3/bin/activate
conda activate circ
conda info --envs

#for file in `ls ${fqdir1}/*.gz`
#do
#	echo "*** $(basename $file ".fastq.gz"):`date` ***"
#	time seqkit sample -s3104 ${file} -n 1000000 -j $thread > ${outdir}/$(basename $file ".fastq.gz").sub.fq
#done

for file in {C1,C2,EC3,Y3,EY3,EY4}
do
	echo "*** $(basename ${fqdir2}/${file}/${file}_R1.fq.gz ".fq.gz"):`date` ***"
	time seqkit sample -s3104 ${fqdir2}/${file}/${file}_R1.fq.gz -n 1000000 -j $thread > ${outdir}/$(basename ${fqdir2}/${file}/${file}_R1.fq.gz ".fq.gz").sub.fq
	echo "*** $(basename ${fqdir2}/${file}/${file}_R2.fq.gz ".fq.gz"):`date` ***"
	time seqkit sample -s3104 ${fqdir2}/${file}/${file}_R2.fq.gz -n 1000000 -j $thread > ${outdir}/$(basename ${fqdir2}/${file}/${file}_R2.fq.gz ".fq.gz").sub.fq
done

