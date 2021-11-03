#!/bin/bash
#File:/media/data7/hxy/PM1_lnc/rRNA.sh 
#For: caculate rRNA residual ratio of samples

fqdir=/media/data7/hxy/PM1_lncRNA/subSam/fq
workdir=/media/data7/hxy/PM1_lncRNA/subSam/rRNA
bowtie2index=/media/data5/hxy/genome/pm1rRNA/pm1rRNAnETS

source ~/miniconda3/etc/profile.d/conda.sh
source ~/miniconda3/bin/activate
conda activate circ
conda info --envs

#build index
#$bowtie2-build /media/data5/hxy/genome/pm1.rRNA.fa /media/data5/hxy/genome/pm1rRNA/pm1rRNA


#rRNA ratio

for file in {C1,C2,EC3,Y3,EY3,EY4}
do
	echo $file "****** Bowtie2 start: `date` ******" && \
	time bowtie2 -p 40 --qc-filter -x $bowtie2index \
	-1 ${fqdir}/${file}_R1.sub.fq -2 ${fqdir}/${file}_R2.sub.fq \
	--no-unal -S ${workdir}/${file}.realrRNA.sam
	echo "******* Bowtie2 End: `date` *******"
done

for file in {SRR5028789,SRR5028790,SRR6435878,SRR6435879}
do
        echo $file "****** Bowtie2 start: `date` ******" && \
        time bowtie2 -p 40 --qc-filter -x $bowtie2index \
        -1 ${fqdir}/${file}_1.sub.fq -2 ${fqdir}/${file}_2.sub.fq \
        --no-unal -S ${workdir}/${file}.realrRNA.sam
        echo "******* Bowtie2 End: `date` *******"
done

