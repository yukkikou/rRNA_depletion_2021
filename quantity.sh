#!/bin/bash
#File: /media/data7/hxy/PM1_lnc/scripts/requantity_total.sh

genome=/media/data5/hxy/genome/pm1.soft.fasta
Hisatindex=/media/data5/hxy/genome/hisatindex/pm1
fqdir=/media/data7/hxy/PM1_lncRNA/subSam/fq
gtf=/media/data7/hxy/PM1_lncRNA/total/rD_tA_hk_smerged.gtf
workdir=/media/data7/hxy/PM1_lncRNA/subSam
thread=65

mkdir -p ${workdir}
source ~/miniconda3/etc/profile.d/conda.sh
source ~/miniconda3/bin/activate
conda activate circ
conda info --envs

for p in {C1,C2,EC3,Y3,EY3,EY4}
do
       echo $p "****** hisat2 start: `date` ******" && \
       time hisat2 \
        --quiet -p ${thread} \
        --dta -x ${Hisatindex} \
        --rna-strandness RF \
        --max-intronlen 2000 \
        --summary-file $workdir/bam/${p}.hisat2.log \
        -1 ${fqdir}/${p}_R1.sub.fq \
        -2 ${fqdir}/${p}_R2.sub.fq |\
        sambamba view -q -t ${thread} --sam-input /dev/stdin --format=bam |\
        sambamba sort -q -t ${thread} --tmpdir=$workdir/tmp -m 100G -o $workdir/bam/${p}.sort.bam /dev/stdin
        echo "******* hisat2 End: `date` *******"

	echo $p "****** stringtie start: `date` ******" && \
	time stringtie $workdir/bam/${p}.sort.bam  \
	-p ${thread} \
	-G ${gtf} -e \
	-b ${workdir}/gtf/${p} \
	-o ${workdir}/gtf/${p}.gtf \
	--rf --conservative
	echo "****** stringtie End: `date` ******"
done

for p in {SRR5028789,SRR5028790,SRR6435878,SRR6435879}
do
        echo $p "****** hisat2 start: `date` ******" && \
        time hisat2 \
         --quiet -p ${thread} \
         --dta -x ${Hisatindex} \
         --rna-strandness RF \
         --max-intronlen 2000 \
         --summary-file $workdir/bam/${p}.hisat2.log \
         -1 ${fqdir}/${p}_1.sub.fq \
         -2 ${fqdir}/${p}_2.sub.fq |\
         sambamba view -q -t ${thread} --sam-input /dev/stdin --format=bam |\
         sambamba sort -q -t ${thread} --tmpdir=$workdir/tmp -m 100G -o $workdir/bam/${p}.sort.bam /dev/stdin
         echo "******* hisat2 End: `date` *******"

        echo $p "****** stringtie start: `date` ******" && \
        time stringtie ${workdir}/bam/${p}.sort.bam  \
        -p ${thread} \
        -G ${gtf} -e \
        -b ${workdir}/gtf/${p} \
        -o ${workdir}/gtf/${p}.gtf \
        --rf --conservative
        echo "****** stringtie End: `date` ******"
done

#conda activate py2
#echo "*** Exp matrix:`date` ***"
#python /media/data7/hxy/PM1_lncRNA/script/stringtie/prepDE.py -i $workdir/rD_prepDE.lst
#python /media/data7/hxy/PM1_lncRNA/script/stringtie/getTPM.py -i $workdir/rD_prepDE.lst
#
#python /media/data7/hxy/PM1_lncRNA/script/stringtie/prepDE.py -i $workdir/hk_prepDE.lst
#python /media/data7/hxy/PM1_lncRNA/script/stringtie/getTPM.py -i $workdir/hk_prepDE.lst
echo "*** END:`date` ***"
