#!/bin/bash
#SBATCH --mail-user=winn.martinluther@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name=AMY22
#SBATCH --partition=cpu-long
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=512G
#SBATCH --output /mforge/research/labs/experpath/maia/m315695/logs/%x.stdout
#SBATCH --error /mforge/research/labs/experpath/maia/m315695/logs/%x.stderr

. ~/.bash_profile
module restore

MASDIR=/mforge/research/labs/experpath/maia

PICARD=/usr/local/biotools/picard/2.27.5/picard.jar
BOWTIE2INDEX=${MASDIR}/shared/Indexes/Bowtie2Index/mm10
CHROMSIZES=${MASDIR}/m237371/public_data/mm10.chrom.sizes

DATADIR=${MASDIR}/shared/230418_VH01077_117_AACJ3TTM5
OUTDIR=${MASDIR}/m315695/CUTandRUN/trimmed_batch2

# mkdir -p "$OUTDIR"


SeqId=(22)
name="AMY"
Pair1=*_L1_R1_*.fastq.gz
Pair2=*_L1_R2_*.fastq.gz

for i in ${SeqId[@]}
do
    READS1+=($name${i}.$Pair1)
    READS2+=($name${i}.$Pair2)
    SAMPLE+=($name${i})
done

#temp file with list of commands to run virtually
trim_file=$(mktemp)
bow_file=$(mktemp)
sam_file=$(mktemp)
sami_file=$(mktemp)
samv_file=$(mktemp)
javan_file=$(mktemp)
javad_file=$(mktemp)
samii_file=$(mktemp)
bedb_file=$(mktemp)
rmx_file=$(mktemp)
samiii_file=$(mktemp)
bc_file=$(mktemp)
sams_file=$(mktemp)
samf_file=$(mktemp)
bedbb_file=$(mktemp)
awkx_file=$(mktemp)
bedt_file=$(mktemp)
bedg_file=$(mktemp)


truncate -s 0 "trim_file"
truncate -s 0 "bow_file"
truncate -s 0 "sam_file"
truncate -s 0 "sami_file"
truncate -s 0 "samv_file"
truncate -s 0 "javan_file"
truncate -s 0 "javad_file"
truncate -s 0 "samii_file"
truncate -s 0 "bedb_file"
truncate -s 0 "rmx_file"
truncate -s 0 "samiii_file"
truncate -s 0 "bc_file"
truncate -s 0 "sams_file"
truncate -s 0 "samf_file"
truncate -s 0 "bedbb_file"
truncate -s 0 "awkx_file"
truncate -s 0 "bedt_file"
truncate -s 0 "bedg_file"

wait
for i in ${!SeqId[@]}
do
    echo "trim_galore --cores 2 --paired ${DATADIR}/${READS1[i]} ${DATADIR}/${READS2[i]} --fastqc --output_dir ${OUTDIR} --nextseq 20 --basename ${SAMPLE[i]}" >> "trim_file"
done

for i in ${!SeqId[@]}
do
    echo "bowtie2 --dovetail --threads 8 -x ${BOWTIE2INDEX} -1 ${OUTDIR}/${SAMPLE[i]}_val_1.fq.gz -2 ${OUTDIR}/${SAMPLE[i]}_val_2.fq.gz -S ${OUTDIR}/${SAMPLE[i]}.sam 2> ${OUTDIR}/${SAMPLE[i]}_bowtie2.log" >> "bow_file"
done

for i in ${!SeqId[@]}
do 
    echo "samtools view -@ 8 -bS -q 30 ${OUTDIR}/${SAMPLE[i]}.sam | samtools sort -O BAM - > ${OUTDIR}/${SAMPLE[i]}_temp.bam" >> "sam_file"
done

for i in ${!SeqId[@]}
do 
    echo "samtools index -@ 16 ${OUTDIR}/${SAMPLE[i]}_temp.bam" >> "sami_file"
done

for i in ${!SeqId[@]}
do 
    echo "samtools view -@ 16 -b ${OUTDIR}/${SAMPLE[i]}_temp.bam \$(seq 1 22 | sed 's/^/chr/') chrX chrY > ${OUTDIR}/${SAMPLE[i]}_q30.bam" >> "samv_file"
done

for i in ${!SeqId[@]}
do 
    echo "java -Xmx32g -jar ${PICARD} MarkDuplicates I=${OUTDIR}/${SAMPLE[i]}_q30.bam O=${OUTDIR}/${SAMPLE[i]}_NoDup.bam M=${OUTDIR}/${SAMPLE[i]}_picard.log REMOVE_DUPLICATES=true" >> "javan_file"
done

#qc
for i in ${!SeqId[@]}
do 
    echo "java -Xmx32g -jar ${PICARD} MarkDuplicates I=${OUTDIR}/${SAMPLE[i]}_q30.bam O=${OUTDIR}/${SAMPLE[i]}_DupMarked.bam M=${OUTDIR}/${SAMPLE[i]}_picard_1.log REMOVE_DUPLICATES=false" >> "javad_file"
done

for i in ${!SeqId[@]}
do 
    echo "samtools index -@ 8 ${OUTDIR}/${SAMPLE[i]}_DupMarked.bam" >> "samii_file"

done
for i in ${!SeqId[@]}
do 
    echo "bedtools bamtobed -i ${OUTDIR}/${SAMPLE[i]}_DupMarked.bam | awk 'BEGIN{OFS="\t"}{print \$1,\$2,\$3,\$6}' | grep --color=auto -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+\$1} END{printf \"%d\t%d\t%d\t%d\t%f\t%f\t%f\n\",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${OUTDIR}/${SAMPLE[i]}_DupMarked_LC" >> "bedb_file"
done

#cleanup
for i in ${!SeqId[@]}
do 
    echo "rm ${OUTDIR}/${SAMPLE[i]}.sam ${OUTDIR}/${SAMPLE[i]}_q30.bam ${OUTDIR}/${SAMPLE[i]}_DupMarked.bam ${OUTDIR}/${SAMPLE[i]}_DupMarked.bam.bai ${OUTDIR}/${SAMPLE[i]}_picard_1.log ${OUTDIR}/${SAMPLE[i]}_temp.bam ${OUTDIR}/${SAMPLE[i]}_temp.bam.bai" >> "rmx_file"
done

for i in ${!SeqId[@]}
do 
    echo "samtools index -@ 8 ${OUTDIR}/${SAMPLE[i]}_NoDup.bam" >> "samiii_file"
done

for i in ${!SeqId[@]}
do 
    echo "bamCoverage --bam ${OUTDIR}/${SAMPLE[i]}_NoDup.bam --outFileName ${OUTDIR}/${SAMPLE[i]}_RPKM.bw --binSize 10 --normalizeUsing RPKM -p 8" >> "bc_file"
done

for i in ${!SeqId[@]}
do 
    echo "samtools sort -n -O BAM ${OUTDIR}/${SAMPLE[i]}_NoDup.bam > ${OUTDIR}/${SAMPLE[i]}_nsorted.bam" >> "sams_file"
done

for i in ${!SeqId[@]}
do 
    echo "samtools fixmate ${OUTDIR}/${SAMPLE[i]}_nsorted.bam ${OUTDIR}/${SAMPLE[i]}_fixed.bam -@ 8" >> "samf_file"
done

for i in ${!SeqId[@]}
do 
    echo "bedtools bamtobed -bedpe -i ${OUTDIR}/${SAMPLE[i]}_fixed.bam > ${OUTDIR}/${SAMPLE[i]}.bedpe" >> "bedbb_file"
done

for i in ${!SeqId[@]}
do 
    echo "awk '\$1==\$4 && \$6-\$2 < 1000 {print \$0}' ${OUTDIR}/${SAMPLE[i]}.bedpe | cut -f 1,2,6 | sort -k1,1 -k2,2n -k3,3n > ${OUTDIR}/${SAMPLE[i]}_fragments.bed" >> "awkx_file"
done

for i in ${!SeqId[@]}
do 
    echo "bedtools genomecov -bg -i ${OUTDIR}/${SAMPLE[i]}_fragments.bed -g ${CHROMSIZES} > ${OUTDIR}/${SAMPLE[i]}.bg" >> "bedt_file"
done

for i in ${!SeqId[@]}
do 
    echo "bedGraphToBigWig ${OUTDIR}/${SAMPLE[i]}.bg ${CHROMSIZES} ${OUTDIR}/${SAMPLE[i]}_rpd.bw" >> "bedg_file"
done

cat "trim_file" | parallel -P 1
wait
cat "bow_file" | parallel -P 1 
wait
cat "sam_file" | parallel -P 1
wait
cat "sami_file" | parallel -P 1
wait
cat "samv_file" | parallel -P 1
wait
cat "javan_file" | parallel -P 1
wait
cat "javad_file" | parallel -P 1
wait
cat "samii_file" | parallel -P 1
wait
cat "bedb_file" | parallel -P 1
wait
cat "rmx_file" | xargs -P 1 -I {} bash -c "{}"
wait
cat "samiii_file" | xargs -P 1 -I {} bash -c "{}"
wait
cat "bc_file" | xargs -P 1 -I {} bash -c "{}"
wait
cat "sams_file" | xargs -P 1 -I {} bash -c "{}"
wait
cat "samf_file" | xargs -P 1 -I {} bash -c "{}"
wait
cat "bedbb_file" | xargs -P 1 -I {} bash -c "{}"
wait
cat "awkx_file" | parallel -P 1
wait
cat "bedt_file" | xargs -P 1 -I {} bash -c "{}"
wait
cat "bedg_file" | xargs -P 1 -I {} bash -c "{}"
wait