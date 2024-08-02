#!/bin/bash
#SBATCH --mail-user=winn.martinluther@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name=BED2BAM
#SBATCH --partition=cpu-long
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --time=4-00:00:00
#SBATCH --mem=256G
#SBATCH --output /mforge/research/labs/experpath/maia/m315695/logs/%x.stdout
#SBATCH --error /mforge/research/labs/experpath/maia/m315695/logs/%x.stderr

. ~/.bash_profile
module restore

MASDIR=/mforge/research/labs/experpath/maia
OUTDIR=${MASDIR}/m315695/CUTandRUN/AVG_AMY
CHROMSIZES=${MASDIR}/m237371/public_data/mm10.chrom.sizes

binding_Ab=(mCh_H3K4me3 mCh_H3k27ac TBX2_H3K4me3 TBX2_H3k27ac)
name="*relaxed.bed"

for i in "${binding_Ab[@]}"
do
    BA_NAM+=("${i}${name}")
done

for i in "${BA_NAM[@]}"
do
    bed_files+=($(find "${OUTDIR}" -type f -name "${i}" -exec basename {} \;) )
done


#actual commands
for ((i=0; i<${#bed_files[@]}; i++))
do
    BTB+=("bedtools bedtobam -i ${OUTDIR}/${bed_files[i]} -g $CHROMSIZES > ${OUTDIR}/${binding_Ab[i]}.sam")
done

for ((i=0; i<${#bed_files[@]}; i++))
do
    SV+=("samtools view -Sb ${OUTDIR}/${binding_Ab[i]}.sam > ${OUTDIR}/${binding_Ab[i]}.bam -@ 8 ")
done

for ((i=0; i<${#bed_files[@]}; i++))
do
    SS+=("samtools sort ${OUTDIR}/${binding_Ab[i]}.bam -o ${OUTDIR}/${binding_Ab[i]}_sorted.bam -@ 8")
done

# for ((i=0; i<${#bed_files[@]}; i++))
# do
#     SI+=("samtools index ${OUTDIR}/${binding_Ab[i]}_sorted.bam -@ 8 ")
# done

for ((i=0; i<${#bed_files[@]}; i++))
do
    RMF+=("rm ${OUTDIR}/${binding_Ab[i]}.sam ${OUTDIR}/${ibinding_Ab[i]}.bam ${OUTDIR}/*.bai")

done

parallel -j 4 ::: "${BTB[@]}"
wait
parallel -j 4 ::: "${SV[@]}"
wait
parallel -j 4 ::: "${SS[@]}"
wait
parallel -j 4 ::: "${SI[@]}"
wait
parallel -j 4 ::: "${RMF[@]}"

unset BTB SV SS SI RMF

#average all of the peaks and the bed files using the IGG
SC+=("SEACR_1.3.sh ${OUTDIR}/${binding_Ab[i]}.bg ${OUTDIR}/IgG.bg norm relaxed ${OUTDIR}/${binding_Ab[i]}_comp")