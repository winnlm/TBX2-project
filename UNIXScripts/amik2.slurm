#!/bin/bash
#SBATCH --mail-user=winn.martinluther@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name=AVG_PEAK
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

Bat4=m315695/CUTandRUN/Trim_Batch4_AMY
Bat2IgG=m237371/CUTandRUN/trimmed_batch2

OUTDIR=${MASDIR}/m315695/CUTandRUN/AVG_AMY
mkdir -p "$OUTDIR"

name="AMY*.bg"
bg_files=$(find "${MASDIR}/${Bat4}" -type f -name "${name}" -exec basename {} \;)
sort=($(printf "%s\n" "${bg_files[@]}" | sort -tA -k2n))

names=("AMY22*.bg" "AMY18*.bg")
igG_files=()

#issue getting the file
for i in "${names[@]}"; do
    found_file=$(find "${MASDIR}/${Bat2IgG}" -type f -name "$i" -exec basename {} \;)
    if [ -n "$found_file" ]; then
        igG_files+=("$found_file")
    fi
done

for ((i=0; i<$((${#sort[@]}/3)); i++ ))
do 
    j=$(($i+6))
    p=$(($j+6))

    trip1+=(${sort[i]})
    trip2+=(${sort[j]})
    trip3+=(${sort[p]})
done

binding_Ab=(mCh_H3K4me3 mCh_Anti_TBX mCh_H3k27ac TBX2_H3K4me3 TBX2_Anti_TBX TBX2_H3k27ac)
 
for ((i=0; i<$((${#sort[@]}/3)); i++ ))
do 
    # Step 1: Averaging multiple replicates
    wiggle+=("wiggletools write_bg ${OUTDIR}/${binding_Ab[i]}.bg mean ${MASDIR}/${Bat4}/${trip1[i]} ${MASDIR}/${Bat4}/${trip2[i]} ${MASDIR}/${Bat4}/${trip3[i]}")
done
wiggle+=("wiggletools write_bg ${OUTDIR}/IgG.bg mean ${MASDIR}/${Bat2IgG}/${igG_files[0]} ${MASDIR}/${Bat2IgG}/${igG_files[1]}")


CHROMSIZES=/mforge/research/labs/experpath/maia/m237371/public_data/mm10.chrom.sizes
for ((i=0; i<$((${#sort[@]}/3)); i++ ))
do 
# Step 2: Convert bedGraphToBigWig
B2BW+=("bedGraphToBigWig ${OUTDIR}/${binding_Ab[i]}.bg $CHROMSIZES ${OUTDIR}/${binding_Ab[i]}.bw")
done
B2BW+=("bedGraphToBigWig ${OUTDIR}/IgG.bg $CHROMSIZES ${OUTDIR}/IgG.bw")

for ((i=0; i<$((${#sort[@]}/3)); i++ ))
do 
# Step 3: Peak calling (Comparing sample vs. control)
SC+=("SEACR_1.3.sh ${OUTDIR}/${binding_Ab[i]}.bg ${OUTDIR}/IgG.bg norm relaxed ${OUTDIR}/${binding_Ab[i]}_comp")
done

for ((i=0; i<$((${#sort[@]}/3)); i++ ))
do 
# Step 4: Get comparative bigwig
BWC+=("bigwigCompare --bigwig1 ${OUTDIR}/${binding_Ab[i]}.bw --bigwig2 ${OUTDIR}/IgG.bw --operation ratio --binSize 10 --outFileName ${OUTDIR}/${binding_Ab[i]}_comp.bw")
done

# echo "Wiggle commands:"
# printf "%s\n" "${wiggle[@]}"
# echo "B2BW commands:"
# printf "%s\n" "${B2BW[@]}"
# echo "SC commands:"
# printf "%s\n" "${SC[@]}"
# echo "BWC commands:"
# printf "%s\n" "${BWC[@]}"

parallel -j 7 ::: "${wiggle[@]}"
wait
parallel -j 7 ::: "${B2BW[@]}"
wait
parallel -j 6 ::: "${SC[@]}"
wait
parallel -j 6 ::: "${BWC[@]}"
wait

unset Wiggle B2BW SC BWC
