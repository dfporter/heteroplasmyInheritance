# Extracts MAPQ >= 20 mtDNA mapping sequences.
# Example: bash extract_unique_mt.sh infolder outfolder

for FULLBAM in $1/*.bam; do
	if [[ ${FULLBAM} == *.mt.bam ]]
	then
		echo "Skipping ${FULLBAM}"
	else
		OUTNAME=$(basename ${FULLBAM} .bam).mt.bam
		#echo "samtools view -h -b ${FULLBAM} MT > $2/${OUTNAME}"
		echo "samtools view -h -b ${FULLBAM} -q 20 chrM > $2/${OUTNAME}"
		echo "*"
		samtools view -h -b ${FULLBAM} chrM > $2/${OUTNAME}
		echo "samtools sort $2/${OUTNAME} $2/$(basename ${OUTNAME} .bam)"
		echo "*"
		samtools sort $2/${OUTNAME} $2/$(basename ${OUTNAME} .bam)
		echo "samtools index $2/${OUTNAME}.bam"
		echo "*"
		samtools index $2/${OUTNAME}
	fi
done

# Old code below:
#echo "****"
#echo "Extracting unique sequences"
#
#for MTBAM in $2/*.mt.bam; do
#	OUTNAME=$(basename ${MTBAM} .bam)
#	#echo "samtools view -h ${MTBAM} | grep '@\|XT:A:U' > $2/${OUTNAME}.sam"
#	#samtools view -h ${MTBAM} | grep '@\|XT:A:U' > $2/${OUTNAME}.sam
#	echo "samtools view -h ${MTBAM} -q 20 > $2/${OUTNAME}.sam"
#	samtools view -h ${MTBAM} -q 20 > $2/${OUTNAME}.sam
#	echo "samtools view -h -Sb $2/${OUTNAME}.sam > $2/${OUTNAME}.bam"
#	samtools view -h -Sb $2/${OUTNAME}.sam > $2/${OUTNAME}.bam
#	echo "samtools sort $2/${OUTNAME}.bam"
#	samtools sort $2/${OUTNAME}.bam $2/${OUTNAME}
#	echo "samtools index $2/${OUTNAME}.bam"
#	samtools index $2/${OUTNAME}.bam
#done
