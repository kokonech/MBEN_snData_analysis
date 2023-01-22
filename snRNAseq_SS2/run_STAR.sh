#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=30g
#PBS -l walltime=10:00:00
#PBS -M k.okonechnikov@dkfz-heidelberg.de
#PBS -j oe
#PBS -o /icgc/dkfzlsdf/analysis/dktk/Ependymoma/external_data/logFolder/STAR_shh

TOOL=/home/okonechn/tools/STAR-STAR_2.4.1d/bin/Linux_x86_64/STAR
GENOME=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/annotations/star_index_gencode19/hg19

echo "##########" 
echo "qsub_runSTAR.sh" 
echo `date`
echo `hostname`
echo "TUMOR: ${TUMOR}"
echo "INPUT: ${INPUT}" 


RESDIR=/b06x-isilon/b06x-m/mbCSF/results/smSeq2Res/${TUMOR}/STAR

if [ ! -d "$RESDIR" ]; then
    mkdir -p $RESDIR
fi


echo "Entering $RESDIR"
cd $RESDIR

SID=$INPUT

mkdir -p $SID
cd $SID

cmd="$TOOL --genomeDir $GENOME --readFilesCommand zcat --runThreadN 8 --outFilterMultimapNmax 1  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${SID}. --readFilesIn $INFILE"
echo $cmd
$cmd

echo "Finished!"

