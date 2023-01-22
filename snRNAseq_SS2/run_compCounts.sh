#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=32g
#PBS -l walltime=5:00:00
#PBS -M k.okonechnikov@dkfz-heidelberg.de
#PBS -j oe
#PBS -o /icgc/dkfzlsdf/analysis/dktk/Ependymoma/external_data/logFolder/compCountsMbSHH/smSeq2

TOOL=/home/okonechn/tools/subread-1.6.4-Linux-x86_64/bin/featureCounts
GENCODE=/icgc/dkfzlsdf/analysis/dktk/Ependymoma/annotations/star_index_gencode19/gencode.v19.annotation.gtf

MAINDIR=/b06x-isilon/b06x-m/mbCSF/results/smSeq2Res


echo "##########"
echo "Run featureCounts"
echo `date`
echo `hostname`
echo "INPUT: ${INPUT}"

SID=$INPUT
PID=$TUMOR

DATADIR=$MAINDIR/$PID/STAR
BAMFILE=$DATADIR/$SID/${SID}.Aligned.sortedByCoord.out.bam
RESFILE=$MAINDIR/$PID/counts/${SID}.counts

echo "Processing $SID"
 
cd $RESDIR   
# option -t gene allows to use full gene annotation with introns
cmd="$TOOL -p -T 8 -t gene --tmpDir $DATADIR -a $GENCODE -o $RESFILE $BAMFILE"
echo $cmd
$cmd

echo "Finished!"


