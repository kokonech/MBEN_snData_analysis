#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l mem=40g
#PBS -l walltime=24:00:00
#PBS -M k.okonechnikov@dkfz-heidelberg.de
#PBS -j oe
#PBS -o /icgc/dkfzlsdf/analysis/dktk/Ependymoma/external_data/logFolder/cellRanger

# human premrna
REFDATA="/icgc/dkfzlsdf/analysis/dktk/Ependymoma/annotations/cell10x/hg19_premrna"
SUFFIX="with_intron"

## MIX

echo "##########" 
echo "runCellRanger.sh" 
echo `date`
echo `hostname`
echo "INPUT: ${INPUT}" 
echo "RESDIR: ${RESDIR}"


SID=$INPUT

echo "Processing $SID"

cd $RESDIR

cmd="/icgc/dkfzlsdf/analysis/dktk/Ependymoma/annotations/cell10x/cellranger-2.1.0/cellranger count --localcores=8 --id=${SID}_${SUFFIX} --fastqs=fq --sample=$SID --transcriptome=$REFDATA"
echo $cmd
$cmd

echo "Finished!"

