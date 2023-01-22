## OVERVIEW ##

The repository contains source code for the single nuclei data analysis for the study

*D.R. Ghasemi, K.Okonechnikov et al "Compartments in medulloblastoma with extensive nodularity are connected through differentiation along the granular precursor lineage"*

[Bioarxiv link](https://www.biorxiv.org/content/10.1101/2022.09.02.506321v1.abstract)

## DIRECTORY CONTENTS ##

### snRNAsesq_10X ###

_10X single nuclei data analysis_

Processing reads with CellRanger:
run_cellranger.sh

Seurat main analysis:
processTumorSample.R

Mark doublets:
runDecontX_perSample.R

SingleR comparison to reference:
runSingleR_perSample.R

Infer CNV copy number profiling:
runInferCNV_perSample.R

Analysis Decoupler and LIANA:
post_analysis subfolder (includes additional documentation)


### snRNAseq_ss2 ###

_Smart-seq2 snRNA-seq data analysis_

Align reads per cell:
run_STAR.sh

Compute counts per cell:
run_compCounts.sh

Merge the cells of a sample into one matrix:
summarizeComputedCounts.py

Convert gene IDs to names:
renameCounts.py

Seurat main analysis:
processTumorSampleSmartSeq2.R

### spatial_RB ###

_Resolve Bioscence spatial data analysis_

Seurat per sample analysis:
resolveSeurat.R

Giotto additional analysis (cell proximity):
giotto_per_sample.R

Combined analysis for merged cohort:
resolveMergedAnalysis.R











