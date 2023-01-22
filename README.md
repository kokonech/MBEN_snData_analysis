## DATA ANALYSIS CODE OVERVIEW ##

** 10X analysis per sample **

Processing with CellRanger:
run_cellranger.sh

Seurat main analysis:
processTumorSample.R

Mark doublets:
runDecontX_perSample.R

SingleR comparison to reference:
runSingleR_perSample.R

Infer CNV calling:
runInferCNV_perSample.R

** Smart-seq2 data analysis per sample **

Align reads per cell:
run_STAR.sh

Compute counts per cell:
run_compCounts.sh

Merge the cells of a sample into one matrix:
summarizeComputedCounts.py

Convert gene IDs to names:
renameCounts.py

Seurat standard analysis:
processTumorSampleSmartSeq2.R

** Resolve Bioscence data analysis **

Seurat per sample analysis:
resolveSeurat.R

Combined analysis for merged samples:
resolveMergedAnalysis.R

Giotto additional analysis (cell proximity):
giotto_per_sample.R









