# README

This folder includes the analysis input, code and output for the analysis with the tools Decoupler-py and LIANA.    
`mben` is a python package with a main script that initializes the analysis and runs the Decoupler functions.   
`seqDataAnalysisSc` is an R package with helper functions for the Liana analysis.  


**Description of the input data files:**   

- *MB_all_t*: rds file with preprocessed tumor samples that is used for the analysis with Decoupler & Liana.  
- *all_t.h5Seurat*: The file was created from *MB_all_t * 
- *all_t.h5ad*: The file was created from *all_t.h5Seurat*  


**To rerun the analysis, follow these steps:**  

- Initialize the analysis by creating an analysis object from the input data. 
  - Install the needed poetry environment by running `poetry install` in the current directory. 
  - Install the local python package with: `poetry run python -m pip install -e .`
  - Run the first code chunks of *./mben/main.ipynb* . 
- To go forward with the Decoupler analysis, run the rest of *./mben/main.ipynb* .
- To go forward with the Liana analysis, use `renv` for the installations and run *./liana.Rmd*

**Decoupler versions**  
  
- Decoupler-py 1.2.0
- Python 3.9
- Omnipath Dorothea (downloaded 03.08.2022)  
  
For more depencencies see the file *pyproject.toml*.  

  
**LIANA versions**  
  
- Liana 0.1.6 (R version)  
- OmnipathR 3.7.0  
- R 4.2.1
  
For more depencencies see the file *renv.lock*.  

