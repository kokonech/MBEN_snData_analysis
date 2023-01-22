# README

The intention of the code in this package is to wrap multiple datasets into one analysis object. 
Both, the analysis object and the datasets have path and analysis parameter definitions as properties. 
Depending on the intended analysis, the datasets can inherit from multiple classes. In this case they inherit from the 'Decoupler' class.
A loop decorator helps with executing functions over multiple datasets and parameter combinations. 
Results are saved in a standardised folder structure. 

## Files

*sc_analysis_baseclass.py:* Holds the analysis class and the basic dataset class.  
*sc_analysis_loops.py:* Wraps the functions into loops that go over all parameter combinations.  
*sc_decoupler_utility.py:* Decoupler class that a dataset can inherit from to receive additional methods, paths and properties for a smooth analysis with decoupler.  
*scfunctions.py:* Any useful function that doesn't make sense in a class context.  
*analysis_params.py:* Config file for all datasets and analysis tools.


## Analysis Object 
  
Initialization:   
  
- The analysis_params from the file are set as property of the analysis.   
- Paths are generated from the basepaths information in analysis_params and then set as 'paths' property of the analysis.   
- Paths and analysis_params are forwarded to the constructor of the datasets. Each dataset holds an extended copy of these properties.   
- A combination of analysis_params and paths is saved as yaml file for each dataset. Therefore, 'paths' gets added as key to the dictionary.  


## Variable Naming   

Plural = Shortcut + 's'  
  
- col     =   column  
- val, v  =   value  
- k       =   key  
- ind, i  =   index (used in loops)
- elem, e =   element (used in loops)
- meta    =   metadata (obs)  
- ds      =   dataset ('Dataset' object)
- param   =   parameter  
- subs (Pl. subsets) = subset  
- cond    =   condition (for comparisons)  
- relAb   =   relative abundance






