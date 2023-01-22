
import yaml
from pathlib import Path
import os

analysis_params = {
    'default': {
        'basepaths':{
                # path from which all needed files can be reached
                "basepath_local": Path(os.path.abspath('../')),
                "storagedir": Path('')
            },
        'priorKnowledge':{
                'dorothea': 
                    {   # levels
                        'levels': [('A', 'B', 'C')]
                    }  
            },
        'decoupler':{
            'methods': [('mlm', ),],
            'meanacts': {
                'groupby': ['vars'],
                'minstd': 0.0
            }
        },
        'liana': {
            "methods": [["natmi", "connectome", "logfc", "sca", "cellphonedb"]],
            "base": "exp(1)",
            "lig_rec": [["all"]]
        }
    }, 
    'other': {
        'all_t': {
            'original_name': 'MB_all_t',
            'decoupler': {},
            'priorKnowledge':{
                    'dorothea': 
                        {}
                    }
                }
         }
}



def update(analysis):
    """ Updates analysis_params values of the analysis obj according to the values above saves a copy of them as a yaml file. """

    with open(analysis_params['default']['basepaths']["basepath_local"] / 'projects' / analysis.proj / 'analysis_params.yaml', 'w+') as file:
        yaml.dump(analysis_params, file)
    from copy import deepcopy
    return deepcopy(analysis_params)