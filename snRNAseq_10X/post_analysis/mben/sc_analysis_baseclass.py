from pathlib import Path
from os.path import exists
from os import makedirs
import dill, copy, yaml
import scanpy as sc
import sc_analysis_loops as scl
import scfunctions as sc_funcs
import memoize


############################
#### Baseanalysis Class ####
############################

class Baseanalysis:
    def __init__(self, name, type, organism, analysis_params, paths):
        self.name = name
        self.type = type
        self.organism = organism

        # set analysis_params to only the dataset relevant param subset combined from default and other
        self.analysis_params = copy.deepcopy(analysis_params['default'])

        if(self.name in analysis_params['other']):
            self.analysis_params = sc_funcs.merge_dicts(self.analysis_params, analysis_params['other'][self.name])

        self._paths = paths
        self._paths.update(paths)
        middlepath = Path(f'{self._paths["analysisversion"]}/{self.organism}/{self.type}')
        middlepath_dataset = copy.deepcopy(middlepath) / self.name 
        middlepath_data = copy.deepcopy(middlepath_dataset) / 'data'
        middlepath_figures = copy.deepcopy(middlepath_dataset) / 'figures'
        middlepath_results = copy.deepcopy(middlepath_dataset) / 'results'
        middlepath_log = copy.deepcopy(middlepath_dataset) / 'log'
        middlepath_subs = copy.deepcopy(middlepath_dataset) / 'subsets'
        self._paths.update({
            'datasetpath': middlepath_dataset,
            # mounted/data/sn/all
            'datapath': self._paths['storagepath_local'] / middlepath_data})
        self._paths.update({
            'datapath_tmp': self._paths['analysispath_local'] / middlepath_data})
        self._paths.update({
            'picklepath': self._paths['datapath_tmp'] / middlepath / f"{self.name}.h5ad",
            'priorknwldgpath': self._paths['analysispath_local'] / middlepath / 'priorKnowledge',
            'figpath': self._paths['analysispath_local'] / middlepath_figures,
            'resultpath': self._paths['analysispath_local'] / middlepath_results,
            'logpath': self._paths['analysispath_local'] / middlepath_log,
            'subsetspath': self._paths['analysispath_local'] / middlepath_subs})   

        self.data = '' # is set in init of analysis obj

        super().__init__()
    
    def get_paths(self):
        self._paths

    def set_paths(self, paths):
        self._paths.update(paths)

    def save_paths(self):
        """ Saves analysis params together with paths to yaml file. """
        self.analysis_params['default']['paths'] = self.get_paths()
        for data in self.datasets:
            data.analysis_params['other'][data.name]['paths'] <- data.get_paths()
            with open(data.get_paths()['datasetpath'] / 'analysis_params.yaml', 'w+') as file:
                yaml.dump(data.analysis_params, file)

    
        
     
########################
#### Analysis Class ####
######################## 

class Analysis: 
    @memoize.Memoize
    def new_dataset(*bases):
        class Dataset(*bases):
            def __init__(self, name, type, organism, analysis_params, paths):
                """ Set up of new dataset. 

                Parameters
                ----------
                name : str
                type: type of sequencing data, sn, sc etc.
                organism : 'human' or 'mouse'
                """
                super().__init__(name, type, organism, analysis_params, paths)
                
        return Dataset

    @staticmethod
    def read_analysis(proj, basepath_local, analysisversion):
        dill.load(Path(basepath_local, '/projects/', proj, '/results/analysis_', analysisversion, '.pickle'))


    ### Init functions ###
    def __init__(self, proj, datasets, version): 
        self.proj = proj
        self.version = version
        ### Init paths based on params ###
        from . import analysis_params
        self.analysis_params = analysis_params.update(self)
        
        import copy
        self.__paths = copy.deepcopy(self.analysis_params['default']['basepaths'])
        self.__paths.update({    
            "analysisversion": self.version,
            "storagepath_hpc": Path('~/../../mnt/sds-hd/sd22b002'),
            "storagepath_local": self.__paths["basepath_local"] / self.__paths["storagedir"] / 'projects' / self.proj,  # mounted
            "analysispath_local":  self.__paths["basepath_local"] / 'projects' / self.proj
        })
        self.__paths.update({
            "datapath": self.__paths["storagepath_local"] / 'data'
        })
        self.__paths.update({
            "scriptpath": self.__paths["analysispath_local"] / 'scripts',
            "figpath": self.__paths["analysispath_local"] / 'figures',
            "logpath": self.__paths["analysispath_local"] / 'log',
            "resultpath": self.__paths["analysispath_local"] / 'results',
            'priorKnowledge': self.__paths['datapath'] / 'priorKnowledge'
        })

        self.datasets = [constructor(name, type, organism, self.analysis_params, copy.deepcopy(self.__paths)) for name, type, organism, constructor in datasets]        
        self.init_datasets()

    def get_paths(self):
        return self.__paths

    
    def init_datasets (self) :
        @scl.loop(self.datasets, True)
        def init(self):
            """ Read hf5ad data if no pickle exists. Reads into data property of dataset. """
            dirpath = self._paths['datapath'] 
            dirpathlocal = self._paths['datapath_tmp']
            filepath_p = dirpathlocal / f'{self.name}.pickle'
            if not exists(dirpathlocal):
                        makedirs(dirpathlocal)
            if(exists(filepath_p)):
                with open(filepath_p, 'rb') as f:
                    self.data = dill.load(f)
                print('Data was read in from local pickle file.')
            else:
                filepath_raw = dirpath/ f"{self.name}.h5ad"
                try: 
                    self.data = sc.read(filepath_raw, cache = True)
                except OSError as e:
                    # probably a new version number
                    if not exists(self._paths['datapath']):
                        # get previous version number
                        import re
                        datapath = str(self._paths['datapath'])
                        r = re.compile('.*v(..).*')
                        numb = int((r.match(datapath)).group(1)) - 1
                        numb = 'v' + str(numb).zfill(2) 
                        datapath_old = re.sub(r'v..', numb, datapath)                        
                        # take data from old version
                        import distutils.dir_util
                        distutils.dir_util.copy_tree(datapath_old, datapath)
                        # retry
                        self.data = sc.read(filepath_raw, cache = True)

                with open(filepath_p, "wb") as dill_file:
                    dill.dump(self.data, dill_file)
                print('Data was read in from h5ad file and saved as pickle file.')
        init()
        self.clean_datasets()

    ### Processing functions ###       
    def clean_datasets (self) :
        """ correct datatypes """
        @scl.loop(self.datasets, True)
        def clean (dataset) : 
            if 'seurat_clusters' in dataset.data.obs.columns : dataset.data.obs['seurat_clusters'] = dataset.data.obs.seurat_clusters.astype("category")
        clean()

    ### Save & Load ###
    # These functions are for saving an intermediate status as pickle file. 
    def save_datasets(self):
        with open(Path(self.paths['datapath'] / 'datasets.pickle'), "wb") as dill_file:
            dill.dump(self.datasets, dill_file)

    def load_datasets(self):
        return dill.load(Path(self.paths['datapath'] / 'datasets.pickle'))


                     



