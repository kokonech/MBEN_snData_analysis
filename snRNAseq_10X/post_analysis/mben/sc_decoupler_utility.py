import dill, re
import sc_analysis_baseclass as sc_classes
import sc_analysis_loops as scl
from copy import deepcopy
from pathlib import Path
from os.path import exists
from os import makedirs
import collections
import scanpy as sc, decoupler as dc, seaborn as sns, pandas as pd #, numpy as np, matplotlib.pyplot as plt, seaborn as sns, matplotlib as mpl


class Decoupler:
    """ This class is a wrapper around the tool Decoupler with some comfort functions. """

    def __init__(self):
        """ A Decoupler object has a list of activity objects and additional paths. """
        super().__init__()
        self.acts = []
        self._paths.update({'dc': 'decoupler'})
        self._paths.update({'actsdir': self._paths['resultpath'] / 'log' / self._paths['dc']})
        self._paths.update({'pseudobulkdir': self._paths['figpath'] / 'counts' / self._paths['dc'] / 'pseudobulk'})
        # liana
        self._paths.update({'liana_resdir': self._paths['resultpath'] / 'log' / 'liana'})
        self._paths.update({'liana_figdir': self._paths['figpath'] / 'log' / 'liana'})
        self._paths.update({'liana_aggdir': self._paths['liana_resdir'] / 'aggregated'})

    def get_all_acts(self, new = False):
        """ Create new Activity object for all parameter combinations. """
        for modeltype in self.analysis_params['priorKnowledge']:
            for modelparams in self.analysis_params['priorKnowledge'][modeltype].values():
                for param in modelparams: 
                    for method in self.analysis_params['decoupler']['methods']:
                        model = self._getmodel(modeltype, param)
                        print(f'Activity calculation starts for modeltype **{modeltype}** with parameter **{param}** and method **{method}**')
                        self._get_acts(model, modeltype, param, method, deepcopy(self._paths), new)   

    # priorKnowledge
    def _getmodel(self, modeltype, param):
        """ Get prior knowledge model with the fitting get method of decoupler. """
        # add to paths
        param_name = param
        if type(param_name) == tuple:
            param_name = str.lower(''.join(param_name))

        dirpath = self._paths['priorknwldgpath'] / modeltype 
        filepath = dirpath / f"{param_name}.pickle"
        if(exists(filepath)):
            with open(filepath, 'rb') as file:
                model = dill.load(file)
            print('Prior Knowledge was read in from pickle files.')
        else:
            model = eval('dc.get_' + modeltype + '(self.organism,' + str(param) + ')')
            if not exists(dirpath):
                makedirs(dirpath)
            with open(filepath, "wb") as dill_file:
                dill.dump(model, dill_file)
            print('Prior Knowledge was read in from Omnipath via Decoupler.')
        return model

    def _get_acts(self, model, modeltype, modelparams, methods, paths, new):
        """ Create new Activity object and add it to dataset. """
        # define activity class
        act = sc_classes.Analysis.new_dataset(Activity)

        # example: sn -> all_t -> raw -> decoupler -> Progeny -> 50_ulmmlm_estimate.pickle
        paths.update({'actsdir': paths['actsdir'] / modeltype})
        dirpath = self._paths['actsdir']
        if not exists(dirpath):
                makedirs(dirpath)

        # Define method names
        methods = list(methods) # Decoupler needs list as input, not tuple
        methods = list(methods) 
        methodnames = deepcopy(methods)
        
        if(len(methods) >= 2):
            is_consensus = True
            # add name for consensus method; it is the concatenation of the methods
            cmethod = ''.join(methods)
            cmethod = cmethod + 'consensus'
            methodnames.insert(0, cmethod) # insert needs the names as a list, not tuple
            methods.insert(0, 'consensus')
        else: 
            is_consensus = False

        # Define model param names
        if type(modelparams) == tuple:
            paramname = str.lower(''.join(modelparams))
        else:
            paramname = str(modelparams)#[0])

        # Assumption: if consensus exists then other results exist as well. 
        # Assumption: Either one method or multiple with consensus
        # Caution: Consensus can get out of date when the results that it's based on get updated. 

        def create_actanndata(base, estimatekey, paths):
                    data = dc.get_acts(base, obsm_key= estimatekey) # with estimates and drop estimates from obsm and save acts
                    data.obsm.pop(estimatekey)
                    # Write anndata
                    data.write(paths['acts_data'])
                    return data

        namedmethods = list(zip(methodnames, methods))
        for namedmethod in namedmethods: 
            modelname = modeltype + paramname
            pvalkey = modelname + '_' + namedmethod[0] + '_pvals'
            estimatekey = modelname + '_' + namedmethod[0] + '_estimate'

            actpaths = deepcopy(paths)
            actpath = paths['actsdir'] / f'{paramname}_{namedmethod[0]}'    
            if not exists(actpath):
                makedirs(actpath) 
            paths.update({'acts_data': actpath / 'data.h5ad'}) 
            paths.update({'mean_acts': actpath / 'mean_acts'}) 
            actpaths.update({'acts_data': actpath / 'data.h5ad'})
            actpaths.update({'mean_acts': actpath / 'mean_acts'}) 

            actpaths.update({'acts_estimate': actpath / 'estimate.csv'})
            actpaths.update({'acts_pvals': actpath / 'pvals.csv'})
            
            actpaths.update({'actdir': actpath})

            if(exists(actpaths['acts_data']) and not new): 
                data = sc.read(actpaths['acts_data'], cache = True)
                print('Activitiy objects were read in from h5ad files.')
            elif(exists(actpaths['acts_estimate']) and not new): # Assumption: When there is an estimate file there is a pvals file, too.
                pvals = pd.read_csv(actpaths['acts_pvals'], index_col=0)
                estimates = pd.read_csv(actpaths['acts_estimate'], index_col=0)
                self.data.obsm[estimatekey] = estimates
                self.data.obsm[pvalkey] = pvals
                data = create_actanndata(self.data, estimatekey, actpaths)
            else:
                if(namedmethod[1] == 'consensus'): 
                    methods.pop(0) # delete the 'consensus' entry
                    new = False # When consensus is newly calculated 'new' is set to false so that the single methods that are part of the consensus calculation don't get recalculated. 
  
                dc.decouple(mat = self.data, net = model, methods=methods, consensus = is_consensus, use_raw = True) 
                # rename new obsm properties (consensus and all the rest)
                self.data.obsm[pvalkey] = self.data.obsm.pop(namedmethod[1] + '_pvals')
                self.data.obsm[estimatekey] = self.data.obsm.pop(namedmethod[1] + '_estimate')

                # Create anndata obj with activity as .X
                data = create_actanndata(self.data, estimatekey, actpaths)

                # save estimates and pvals separately, delete results from dataset obj
                pvals = self.data.obsm.pop(pvalkey)
                estimates = self.data.obsm.pop(estimatekey)
                pvals.to_csv(actpaths['acts_pvals'])
                estimates.to_csv(actpaths['acts_estimate'])
            # Add activity                
            self.acts.append(act(data, modeltype, modelparams, namedmethod, deepcopy(actpaths)))
                

class Activity(): 
    """ An Activity is a Dataset that holds an Anndata object that stores activity estimations in its X property. These estimations are derived from the tool Decoupler. """
    def __init__(self, data, modeltype, modelparams, namedmethod, paths):
        """ Adds paths and reads in anndata object if it already exists or creates it newly via Decoupler. """
        super().__init__()
        self.data = data
        self.modeltype = modeltype
        self.modelparams = modelparams
        self.namedmethod = namedmethod
        self.paths = paths
        self.mean_acts = {}
        
    
    def get_mean_acts(self, minstd, groupby):
        meanactname = f'{groupby}_{minstd}'
        dirpath = self.paths['mean_acts']
        meanactpath = self.paths['mean_acts'] / f'{meanactname}.csv'
        if(exists(meanactpath)):
            self.mean_acts.update({meanactname: pd.read_csv(meanactpath, index_col=0)})
        else:
            obsvalues = self.data.obs[groupby]
            counts = collections.Counter(obsvalues)
            if(len(counts) is not sum(counts.values())):
                if(len(counts) <= 20): 
                    if not exists(dirpath):
                        makedirs(dirpath)
                    dcresult = dc.summarize_acts(self.data, groupby=groupby, min_std=minstd)
                    if(not dcresult.empty):
                        self.mean_acts.update({meanactname: dcresult})
                        self.mean_acts[meanactname].to_csv(meanactpath)

    # The first column is called Unnamed: 0  because these are the row names that were saved beforehand
    def plot_mean_acts(self):
        dirpath = self.paths['mean_acts']
        dirpath = str(dirpath).replace('results', 'figures')
        if not exists(dirpath):
            makedirs(dirpath)

        @scl.loop(self.mean_acts, True)
        def __plot_mean_acts(meanactname, self): 
            meanactpath = str(Path(dirpath,f'{meanactname}.pdf')) # str / str doesn't work, when using / one of the elements must be of type Path
            mean_acts = self.mean_acts[meanactname]

            # plot factor levels
            obsname = re.search('(.*)_.*', meanactname).group(1)
            obsvalues = self.data.obs[obsname]
            counts = collections.Counter(obsvalues)

            if(len(counts) <= 20): 
                if(len(counts.keys()) >= 2): 
                    if(exists(meanactpath)):
                        import glob
                        import numpy
                        #import matplotlib.pyplot

                        filenames = sorted(glob.glob(meanactpath))
                        s = ''
                        for filename in filenames: 
                            s = s + str('\n###' + obsname + '  \n')
                            s = s + str('<p><img src="'+ filename+ '" style="width:500px;height:600px;"></p>')
                    else:
                        sns.set(font_scale=0.5)
                        sns.set(rc={'figure.figsize':(15,20)})
                        if(len(mean_acts.columns) >= 2):
                            fig = sns.clustermap(mean_acts, xticklabels=(mean_acts).columns, vmin=-2, vmax=2, cmap='coolwarm')
                            fig = fig.fig
                            if not len(mean_acts.columns) >= 50: 
                                fig.set_figwidth(len(mean_acts.columns))
                            else:
                                fig.set_figwidth(100)
                            fig.suptitle(obsname) 
                            fig.savefig(meanactpath)
                        else: 
                            fig = sns.heatmap(mean_acts, xticklabels=(mean_acts).columns, vmin=-2, vmax=2, cmap='coolwarm')
        __plot_mean_acts(self)

