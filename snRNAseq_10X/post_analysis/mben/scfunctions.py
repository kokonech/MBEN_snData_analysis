"""Docstring for the scfunctions.py module.

v1: matches to decoupler_v1, putting the functions into classes

Modules names should have short, all-lowercase names.  The module name may
have underscores if this improves readability.

"""
# for markdown magic
from __future__ import absolute_import
from functools import reduce
from copy import deepcopy
from mergedeep import merge

def deep_get(dictionary, *keys):
    return reduce(lambda d, key: d.get(key) if d else None, *keys, dictionary)

def getpath(nested_dict, value, prepath=()):
    """ Get a tuple of keys that lead to a given value in a given nested dictionary.

    Code taken from:
    [Stackoverflow - answer from jsf](https://stackoverflow.com/questions/22162321/search-for-a-value-in-a-nested-dictionary-python)
    """
    for k, v in nested_dict.items():
        path = prepath + (k,)
        if value in v: # found value
            return path
        elif hasattr(v, 'items'): # v is a dict
            p = getpath(v, value, path) # recursive call
            if p is not None:
                return p

def dict_replace(dict, v, path) -> dict:
    pathlen = len(path)
    r = range(0,pathlen-1)
    (reduce(lambda d,i: d[path[i]], r, dict))[path[pathlen-1]] = v
    return dict



def merge_dicts(dict_1: dict, dict_2: dict) -> dict:
    """First updates dict_2 with dict_1. 
    Then, inserts the contents of dict_1 at the positions of the 'add' entries in the values of dict_2. 

    UseCase
    -------
    dict_1 with default values for an analysis. dict_2 with intended changes to these default values. 
    These can be replacements or additions. Values that should be extended must be placed in a list in both dicts. 
    dict_2 has an 'add' entry wherever the default values from dict_1 should be inserted. 
        
    Parameters
    ----------
    dict_1 : dict
    dict_2 : dict

    Variables & Functions
    ---------------------
    placeholder : 'add'
    keys: all keys till one key value contains the placeholder
    value: the list that contains the placeholder
    ind: index of the placeholder in the list
    getpath() 
    deep_get()
    copy:deepcopy()

    long_var_name : {'hi', 'ho'}, optional
        Choices in brackets, default first when optional.

    Returns
    -------
    merged : dict

    Examples
    --------
    >>> dict_1 = {'pets': {'dog': {'name': 'Bello', 'sound': ['wuff'], 'isplayful': True}, 'cat': {'name': 'Kitty', 'food': 'fish'}}}
    >>> dict_2 = {'pets': {'dog': {'name': 'Wauwau', 'sound': ['wau', 'add', 'grrr']}, 'cat': {'name': 'Miezy', 'sound': 'maunz', 'food': ['add', 'mice']}}}
    >>> merge_dicts(dict_1, dict_2)
    {'pets': {'dog': {'name': 'Wauwau', 'sound': ['wau', 'wuff', 'grrr'], 'isplayful': True}, 'cat': {'name': 'Miezy', 'food': ['fish', 'mice'], 'sound': 'maunz'}}}
    """
    placeholder = 'add'
    keys = getpath(dict_2, placeholder)
    if(keys == None):
        return merge(deepcopy(dict_1), deepcopy(dict_2)) # to get the missing keys from dict_1
    else: 
        value = list(deep_get(dict_2, keys))
        ind = value.index(placeholder)
        value = value[:ind]+ deep_get(dict_1, keys) + value[ind+1:]
        merged = dict_replace(deepcopy(dict_2), value, keys)
        return merge_dicts(dict_1, merged)

