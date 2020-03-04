import warnings
import pandas as pd
from pathlib import Path
from pybiomart import Dataset


GENE_LOOKUP_PATH = Path("data/ens2name_lookup.csv")


def download_lookup(path = GENE_LOOKUP_PATH) -> pd.DataFrame:
    '''
    Download via biomart lookup table for converting
     gene ids to external ids and names.
     
    :param path
    '''
    
    dataset = Dataset(name='hsapiens_gene_ensembl',
                      host='http://www.ensembl.org')
        
    attrs_selected = ['ensembl_gene_id',
                      'external_gene_name']
    
    q = dataset.query(attributes=attrs_selected)
    
    # short colnames
    q.rename(axis=1,
             inplace=True,
             mapper={
                 "Gene name": "name",
                 "Gene stable ID": "ensembl_id"
             })
    
    q.to_csv(path, index=True, header=True)
    
    return q


def ensembl_to_names(genes, lookup_path: str = GENE_LOOKUP_PATH) -> list:
    '''
    Converts gene ensembl ids to names.
    
    :param genes - list-like of ensembl ids
    :param lookup_path - path to lookup file
    '''
    
    fp = Path(lookup_path)
    
    if fp.is_file():
        lk = pd.read_csv(fp)
    else:
        lk = download_lookup(path=fp)
    
    return lk.set_index('ensembl_id').loc[genes]['name'].to_list()


def get_genesets(gmts: list, terms: list) -> dict:
    '''
    Gets gene sets from .gmt files.
    
    :param gmts - list of paths to .gmt files to be searched
    :param setnames - list of names of gene sets to extract.
    
    :return gsets, e.g. {'adipokine': ['gene1', 'gene2', ...], 'mitochondria': ['gene1', ...], ...}
    '''
    
    gsets = {}
    
    for gmt in gmts:
        
        with open(gmt, 'r') as f:    
            for line in f:
                process = line.split("\t\t")[0]
                if process in terms:
                    gsets[process] = line.split("\t\t")[1].split()
                    
#                     assert gsets[process][-1] != '\n'
#                     warnings.warn("Will remove last gene from last term!", UserWarning)
    
    return gsets