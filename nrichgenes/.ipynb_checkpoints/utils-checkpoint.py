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