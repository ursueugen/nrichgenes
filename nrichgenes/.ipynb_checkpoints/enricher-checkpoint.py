import requests
from pathlib import Path
import warnings

import pandas as pd
import gseapy as gp
from .utils import ensembl_to_names


GENESETS_DIR  = Path('data/genesets')
GENESETS_URLS = {'KEGG_Human_2019': "http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Human",
                 'GO_Process_2018': 'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018'}
GENESETS_ENRICHR = ['KEGG_2019_Human', 'GO_Biological_Process_2018', 'GO_Molecular_Function_2018'] 


class Enricher(object):
    '''
    Class for performing enrichment analysis.
    
    :attributes
        gene_list;
        gene_rnk;
        description;
        outdir;
    :methods
        enrichr();
        prerank();
        get_results();
    ____________________________________________________________
    https://gseapy.readthedocs.io/en/master/gseapy_example.html
    https://gseapy.readthedocs.io/en/master/gseapy_tutorial.html    
    '''
    
    
    def __init__(self, genes: list or pd.Series or pd.DataFrame, description: str, outdir: str):
        '''
        :param genes
        :param description
        :param outdir
        '''
        
        # Resolve type of genes argument
        if type(genes) == list:
            self.gene_list = genes
            self.gene_rnk  = None
            self.gene_df   = None
        elif type(genes) == pd.Series:
            self.gene_list = genes.index.to_list()
            self.gene_rnk  = genes
            self.gene_df   = None
        elif type(genes) == pd.DataFrame:
            self.gene_list = genes.columns.to_list()
            self.gene_rnk  = None
            self.gene_df   = genes
        else:
            raise TypeError("Type of 'genes' is not recognized.")
            
        if self.gene_list[0][:3] == 'ENS':
            # assume ensembl ids as gene labels
            self.gene_names = self._ens_to_names(self.gene_list)
            if self.gene_rnk:
                self.gene_rnk.index = self._ens_to_names(self.gene_rnk.index)
        
        self.description = description
        self.outdir = str(outdir)
    
        self.enrichr_gene_sets = GENESETS_ENRICHR
        self.local_gene_sets = self._get_genesets(path=GENESETS_DIR, urls=GENESETS_URLS)
        
        self.enrichr_flag  = False
        self.prerank_flag  = False
        self.activity_flag = False
        

    def get_results(self) -> dict:
        '''Returns dict of results of enrichr and prerank analysis.'''
        
        if not (self.activity_flag & self.prerank_flag & self.enrich_flag):
            raise LookupError("No results are available yet.")
        
        res = {}
        
        if self.enrichr_flag:
            res['enrichr'] = self.enrichr_res
        if self.prerank_flag:
            res['prerank'] = self.prerank_res
        if self.activity_flag:
            res['activity'] = self.measured_activity
            
        return res
    
    
    def enrichr(self, enrichr_kws: dict = {}):
        '''
        Performs enrichment analysis. Returns enrichr object.
        '''
        
        gene_names = self.gene_list
            
        enr = gp.enrichr(gene_list   = gene_names,
                         description = self.description,
                         gene_sets   = self.enrichr_gene_sets,
                         outdir      = str(self.outdir),
                         **enrichr_kws)
        
        # enr.results
        self.enrichr_res = enr
        self.enrichr_flag = True
#         self.enrichr_paths = {gs: self.outdir + "/" + ".".join([gs, self.description, 'enrichr.reports.txt'])
#                              for gs in self.gene_sets}
        
        return enr
    
    
    def prerank(self, gene_sets: None or dict = None, prerank_kws: dict = {}):
        '''
        Performs prerank analysis on the input genes. 
        
        :param gene_sets   - dict with gene sets, e.g. {'name': {'set1': ['gene1', ...], 'set2': ['gene1', ...]}}
        :param prerank_kws - keyword args for prerank function of gseapy
        '''    
        
        rnk = self.gene_rnk
        res = {}
        
        if gene_sets:
            gsets = gene_sets
        else:
            gsets = self.local_gene_sets
        
        # change variable name for v, misleading!
        for gset, v in gsets.items():
            
            # need to put the results in different directories to avoid overwriting
            warnings.warn("For multiple processes, the directory files are overwritten!", UserWarning)
            
            prerank   = gp.prerank(rnk = rnk, 
                                   gene_sets=v,
                                   outdir = str(self.outdir),
                                   format = 'png',
                                   **prerank_kws)
            res[gset] = prerank
            
        self.prerank_res = res
        self.prerank_flag = True
        return res
    
    
    def measure_activity(self, processes: list = None,
                         gmt: str = None,
                         genesets: dict = None,
                         prerank_kws: dict = {}):
        '''
        Measure the enrichment score for a process.
        
        :param processes - has to be valid process in selected gmt
        :param gmt - in {'GO_Process_2018', 'KEGG_Human_2019'}
        :param genesets
        
        :return prerank object?
        '''
        
        rnk = self.gene_rnk
        df  = self.gene_df
        
        if processes:
            
            if not gmt:
                raise TypeError("If input processes, need to indicate desired gmt to be searched.")
            elif gmt not in ['GO_Process_2018', 'KEGG_Human_2019']:
                raise ValueError("gmt has to be in 'GO_Process_2018', 'KEGG_Human_2019'")
            
            fp  = self.local_gene_sets[gmt]
            name = "-".join([gmt, 'selectedterms'])
            gsets = {name: {}}
            with open(fp, 'r') as f:
                for line in f:
                    process = line.split("\t\t")[0]
                    if process in processes:
                        gsets[name][process] = line.split("\t\t")[1].split("\t")[:-1]  # exclude ['\n'], will bias last line
                        warnings.warn("Will remove last gene from last term!", UserWarning)
            
        elif genesets:        
            if len(genesets.keys()) != 1:
                raise ValueError("Genesets dict has to be of format: {'name': {'set1': ['gene1', ...], 'set2': ['gene1', ...]}}.")
            gsets = genesets
        
        else:
            raise ValueError("Require some input. Either genesets or (processes and gmt).")
        
        raise NotImplementedError("ADD HANDLING OF DATAFRAMES FOR GENE RANKINGS!")
        prerank = self.prerank(gene_sets = gsets, prerank_kws = prerank_kws)
        
        self.measured_activity = prerank
        self.activity_flag = True
        return prerank
        
    
    def _get_genesets(self, path, urls):
        '''Downloads .gmt gene sets. Returns dict of paths.'''
    
        gmt_paths = {}  # store paths to gmt files

        if not path.is_dir():
            path.mkdir()

        for gset, url in urls.items():

            fp = path / ".".join([gset, 'gmt'])

            if fp.is_file():
                gmt_paths[gset] = str(fp)
                continue

            r = requests.request(url=url, method='GET')
            if r:    
                with open(fp, 'w+') as f:
                    f.write(r.text)
                gmt_paths[gset] = (fp)

        return gmt_paths
    
    
    def _ens_to_names(self, genelist: list):
        return ensembl_to_names(genelist)