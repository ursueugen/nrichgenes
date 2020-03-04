from pathlib import Path
import pandas as pd
import json
import unittest
import shutil

from nrichgenes import utils, enricher


TEST_DIR = Path('tests')
TEST_DATA = TEST_DIR / 'test_data.json'
TEST_ENRICHMENT = TEST_DIR / 'test_enrichment/'


class TestPackage(unittest.TestCase):
    '''Simple tests to evaluate basic functionality. Not in depth.'''
    
    
    def setUp(self):
        with open(TEST_DATA, 'r') as f:
            data = json.load(f)
        self.genes  = data['genes']
        self.values = data['values']
        self.prerank_kws = data['prerank_kws']
        self.custom_genesets = data["geneset"]
        self.s = pd.Series(self.values, index=self.genes)
        self.df = pd.DataFrame(data=data['array'], columns=data['genes'])
        
        if not TEST_ENRICHMENT.is_dir():
            TEST_ENRICHMENT.mkdir()
        
        
    def test1_download_lookup(self):
        lk = utils.download_lookup()
        self.assertEqual(lk.columns.to_list(), ['ensembl_id', 'name'])
    
    
    def test1_ensembl_to_names(self):
        genes = ["ENSG00000139618", "ENSG00000012048"]
        names = utils.ensembl_to_names(genes)
        self.assertEqual(names, ["BRCA2", "BRCA1"])
    
    
    def test1_prerank_enricher(self):
        genes  = self.genes
        values = self.values
        s      = self.s
        
        enr = enricher.Enricher(genes = s,
                               description = 'test1_prerank_enrichr',
                               outdir = TEST_ENRICHMENT)
        
        enr.prerank(prerank_kws = self.prerank_kws)
        enr.enrichr()
        res = enr.get_results()
        shutil.rmtree(TEST_ENRICHMENT)
        self.assertTrue(set(res.keys()) == {'prerank', 'enrichr'})

    
    def test1_measure_activity(self):
        '''Measure the activity of single process from gmt.'''
        s = self.s
        prerank_kws = self.prerank_kws
        
        enr = enricher.Enricher(genes = s,
                               description = 'test1_measureactivity',
                               outdir = TEST_ENRICHMENT)
        
        enr.measure_activity(processes = ['Adipocytokine signaling pathway'],
                            gmt = 'KEGG_Human_2019',
                            prerank_kws = prerank_kws)
        shutil.rmtree(TEST_ENRICHMENT)
        
    
    def test2_measure_activity(self):
        '''Measure the activity of multiple processes from gmt.'''
        s = self.s
        prerank_kws = self.prerank_kws
        
        enr = enricher.Enricher(genes = s,
                               description = 'test2_measureactivity',
                               outdir = TEST_ENRICHMENT)
        
        enr.measure_activity(processes = ['Adipocytokine signaling pathway',
                                          "AGE-RAGE signaling pathway in diabetic complications",
                                          "Cell adhesion molecules (CAMs)"],
                            gmt = 'KEGG_Human_2019',
                            prerank_kws = prerank_kws)
        shutil.rmtree(TEST_ENRICHMENT)
        
    
    def test3_measure_activity(self):
        '''Measure the activity from predefined gene set.'''
        s = self.s
        prerank_kws = self.prerank_kws
        custom_genesets = self.custom_genesets
        
        enr = enricher.Enricher(genes = s,
                               description = 'test3_measureactivity',
                               outdir = TEST_ENRICHMENT)
        
        enr.measure_activity(genesets = custom_genesets,
                            prerank_kws = prerank_kws)
        
        shutil.rmtree(TEST_ENRICHMENT)
        
    
    def test4_measure_activity(self):
        '''Measure the activity for pd.DataFrame of gene expression from predefined gene set.'''
        
        df = self.df
        
        enr = enricher.Enricher(genes = df,
                               description = 'test4_measureactivity',
                               outdir = TEST_ENRICHMENT)
        
        enr.measure_activity(genesets = custom_genesets,
                            prerank_kws = prerank_kws)
        
#         shutil.rmtree(TEST_ENRICHMENT)
    
    
    def test5_measure_activity(self):
        '''Measure the activity for pd.DataFrame of gene expression of multiple processes from gmt.'''
        pass
    
    
    def test6_measure_activity(self):
        '''Measure the activity for pd.DataFrame of gene expression of all processes from multiple gmt.'''
        pass
    
    
    def test6_measure_activity(self):
        '''Measure the activity for pd.DataFrame of gene expression of multiple processes from multiple gmt.'''
        pass
    
    
if __name__ == '__main__':
    unittest.main()