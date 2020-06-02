#!/usr/bin/env python


__author__ = 'Monika Krzak'
__date__ = '2020-05-29'
__version__ = '0.0.1'

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import os



def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object and markers from database. Plot expression \
            expression of markers in each cluster.
            """
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '--markers_database',
        action='store',
        dest='markers_database',
        default='none',
        help='Markers from database to be dotplotted.'
    )

    options = parser.parse_args()

    adata = sc.read_h5ad(filename=options.h5)

    markers=pd.read_table(options.markers_database, sep='\t', header=None)[0]

    
    for i in markers:
    
        base=os.path.basename(os.path.dirname(i))
        df = pd.read_table(i).copy()

        if 'p_value_adj' in df.columns: 
            df = df[['cell_type','hgnc_symbol', 'p_value_adj']].copy()
            df = df.sort_values(by=['cell_type', 'p_value_adj'], ascending=[True, True]).copy()
        else:
            df = df[['cell_type','hgnc_symbol']].copy()
      
        
        cell_type=np.unique(df['cell_type']).copy()
    
   
        for j in cell_type:
                
            marker_genes = df.loc[df['cell_type'] == j]['hgnc_symbol'].copy()
            marker_genes_found = adata.var['gene_symbols'][adata.var['gene_symbols'].isin(marker_genes)].copy()
        
    
            if marker_genes_found.size>100: 
                marker_genes_found = marker_genes_found[0:100]

            
            j_out = j.replace(' ', '_').replace('.', '_').replace('/', '-')

            try:   
                sc.pl.dotplot(adata, marker_genes_found,
                        groupby='leiden', 
                        gene_symbols='gene_symbols',
                        dendrogram=True,
                        show=False,
                        save=base + ('-{}.png'.format(j_out)))
            except:
                pass
        
        
    
    

if __name__ == '__main__':
    main()
