#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-03-13'
__version__ = '0.0.1'

import argparse
import csv
import os
import random
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import plotnine as plt9
from kneed import KneeLocator

# import matplotlib
# import matplotlib.pyplot as plt

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)

# sc verbosity: errors (0), warnings (1), info (2), hints (3)
# sc.settings.verbosity = 3
# sc.logging.print_versions()
# sc.settings.set_figure_params(dpi=80)


def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Filter and merge 10x data. Save to AnnData object.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-txd', '--tenxdata_path',
        action='store',
        dest='txd',
        required=True,
        help='Path to 10x data.'
    )

    parser.add_argument(
        '-ncpu', '--number_cpu',
        action='store',
        dest='ncpu',
        default=2,
        type=int,
        help='Number of CPUs to use.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='adata',
        help='Basename of output anndata file, assuming output in current \
            working directory. Will have .h5ad appended.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save = 300)

    output_file = 'test'
    verbose = True

    # Load a file of the samples to analyse
    file = options.txd
    file = "/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/results/iget_cellranger/full_data/blood/cd/Crohns_Disease_Collection_Study8727394/raw_feature_bc_matrix"
    adata = sc.read_10x_mtx(
        path=file,
        # var_names='gene_symbols',
        var_names='gene_ids',
        make_unique=False
    )

    df = pd.DataFrame({
        'umi_counts': np.sort(
            (np.array(adata.X.sum(axis=1))).flatten()
        )[::-1] + 1
    })
    df['barcode'] = df.index + 1

    kneedle_dict = {}
    output_dict = {}

    # Get kneedle estimate using raw data OR fit smoothing spline. Also
    # maximize sensitivity.
    # NOTE: don't fit polynomial makes no sense for our data
    # for fit in [None, 'interp1d', 'polynomial']:
    for fit in [None, 'interp1d']:
        key = 'spline={}'.format(fit)
        sensitivity = 1000
        while True:
            if fit:
                kneedle_dict[key] = KneeLocator(
                    df['barcode'],
                    df['umi_counts'],
                    curve='convex',
                    direction='decreasing',
                    S=sensitivity,
                    interp_method=fit
                )
            else:
                kneedle_dict[key] = KneeLocator(
                    df['barcode'],
                    df['umi_counts'],
                    curve='convex',
                    direction='decreasing',
                    S=sensitivity
                )
            if kneedle_dict[key].knee is None:
                sensitivity -= 5
            else:
                if verbose:
                    print('S:\t{}\nknee:\t{}\nelbow:\t{}'.format(
                        sensitivity,
                        round(kneedle_dict[key].knee, 3),
                        round(kneedle_dict[key].elbow, 3)
                    ))
                output_dict[key] = {
                    'sensitivity': sensitivity,
                    'knee': kneedle_dict[key].knee,
                    'elbow': kneedle_dict[key].elbow,
                    # 'norm_knee': kneedle_dict[key].norm_knee,
                    # 'norm_elbow': kneedle_dict[key].norm_elbow,
                    'interp_method': fit
                }
                break
    output_df = pd.DataFrame(output_dict).transpose().reset_index(drop=True)

    gplt = plt9.ggplot(df, plt9.aes(x='barcode', y='umi_counts'))
    gplt = gplt + plt9.theme_bw()
    gplt = gplt + plt9.geom_line()
    gplt = gplt + plt9.geom_vline(
        xintercept=output_dict['spline=None']['knee'],
        linetype='dotted',
        color='grey',
        size=0.5
    )
    gplt = gplt + plt9.scale_x_continuous(
        trans='log10',
        labels=comma_labels,
        minor_breaks=0
    )
    gplt = gplt + plt9.scale_y_continuous(
        trans='log10',
        labels=comma_labels,
        minor_breaks=0
    )
    gplt = gplt + plt9.labs(
        title='',
        y='UMI counts',
        x='Barcode index, sorted by UMI count'
    )
    gplt.save(
        '{}.png'.format(output_file),
        dpi=300,
        width=5,
        height=4
    )
    # fig, ax = plt.subplots(figsize=(10, 7))
    #
    # ax.loglog(knee, range(len(knee)), linewidth=5, color="g")
    # ax.axvline(x=knee[expected_num_cells], linewidth=3, color="k")
    # ax.axhline(y=expected_num_cells, linewidth=3, color="k")
    #
    # ax.set_xlabel("UMI Counts")
    # ax.set_ylabel("Set of Barcodes")
    #
    # plt.grid(True, which="both")
    # plt.show()


if __name__ == '__main__':
    main()
