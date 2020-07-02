#!/usr/bin/env python


__author__ = 'Leland Taylor'
__date__ = '2020-07-03'
__version__ = '0.0.1'

import argparse
#import csv
import os
import random
import numpy as np
import pandas as pd
import scanpy as sc
import plotnine as plt9
from kneed import KneeLocator
from distutils.version import LooseVersion
from scipy.interpolate import UnivariateSpline


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


def estimate_cutoffs_plot(
    output_file,
    df_plt,
    df_knee_output,
    df_fit=None,
    scale_x_log10=False,
    save_plot=True
):
    """Plot UMI counts by sorted cell barcodes."""
    gplt = plt9.ggplot()
    gplt = gplt + plt9.theme_bw()
    if len(df_plt) <= 50000:
        gplt = gplt + plt9.geom_point(
            mapping=plt9.aes(x='barcode', y='umi_counts'),
            data=df_plt,
            alpha=0.1,
            size=0.25
        )
    else:
        gplt = gplt + plt9.geom_line(
            mapping=plt9.aes(x='barcode', y='umi_counts'),
            data=df_plt,
            alpha=0.75,
            size=0.5
        )
    gplt = gplt + plt9.geom_vline(
        mapping=plt9.aes(xintercept='knee', color='method'),
        data=df_knee_output,
        alpha=0.75,
        linetype='dashdot'
    )
    gplt = gplt + plt9.scale_color_brewer(
        palette='Dark2',
        type='qual'
    )
    if scale_x_log10:
        gplt = gplt + plt9.scale_x_continuous(
            trans='log10',
            labels=comma_labels,
            minor_breaks=0
        )
    else:
        gplt = gplt + plt9.scale_x_continuous(
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
        x='Barcode index, sorted by UMI count',
        color='Cutoff method'
    )
    # Add the fit of the droplet utils model
    if df_fit:
        gplt = gplt + plt9.geom_line(
            mapping=plt9.aes(x='x', y='y'),
            data=df_fit,
            alpha=1,
            color='yellow'
        )
    if save_plot:
        gplt.save(
            '{}.png'.format(output_file),
            dpi=300,
            width=5,
            height=4
        )
    return gplt


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

    # parser.add_argument(
    #     '-ncpu', '--number_cpu',
    #     action='store',
    #     dest='ncpu',
    #     default=2,
    #     type=int,
    #     help='Number of CPUs to use.\
    #         (default: %(default)s)'
    # )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='adata',
        help='Basename of output files, assuming output in current \
            working directory. \
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    output_file = 'test'
    verbose = True

    # Load a file of the samples to analyse
    file = options.txd
    #file = "/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/scrna_cellranger/results/iget_cellranger/full_data/blood/cd/Crohns_Disease_Collection_Study8727394/raw_feature_bc_matrix"
    adata = sc.read_10x_mtx(
        path=file,
        # var_names='gene_symbols',
        var_names='gene_ids',
        make_unique=False
    )

    # Make the umi dataframe
    df = pd.DataFrame({
        'umi_counts': np.sort(
            (np.array(adata.X.sum(axis=1))).flatten()
        )[::-1] + 1
    })
    df['barcode'] = df.index + 1

    # Calculate the knee and inflection point as in DropletUtils::barcodeRanks
    # https://github.com/MarioniLab/DropletUtils/blob/master/R/barcodeRanks.R
    # Numerical differentiation to identify bounds for spline fitting.
    # The upper/lower bounds are defined at the plateau and inflection.
    #
    # The lower bound on the total UMI count, at or below which all barcodes
    # are assumed to correspond to empty droplets
    lower_bound = 100
    df_analysis = df.loc[df['umi_counts'] >= lower_bound, :]
    d1n = np.diff(
        np.log10(df_analysis['umi_counts'])
    )/np.diff(np.log10(df_analysis['barcode']))
    right_edge = np.argmin(d1n)
    left_edge = np.argmax(d1n[0:right_edge])

    # We restrict to this region, thereby simplifying the shape of the curve.
    # This allows us to get a decent fit with low df for stable differentiation
    # Smoothing to avoid error multiplication upon differentiation.
    # Minimizing the signed curvature and returning the total for the knee
    # point
    x = np.log10(df_analysis['barcode'][left_edge:right_edge])
    y = np.log10(df_analysis['umi_counts'][left_edge:right_edge])
    dropletutils_fit = UnivariateSpline(x, y)
    fit1 = dropletutils_fit.derivative(n=1)
    fit2 = dropletutils_fit.derivative(n=2)
    curvature = fit1(x) / np.power(1 + np.power(fit2(x), 2), 1.5)
    dropletutils_knee = 10 ** x[np.argmin(curvature)]
    dropletutils_inflection = 10 ** x[right_edge-1]

    kneedle_dict = {}
    cell_estimate_outdict = {
        'dropletutils::knee': {
            'knee': dropletutils_knee,
            'method': 'dropletutils::barcoderanks::knee'
        },
        'dropletutils::inflection': {
            'knee': dropletutils_inflection,
            'method': 'dropletutils::barcoderanks::inflection'
        }
    }

    # Get kneedle estimate for what is a cell
    # Get kneedle estimate using raw data OR fit smoothing spline. Also
    # maximize sensitivity.
    for fit in [None, 'interp1d', 'polynomial']:
        key = 'kneedle::spline={}'.format(fit)
        sensitivity = 1000
        while True:
            if fit:
                kneedle_dict[key] = KneeLocator(
                    x,
                    y,
                    curve='convex',
                    direction='decreasing',
                    S=sensitivity,
                    interp_method=fit
                )
            else:
                kneedle_dict[key] = KneeLocator(
                    x,
                    y,
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
                        round(10 ** kneedle_dict[key].knee, 3),
                        round(10 ** kneedle_dict[key].elbow, 3)
                    ))
                cell_estimate_outdict[key] = {
                    'sensitivity': sensitivity,
                    'knee': 10 ** kneedle_dict[key].knee,
                    'elbow': 10 ** kneedle_dict[key].elbow,
                    'method': key
                }
                break

    # Below is another way to esimate cells based on expected number of cells
    # ...this is what CellRanger v2? uses.
    # https://scrnaseq-course.cog.sanger.ac.uk/website/processing-raw-scrna-seq-data.html
    # n_cells <- length(truth[,1])
    # totals <- umi_per_barcode[,2]
    # totals <- sort(totals, decreasing = TRUE)
    # # 99th percentile of top n_cells divided by 10
    # thresh = totals[round(0.01*n_cells)]/10
    # plot(totals, xlim=c(1,8000))
    # abline(h=thresh, col="red", lwd=2)

    # Make a dataframe of all of the different knee calculations
    # NOTE: here knee is the number of cells that would be kept at that
    # threshold
    df_knee_output = pd.DataFrame(
        cell_estimate_outdict
    ).transpose().reset_index(drop=True)
    # Add the fit of the droplet utils model
    df_fit = pd.DataFrame({
        'y': 10 ** dropletutils_fit(x),
        'x': 10 ** x
    })

    # Make a plot of the different cutoffs
    _ = estimate_cutoffs_plot(
        '{}-cell_estimate_cutoffs-zoomed'.format(output_file),
        df_analysis,
        df_knee_output,
        # df_fit=df_fit,
        scale_x_log10=False
    )
    _ = estimate_cutoffs_plot(
        '{}-cell_estimate_cutoffs'.format(output_file),
        df,
        df_knee_output,
        # df_fit=df_fit,
        scale_x_log10=True
    )

    # Save all of our estimates
    df_knee_output.to_csv(
        '{}-knee_esimates.tsv.gz'.format(output_file),
        sep='\t',
        compression=compression_opts,
        index=False,
        header=True
    )
    with open('{}-expected_cells.txt'.format(output_file), 'w') as f:
        f.write(str(round(
            cell_estimate_outdict['dropletutils::inflection']['knee']
        )))
    # with open('{}-knee_esimate.txt'.format(output_file), 'w') as f:
    #     f.write(str(knee))
    # with open('{}-total_droplets_included.txt'.format(output_file), 'w') as f:
    #     f.write(str(total_droplets_included))

    # Testing out how to set passed via the --total-droplets-included command.
    #
    #
    # Now get an estimate of the knee to identify likely empty cells to pass
    # to cellbender to learn an ambient RNA signature.
    kneedle_dict = {}
    total_drops_estimate_outdict = {}

    # Get kneedle estimate using raw data OR fit smoothing spline. Also
    # maximize sensitivity.
    for fit in [None, 'interp1d', 'polynomial']:
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
                total_drops_estimate_outdict[key] = {
                    'sensitivity': sensitivity,
                    'knee': kneedle_dict[key].knee,
                    'elbow': kneedle_dict[key].elbow,
                    'norm_knee': kneedle_dict[key].norm_knee,
                    'norm_elbow': kneedle_dict[key].norm_elbow,
                    'interp_method': fit
                }
                break

    # Make a dataframe of the knee estimates.
    df_knee_output_total_drops_estimate = pd.DataFrame(
        total_drops_estimate_outdict
    ).transpose().reset_index(drop=True)

    # Save all of our estimates
    df_knee_output_total_drops_estimate.to_csv(
        '{}-knee_total_drops_estimate.tsv.gz'.format(output_file),
        sep='\t',
        compression=compression_opts,
        index=False,
        header=True
    )

    # Select our knee cutoff.
    knee = df_knee_output_total_drops_estimate['spline=None']['knee']
    # Based on the knee, add in a few thousand cells that we think are empty
    # these will be used to build the background predictor in cellbdender
    # and passed via the --total-droplets-included command.
    total_droplets_included = knee + 20000

    with open('{}-knee_esimate.txt'.format(output_file), 'w') as f:
        f.write(str(knee))
    with open('{}-total_droplets_included.txt'.format(output_file), 'w') as f:
        f.write(str(total_droplets_included))

    # Make a plot of the different cutoffs
    gplt = estimate_cutoffs_plot(
        '{}-total_drops_estimate_cutoffs-zoomed'.format(output_file),
        df.loc[df['umi_counts'] >= 10, :],
        df_knee_output_total_drops_estimate,
        # df_fit=df_fit,
        scale_x_log10=False,
        save_plot=False
    )
    gplt = gplt + plt9.geom_vline(
        xintercept=total_droplets_included,
        alpha=1.0,
        linetype='solid'
    )
    gplt.save(
        '{}.png'.format(output_file),
        dpi=300,
        width=5,
        height=4
    )

    gplt = estimate_cutoffs_plot(
        '{}-total_drops_estimate_cutoffs'.format(output_file),
        df,
        df_knee_output_total_drops_estimate,
        # df_fit=df_fit,
        scale_x_log10=True
    )


if __name__ == '__main__':
    main()
