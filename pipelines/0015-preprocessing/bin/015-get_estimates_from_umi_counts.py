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

# To resolve strange TclError for interactive job
import matplotlib
matplotlib.use('Agg')  # Agg for png and pdf for pdf

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
os.environ['PYTHONHASHSEED'] = str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)

# Set valid methods for estimators
valid_methods = [
    'dropletutils::barcoderanks::inflection',
    'dropletutils::barcoderanks::knee',
    'kneedle::spline=None',
    'kneedle::spline=interp1d',
    'kneedle::spline=polynomial'
]



def comma_labels(x_list):
    """Change list of int to comma format."""
    result = []
    for x in x_list:
        result.append(format(int(x), ','))
    return(result)


def estimate_cutoffs_plot(
    output_file,
    df_plt,
    df_cell_estimate_cutoff,
    df_fit=None,
    scale_x_log10=False,
    save_plot=True
):
    """Plot UMI counts by sorted cell barcodes."""
    if min(df_plt['umi_counts']) <= 0:
        fix_log_scale = min(df_plt['umi_counts']) + 1
        df_plt['umi_counts'] = df_plt['umi_counts'] + fix_log_scale
    gplt = plt9.ggplot()
    gplt = gplt + plt9.theme_bw()
    if len(df_plt) <= 50000:
        gplt = gplt + plt9.geom_point(
            mapping=plt9.aes(x='barcode', y='umi_counts'),
            data=df_plt,
            alpha=0.05,
            size=0.1
        )
    else:
        gplt = gplt + plt9.geom_line(
            mapping=plt9.aes(x='barcode', y='umi_counts'),
            data=df_plt,
            alpha=0.25,
            size=0.75,
            color='black'
        )
    gplt = gplt + plt9.geom_vline(
        mapping=plt9.aes(xintercept='n_cells', color='method'),
        data=df_cell_estimate_cutoff,
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
        color='Cutoff'
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


def dropletutils_cutoff(df_analysis):
    """Calculate the knee and inflection point as in DropletUtils."""
    # Calculate the knee and inflection point as in DropletUtils::barcodeRanks
    # https://github.com/MarioniLab/DropletUtils/blob/master/R/barcodeRanks.R
    # Numerical differentiation to identify bounds for spline fitting.
    # The upper/lower bounds are defined at the plateau and inflection.
    #
    # The lower bound on the total UMI count, at or below which all barcodes
    # are assumed to correspond to empty droplets
    # df_analysis = df.loc[df['umi_counts'] >= lower_bound, :]
    if min(df_analysis['umi_counts']) <= 0:
        fix_log_scale = min(df_analysis['umi_counts']) + 1
        df_analysis['umi_counts'] = df_analysis['umi_counts'] + fix_log_scale
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
    x = np.log10(df_analysis['barcode'][left_edge:right_edge].values)
    y = np.log10(df_analysis['umi_counts'][left_edge:right_edge].values)
    dropletutils_fit = UnivariateSpline(x, y)
    fit1 = dropletutils_fit.derivative(n=1)
    fit2 = dropletutils_fit.derivative(n=2)
    curvature = fit1(x) / np.power(1 + np.power(fit2(x), 2), 1.5)
    dropletutils_knee = 10 ** y[np.argmin(curvature)]
    dropletutils_inflection = df_analysis['umi_counts'][right_edge]
    dropletutils_knee_ncells = 10 ** x[np.argmin(curvature)]
    dropletutils_inflection_ncells = df_analysis['barcode'][right_edge]

    result = [
        {
            'method': 'dropletutils::barcoderanks::knee',
            'umi_counts_cutoff': dropletutils_knee,
            'n_cells': dropletutils_knee_ncells
        },
        {
            'method': 'dropletutils::barcoderanks::inflection',
            'umi_counts_cutoff': dropletutils_inflection,
            'n_cells': dropletutils_inflection_ncells
        }
    ]

    return right_edge, left_edge, result


def kneedle_cutoff(df_analysis, verbose=True):
    """Get kneedle estimate for what is a cell."""
    results = []
    # Get kneedle estimate using raw data OR fit smoothing spline. Also
    # maximize sensitivity.
    # NOTE: Much better estimates with kneedle (that is similar to cellranger3)
    # if this is done the original value space - not log space.
    kneedle_dict = {}
    for fit in [None, 'interp1d', 'polynomial']:
        key = 'kneedle::spline={}'.format(fit)
        sensitivity = 5000
        while True:
            if fit:
                kneedle_dict[key] = KneeLocator(
                    df_analysis['barcode'].values,
                    df_analysis['umi_counts'].values,
                    curve='convex',
                    direction='decreasing',
                    S=sensitivity,
                    interp_method=fit
                )
            else:
                kneedle_dict[key] = KneeLocator(
                    df_analysis['barcode'].values,
                    df_analysis['umi_counts'].values,
                    curve='convex',
                    direction='decreasing',
                    S=sensitivity
                )
            if kneedle_dict[key].knee is None:
                sensitivity -= 100
            else:
                if verbose:
                    print('S:\t{}\nknee:\t{}\nelbow:\t{}'.format(
                        sensitivity,
                        round(kneedle_dict[key].knee, 3),
                        round(kneedle_dict[key].elbow, 3)
                    ))
                results.append({
                    'method': key,
                    'umi_counts_cutoff': df_analysis['umi_counts'][
                        int(kneedle_dict[key].knee)
                    ],
                    'n_cells': kneedle_dict[key].knee,
                    'sensitivity': sensitivity
                    # 'elbow': 10 ** kneedle_dict[key].elbow  # same as knee
                })
                break
    return results


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Estimates cutoffs for cellbender remove background from total umi
            counts. NOTE: different cutoff estimators work better on different
            tissues/protocols that have slightly different UMI count plots.
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
        '--expected_ncells',
        action='store',
        dest='expected_ncells',
        default=0,
        type=int,
        help='Expected number of cells. If != 0, this number of cells \
            will be used for expected_cells.txt. (default: %(default)s).'
    )

    parser.add_argument(
        '--total_ndroplets',
        action='store',
        dest='total_ndroplets',
        default=0,
        type=int,
        help='Expected total number of droplets. If != 0, this number of \
            will be used for total_droplets_included.txt. \
            (default: %(default)s).'
    )

    parser.add_argument(
        '--method_expected_ncells',
        action='store',
        dest='method_expected_ncells',
        default='dropletutils::barcoderanks::inflection',
        help='Method to use to estimate number of cells if \
            expected_ncells == 0. Valid options = [{}]. \
            (default: %(default)s)'.format(','.join(valid_methods))
    )

    parser.add_argument(
        '--method_total_ndroplets',
        action='store',
        dest='method_total_ndroplets',
        default='dropletutils::barcoderanks::inflection',
        help='Method to use to estimate total number of droplets if \
            total_ndroplets == 0. Valid options = [{}]. \
            (default: %(default)s)'.format(','.join(valid_methods))
    )

    parser.add_argument(
        '--lower_bound_expected_ncells',
        action='store',
        dest='lb_expected_ncells',
        default=100,
        type=int,
        help='The lower bound on the total UMI count at or below which all \
            barcodes are assumed to correspond to empty droplets (for \
            estimating expected ncells). \
            (default: %(default)s)'
    )

    parser.add_argument(
        '--lower_bound_total_droplets_included',
        action='store',
        dest='lb_total_droplets_included',
        default=10,
        type=int,
        help='The lower bound on the total UMI count used to estimate cells \
            with ambient RNA that will be included in cellbender to \
            estimate the background. \
            (default: %(default)s)'
    )

    parser.add_argument(
        '--total_ndroplets_subtract_factor',
        action='store',
        dest='total_ndroplets_subtract_factor',
        default=2500,
        type=int,
        help='Subtract this number from total number of droplets \
            (either estimated or passed directly) to get the final \
            total_ndroplets. (default: %(default)s)'
    )

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

    # Get clean names
    expected_ncells = options.expected_ncells
    total_ndroplets = options.total_ndroplets
    lb_cell_estimate = options.lb_expected_ncells
    lb_total_droplets_included = options.lb_total_droplets_included
    output_file = options.of
    verbose = True

    # Check the methods
    if options.method_expected_ncells not in valid_methods:
        raise Exception('Invalid method_expected_ncells.')
    if options.method_total_ndroplets not in valid_methods:
        raise Exception('Invalid method_total_ndroplets.')

    # Load a file of the samples to analyse
    file = options.txd
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
        )[::-1]
    })
    df['barcode'] = df.index + 1

    # Make an analysis df using lower bound
    df_analysis = df.loc[df['umi_counts'] >= lb_cell_estimate, :]
    # Get dropletutils estimates
    right_edge, left_edge, result_dropletutils = dropletutils_cutoff(
        df_analysis
    )
    # Get kneedle estimates
    result_kneedle = kneedle_cutoff(df_analysis)
    # Save the results to a list we will turn into a pandas matrix
    cell_estimate_outdict = []
    for i in result_dropletutils:
        cell_estimate_outdict.append(i)
    for i in result_kneedle:
        cell_estimate_outdict.append(i)

    # If we have expected number of cells, then estimate like CellRanger2
    # assumes a ~10-fold range of library sizes for real cells and estimates
    # this range using the expected number of cells.
    # https://scrnaseq-course.cog.sanger.ac.uk/website/processing-raw-scrna-seq-data.html
    if expected_ncells != 0:
        # 99th percentile of top n_cells divided by 10
        cellranger_expected_knee = df['umi_counts'][
            round(0.01*expected_ncells)
        ]/10
        cell_estimate_outdict.append({
            'method': 'cellrangerv2::expected_ncells',
            'umi_counts_cutoff': cellranger_expected_knee,
            'n_cells': (df['umi_counts'] >= cellranger_expected_knee).sum()
        })
        cell_estimate_outdict.append({
            'method': 'expected_ncells',
            'umi_counts_cutoff': df['umi_counts'][expected_ncells],
            'n_cells': expected_ncells
        })

    # Make a dataframe of all of the different knee calculations
    # NOTE: here knee is the number of cells that would be kept at that
    # threshold
    df_cell_estimate_cutoff = pd.DataFrame(
        cell_estimate_outdict
    )
    # Add the fit of the droplet utils model
    # df_fit = pd.DataFrame({
    #     'y': 10 ** dropletutils_fit(x),
    #     'x': 10 ** x
    # })

    # Make a plot of the different cutoffs
    _ = estimate_cutoffs_plot(
        '{}-cell_estimate_cutoffs-zoomed'.format(output_file),
        df_analysis,
        df_cell_estimate_cutoff,
        # df_fit=df_fit,
        scale_x_log10=False
    )
    _ = estimate_cutoffs_plot(
        '{}-cell_estimate_cutoffs'.format(output_file),
        df,
        df_cell_estimate_cutoff,
        # df_fit=df_fit,
        scale_x_log10=True
    )

    # Save all of our estimates
    df_cell_estimate_cutoff.to_csv(
        '{}-cell_estimate_cutoff.tsv.gz'.format(output_file),
        sep='\t',
        compression=compression_opts,
        index=False,
        header=True
    )
    # Get the expected number of cells based from the umi plot.
    # By default this is 'dropletutils::barcoderanks::inflection'
    method_expected_ncells = options.method_expected_ncells
    if expected_ncells != 0:
        method_expected_ncells = 'expected_ncells'
    filt = df_cell_estimate_cutoff['method'] == method_expected_ncells
    expected_cells = df_cell_estimate_cutoff.loc[filt, 'n_cells'].values[0]
    with open('{}-expected_cells.txt'.format(output_file), 'w') as f:
        f.write(str(int(expected_cells)))

    # Run a similar proceedure but for total-droplets-included.
    # ...these will be the lower bound of the droplets used to estimate the
    # background signals
    # Make an analysis df using lower bound
    df_analysis = df.loc[df['umi_counts'] >= lb_total_droplets_included, :]
    # Get dropletutils estimates
    right_edge, left_edge, result_dropletutils = dropletutils_cutoff(
        df_analysis
    )
    # Get kneedle estimates
    result_kneedle = kneedle_cutoff(df_analysis)
    # Save the results to a list we will turn into a pandas matrix
    total_droplets_estimate_outdict = []
    for i in result_dropletutils:
        total_droplets_estimate_outdict.append(i)
    for i in result_kneedle:
        total_droplets_estimate_outdict.append(i)

    # If we have an expected total number of droplets use that.
    if total_ndroplets != 0:
        total_droplets_estimate_outdict.append({
            'method': 'expected_total_ndroplets',
            'umi_counts_cutoff': df['umi_counts'][total_ndroplets],
            'n_cells': total_ndroplets
        })

    df_total_droplets_cutoff = pd.DataFrame(
        total_droplets_estimate_outdict
    )

    # Save all of our estimates
    df_total_droplets_cutoff.to_csv(
        '{}-total_droplets_cutoff.tsv.gz'.format(output_file),
        sep='\t',
        compression=compression_opts,
        index=False,
        header=True
    )
    # Get the expected number of cells based from the umi plot.
    # By default this is 'dropletutils::barcoderanks::inflection'
    method_total_ndroplets = options.method_total_ndroplets
    if total_ndroplets != 0:
        method_expected_ncells = 'expected_total_ndroplets'
    filt = df_total_droplets_cutoff['method'] == method_total_ndroplets
    total_droplets_cutoff_knee = df_total_droplets_cutoff.loc[
        filt, 'n_cells'
    ].values[0]
    # Based on the knee, add in a few thousand cells that we think are empty
    # these will be used to build the background predictor in cellbdender
    # and passed via the --total-droplets-included command.
    total_droplets_included = (
        total_droplets_cutoff_knee - options.total_ndroplets_subtract_factor
    )
    # Save total_droplets_included
    with open('{}-total_droplets_included.txt'.format(output_file), 'w') as f:
        f.write(str(int(total_droplets_included)))

    total_droplets_estimate_outdict.append({
        'method': 'total_droplets_included',
        'umi_counts_cutoff': 0,
        'n_cells': total_droplets_included
    })

    # Make a plot of the different cutoffs
    _ = estimate_cutoffs_plot(
        '{}-total_drops_estimate_cutoffs-zoomed'.format(output_file),
        df_analysis,
        pd.DataFrame(total_droplets_estimate_outdict),
        # df_fit=df_fit,
        scale_x_log10=False,
        save_plot=True
    )
    _ = estimate_cutoffs_plot(
        '{}-total_drops_estimate_cutoffs'.format(output_file),
        df,
        pd.DataFrame(total_droplets_estimate_outdict),
        # df_fit=df_fit,
        scale_x_log10=True
    )

    # Plot the final estimates togher
    final_cutoffs = []
    final_cutoffs.append({
        'method': 'Estimated expected # cells',
        'umi_counts_cutoff': 0,
        'n_cells': expected_cells
    })
    final_cutoffs.append({
        'method': 'Estimated total # droplets',
        'umi_counts_cutoff': 0,
        'n_cells': total_droplets_included
    })
    df_plt = df.loc[df['barcode'] <= total_droplets_cutoff_knee + 100, :]
    _ = estimate_cutoffs_plot(
        '{}-final_estimates'.format(output_file),
        df_plt,
        pd.DataFrame(final_cutoffs),
        # df_fit=df_fit,
        scale_x_log10=False,
        save_plot=True
    )
    _ = estimate_cutoffs_plot(
        '{}-final_estimates-scale_x_log10'.format(output_file),
        df_plt,
        pd.DataFrame(final_cutoffs),
        # df_fit=df_fit,
        scale_x_log10=True,
        save_plot=True
    )


if __name__ == '__main__':
    main()
