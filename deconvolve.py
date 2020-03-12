#!/usr/bin/python3 -u
import numpy as np
import pandas as pd
from scipy import optimize
import argparse
import os.path as op
import sys
from multiprocessing import Pool
import math
import matplotlib.pylab as plt
import matplotlib.cm
import matplotlib.colors

ATLAS_FILE = './reference_atlas.csv'
OUT_PATH = '.'

# Plotting parameters:
NR_CHRS_XTICKS = 30         # number of characters to be printed of the xticks
FIG_SIZE = (15, 7)          # figure size
COLOR_MAP = 'tab10'         # color map. See https://matplotlib.org/users/colormaps.html
#COLOR_MAP = 'Vega10'
# tissues with less than OTHERS_THRESH contribution will be clustered to 'other' (black):
OTHERS_THRESH = 0.01


####################################
#       Plotting methods           #
####################################

def hide_small_tissues(df):
    """
    tissues with very small contribution are grouped to the 'other' category.
    :return: The DataFrame with the new category ('other'),
             where the low-contribution tissues are set to 0.
    """
    others = df[df < OTHERS_THRESH].sum()
    df[df < OTHERS_THRESH] = 0.0
    df = df.append(others.rename('other'))
    return df


def gen_bars_colors_hatches(nr_tissues):
    """
    Generate combinations of colors and hatches for the tissues bars
    Every tissue will get a tuple of (color, hatch)
    the last tuple is for the 'other' category, and is always black with no hatch.
    :return: a list of tuples, with length == nr_tissues
    """
    matplotlib.rcParams['hatch.linewidth'] = 0.3
    hatches = [None, 'xxx', '...', 'O', '++'][:nr_tissues // 7]

    nr_colors = int(math.ceil(nr_tissues / len(hatches)) + 1)

    # generate bars colors:
    cmap = matplotlib.cm.get_cmap(COLOR_MAP)
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=float(nr_colors))
    colors = [cmap(norm(k)) for k in range(nr_colors)]

    def get_i_bar_tuple(i):
        color_ind = i % nr_colors
        hatch_ind = int(i // math.ceil(nr_tissues / len(hatches)))
        return colors[color_ind], hatches[hatch_ind]

    colors_hatches_list = [get_i_bar_tuple(i) for i in range(nr_tissues - 1)]
    return colors_hatches_list + [((0, 0, 0, 1), None)]


def plot_res(df, outpath, show=False):

    df = hide_small_tissues(df)
    nr_tissues, nr_samples = df.shape

    # generate bars colors and hatches:
    colors_hatches = gen_bars_colors_hatches(nr_tissues)

    plt.figure(figsize=FIG_SIZE)
    r = [i for i in range(nr_samples)]
    bottom = np.zeros(nr_samples)
    for i in range(nr_tissues):
        plt.bar(r, list(df.iloc[i, :]),
                edgecolor='white',
                width=0.85,
                label=df.index[i],
                bottom=bottom,
                color=colors_hatches[i][0],
                hatch=colors_hatches[i][1])
        bottom += np.array(df.iloc[i, :])

    # Custom x axis
    plt.xticks(r, [w[:NR_CHRS_XTICKS] for w in df.columns], rotation='vertical', fontsize=9)
    plt.xlabel("sample")
    plt.xlim(-.6, nr_samples - .4)

    # Add a legend and a title
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
    plt.title('Deconvolution Results\n' + op.basename(outpath))

    # adjust layout, save and show
    plt.tight_layout(rect=[0, 0, .83, 1])
    plt.savefig(outpath + '_deconv_plot.png')
    if show:
        plt.show()


####################################
#     Deconvolve class             #
####################################


class Deconvolve:
    def __init__(self, atlas_path, samp_path, out_dir, resid, slim=False, plot=False):
        self.out_dir = out_dir                      # Output dir to save mixture results and plot
        self.slim = slim                            # Write results table w\o indexes and header (bool)
        self.plot = plot                            # Plot results (bool)
        self.resid = resid                          # Output residuals as well
        self.out_bname = self.get_bname(samp_path)  # output files path w/o extension

        # Load input files:
        self.atlas = self.load_atlas(atlas_path)    # Atlas
        self.samples = self.load_sample(samp_path)  # Samples to deconvolve

    def get_bname(self, samp_path):
        """
        Compose output files path:
        join the out_dir path with the basename of the samples file
        remove csv and gz extensions.
        """
        base_fname = op.basename(samp_path)

        if base_fname.endswith('.gz'):
            base_fname = op.splitext(base_fname)[0]
        base_fname = op.splitext(base_fname)[0]
        return op.join(self.out_dir, base_fname)

    @staticmethod
    def load_atlas(atlas_path):
        """
        Read the atlas csv file, save data in self.atlas
        :param atlas_path: Path to the atlas csv file
        """
        # validate path:
        Deconvolve._validate_csv_file(atlas_path)

        # Read atlas, sort it and drop duplicates
        # print('atlas_path', atlas_path)
        df = pd.read_csv(atlas_path)
        df.rename(columns={list(df)[0]: 'acc'}, inplace=True)
        df = df.sort_values(by='acc').drop_duplicates(subset='acc').reset_index(drop=True)
        return df

    @staticmethod
    def _validate_csv_file(csv_path):
        """
        Validate an input csv file. Raise an exception or print warning if necessary.
        :param csv_path: input csv path
        """
        err_msg = ''

        # check if file exists and ends with 'csv':
        if not op.isfile(csv_path):
            err_msg = 'no such file:\n%s' % csv_path
        elif not (csv_path.endswith('csv') or csv_path.endswith('csv.gz')):
            err_msg = 'file must end with ".csv[.gz]":\n%s' % csv_path

        # take a peek and validate the file format
        else:
            input_head = pd.read_csv(csv_path, nrows=4)

            # at least two columns:
            if input_head.shape[1] < 2:
                err_msg = 'file must contain at least 2 columns (accessions and a values). '

            # first column must be Illumina IDs column
            elif not str(input_head.iloc[0, 0]).startswith('cg'): 
                    err_msg = 'invalid Illumina ID column'

            # print a warning if the second column in the csv file has a numeric header
            # (this probably means there is no header)
            if input_head.columns[1].replace('.', '', 1).isdigit():
                print('Warning: input files should have headers', file=sys.stderr)

        if err_msg:
            err_msg = op.basename(csv_path) + ': ' + err_msg
            raise ValueError(err_msg)

    @staticmethod
    def decon_single_samp(samp, atlas):
        """
        Deconvolve a single sample, using NNLS, to get the mixture coefficients.
        :param samp: a vector of a single sample
        :param atlas: the atlas DadtaFrame
        :return: the mixture coefficients (of size 25)
        """

        name = samp.columns[1]

        # remove missing sites from both sample and atlas:
        data = samp.merge(atlas, on='acc', how='inner').copy().dropna(axis=0)
        if data.empty:
            print('Warning: skipping an empty sample {}'.format(name), file=sys.stderr)
            # print('Dropped {} missing sites'.format(self.atlas.shape[0] - red_atlas.shape[0]))
            return np.nan
        print('{}: {} sites'.format(name, data.shape[0]), file=sys.stderr)
        del data['acc']

        samp = data.iloc[:, 0]
        red_atlas = data.iloc[:, 1:]

        # get the mixture coefficients by deconvolution (non-negative least squares)
        mixture, residual = optimize.nnls(red_atlas, samp)
        mixture /= np.sum(mixture)
        return mixture, residual

    def load_sample(self, samp_path):
        """
        Read samples csv file. Reduce it to the atlas sites, and save data in self.samples
        Note: samples file must contain a header line.
        """

        # validate path:
        Deconvolve._validate_csv_file(samp_path)

        samples = pd.read_csv(samp_path)
        samples.rename(columns={list(samples)[0]: 'acc'}, inplace=True)
        samples = samples.sort_values(by='acc').drop_duplicates(subset='acc').reset_index(drop=True)
        samples = samples.merge(self.atlas['acc'].to_frame(), how='inner', on='acc')
        return samples

    def run(self):

        # run deconvolution on all samples in parallel
        processes = []
        with Pool() as p:
            for i, smp_name in enumerate(list(self.samples)[1:]):
                params = (self.samples[['acc', smp_name]], self.atlas)
                processes.append(p.apply_async(Deconvolve.decon_single_samp, params))
            p.close()
            p.join()

        self.samples = self.samples.iloc[:, 1:]

        # collect the results to 'res_table':
        arr = [pr.get() for pr in processes]
        res_table = np.empty((self.atlas.shape[1] - 1, self.samples.shape[1]))
        resids_table = np.empty((self.samples.shape[1], 1))
        for i in range(len(arr)):
            res_table[:, i], resids_table[i] = arr[i]
        df = pd.DataFrame(res_table, columns=self.samples.columns, index=list(self.atlas.columns)[1:])

        # Dump results
        out_path = self.out_bname + '_deconv_output.csv'

        if self.slim:   # without indexes and header line
            df.to_csv(out_path, index=None, header=None, float_format='%.3f')
        else:
            df.to_csv(out_path, float_format='%.3f')

        if self.resid:
            rf = pd.DataFrame(resids_table, columns=['Residuals'], index=self.samples.columns)
            rf.to_csv(self.out_bname + '_residuals.csv', float_format='%.3f')

        # Plot pie charts
        plot_res(df, self.out_bname, self.plot)


####################################
#            main                  #
####################################


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--atlas_path', '-a', default=ATLAS_FILE,
                        help='Path to Atlas csv file.\nThe first column must be'
                             ' Illumina IDs (e.g cg00000029)')

    parser.add_argument('samples_path',
                        help='Path to samples csv file. It must have a header line.\n'
                             'The first column must be Illumina IDs (e.g cg00000029)')

    parser.add_argument('--slim', action='store_true',
                        help='Write the results table *without indexes and header line*')

    parser.add_argument('--residuals', '-r', action='store_true',
                        help='Output residuals to a separate file')

    parser.add_argument('--plot', action='store_true',
                        help='Plot stacked bars of the results')

    parser.add_argument('--out_dir', '-o', default=OUT_PATH, help='Output directory')

    args = parser.parse_args()

    Deconvolve(atlas_path=args.atlas_path,
               samp_path=args.samples_path,
               out_dir=args.out_dir,
               resid=args.residuals,
               slim=args.slim,
               plot=args.plot).run()


if __name__ == '__main__':
    main()
