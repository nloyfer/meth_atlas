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

ATLAS_FILE = 'Atlas.LUMP0.7.100bpblock.100CpGs+ECC500.csv.gz'
OUT_PATH = '.'
NR_CHRS_XTICKS = 30         # number of characters to be printed of the xticks
FIG_SIZE = (15, 7)          # figure size
COLOR_MAP = 'Vega20'        # color map. e.g gist_ncar, nipy_spectral, etc.
                            # See https://matplotlib.org/users/colormaps.html


###

def hide_small_tissues(df):
    others = df[df < 0.01].sum()
    df[df < 0.01] = 0
    df = df.append(others.rename('other'))
    return df


def gen_bars_colors_hatches(nr_tissues):
    matplotlib.rcParams['hatch.linewidth'] = 0.3
    hatches = [None, 'xxx', '...']

    nr_colors = int(math.ceil(nr_tissues / len(hatches)) + 1)

    # generate bars colors:
    cmap = matplotlib.cm.get_cmap(COLOR_MAP)
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=float(nr_colors))
    my_colors = [cmap(norm(k)) for k in range(nr_colors)]

    def get_bar_color(i):
        return my_colors[i % len(my_colors)], \
               hatches[int(i // math.ceil(nr_tissues / len(hatches)))]

    return [get_bar_color(x) for x in range(nr_tissues - 1)] + \
           [((.5, .5, .5, 1.), None)]


def plot_res(df, outpath):

    df = hide_small_tissues(df)
    nr_tissues = df.shape[0]

    # generate bars colors and hatches:
    colors_hatches = gen_bars_colors_hatches(nr_tissues)

    plt.figure(figsize=FIG_SIZE)
    r = [i for i in range(df.shape[1])]
    pre = np.zeros((1, df.shape[1]))
    for i in range(nr_tissues):
        adj_pre = [x for x in pre.flatten()]
        plt.bar(r,
                list(df.iloc[i, :]),
                edgecolor='white',
                width=0.85,
                label=df.index[i],
                bottom=adj_pre,
                color=colors_hatches[i][0],
                hatch=colors_hatches[i][1])
        pre += np.array(df.iloc[i, :])

    # Custom x axis
    plt.xticks(r, [w[:NR_CHRS_XTICKS] for w in df.columns], rotation='vertical')
    plt.xlabel("sample")
    plt.xlim(-.6, len(r) - .4)

    # Add a legend and a title
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
    plt.title('Deconvolution Results\n' + op.basename(outpath))

    # adjust layout, save and show
    plt.tight_layout(rect=[0, 0, .83, 1])
    plt.savefig(outpath + '_deconv_plot.png')
    plt.show()


def decon_single_samp(samp, atlas):
    """
    Deconvolve a single sample, using NNLS, to get the mixture coefficients.
    :param samp: a vector of a single sample
    :param atlas: the atlas DadtaFrame
    :return: the mixture coefficients (of size 25)
    """
    # remove missing sites from both sample and atlas:
    data = pd.concat([samp, atlas], axis=1).dropna(axis=0)
    if data.empty:
        print('Warning: skipping an empty sample ', file=sys.stderr)
        # print('Dropped {} missing sites'.format(self.atlas.shape[0] - red_atlas.shape[0]))
        return np.nan
    samp = data.iloc[:, 0]
    red_atlas = data.iloc[:, 1:]

    # get the mixture coefficients by deconvolution (non-negative least squares)
    mixture, residual = optimize.nnls(red_atlas, samp)
    mixture /= np.sum(mixture)
    return mixture


class Deconvolve:
    def __init__(self, atlas_path, samp_path, out_dir, slim=False, plot=False):
        self.out_dir = out_dir                      # Output dir to save mixture results and plot
        self.slim = slim                            # Write results table w\o indexes and header (bool)
        self.plot = plot                            # Plot results (bool)
        self.out_bname = self.get_bname(samp_path)

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
        return {'table': df.iloc[:, 1:], 'acc': df['acc'].to_frame(), 'tissues': list(df.columns)[1:]}

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
        samples = samples.merge(self.atlas['acc'], how='inner', on='acc')
        del samples['acc']
        return samples

    def run(self):

        # run deconvolution on all samples in parallel
        processes = []
        with Pool() as p:
            for i, smp_name in enumerate(list(self.samples)):
                params = (self.samples[smp_name], self.atlas['table'])
                processes.append(p.apply_async(decon_single_samp, params))
            p.close()
            p.join()

        # collect the results to 'res_table':
        arr = [pr.get() for pr in processes]
        res_table = np.empty((self.atlas['table'].shape[1], self.samples.shape[1]))
        for i in range(len(arr)):
            res_table[:, i] = arr[i]
        df = pd.DataFrame(res_table, columns=self.samples.columns, index=self.atlas['tissues'])

        # Dump results
        out_path = self.out_bname + '_deconv_output.csv'

        if self.slim:   # without indexes and header line
            df.to_csv(out_path, index=None, header=None, float_format='%.3f')
        else:
            df.to_csv(out_path, float_format='%.3f')

        # Plot pie charts
        if self.plot:
            plot_res(df, self.out_bname)


def main(args):
    try:
        Deconvolve(atlas_path=args.atlas_path,
                   samp_path=args.samples_path,
                   out_dir=args.out_dir,
                   slim=args.slim,
                   plot=args.plot).run()
    except Exception as e:
        print(e)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--atlas_path', '-a', default=ATLAS_FILE,
                        help='Path to Atlas csv file.\n'
                             'The first column must be Illumina IDs '
                           '(e.g cg00000029)')
    parser.add_argument('samples_path',
                        help='Path to samples csv file. It must have a header line.\n'
                             'The first column must be Illumina IDs (e.g cg00000029)')
    parser.add_argument('--slim', action='store_true',
                        help='Write the results table *without indexes and header line*')
    parser.add_argument('--plot', action='store_true',
                        help='Plot pie charts of the results')
    parser.add_argument('--out_dir', '-o', default=OUT_PATH, help='Ouput directory')

    args = parser.parse_args()
    main(args)
