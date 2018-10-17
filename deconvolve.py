#!/usr/bin/python3 -u
import numpy as np
import pandas as pd
from scipy import optimize
import argparse
import os.path as op
import sys
from multiprocessing import Pool


ATLAS_FILE = 'atlas.csv'
ATLAS_FILE = 'Atlas.LUMP0.7.100bpblock.100CpGs+ECC500.csv.gz'
OUT_PATH = '.'
NR_CHRS_XTICKS = 30         # number of characters to be printed of the xticks
FIG_SIZE = (15, 7)          # figure size
COLOR_MAP = 'gist_ncar'     # color map. e.g gist_ncar, nipy_spectral, etc.
                            # See https://matplotlib.org/users/colormaps.html


def decon_single_samp(samp, samp_name, atlas):
    """
    Deconvolve a single sample, using NNLS, to get the mixture coefficients.
    :param samp: a vector of a single sample
    :param atlas: the atlas DadtaFrame
    :return: the mixture coefficients (of size 25)
    """
    # remove missing sites from both sample and atlas:
    data = pd.concat([samp, atlas], axis=1).dropna(axis=0)
    if data.empty:
        print('Warning: sample {} is empty'.format(samp_name), file=sys.stderr)
        # print('Dropped {} missing sites'.format(self.atlas.shape[0] - red_atlas.shape[0]))
        return np.nan
    samp = data.iloc[:, 0]
    red_atlas = data.iloc[:, 1:]

    # get the mixture coefficients by deconvolution (non-negative least squares)
    mixture, residual = optimize.nnls(red_atlas, samp)
    mixture /= np.sum(mixture)
    return (mixture, samp_name)


class Deconvolve:
    def __init__(self, atlas_path, samp_path, out_dir, slim=False, plot=False):
        self.out_dir = out_dir      # Output dir to save mixture results and plot
        self.slim = slim            # Write results table without indexes and header line (bool)
        self.plot = plot            # Plot results (bool)
        self.samp_path = samp_path
        self.out_bname = self.get_bname()

        # Load input files:
        self.atlas = self.load_atlas(atlas_path)        # Atlas
        self.samples = self.load_sample()               # Samples to deconvolve

    def get_bname(self):
        """
        Compose output files path:
        join the out_dir path with the basename of the samples file
        remove csv and gz extensions.
        """
        base_fname = op.basename(self.samp_path)

        if base_fname.endswith('.gz'):
            base_fname = op.splitext(base_fname)[0]
        base_fname = op.splitext(base_fname)[0]
        out_path = op.join(self.out_dir, base_fname)
        return out_path

    def load_atlas(self, atlas_path):
        """
        Read the atlas csv file, save data in self.atlas
        :param atlas_path: Path to the atlas csv file
        """
        # validate path:
        self._validate_csv_file(atlas_path)

        # Read atlas, sort it and drop duplicates
        # print('atlas_path', atlas_path)
        df = pd.read_csv(atlas_path)
        df.rename(columns={list(df)[0]: 'acc'}, inplace=True)
        df = df.sort_values(by='acc').drop_duplicates(subset='acc').reset_index(drop=True)
        return {'table': df.iloc[:, 1:], 'acc': df['acc'].to_frame(), 'tissues': list(df.columns)[1:]}

    def _validate_csv_file(self, csv_path):
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

    def load_sample(self):
        """
        Read samples csv file. Reduce it to the atlas sites, and save data in self.samples
        Note: samples file must contain a header line.
        """

        # validate path:
        self._validate_csv_file(self.samp_path)

        samples = pd.read_csv(self.samp_path)
        samples.rename(columns={list(samples)[0]: 'acc'}, inplace=True)
        samples = samples.sort_values(by='acc').drop_duplicates(subset='acc').reset_index(drop=True)
        samples = samples.merge(self.atlas['acc'], how='inner', on='acc')
        del samples['acc']
        return samples

    def plot_res2(self, df):
        import matplotlib.pylab as plt
        import matplotlib.cm

        import seaborn as sns

        nr_tissues = len(self.atlas['tissues'])
        # plt.figure(figsize=FIG_SIZE)
        sns.set_style('ticks')
        # sns.set()
        current_palette_4 = sns.color_palette("hls", 25)
        sns.set_palette(current_palette_4)
        df.T.plot(kind='bar', stacked=True)

        # barWidth = 0.85
        r = [i for i in range(self.samples.shape[1])]
        # pre = np.zeros((1, df.shape[1]))
        # cmap = matplotlib.cm.get_cmap(COLOR_MAP)
        # norm = matplotlib.colors.Normalize(vmin=0.0, vmax=float(nr_tissues))
        # my_colors = [cmap(norm(k)) for k in range(nr_tissues)]
        #
        # # print(my_colors)
        #
        # for i in range(nr_tissues):
        #     adj_pre = [x for x in pre.flatten()]
        #     plt.bar(r,
        #             list(df.iloc[i, :]),
        #             edgecolor='white',
        #             width=barWidth,
        #             label=self.atlas['tissues'][i],
        #             bottom=adj_pre,
        #             color=my_colors[i])
        #     pre += np.array(df.iloc[i, :])

        # Custom x axis
        plt.xticks(r, [w[:NR_CHRS_XTICKS] for w in self.samples.columns], rotation='vertical')
        plt.xlabel("sample")
        # Add a legend
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
        plt.title('Deconvolution Results\n' + self.out_bname)
        plt.tight_layout(rect=[0, 0, .83, 1])
        # plt.savefig(self.out_bname + '_deconv_plot.png')
        plt.show()

    def _remove_small_tissues(self, df):
        df['maxVal'] = np.max(df, axis=1)
        good = df.loc[df['maxVal'] >= 0.01]
        others = df.loc[df['maxVal'] < 0.01]
        others = others.sum(axis=0)
        good = good.append(others.rename('other'))
        del good['maxVal']
        return good

    def plot_res(self, df):
        import matplotlib.pylab as plt
        import matplotlib.cm
        import math
        import matplotlib as mpl

        # df = self._remove_small_tissues(df)

        tissues = self.atlas['tissues']
        # tissues = list(df.index)
        mpl.rcParams['hatch.linewidth'] = 0.3

        hatches = [None, 'xxx', '...']

        nr_tissues = len(tissues)
        nr_colors = int(math.ceil(nr_tissues / len(hatches)) + 1)

        plt.figure(figsize=FIG_SIZE)

        r = [i for i in range(self.samples.shape[1])]
        pre = np.zeros((1, df.shape[1]))

        # generate colors:
        cmap = matplotlib.cm.get_cmap(COLOR_MAP)
        norm = matplotlib.colors.Normalize(vmin=0.0, vmax=float(nr_colors))
        my_colors = [cmap(norm(k)) for k in range(1, nr_colors)]

        for i in range(nr_tissues):
            adj_pre = [x for x in pre.flatten()]
            plt.bar(r,
                    list(df.iloc[i, :]),
                    edgecolor='white',
                    width=0.85,
                    label=tissues[i],
                    bottom=adj_pre,
                    color=my_colors[i % len(my_colors)],
                    # hatch=hatches[i % len(hatches)])
                    hatch=hatches[int(i // math.ceil(nr_tissues / len(hatches)))])
            pre += np.array(df.iloc[i, :])

        # Custom x axis
        plt.xticks(r, [w[:NR_CHRS_XTICKS] for w in self.samples.columns], rotation='vertical')
        plt.xlabel("sample")
        # Add a legend
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
        plt.title('Deconvolution Results\n' + op.basename(self.out_bname))
        plt.tight_layout(rect=[0, 0, .83, 1])
        plt.savefig(self.out_bname + '_deconv_plot.png')
        plt.show()


    def run(self):

        # run deconvolution on all samples in parallel
        processes = []
        with Pool() as p:
            for i, samp_name in enumerate(list(self.samples)):
                processes.append(p.apply_async(decon_single_samp,
                                               (self.samples[samp_name],
                                                samp_name,
                                                self.atlas['table'])))
            p.close()
            p.join()

        # collect the results to 'res_table':
        arr = [pr.get() for pr in processes]
        res_table = np.empty((self.atlas['table'].shape[1], self.samples.shape[1]))
        for i in range(len(arr)):
            res_table[:, i] = arr[i][0]
            # print(arr[i][1], end=', ')
        print('')
        df = pd.DataFrame(res_table, columns=self.samples.columns, index=self.atlas['tissues'])

        # Dump results
        out_path = self.out_bname + '_deconv_output.csv'

        if self.slim:   # without indexes and header line
            df.to_csv(out_path, index=None, header=None, float_format='%.3f')
        else:
            df.to_csv(out_path, float_format='%.3f')

        # Plot pie charts
        if self.plot:
            self.plot_res(df)


def main(args):
    try:
        Deconvolve(atlas_path=args.atlas_path, samp_path=args.samples_path, out_dir=args.out_dir, slim=args.slim,
                   plot=args.plot).run()
    except Exception as e:
        print(e)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--atlas_path', '-a', default=ATLAS_FILE,
                        help='Path to Atlas csv file.\n'
                             'The first column must be Illumina IDs (e.g cg00000029)')
    parser.add_argument('samples_path', help='Path to samples csv file. It must have a header line.\n'
                                             'The first column must be Illumina IDs (e.g cg00000029)')
    parser.add_argument('--slim', action='store_true', help='Write the results table *without indexes and header line*')
    parser.add_argument('--plot', action='store_true', help='Plot pie charts of the results')
    parser.add_argument('--out_dir', '-o', default=OUT_PATH, help='Ouput directory')

    args = parser.parse_args()
    main(args)
