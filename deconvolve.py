#!/usr/bin/python3 -u
import numpy as np
import pandas as pd
from scipy import optimize
import argparse
import os.path as op
import sys

ATLAS_FILE = 'atlas.csv'
OUT_PATH = 'out.csv'
NR_CHRS_XTICKS = 30         # number of characters to be printed of the xticks
FIG_SIZE = (15, 7)         # figure size
COLOR_MAP = 'gist_ncar'     # color map. e.g gist_ncar, nipy_spectral, etc.
                            # See https://matplotlib.org/users/colormaps.html



def decon_single_samp2(samp, atlas):
    """
    Deconvolve a single sample, using NNLS, to get the mixture coefficients.
    :param samp: a vector of a single sample
    :return: the mixture coefficients (of size 25)
    """
    # remove missing sites from both sample and atlas:
    data = pd.concat([samp, atlas], axis=1).dropna(axis=0)
    samp = data.iloc[:, 0]
    red_atlas = data.iloc[:, 1:]
    # print('Dropped {} missing sites'.format(self.atlas.shape[0] - red_atlas.shape[0]))

    # get the mixture coefficients by deconvolution (non-negative least squares)
    mixture, residual = optimize.nnls(red_atlas, samp)
    mixture /= np.sum(mixture)
    return mixture.round(3)

class Deconvolve:
    def __init__(self, atlas_path, samp_path, out_path=None, slim=False, plot=False):
        self.out_path = out_path    # Output path to save mixture results
        self.slim = slim            # Write results table without indexes and header line
        self.plot = plot            # Plot results

        # Load input files:
        self.atlas = self.load_atlas(atlas_path)        # Atlas
        self.samples = self.load_sample(samp_path)      # Samples to deconvolve

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
        elif not csv_path.endswith('csv'):
            err_msg = 'file must end with ".csv":\n%s' % csv_path

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
        :param samp_path: Path to samples csv file. Each column is a sample. Note: file must contain a header line.
        """

        # validate path:
        self._validate_csv_file(samp_path)

        samples = pd.read_csv(samp_path)
        samples.rename(columns={list(samples)[0]: 'acc'}, inplace=True)
        samples = samples.sort_values(by='acc').drop_duplicates(subset='acc').reset_index(drop=True)
        samples = samples.merge(self.atlas['acc'], how='inner', on='acc')
        del samples['acc']
        return samples

    def decon_single_samp(self, samp):
        """
        Deconvolve a single sample, using NNLS, to get the mixture coefficients.
        :param samp: a vector of a single sample
        :return: the mixture coefficients (of size 25)
        """
        # remove missing sites from both sample and atlas:
        data = pd.concat([samp, self.atlas['table']], axis=1).dropna(axis=0)
        samp = data.iloc[:, 0]
        red_atlas = data.iloc[:, 1:]
        # print('Dropped {} missing sites'.format(self.atlas.shape[0] - red_atlas.shape[0]))

        # get the mixture coefficients by deconvolution (non-negative least squares)
        mixture, residual = optimize.nnls(red_atlas, samp)
        mixture /= np.sum(mixture)
        return mixture.round(3)

    def plot_res(self, df):
        import matplotlib.pylab as plt
        plt.figure(figsize=FIG_SIZE)

        barWidth = 0.85
        r = [i for i in range(self.samples.shape[1])]
        pre = np.zeros((1, df.shape[1]))
        my_colors = np.random.rand(25, 3)
        my_colors = [(my_colors[x, 0], my_colors[x, 1], my_colors[x, 2]) for x in range(25)]
        import matplotlib
        cmap = matplotlib.cm.get_cmap(COLOR_MAP)
        norm = matplotlib.colors.Normalize(vmin=0.0, vmax=25.0)
        my_colors = [cmap(norm(k)) for k in range(25)]

        print(my_colors)
        for i in range(len(self.atlas['tissues'])):
            adj_pre = [x for x in pre.flatten()]
            plt.bar(r,
                    list(df.iloc[i, :]),
                    edgecolor='white',
                    width=barWidth,
                    label=self.atlas['tissues'][i],
                    bottom=adj_pre,
                    color=my_colors[i])
            pre += np.array(df.iloc[i, :])

        # Custom x axis
        plt.xticks(r,
                   [w[:NR_CHRS_XTICKS] for w in self.samples.columns],
                   rotation='vertical')
        # plt.subplots_adjust(bottom=.25)
        plt.xlabel("sample")
        # Add a legend
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
        plt.title('Deconvolution Results')
        plt.tight_layout(rect=[0, 0, .83, 1])
        plt.show()

    def run(self):

        res_table = np.empty((self.atlas['table'].shape[1], self.samples.shape[1]))
        for i, samp_name in enumerate(list(self.samples)):
            res_table[:, i] = self.decon_single_samp(self.samples[samp_name])
        # print(res_table)

        df = pd.DataFrame(res_table, columns=self.samples.columns, index=self.atlas['tissues'])

        # Dump results
        if self.slim:   # without indexes and header line
            df.to_csv(self.out_path, index=None, header=None)
        else:
            df.to_csv(self.out_path)

        # Plot pie charts
        if self.plot:
            self.plot_res(df)


    def run_p(self):
        from multiprocessing import Pool

        res_table = np.empty((self.atlas['table'].shape[1], self.samples.shape[1]))
        # print(res_table)

        processes = []
        with Pool() as p:
            for i, samp_name in enumerate(list(self.samples)):
                processes.append(p.apply_async(decon_single_samp2, (self.samples[samp_name], self.atlas['table'])))
            p.close()
            p.join()

        arr = [pr.get() for pr in processes]
        for i in range(len(arr)):
            res_table[:, i] = arr[i]

        df = pd.DataFrame(res_table, columns=self.samples.columns, index=self.atlas['tissues'])

        # Dump results
        if self.slim:   # without indexes and header line
            df.to_csv(self.out_path, index=None, header=None)
        else:
            df.to_csv(self.out_path)

        # Plot pie charts
        if self.plot:
            self.plot_res(df)


def main(args):
    try:
        Deconvolve(atlas_path=args.atlas_path, samp_path=args.samples_path, out_path=args.out_path, slim=args.slim,
                   plot=args.plot).run_p()
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
    parser.add_argument('--out_path', '-o', default=OUT_PATH, help='Ouput path')
    args = parser.parse_args()
    main(args)
