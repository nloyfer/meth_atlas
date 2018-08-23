#!/usr/bin/python3 -u
import numpy as np
import pandas as pd
from scipy import optimize
import argparse
import os.path as op
import math

ATLAS_FILE = 'atlas.csv'
# TEST_FILE = '/cs/cbio/tommy/Meth/DATA/Atlas/Whole_blood_EPIC.csv'
OUT_PATH = 'out.csv'

pd.set_option('display.width', 160)
np.set_printoptions(linewidth=160)

# todo: deal with header-less input files (self.load_sample)

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
        if not (op.isfile(atlas_path) and atlas_path.endswith('.csv')):
            raise IOError('Invalid atlas path:\n%s' % atlas_path)

        # Read atlas, sort it and drop duplicates
        # print('atlas_path', atlas_path)
        df = pd.read_csv(atlas_path)
        df.rename(columns={list(df)[0]: 'acc'}, inplace=True)
        df = df.sort_values(by='acc').drop_duplicates(subset='acc').reset_index(drop=True)
        return {'table': df.iloc[:, 1:], 'acc': df['acc'].to_frame(), 'tissues': list(df.columns)[1:]}

    def load_sample(self, samp_path):
        """
        Read samples csv file. Reduce it to the atlas sites, and save data in self.samples
        :param samp_path: Path to samples csv file. Each column is a sample. Note: file must contain a header line.
        """

        # validate path:
        if not (op.isfile(samp_path) and samp_path.endswith('csv')):
            raise IOError('Invalid samples path:\n%s' % samp_path)

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

    def _plot_dims(self):
        N = self.samples.shape[1]
        nrows = int(math.sqrt(N))
        ncols = N // nrows
        if nrows * ncols < N:
            ncols += 1
        return nrows, ncols

    def plot_res(self, df):
        import matplotlib.pylab as plt

        fig, axes = plt.subplots(*self._plot_dims())
        for i, ax in enumerate(axes.flatten()):
            if i >= df.shape[1]:
                break
            ax.pie(df.iloc[:, i], autopct="%.1f%%")
            ax.set_title(df.columns.values[i])

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


def main(args):
    try:
        Deconvolve(atlas_path=args.atlas_path, samp_path=args.samples_path, out_path=args.out_path, slim=args.slim,
                   plot=args.plot).run()
    except Exception as e:
        print(e)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--atlas_path', '-a', default=ATLAS_FILE,
                        help='Path to Atlas csv file')
    parser.add_argument('samples_path', help='Path to samples csv file. It must have a header line.')
    parser.add_argument('--slim', action='store_true', help='Write the results table *without indexes and header line*')
    parser.add_argument('--plot', action='store_true', help='Plot pie charts of the results')
    parser.add_argument('--out_path', '-o', default=OUT_PATH, help='Ouput path')
    args = parser.parse_args()
    main(args)
