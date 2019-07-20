# Methylation Atlas Deconvolution

This is a standalone program for [deconvolution of array methylome](https://www.nature.com/articles/s41467-018-07466-6)[1].
It uses an input **reference atlas** file to deconvolve a given **sample**, or multiple samples.
Outputs a csv file, and plots a stacked bars figure.

### atlas
A reference atlas file. 
- csv file
- Contains a header (columns names).
- The first column must be Illumina IDs.

The reference atlas used on the paper is supplied in this repository - *reference_atlas.csv*.
The full reference atlas, with ~390K sites (before the feature selection process) is also supplied - *full_reference_atlas.csv.gz*. For better deconvolurions results, it's recommended not to use the full atlas, but the smaller one.

### samples
A file containing one or more samples, with similar requirements as the atlas file (csv, header, index column).
The CpG (Illumina ID) column may contain different CpG sites than the ones in the atlas files, as long as they share some sites.

An example dummy file is supplied, *examples.csv*.

---

### Usage

```
usage: deconvolve.py [-h] [--atlas_path ATLAS_PATH] [--slim] [--plot]
                     [--out_dir OUT_DIR]
                     samples_path

positional arguments:
  samples_path          Path to samples csv file. It must have a header line.
                        The first column must be Illumina IDs (e.g cg00000029)

optional arguments:
  -h, --help            show this help message and exit
  --atlas_path ATLAS_PATH, -a ATLAS_PATH
                        Path to Atlas csv file. The first column must be
                        Illumina IDs (e.g cg00000029)
  --slim                Write the results table *without indexes and header
                        line*
  --plot                Display the stacked bars figure
  --out_dir OUT_DIR, -o OUT_DIR
                        Output directory
```

---
### Example
```
deconvolve.py -a reference_atlas.csv examples.csv
```
will deconvolve all samples given as columns in *examples.csv*, dump the resulting coefficients to a csv file named *examples_deconv_output.csv*, plot them, and dump the figure to *examples_deconv_plot.png*.
![Image of bar plot](https://github.com/nloyfer/meth_atlas/blob/master/examples_deconv_plot.png)


### Pre-processing
We recommend pre-processing the IDAT files using the script and reference sample listed [here](https://github.com/nloyfer/meth_atlas/tree/master/pre_process#pre-processing).

------------------
[1]: Moss, J.; Magenheim, J.; Neiman, D.; Zemmour, H.; Loyfer, N.; Korach, A.; Samet, Y.; Maoz, M.; Druid, H.; Arner, P.; et al. Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease. Nat. Commun. 2018, 9, 5068.
