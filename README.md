Methylation Atlas Deconvolution




usage: deconvolve.py [-h] [--atlas_path ATLAS_PATH] [--slim] [--plot]
                     [--out_path OUT_PATH]
                     samples_path

positional arguments:
  samples_path          Path to samples csv file. It must have a header line.

optional arguments:
  -h, --help            show this help message and exit
  --atlas_path ATLAS_PATH, -a ATLAS_PATH
                        Path to Atlas csv file
  --slim                Write the results table *without indexes and header
                        line*
  --plot                Plot pie charts of the results
  --out_path OUT_PATH, -o OUT_PATH
                        Ouput path
