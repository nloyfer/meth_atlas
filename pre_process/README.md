# Pre-Processing
For optimal results, we recommend using `process_array.R` to generate the samples csv file from a directory containing `IDAT` files.
This script normalises the data using an arbitrary reference sample (`ref_sample.RData`), and filters by p-value, sex chromosomes and bead number. It also drops SNPs and non-CpG sites.

#### Requirenments:
- R [minfi package](https://bioconductor.org/packages/release/bioc/html/minfi.html).

To install minfi:

```R
# try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("minfi")
```

### Usage:
```
Rscript process_array.R {idat directory} {output csv path} {path to ref_sample.RData}
```

