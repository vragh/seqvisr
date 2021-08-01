seqvisr
================

[![Build
Status](https://travis-ci.com/vragh/seqvisr.svg?branch=main)](https://travis-ci.com/vragh/seqvisr)

# seqvisr

Biological sequence visualization functions in `R`.

Consists of two functions: `msavisr()` and `pdomvisr()`.

`msavisr()` takes `FASTA`-formatted multiple sequence alignments (MSAs)
as inputs and produces visualizations like the example below:
<img src="./docs/articles/seqvisr_files/figure-html/msavisr_examples-1.png" alt="msavisr" width="60%" height="60%"/>

`pdomvisr()` takes appropriately formatted tabular data and produces
domain structure diagrams like the one below:

<img src="./docs/articles/seqvisr_files/figure-html/label_small-1.png" alt="pdomvisr" width="60%" height="60%"/>

## Installation

Ensure that the package [`devtools`](https://github.com/r-lib/devtools)
is installed.

Then execute
`devtools::install_github("vragh/seqvisr", build_manual = TRUE, build_vignettes = TRUE)`
from within `R` to install the package. If the manual and package
vignette are not necessary,
`build_manual = TRUE, build_vignettes = TRUE` can be omitted.

## Usage

See `?seqvisr`, `?msavisr`, and `?pdomvisr` (from within `R`) to access
documentation and examples. Run `browseVignettes(seqvisr)` to access the
vignette.

## License

Released under GPL-3. A copy of the license file can be found in the
file `LICENSE` (`LICENCE` is just for `R` purposes).
