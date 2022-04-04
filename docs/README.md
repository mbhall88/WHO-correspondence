Description of the files in this directory.

## `gentb-samplesheet.csv`

This was taken
from [Additional File 1](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00953-4#MOESM1)
from
the [GenTB paper](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00953-4)
.

It aggregates accessions and phenotypes from across multiple datasets.

## `who-samplesheet.csv`

This was extracted from sheet S1 of Supplementary Appendix 2 in
the [WHO catalogue paper](https://doi.org/10.1016/S2666-5247(21)00301-3).

Specifically, the following columns were extracted:

- `ena_project`
- `ena_sample`
- `ena_experiment`
- `ena_run`

## `samplesheet.csv`

This is the "master" sheet for this project. See
the [data cleaning notebook](../workflow/notebooks/data-cleaning.ipynb) for the full
details of how this file was generated.

Essentially, all samples used in the WHO paper were removed from the non-WHO datasets in
the samplesheets in this directory and the remaining data were collated into this file.