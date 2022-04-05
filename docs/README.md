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

## `dlm-mics.csv`

Gathering Delamanid (DLM) phenotypes from https://journals.asm.org/doi/full/10.1128/JCM.01304-20

They used a critical concentraction of >0.06μg/ml - based on the [WHO technical guidelines](https://apps.who.int/iris/bitstream/handle/10665/260470/WHO-CDS-TB-2018.5-eng.pdf).

MICs and accessions where taken from the supplementary Excel spreadsheet https://journals.asm.org/doi/suppl/10.1128/JCM.01304-20/suppl_file/jcm.01304-20-sd004.xlsx

This file was taken from the following columns in the "Mutations list for DLM" sheet of the above-mentioned Excel file.
- ISOLATE 1 (renamed to ISOLATE 1 DLM MIC (ug/ml) in this file)
- ISOLATE 2 (renamed to ISOLATE 2 DLM MIC (ug/ml) in this file)
- DRS Sample selected (Isolate 1)
- DRS Sample selected (Isolate 2)

The purpose of this file was to match up the isolate names and phenotypes with the accessions in `dlm-accessions.csv`.

## `dlm-accessions.csv`

Gathering Delamanid (DLM) phenotypes from https://journals.asm.org/doi/full/10.1128/JCM.01304-20

They used a critical concentraction of >0.06μg/ml - based on the [WHO technical guidelines](https://apps.who.int/iris/bitstream/handle/10665/260470/WHO-CDS-TB-2018.5-eng.pdf).

MICs and accessions where taken from the supplementary Excel spreadsheet https://journals.asm.org/doi/suppl/10.1128/JCM.01304-20/suppl_file/jcm.01304-20-sd004.xlsx

This data was extracted from the "WHO database column" of that Excel file.

Phenotypic DST of first and second line dugs. In brackets are indicated the used critical concentration to defined Suceptibility (S) or Resistance (R) to each drug. N/A means "not available" (reference: Zignol M et al., 2018. Lancet Infect Dis 18:675-683. doi:10.1016/S1473-3099(18)30073-2).

## `huang-samplesheet.csv`

This CSV file contains the phenotypes for the isolates sequenced in [Huang *et al.*](https://doi.org/10.1093/cid/ciy883).

In particular, this dataset provides linezolid phenotypes. 

This file was extracted from [Supplementary Table S6](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/cid/69/3/10.1093_cid_ciy883/1/ciy883_suppl_supplementary_table_s6.docx?Expires=1652148910&Signature=D7iy~iS4ZXFMCAifGgjnpD-OslN6pINjGYhFbqx1RI2unpvW9gaZ2CkiwXLd3cBagIABut8U4qKOXY11mOVw9LMZohqZNtkibKuu7SFgBJ-c2vMz9h10GNKHj5Ya98dg6AT7IPVsTk77OHopWoFsE6JlgbeCPIlH4i-kqBpPVQi~0fqr~hPmvO8Q6PolApyyi6W9hpuzejXT7cRipy1UH69AY7E1upLmVcsWqBnxptVDon3SwiuAoiGq3KrPIvX54vOJUfkTJ9S6uhg6Q-wJDnEb1DhoNtNRfTXmtuNqFrMBUwieDNhS9e0Jt985tMGFzrBQupjFZ7EIOa88FyZ3bA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA).

The sampleID can be mapped to run accessions in `PRJNA436454-huang.tsv`.

## `PRJNA436454-huang.tsv`

This is the report downloaded for bioproject [PRJNA436454](https://www.ebi.ac.uk/ena/browser/view/PRJNA436454).

It can be used to map the isolates in `huang-samplesheet.csv` to ENA runs to download the fastq files.

The `sample_title` column can be used to match to `sampleID` in `huang-samplesheet.csv`.

## `samplesheet.csv`

This is the "master" sheet for this project. See
the [data cleaning notebook](../workflow/notebooks/data-cleaning.ipynb) for the full
details of how this file was generated.

Essentially, all samples used in the WHO paper were removed from the non-WHO datasets in
the samplesheets in this directory and the remaining data were collated into this file.