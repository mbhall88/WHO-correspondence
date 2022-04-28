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

Gathering Delamanid (DLM) phenotypes
from https://journals.asm.org/doi/full/10.1128/JCM.01304-20

They used a critical concentraction of >0.06μg/ml - based on
the [WHO technical guidelines](https://apps.who.int/iris/bitstream/handle/10665/260470/WHO-CDS-TB-2018.5-eng.pdf)
.

MICs and accessions where taken from the supplementary Excel
spreadsheet https://journals.asm.org/doi/suppl/10.1128/JCM.01304-20/suppl_file/jcm.01304-20-sd004.xlsx

This file was taken from the following columns in the "Mutations list for DLM" sheet of
the above-mentioned Excel file.

- ISOLATE 1 (renamed to ISOLATE 1 DLM MIC (ug/ml) in this file)
- ISOLATE 2 (renamed to ISOLATE 2 DLM MIC (ug/ml) in this file)
- DRS Sample selected (Isolate 1)
- DRS Sample selected (Isolate 2)

The purpose of this file was to match up the isolate names and phenotypes with the
accessions in `dlm-accessions.csv`.

## `dlm-accessions.csv`

Gathering Delamanid (DLM) phenotypes
from https://journals.asm.org/doi/full/10.1128/JCM.01304-20

They used a critical concentraction of >0.06μg/ml - based on
the [WHO technical guidelines](https://apps.who.int/iris/bitstream/handle/10665/260470/WHO-CDS-TB-2018.5-eng.pdf)
.

MICs and accessions where taken from the supplementary Excel
spreadsheet https://journals.asm.org/doi/suppl/10.1128/JCM.01304-20/suppl_file/jcm.01304-20-sd004.xlsx

This data was extracted from the "WHO database column" of that Excel file.

Phenotypic DST of first and second line dugs. In brackets are indicated the used
critical concentration to defined Suceptibility (S) or Resistance (R) to each drug. N/A
means "not available" (reference: Zignol M et al., 2018. Lancet Infect Dis 18:675-683.
doi:10.1016/S1473-3099(18)30073-2).

## `huang-samplesheet.csv`

This CSV file contains the phenotypes for the isolates sequenced in [Huang *et
al.*](https://doi.org/10.1093/cid/ciy883).

In particular, this dataset provides linezolid phenotypes.

This file was extracted
from [Supplementary Table S6](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/cid/69/3/10.1093_cid_ciy883/1/ciy883_suppl_supplementary_table_s6.docx?Expires=1652148910&Signature=D7iy~iS4ZXFMCAifGgjnpD-OslN6pINjGYhFbqx1RI2unpvW9gaZ2CkiwXLd3cBagIABut8U4qKOXY11mOVw9LMZohqZNtkibKuu7SFgBJ-c2vMz9h10GNKHj5Ya98dg6AT7IPVsTk77OHopWoFsE6JlgbeCPIlH4i-kqBpPVQi~0fqr~hPmvO8Q6PolApyyi6W9hpuzejXT7cRipy1UH69AY7E1upLmVcsWqBnxptVDon3SwiuAoiGq3KrPIvX54vOJUfkTJ9S6uhg6Q-wJDnEb1DhoNtNRfTXmtuNqFrMBUwieDNhS9e0Jt985tMGFzrBQupjFZ7EIOa88FyZ3bA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)
.

The sampleID can be mapped to run accessions in `PRJNA436454-huang.tsv`.

## `PRJNA436454-huang.tsv`

This is the report downloaded for
bioproject [PRJNA436454](https://www.ebi.ac.uk/ena/browser/view/PRJNA436454).

It can be used to map the isolates in `huang-samplesheet.csv` to ENA runs to download
the fastq files.

The `sample_title` column can be used to match to `sampleID` in `huang-samplesheet.csv`.

## `bainomugisa-pza-phenotypes.csv`

Biosamples and pyrazinamide phenotypes for 100 isolates from [Arnold's MGen paper](https://doi.org/10.1099/mgen.0.000147).

## `bainomugisa-PRJNA385247.tsv`

This is the report downloaded for bioproject [PRJNA385247](https://www.ebi.ac.uk/ena/browser/view/PRJNA385247).

It can be used to map the isolates from `bainomugisa-pza-phenotypes.csv` to ENA run accessions.

## `smith-PRJNA650381.tsv`

This is the accessions for the Illumina data from https://doi.org/10.1128/JCM.00583-20

It was provided by Pascal Lapierre via email as the Illumina data was not part of the original submission, but has now been added to the same bioproject as the nanopore data.

## `smith-phenotypes.csv`

These phenotypes were taken from [Supplementary table 6](https://journals.asm.org/doi/suppl/10.1128/JCM.00583-20/suppl_file/jcm.00583-20-s0006.xlsx) of https://doi.org/10.1128/JCM.00583-20

I extracted the DST phenotypes only and removed all other gDST columns. The sample names can be used to match the phenotypes to the accessions in `smith-PRJNA650381.tsv`. I also changed the FLQ column to Ofloxacin as they state in the paper that they only did CFX DST and that it was considered a proxy for all fluoroquinolones.

## `samplesheet.csv`

This is the "master" sheet for this project. See
the [data cleaning notebook](../workflow/notebooks/data-cleaning.ipynb) for the full
details of how this file was generated.

Essentially, all samples used in the WHO paper were removed from the non-WHO datasets in
the samplesheets in this directory and the remaining data were collated into this file.

## `who-catalogue-trimmed.tsv.gz`

A slimmed down version of the WHO catalogue. We only keep the following columns

    - drug
    - variant
    - genome_position
    - Present_SOLO_R
    - Present_SOLO_SR
    - Present_S
    - Present_R
    - Absent_S
    - Absent_R
    - grading

The catalogue in the WHO paper was downloaded
from [Supplementary appendix 2](https://www.thelancet.com/cms/10.1016/S2666-5247(21)00301-3/attachment/77fc876a-afad-4c17-856e-6cc0d5951c29/mmc2.xlsx)
on 04/04/2022.

## `who-catalogue-errors.tsv`

These are variants where the insertion/deletion contradicts itself.

For example, `whiB7_192_del_1_gc_g` indicates that we have a deletion (`del`) of 1 base
- `gc->g`.

Any variant in this file indicates a length that is not consistent with the actual
difference in length between the reference and alternate allele.

For example, `whiB7_192_del_2_gc_g` would be considered an error.

## `who-panel.tsv`

This is the cleaned catalogue, minus the variants in `who-catalogue-errors.tsv`. The
columns in the file describe:

- **gene**: gene of the variant
- **mutation**: the gene coordinates of the variant in the respective gene. i.e. `A2T`
  is an `A` changed to a `T` at position 2. Whether the position is in nucleic acid or
  amino acid space depends on the value of **alphabet**.
- **alphabet**: What residue-space the variant is in - DNA or PROT (protein)
- **grading**: The WHO grading of the variant:
    1. Associated with resistance
    2. Associated with resistance - interim
    3. Uncertain significance
    4. Not associated with resistance - interim
    5. Not associated with resistance

This file is intended to be used for constructing a Mykrobe custom panel, with the
removal of the header, grading column, and drugs column. In addition, it is used to
construct the `var2res.json` file Mykrobe uses to map variants to the drug(s) they are
associated with.

## `who-results.csv`

These are the raw sensitivity and specificity numbers used to make figure 3 in the WHO catalogue paper.

It was taken from pages 7-9 of the [WHO catalogue report](https://www.who.int/publications/i/item/9789240028173).