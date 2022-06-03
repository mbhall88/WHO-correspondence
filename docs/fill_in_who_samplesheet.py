"""The accession from the WHO spreadsheet are somewhat incomplete.
For example, some rows have an ENA biosample accession in the run_accession column.
To make downstream processing easier, we want to try and get run accessions for as
many rows as possible.
Where this is not possible, we leave the existing accession as they are, but put them
the correct columns.
The final output sheet has a column for each accession type as listed in
https://ena-docs.readthedocs.io/en/latest/submit/general-guide/accessions.html
"""
import re
import sys
from dataclasses import dataclass
from enum import Enum
from typing import Union

from pysradb.search import EnaSearch, SraSearch

ENA_KEEP_METADATA = [
    "study_accession",
    "sample_accession",
    "run_accession",
    "experiment_accession",
]
SRA_KEEP_METADATA = [
    "study_accession",
    "sample_accession",
    "run_1_accession",
    "experiment_accession",
]


class AccessionType(Enum):
    BIOPROJECT = re.compile(r"PRJ([EDN])[A-Z]\d+")
    STUDY = re.compile(r"([EDS])RP\d{6,}")
    BIOSAMPLE = re.compile(r"SAM([EDN])[A-Z]?\d+")
    SAMPLE = re.compile(r"([EDS])RS\d{6,}")
    EXPERIMENT = re.compile(r"([EDS])RX\d{6,}")
    RUN = re.compile(r"([EDS])RR\d{6,}")

    @staticmethod
    def from_str(s: str) -> "AccessionType":
        for member in AccessionType:
            regex = member.value
            if regex.search(s):
                return member
        raise ValueError(f"{s} does not match a known type")


@dataclass
class Accession:
    run: str = ""
    sample: str = ""
    biosample: str = ""
    experiment: str = ""
    bioproject: str = ""
    study: str = ""

    @staticmethod
    def from_line(s: str, delim: str = ",") -> "Accession":
        acc = Accession()
        fields = [f for f in s.rstrip().split(delim) if f]
        for f in fields:
            acc_type = AccessionType.from_str(f)
            if acc_type is AccessionType.BIOSAMPLE:
                acc.biosample = f
            elif acc_type is AccessionType.RUN:
                acc.run = f
            elif acc_type is AccessionType.EXPERIMENT:
                acc.experiment = f
            elif acc_type is AccessionType.SAMPLE:
                acc.sample = f
            elif acc_type is AccessionType.STUDY:
                acc.study = f
            elif acc_type is AccessionType.BIOPROJECT:
                acc.bioproject = f
            else:
                raise NotImplementedError(
                    f"Don't know how to deal with AccessionType {acc_type}"
                )

        return acc

    def most_specific(self) -> Union[tuple[str, AccessionType], tuple[str, None]]:
        """Returns the most specific accession"""
        if self.run:
            return self.run, AccessionType.RUN
        elif self.experiment:
            return self.experiment, AccessionType.EXPERIMENT
        elif self.biosample:
            return self.biosample, AccessionType.BIOSAMPLE
        elif self.sample:
            return self.sample, AccessionType.SAMPLE
        elif self.bioproject:
            return self.bioproject, AccessionType.BIOPROJECT
        elif self.study:
            return self.study, AccessionType.STUDY
        else:
            return "", None

    def to_row(self, delim: str = ",") -> str:
        return delim.join(
            [
                self.bioproject,
                self.study,
                self.biosample,
                self.sample,
                self.experiment,
                self.run,
            ]
        )


def fetch_metadata(query: str):
    instance = EnaSearch(verbosity=2, accession=query)
    instance.search()
    df = instance.get_df()
    if df.empty:
        instance = SraSearch(verbosity=2, accession=query)
        instance.search()
        df = instance.get_df()
        if df.empty:
            return None
        keep = SRA_KEEP_METADATA
    else:
        keep = ENA_KEEP_METADATA
    return df.loc[:, keep]


def eprint(msg: str):
    print(msg, file=sys.stderr)


def main():
    input_file = sys.argv[1]
    if len(sys.argv) == 3:
        out_fp = open(sys.argv[2], "w")
    else:
        out_fp = sys.stdout

    print(
        ",".join(["bioproject", "study", "biosample", "sample", "experiment", "run"]),
        file=out_fp,
    )

    n_rows = 0
    seen_runs = set()
    rows_written = 0
    no_accs = 0
    only_proj = 0

    with open(input_file) as fp:
        _ = next(fp)  # skip header
        for line in fp:
            rows_to_write = []
            n_rows += 1
            if not line.strip(",").strip():
                no_accs += 1
                continue

            acc = Accession.from_line(line.replace("FAIL", ""))
            query, acc_type = acc.most_specific()

            if acc_type in (AccessionType.BIOPROJECT, AccessionType.STUDY):
                # we don't try and fetch data for project/study types
                only_proj += 1
                rows_to_write.append(acc)
            elif acc_type is AccessionType.RUN:
                if acc.run in seen_runs:
                    continue
                else:
                    seen_runs.add(acc.run)
                    rows_to_write.append(acc)
            else:
                metadata = fetch_metadata(query)
                if metadata is not None and not metadata.empty:
                    for row in metadata.values:
                        new_acc = Accession.from_line(",".join(row))
                        if new_acc.run:
                            if new_acc.run in seen_runs:
                                continue
                            else:
                                seen_runs.add(new_acc.run)
                        rows_to_write.append(new_acc)
                else:
                    if acc.run and acc.run not in seen_runs:
                        seen_runs.add(acc.run)
                    rows_to_write.append(acc)

            for row in rows_to_write:
                rows_written += 1
                print(row.to_row(","), file=out_fp)

    eprint(f"{n_rows} rows in original file")
    eprint(f"{rows_written} rows written to new file")
    eprint(f"{no_accs} rows had no accession of any sort")
    eprint(f"{only_proj} rows only had project/study accessions")


if __name__ == "__main__":
    main()
