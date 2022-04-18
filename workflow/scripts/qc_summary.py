import sys

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path
import pandas as pd

stats_paths = list(map(Path, sorted(snakemake.input.stats_paths)))
keep_ids_paths = list(map(Path, sorted(snakemake.input.keep_ids_paths)))
contam_ids_paths = list(map(Path, sorted(snakemake.input.contam_ids_paths)))
unmapped_ids_paths = list(map(Path, sorted(snakemake.input.unmapped_ids_paths)))

assert (
    len(stats_paths)
    == len(keep_ids_paths)
    == len(contam_ids_paths)
    == len(unmapped_ids_paths)
)

it = zip(stats_paths, keep_ids_paths, contam_ids_paths, unmapped_ids_paths)

data = [["run", "coverage", "f_keep", "f_contam", "f_unmapped"]]

for s, k, c, u in it:
    r1 = s.parts[-2]
    r2 = k.parts[-2]
    r3 = c.parts[-2]
    r4 = u.parts[-2]
    assert r1 == r2 == r3 == r4
    run = r1

    stats = pd.read_csv(s, sep="\t")
    sum_lens = list(stats["sum_len"])
    assert len(sum_lens) == 1
    sum_len = sum_lens[0]

    cov = sum_len / snakemake.params.genome_size

    with open(k) as fp:
        n_keep_ids = sum(1 for line in fp)

    with open(c) as fp:
        n_contam_ids = sum(1 for line in fp)

    with open(u) as fp:
        n_unmapped_ids = sum(1 for line in fp)

    N = n_unmapped_ids + n_contam_ids + n_keep_ids

    data.append([run, cov, n_keep_ids / N, n_contam_ids / N, n_unmapped_ids / N])

with open(snakemake.output.summary, "w") as fp:
    for row in data:
        print(",".join(map(str, row)), file=fp)
