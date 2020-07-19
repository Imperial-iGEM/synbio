"""Example of a Combinatorial Golden Gate assembly with steps and output."""

import os

from Bio.SeqIO import parse

from synbio.designs import Combinatorial
from synbio.protocols import GoldenGate

DIR_NAME = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.path.join(DIR_NAME, "..", "..", "data", "goldengate")

def read(filename):
    """Read in a single Genbank file from the test directory."""

    return next(parse(os.path.join(DATA_DIR, filename), "genbank"))

records = []
for (_, _, filenames) in os.walk(DATA_DIR):
    for file in filenames:
        gb = os.path.join(DATA_DIR, file)
        if not gb.endswith(".gb"):
            continue
        for record in parse(gb, "genbank"):
            records.append(record)

record_sets = []
for f_type in ["promoter", "RBS", "CDS", "terminator"]:

    def test(r):
        return any(
            f.type == f_type and f.location.start < 50 for f in r.features
        )

    new_bin = [r for r in records if test(r)][:5]
    record_sets.append(new_bin)  # add a new bin

records = [r for record_set in record_sets for r in record_set] + [
    read("DVK_AE.gb")
]

# create a combinatorial library design from multiple "bins"
design = Combinatorial(records, linear = False)

# create a protocol using Golden Gate as the sole composite step and run
protocol = GoldenGate(
    name="Combinatorial Golden Gate", design=design, include=["KanR"], min_count=5
)

# export all the output plasmids to a multi-FASTA
protocol.to_fasta("composite_parts.fasta")

# export plate layouts
protocol.to_csv("plate_layouts.csv")

# export human protocol
protocol.to_txt("protocol.txt")

# export a hamilton picklist
protocol.to_picklists("robotic_picklist.gwl", platform="tecan")
