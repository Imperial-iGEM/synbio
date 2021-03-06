"""Example of a Golden Gate Assembly from a single list of SeqRecords into a plasmid.

The SeqRecords are combined into new plasmid SeqRecords:
all new combinations of SeqRecords. An error is thrown if no new plasmids are possible.
"""

import os

from Bio.SeqIO import parse

from synbio.designs import Plasmid
from synbio.protocols import GoldenGate


def read(filename):
    """read a SeqRecord from a file from a Genbank in the "goldengate" dir"""

    filename = os.path.join("..", "data", "goldengate", filename)
    return next(parse(filename, "genbank"))


# read in each SeqRecord: a promoter, RBS, CDS, terminator and backbone
design = Plasmid(
    read(r)
    for r in [
        "J23100_AB.gb",
        "B0032m_BC.gb",
        "C0012m_CD.gb",
        "B0015_DE.gb",
        "DVK_AE.gb",
    ]
)

# create a protocol using GoldenGate as the sole composite step and run
# filter on KanR in a backbone
protocol = GoldenGate(name="Golden Gate Plasmid", design=design, min_count=5)

# export the composite plasmid
protocol.to_fasta("plasmid.fasta")

# export plate layouts
protocol.to_csv("plates.csv")

# export human readable protocol
protocol.to_txt("protocol.txt")
