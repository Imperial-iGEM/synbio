"""Subcloning for digestion and ligation of SeqRecords together."""

from collections import defaultdict
import logging
from typing import Dict, List, Set, Tuple, Iterable, Union

from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
from Bio.Restriction import RestrictionBatch, BsaI, BpiI
from Bio.Restriction.Restriction import RestrictionType
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import networkx as nx
from networkx.algorithms.cycles import simple_cycles
from networkx.exception import NetworkXNoCycle

from ..containers import content_id
from ..designs import Plasmid, CombinatorialBins


CATALYZE_CACHE: Dict[str, List[Tuple[str, SeqRecord, str]]] = {}
"""Store the catalyze results of each SeqRecord. Avoid lots of string searches."""


def goldengate(
    record_set: Iterable[List[SeqRecord]],
    include: List[str] = None,
    min_count: int = -1,
) -> List[Tuple[List[SeqRecord], List[SeqRecord]]]:
    """Simulate a digestion and ligation using BsaI and BpiI.

    Accepts lists of SeqRecords that may combine and returns the lists
    of SeqRecords that may circularize into new vectors.

    Arguments:
        record_set {Iterable[List[SeqRecord]]} -- possible combinations of fragments

    Keyword Arguments:
        include {List[str]} -- the feature to filter assemblies on (default: {""})
        min_count {int} -- minimum number of SeqRecords for an assembly to be considered

    Returns:
        List[Tuple[List[SeqRecord], List[SeqRecord]]] -- list of tuples with:
            1. plasmids that will form
            2. SeqRecords that went into each formed plasmid
    """

    return subclone_many(record_set, [BsaI, BpiI], include, min_count)


def subclone_many(
    design: Iterable[List[SeqRecord]],
    enzymes: List[RestrictionType],
    include: List[str] = None,
    min_count: int = -1,
) -> List[Tuple[List[SeqRecord], List[SeqRecord]]]:
    """Parse a single list of SeqRecords to find all circularizable plasmids.

    Turn each SeqRecord's post-digest seqs into a graph where the nodes are
    the overhangs and the edges are the linear fragments
    post-digest/catalyzing with BsaI/BpiI.

    Arguments:
        record_set {List[SeqRecord]} -- single record set that might circularize
        enzymes {List[Enzyme]} -- list of enzymes to digest the input records with
        include {str} -- the include to filter assemblies
        min_count {int} -- mininum number of SeqRecords for an assembly to be considered

    Returns:
        List[Tuple[List[SeqRecord], List[SeqRecord]]] -- list of tuples with:
            1. plasmids that will form
            2. SeqRecords that went into each formed plasmid
    """

    seen_fragment_ids: Set[str] = set()
    all_plasmids_and_fragments: List[Tuple[List[SeqRecord], List[SeqRecord]]] = []
    for record_set in design:
        plasmids_and_fragments = subclone(record_set, enzymes, include, min_count)
        for plasmids, fragments in plasmids_and_fragments:

            # we don't want to re-use the fragment combination more than once
            fragment_ids = _hash_fragments(fragments)
            if fragment_ids in seen_fragment_ids:
                continue
            seen_fragment_ids.add(fragment_ids)

            all_plasmids_and_fragments.append((plasmids, fragments))
    return all_plasmids_and_fragments


def subclone(
    record_set: List[SeqRecord],
    enzymes: List[RestrictionType],
    include: List[str] = None,
    min_count: int = -1,
) -> List[Tuple[List[SeqRecord], List[SeqRecord]]]:
    """Parse a single list of SeqRecords to find all circularizable plasmids.

    Turn each SeqRecord's post-digest seqs into a graph where the nodes are
    the overhangs and the edges are the linear fragments
    post-digest/catalyzing with BsaI/BpiI.

    Arguments:
        record_set {List[SeqRecord]} -- single record set that might circularize
        enzymes {List[Enzyme]} -- list of enzymes to digest the input records with
        include {str} -- the include to filter assemblies
        min_count {int} -- mininum number of SeqRecords for an assembly to be considered

    Returns:
        List[Tuple[List[SeqRecord], List[SeqRecord]]] -- list of tuples with:
            1. plasmids that will form
            2. SeqRecords that went into each formed plasmid
    """

    graph = nx.MultiDiGraph()

    include = include or []
    seen_seqs: Set[str] = set()  # stored list of input seqs (not new combinations)
    for record in record_set:
        seen_seqs.add(str(record.seq + record.seq))
        seen_seqs.add(str((record.seq + record.seq).reverse_complement()))

        for left, frag, right in _catalyze(record, enzymes):
            # print(left, right, record.id)
            graph.add_node(left)
            graph.add_node(right)
            graph.add_edge(left, right, frag=frag)

    try:  # find all circularizable cycles
        # print(enzymes, graph.nodes)
        cycles = simple_cycles(graph)
    except NetworkXNoCycle:
        return []

    # get the fragments, enzymes back out of the cycle
    ids_to_fragments: Dict[str, List[SeqRecord]] = defaultdict(list)
    ids_to_plasmids: Dict[str, List[SeqRecord]] = defaultdict(list)
    for cycle in cycles:
        if min_count > 0 and len(cycle) < min_count:
            continue

        combinations = CombinatorialBins()
        for i, overhang in enumerate(cycle):
            next_overhang = cycle[(i + 1) % len(cycle)]
            record_bin = []
            for out_edge in graph.out_edges:
                src, dest, index = out_edge
                if src != overhang or dest != next_overhang:
                    continue
                record_bin.append(graph.edges[src, dest, index]["frag"])
            combinations.append(record_bin)

        for fragments in combinations:
            # make sure it's not just a re-ligation of insert + backbone
            new_plasmid = "".join([str(f.seq) for f in fragments])
            if any(new_plasmid in seq for seq in seen_seqs):
                continue

            # filter for plasmids that have the include feature
            if include:
                has_include = any(_has_feature(f, include) for f in fragments)
                if not has_include:
                    continue

            # try and re-order the fragments to match the input order
            fragment_ids = [content_id(r) for r in record_set]
            fragment_indexes = [fragment_ids.index(f.id) for f in fragments]
            fragment_min_index = min(fragment_indexes)
            fragment_first = fragment_indexes.index(fragment_min_index)
            fragments = fragments[fragment_first:] + fragments[:fragment_first]

            # create the composite plasmid
            plasmid = SeqRecord(Seq("", IUPACUnambiguousDNA()))
            for fragment in fragments:
                plasmid += fragment.upper()
            plasmid.id = "+".join(f.id for f in fragments if f.id != "<unknown id>")
            seen_seqs.add(str(plasmid.seq + plasmid.seq))
            seen_seqs.add(str((plasmid.seq + plasmid.seq).reverse_complement()))

            # make a unique id for the fragments
            fragments_id = _hash_fragments(fragments)
            ids_to_fragments[fragments_id] = fragments
            ids_to_plasmids[fragments_id].append(plasmid)

    plasmids_and_fragments: List[Tuple[List[SeqRecord], List[SeqRecord]]] = []
    for ids, fragments in ids_to_fragments.items():
        plasmids = ids_to_plasmids[ids]
        plasmids_and_fragments.append((plasmids, fragments))
    return plasmids_and_fragments


def _hash_fragments(record_set: List[SeqRecord]) -> str:
    """Create a unique ID for a list of records
    
    Arguments:
        record_set {List[SeqRecord]} -- set of Records to make unique ID for
    
    Returns:
        str -- unique ID concatenating records IDs/Seqs
    """

    fragment_ids = [content_id(f) for f in record_set]
    fragment_ids = sorted(fragment_ids)
    return "".join(fragment_ids)


def _catalyze(
    record: SeqRecord, enzymes: List[RestrictionType]
) -> List[Tuple[str, SeqRecord, str]]:
    """Catalyze a SeqRecord and return all post-digest SeqRecords with overhangs.

    Overhangs are returned as the overhang plus the position of the cut
    in the 5' end (^) and 3' end (_). So a 5' overhang may be:
    ^AAAA_. But a 3' overhang may be: _AAAA^.

    Arguments:
        record {SeqRecord} -- the SeqRecord to digest with enzymes
        enzymes {List[Enzyme]} -- list of enzymes to digest the input records with

    Returns:
        List[Tuple[str, SeqRecord, str]] --
            left overhang, cut fragment, right overhang
    """

    record_id = _record_id(record)
    if record_id != "<unknown id>" and record_id in CATALYZE_CACHE:
        return CATALYZE_CACHE[record_id]

    record = record.upper()
    batch = RestrictionBatch(enzymes)
    batch_sites = batch.search(record.seq, linear=False)

    # order all cuts with enzymes based on index
    cuts_seen: Set[int] = set()
    enzyme_cuts: List[Tuple[RestrictionType, int]] = []
    for enzyme, cuts in batch_sites.items():
        for cut in cuts:
            if cut in cuts_seen:
                continue
            cuts_seen.add(cut)
            enzyme_cuts.append((enzyme, cut - 1))  # revert to 0-based
    enzyme_cuts = sorted(enzyme_cuts, key=lambda x: x[1])

    # list of left/right overhangs for each fragment
    frag_w_overhangs: List[Tuple[str, SeqRecord, str]] = []
    for i, (enzyme, cut) in enumerate(enzyme_cuts):
        next_enzyme, next_cut = enzyme_cuts[(i + 1) % len(enzyme_cuts)]

        enzyme_len = len(enzyme.ovhgseq)
        next_enzyme_len = len(next_enzyme.ovhgseq)

        cut = cut + enzyme_len if enzyme.is_3overhang() else cut
        next_cut = (
            next_cut + next_enzyme_len if next_enzyme.is_3overhang() else next_cut
        )

        cut_rc = cut if enzyme.is_3overhang() else cut + enzyme_len
        next_cut_rc = (
            next_cut if next_enzyme.is_3overhang() else next_cut + next_enzyme_len
        )

        # find the cutsite sequences
        left = (
            record[cut : cut - enzyme_len]
            if enzyme.is_3overhang()
            else record[cut : cut + enzyme_len]
        )
        right = (
            record[next_cut : next_cut - next_enzyme_len]
            if next_enzyme.is_3overhang()
            else record[next_cut : next_cut + next_enzyme_len]
        )
        left_rc = right.reverse_complement()
        right_rc = left.reverse_complement()

        left = str(left.seq)
        right = str(right.seq)
        left_rc = str(left_rc.seq)
        right_rc = str(right_rc.seq)

        if next_enzyme.is_3overhang():
            left += "^"
            right += "^"
            left_rc += "^"
            right_rc += "^"
        else:
            left = "^" + left
            right = "^" + right
            left_rc = "^" + left_rc
            right_rc = "^" + right_rc

        frag = record[cut:next_cut]
        frag_rc = record[cut_rc:next_cut_rc].reverse_complement()
        frag_rc.id = record.id

        if next_cut < cut:  # wraps around the zero-index
            frag = (record + record)[cut : next_cut + len(record)]
            frag.id = record.id
            frag_rc = (record + record)[
                cut_rc : next_cut_rc + len(record)
            ].reverse_complement()
            frag_rc.id = record.id

        frag_w_overhangs.append((left, frag, right))
        frag_w_overhangs.append((left_rc, frag_rc, right_rc))

    CATALYZE_CACHE[record_id] = frag_w_overhangs  # store for future look-ups

    return frag_w_overhangs


def _has_feature(record: SeqRecord, include: List[str]) -> bool:
    """Return whether any of a record's features/qualifiers match the include specified.

    Arguments:
        record {SeqRecord} -- the record being checked for include
        include {List[str]} -- the include to filter for

    Returns:
        bool -- whether the record has any features or qualifiers with specified include
    """

    features: Set[str] = set()
    for feature in record.features:
        features.add(feature.id.lower())
        for _, value in feature.qualifiers.items():
            for v in value:
                features.add(v.lower())

    for feature in features:
        for keyword in include:
            if keyword.lower() in feature:
                return True

    return False


def _record_id(record: SeqRecord) -> str:
    """Get the record id, unique, for a seqRecord."""

    return record.id if record.id != "<unknown id>" else str(record.seq)
