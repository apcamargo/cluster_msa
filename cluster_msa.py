#!/usr/bin/env python

import argparse
import bz2
import csv
import gzip
import lzma
import pathlib
import sys
import textwrap
from collections import Counter
from contextlib import contextmanager
from enum import Enum, auto
from pathlib import Path

import numpy as np
from polyleven import levenshtein
from sklearn.cluster import DBSCAN, HDBSCAN


class Compression(Enum):
    bzip2 = auto()
    gzip = auto()
    xz = auto()
    zstd = auto()
    uncompressed = auto()


def is_compressed(filepath: Path):
    with open(filepath, "rb") as fin:
        signature = fin.peek(8)[:8]
        if tuple(signature[:2]) == (0x1F, 0x8B):
            return Compression.gzip
        elif tuple(signature[:3]) == (0x42, 0x5A, 0x68):
            return Compression.bzip2
        elif tuple(signature[:7]) == (0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00, 0x00):
            return Compression.xz
        elif tuple(signature[:4]) == (0x28, 0xB5, 0x2F, 0xFD):
            return Compression.zstd
        else:
            return Compression.uncompressed


@contextmanager
def open_file(filepath):
    filepath_compression = is_compressed(filepath)
    if filepath_compression == Compression.gzip:
        fin = gzip.open(filepath, "rt")
    elif filepath_compression == Compression.bzip2:
        fin = bz2.open(filepath, "rt")
    elif filepath_compression == Compression.xz:
        fin = lzma.open(filepath, "rt")
    else:
        fin = open(filepath, "r")
    try:
        yield fin
    finally:
        fin.close()


class Sequence:
    def __init__(self, header: str, seq: str):
        self._header = header
        self._seq = seq.encode("ascii")

    @property
    def header(self):
        return self._header

    @property
    def accession(self):
        return self._header.split()[0]

    @property
    def seq(self):
        return self._seq.decode()

    def __str__(self):
        return (
            f">{self.header}\n{textwrap.fill(self.seq, 60, break_on_hyphens=False)}\n"
        )


def read_fasta(filepath, uppercase=False, strip_n=False, compress=False):
    with open_file(filepath) as fin:
        last = None
        while True:
            if not last:
                for line in fin:
                    if line[0] == ">":
                        last = line[:-1]
                        break
            if not last:
                break
            name, seqs, last = last[1:], [], None
            for line in fin:
                if line[0] == ">":
                    last = line[:-1]
                    break
                seqs.append(line[:-1])
            seqs = "".join(seqs)
            if uppercase:
                seqs = seqs.upper()
            if strip_n:
                seqs = seqs.strip("nN")
            if len(seqs):
                yield Sequence(name, seqs)
            if not last:
                break


def remove_lowercase_columns(seqs):
    alignment_length = len(seqs[0])
    mask = [True] * alignment_length
    for seq in seqs:
        for i, char in enumerate(seq):
            if char.islower():
                mask[i] = False
    return list(
        map(lambda x: "".join([char for i, char in enumerate(x) if mask[i]]), seqs)
    )


def write_fasta(accessions, seqs, outfile, max_gap_column=0.5):
    num_seqs = len(seqs)
    non_empty_columns = [
        col
        for col in zip(*seqs)
        if sum(c == "-" for c in col) / num_seqs <= max_gap_column
    ]
    processed_seqs = ["".join(seq) for seq in zip(*non_empty_columns)]
    with open(outfile, "w") as fo:
        for a, s in zip(accessions, processed_seqs):
            record = Sequence(a, s)
            fo.write(str(record))


def get_consensus_seq(seqs):
    consensus = ""
    residues = "ACDEFGHIKLMNPQRSTVWY-"
    for i in range(len(seqs[0])):
        res_array = [x[i] for x in seqs]
        res_count = np.array([res_array.count(a) for a in list(residues)])
        vote = np.argmax(res_count)
        consensus += residues[vote]
    return consensus


def encode_seqs_bl62(seqs):
    max_len = max(len(seq) for seq in seqs)
    arr = [
        [BL62NP.get(c, BL62NP["-"]) for c in seq.upper().ljust(max_len, "-")]
        for seq in seqs
    ]
    return np.array(arr).reshape(len(seqs), max_len * len(BL62NP["-"]))


def encode_seqs_ohe(seqs):
    alphabet = "ACDEFGHIKLMNPQRSTVWY-"
    char_to_index = {c: i for i, c in enumerate(alphabet)}
    max_len = max(len(seq) for seq in seqs)
    arr = np.zeros((len(seqs), max_len, len(alphabet)), dtype=int)
    for j, seq in enumerate(seqs):
        for i, c in enumerate(seq.upper()):
            arr[j, i, char_to_index.get(c, len(alphabet) - 1)] = 1
    return arr.reshape(len(seqs), max_len * len(alphabet))


def pick_dbscan_epsilon(data, min_samples, min_eps, max_eps, eps_step):
    n_clusters = []
    eps_test_vals = np.arange(min_eps, max_eps + 1e-6, eps_step)
    for eps in eps_test_vals:
        clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(data)
        n = len(set(clustering.labels_))
        if n == 1 and len(n_clusters) and n_clusters[-1] > 1:
            break
        n_clusters.append(n)
    return eps_test_vals[np.argmax(n_clusters)]


BL62NP = {
    "A": [
        -0.31230882,
        -0.53572156,
        -0.01949946,
        -0.12211268,
        -0.70947917,
        -0.42211092,
        0.02783931,
        0.02637933,
        -0.41760305,
        0.21809875,
        0.53532768,
        0.04833016,
        0.07877711,
        0.50464914,
        -0.26972087,
        -0.52416842,
    ],
    "R": [
        0.29672002,
        0.29005364,
        0.18176298,
        -0.05103382,
        -0.34686519,
        0.58024228,
        -0.49282931,
        0.62304281,
        -0.09575202,
        0.30115555,
        0.09913529,
        0.1577466,
        -0.94391939,
        -0.10505925,
        0.05482389,
        0.38409897,
    ],
    "N": [
        -0.42212537,
        0.12225749,
        0.16279646,
        0.60099009,
        0.19734216,
        0.42819919,
        -0.33562418,
        0.17036334,
        0.4234109,
        0.46681561,
        -0.50347222,
        -0.37936876,
        0.1494825,
        0.32176759,
        0.28584684,
        0.68469861,
    ],
    "D": [
        0.18599294,
        -0.44017825,
        -0.4476952,
        0.34340976,
        0.44603553,
        0.40974629,
        -0.60045935,
        -0.09056728,
        0.22147919,
        -0.33029418,
        0.55635594,
        -0.54149972,
        0.05459062,
        0.57334159,
        -0.06227118,
        0.65299872,
    ],
    "C": [
        -0.19010428,
        0.64418792,
        -0.85286762,
        0.21380295,
        0.37639516,
        -0.67753593,
        0.38751609,
        0.55746524,
        0.01443766,
        0.1776535,
        0.62853954,
        -0.15048523,
        0.55100206,
        -0.21426656,
        0.3644061,
        -0.0018255,
    ],
    "Q": [
        0.7350723,
        0.10111267,
        0.55640019,
        -0.18226966,
        0.51658102,
        -0.19321508,
        -0.46599027,
        -0.02989911,
        0.4036196,
        -0.11978213,
        -0.29837524,
        -0.30232765,
        -0.36738065,
        -0.1379793,
        0.04362871,
        0.33553714,
    ],
    "E": [
        0.41134047,
        0.13512443,
        0.62492322,
        -0.10120261,
        -0.03093491,
        0.23751917,
        -0.68338694,
        0.05124762,
        0.41533821,
        0.46669353,
        0.31467277,
        -0.02427587,
        0.15361135,
        0.70595112,
        -0.27952632,
        0.32408931,
    ],
    "G": [
        -0.33041265,
        -0.43860065,
        -0.5509376,
        -0.04380843,
        -0.35160935,
        0.25134855,
        0.53409314,
        0.54850824,
        0.59490287,
        0.32669345,
        -0.45355268,
        -0.56317041,
        -0.55416297,
        0.18117841,
        -0.71600849,
        -0.08989825,
    ],
    "H": [
        -0.40366849,
        0.10978974,
        0.0280101,
        -0.46667987,
        -0.45607028,
        0.54114052,
        -0.77552923,
        -0.10720425,
        0.55252091,
        -0.34397153,
        -0.59813694,
        0.15567728,
        0.03071009,
        -0.02176143,
        0.34442719,
        0.14681541,
    ],
    "I": [
        0.19280422,
        0.35777863,
        0.06139255,
        0.20081699,
        -0.30546596,
        -0.56901549,
        -0.15290953,
        -0.31181573,
        -0.74523217,
        0.22296016,
        -0.39143832,
        -0.16474685,
        0.58064427,
        -0.77386654,
        0.19713107,
        -0.49477418,
    ],
    "L": [
        -0.16133903,
        0.22112761,
        -0.53162136,
        0.34764073,
        -0.08522381,
        -0.2510216,
        0.04699411,
        -0.25702389,
        -0.8739765,
        -0.24171728,
        -0.24370533,
        0.42193635,
        0.41056913,
        -0.60378211,
        -0.65756832,
        0.0845203,
    ],
    "K": [
        -0.34792144,
        0.18450939,
        0.77038332,
        0.63868511,
        -0.06221681,
        0.11930421,
        0.04895523,
        -0.22463059,
        -0.03268844,
        -0.58941354,
        0.11640045,
        0.32384901,
        -0.42952779,
        0.58119471,
        0.07288662,
        0.26669673,
    ],
    "M": [
        0.01834555,
        -0.16367754,
        0.34900298,
        0.45087949,
        0.47073855,
        -0.37377404,
        0.0606911,
        0.2455703,
        -0.55182937,
        -0.20261009,
        0.28325423,
        -0.04741146,
        0.30565238,
        -0.62090653,
        0.17528413,
        -0.60434975,
    ],
    "F": [
        -0.55464981,
        0.50918784,
        -0.21371646,
        -0.63996967,
        -0.37656862,
        0.27852662,
        0.3287838,
        -0.56800869,
        0.23260763,
        -0.20653106,
        0.63261439,
        -0.22666691,
        0.00726302,
        -0.60125196,
        0.07139961,
        -0.35086639,
    ],
    "P": [
        0.94039731,
        -0.25999326,
        0.43922549,
        -0.485738,
        -0.20492235,
        -0.26005626,
        0.68776626,
        0.57826888,
        -0.05973995,
        -0.1193658,
        -0.12102433,
        -0.22091354,
        0.43427913,
        0.71447886,
        0.32745991,
        0.03466398,
    ],
    "S": [
        -0.13194625,
        -0.12262688,
        0.18029209,
        0.16555524,
        0.39594125,
        -0.58110665,
        0.16161717,
        0.0839783,
        0.0911945,
        0.34546976,
        -0.29415349,
        0.29891936,
        -0.60834721,
        0.5943593,
        -0.29473819,
        0.4864154,
    ],
    "T": [
        0.40850093,
        -0.4638894,
        -0.39732987,
        -0.01972861,
        0.51189582,
        0.10176704,
        0.37528519,
        -0.41479418,
        -0.1932531,
        0.54732221,
        -0.11876511,
        0.32843973,
        -0.259283,
        0.59500132,
        0.35168375,
        -0.21733727,
    ],
    "W": [
        -0.50627723,
        -0.1973602,
        -0.02339884,
        -0.66846048,
        0.62696606,
        0.60049717,
        0.69143364,
        -0.48053591,
        0.17812208,
        -0.58481821,
        -0.23551415,
        -0.06229112,
        0.20993116,
        -0.72485884,
        0.34375662,
        -0.23539168,
    ],
    "Y": [
        -0.51388312,
        -0.2788953,
        0.00859533,
        -0.5247195,
        -0.18021544,
        0.28372911,
        0.10791359,
        0.13033494,
        0.34294013,
        -0.70310089,
        -0.13245433,
        0.48661081,
        0.08451644,
        -0.69990992,
        0.0408274,
        -0.47204888,
    ],
    "V": [
        0.68546275,
        0.22581365,
        -0.32571833,
        0.34394298,
        -0.43232367,
        -0.5041842,
        0.04784017,
        -0.53067936,
        -0.50049908,
        0.36874221,
        0.22429186,
        0.4616482,
        0.11159174,
        -0.26827959,
        -0.39372848,
        -0.40987423,
    ],
}
BL62NP["-"] = np.mean(list(BL62NP.values()), axis=0).tolist()


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="""
    Cluster sequences in a MSA using HDBSCAN algorithm and write FASTA file for each cluster.
    """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=False,
    )

    p1 = p.add_argument_group("output options")
    p2 = p.add_argument_group("MSA filtering options")
    p3 = p.add_argument_group("clustering options")
    p.add_argument(
        "input", action="store", help="Multiple sequence alignment (FASTA format)."
    )
    p.add_argument("output", action="store", help="Output directory.")
    p.add_argument("prefix", action="store", help="Prefix for the output files.")
    p1.add_argument(
        "--write-tabular-outputs",
        action="store_true",
        help="Write tabulars outputs that describe cluster assignments and cluster statistics.",
    )
    p2.add_argument(
        "--remove-lowercase-columns",
        action="store_true",
        help="Remove columns with at least one lowercase residue (before distance computation).",
    )
    p2.add_argument(
        "--max-gap-sequence",
        action="store",
        default=0.25,
        type=float,
        help="Remove sequences with gaps representing more than this fraction of the sequence (before distance computation).",
    )
    p2.add_argument(
        "--max-gap-column",
        action="store",
        default=0.75,
        type=float,
        help="Remove columns with gaps representing more than this fraction of the column (after distance computation).",
    )
    p3.add_argument(
        "--use-one-hot-encoding",
        action="store_true",
        help="Use one-hot encoding instead of the BLOSUM62-like embedding.",
    )
    p3.add_argument(
        "--use-hdbscan",
        action="store_true",
        help="Use HDBSCAN instead of DBSCAN.",
    )
    p3.add_argument(
        "--min-cluster-size",
        action="store",
        default=3,
        type=int,
        help="Minimum cluster size to consider.",
    )
    p3.add_argument(
        "--min-samples",
        action="store",
        default=3,
        type=int,
        help="The minimum number of neighbors required for a sequence to be considered a core data point.",
    )
    p3.add_argument(
        "--min-eps",
        action="store",
        default=0.1,
        type=float,
        help="Minimum epsilon value to be considered in the evaluation range (only for DBSCAN).",
    )
    p3.add_argument(
        "--max-eps",
        action="store",
        default=45,
        type=float,
        help="Maximum epsilon value to be considered in the evaluation range (only for DBSCAN).",
    )
    p3.add_argument(
        "--eps-step",
        action="store",
        default=0.1,
        type=float,
        help="Step size for incrementing epsilon values within the evaluation range (only for DBSCAN).",
    )
    p3.add_argument(
        "--set-epsilon-value",
        action="store",
        type=float,
        help="Set the maximum distance between two sequences for them to be considered neighbors, disabling the automatic epsilon selection (only for DBSCAN).",
    )
    p3.add_argument(
        "--cluster-selection-epsilon",
        action="store",
        default=0.0,
        type=float,
        help="Clusters with a distance below this threshold will be merged (only for HDBSCAN).",
    )

    # Print help if no arguments are provided
    if len(sys.argv) == 1:
        p.print_help(sys.stderr)
        sys.exit(1)

    args = p.parse_args()

    output_dir = pathlib.Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    accs, seqs = [], []
    for seq in read_fasta(args.input):
        accs.append(seq.accession)
        seqs.append(seq.seq)

    if args.remove_lowercase_columns:
        seqs = remove_lowercase_columns(seqs)

    if not len(seqs[0]):
        print(
            "Error: No columns left after filtering out lowercase columns.",
            file=sys.stderr,
        )
        sys.exit(1)

    data = [{"sequence_name": acc, "sequence": seq} for acc, seq in zip(accs, seqs)]
    aln_len = len(data[0]["sequence"])

    # Calculate frac_gaps and filter sequences
    data = [
        {**row, "frac_gaps": row["sequence"].count("-") / aln_len}
        for row in data
        if row["sequence"].count("-") / aln_len < args.max_gap_sequence
    ]

    if not data:
        print(
            "Error: No sequences left after filtering by gap fraction.", file=sys.stderr
        )
        sys.exit(1)

    if args.use_one_hot_encoding:
        encoded_seqs = encode_seqs_ohe([row["sequence"] for row in data])
    else:
        encoded_seqs = encode_seqs_bl62([row["sequence"] for row in data])

    if args.use_hdbscan:
        clustering = HDBSCAN(
            min_samples=args.min_samples,
            cluster_selection_epsilon=args.cluster_selection_epsilon,
        ).fit(encoded_seqs)
    else:
        if args.set_epsilon_value:
            clustering = DBSCAN(
                min_samples=args.min_samples, eps=args.set_epsilon_value
            ).fit(encoded_seqs)
        else:
            eps = pick_dbscan_epsilon(
                encoded_seqs,
                args.min_samples,
                args.min_eps,
                args.max_eps,
                args.eps_step,
            )
            clustering = DBSCAN(eps=eps, min_samples=args.min_samples).fit(encoded_seqs)

    cluster_sizes = Counter(clustering.labels_)
    for row, label in zip(data, clustering.labels_):
        if cluster_sizes[label] < args.min_cluster_size:
            label = -1
        row["clustering_label"] = label
    clusters = sorted(
        {row["clustering_label"] for row in data if row["clustering_label"] >= 0}
    )

    cluster_metadata = []
    for c in clusters:
        cluster_data = [row for row in data if row["clustering_label"] == c]
        write_fasta(
            [row["sequence_name"] for row in cluster_data],
            [row["sequence"] for row in cluster_data],
            f"{pathlib.Path(args.output) / args.prefix}_cluster_{c}.afa",
            args.max_gap_column,
        )

        if args.write_tabular_outputs:
            cs = get_consensus_seq([row["sequence"] for row in cluster_data])
            avg_dist_to_cs = np.mean(
                [1 - levenshtein(row["sequence"], cs) / aln_len for row in cluster_data]
            )
            cluster_metadata.append(
                {
                    "cluster_index": c,
                    "cluster_size": len(cluster_data),
                    "avg_identity_to_consensus": f"{avg_dist_to_cs:.4f}",
                    "consensus_sequence": cs,
                }
            )

    if args.write_tabular_outputs:
        assignments_outfile = (
            pathlib.Path(args.output) / f"{args.prefix}_cluster_assignments.tsv"
        )
        with open(assignments_outfile, "w", newline="") as csvfile:
            writer = csv.DictWriter(
                csvfile,
                fieldnames=["sequence_accession", "cluster_index"],
                delimiter="\t",
            )
            writer.writeheader()
            for row in data:
                writer.writerow(
                    {
                        "sequence_accession": row["sequence_name"],
                        "cluster_index": row["clustering_label"],
                    }
                )

        data_outfile = pathlib.Path(args.output) / f"{args.prefix}_cluster_data.tsv"
        with open(data_outfile, "w", newline="") as csvfile:
            writer = csv.DictWriter(
                csvfile,
                fieldnames=[
                    "cluster_index",
                    "cluster_size",
                    "avg_identity_to_consensus",
                    "consensus_sequence",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            for row in cluster_metadata:
                writer.writerow(row)
