import argparse
import re
import gzip
import sys
from pathlib import Path
from typing import FrozenSet, Dict, List, Tuple
from Bio.SeqIO.FastaIO import SimpleFastaParser

FLYBASE_ID_REGEX = re.compile(r"^\s*FB\w\w\d+\s*$")
PARENT_REGEX = re.compile(r"parent=([\w,]+);")

# Types
FlyBaseIdList = FrozenSet[str]
ParsedFastaRecords = Tuple[Tuple[str, str], ...]
FastaIdx = Dict[str, List[int]]


def get_parent_ids(fasta_header: str) -> FlyBaseIdList or None:
    """
    Returns a frozen Set of FlyBase IDs from the parent= attribute of the FASTA header.
    :param fasta_header: FASTA header
    :return: Frozen Set of FlyBase IDs
    """
    parent_tags = PARENT_REGEX.findall(fasta_header)
    parent_ids = []
    for tag in parent_tags:
        parent_ids.extend(tag.split(","))
    return frozenset(parent_ids)


def read_fasta(
    fasta_file: Path,
) -> Tuple[ParsedFastaRecords, FastaIdx] or None:
    """
    Reads a fasta file and returns a tuple of tuples (header, sequence) and a dictionary of FlyBase IDs
    as the key and a list of FASTA record indices as the value. The indices correspond to the
    position of the FASTA record in the first tuple.

    :param fasta_file: Path to the FASTA file
    :return: Tuple of ParsedFastaRecords (tuple of header and sequence) and a dictionary of FlyBase IDs as the key
             and a list of FASTA record indices as the value.
    """
    records = []
    idx = {}
    try:
        # Read in gzipped FASTA file
        if ".gz" in fasta_file.suffixes:
            with gzip.open(fasta_file, "rt") as gzfh:
                # Use the SimpleFastaParser for speed since we don't need to parse the FASTA header.
                for i, record in enumerate(SimpleFastaParser(gzfh)):
                    records.append(record)
                    parent_ids = get_parent_ids(record[0])
                    for parent_id in parent_ids:
                        # Set the index for the parent ID, adding it or creating a new list.
                        idx[parent_id] = idx.get(parent_id, []) + [i]
        else:
            print("Only gzipped FASTA files are supported currently", file=sys.stderr)
        return tuple(records), idx
    except FileNotFoundError:
        print(f"{fasta_file} not found", file=sys.stderr)
        return None
    except IOError:
        print(f"IO Error while reading {fasta_file}", file=sys.stderr)
        return None


def fbid_to_fasta(
    ids: FlyBaseIdList, fasta: ParsedFastaRecords, index: FastaIdx, output: Path
) -> None:
    """
    Writes the FASTA records corresponding to the supplied FlyBase IDs to the output file.

    :param ids: User supplied FlyBase IDs
    :param fasta: Tuple of tuples (header, sequence)
    :param index: Dictionary of FlyBase IDs as the key and a list of FASTA record indices as the value
    :param output:  Path to the output file default: <input_file>.fasta
    :return: None
    """
    with output.open("w") as fh:
        for fbid in ids:
            for i in index[fbid]:
                fh.write(f">{fasta[i][0]}\n{fasta[i][1]}\n")
    return None


def read_id_file(id_file: Path) -> FlyBaseIdList:
    """
    Reads a file containing user supplied FlyBase IDs and returns them as a frozen set.
    :param id_file: The ID file to read in.
    :return: Frozen set of FlyBase IDs
    """
    try:
        # Read in the IDs, ignoring non FlyBase IDs
        with id_file.open() as fh:
            return frozenset(
                [line.strip() for line in fh if FLYBASE_ID_REGEX.match(line)]
            )
    except FileNotFoundError:
        print(f"{id_file} not found", file=sys.stderr)
        return frozenset()
    except IOError:
        print(f"IO Error while reading {id_file}", file=sys.stderr)
        return frozenset()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extracts FlyBase FASTA records based on ID or parent IDs.",
        epilog="e.g.: python3 fasta/flybase_id_to_fasta.py --fasta dmel-all-CDS.fasta.gz ids1.txt ids2.txt",
    )
    parser.add_argument(
        "idfile",
        type=Path,
        nargs="+",
        help="File(s) containing FlyBase IDs to extract.",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="The FASTA file to extract records from.",
    )
    args = parser.parse_args()

    # Read fasta file into memory.
    # Note: this will be slow or not work for large files.
    # This currently works for all FlyBase gzipped FASTA files.
    fasta_records, fasta_idx = read_fasta(args.fasta)

    # Loop over all user ID lists.
    for file in args.idfile:
        # Read in the user supplied FlyBase IDs.
        ids = read_id_file(file)
        # Set up outout file.
        output = Path(file.stem + ".fasta")
        if ids:
            # Write requested FASTA records to the output file.
            fbid_to_fasta(ids, fasta_records, fasta_idx, output=output)
