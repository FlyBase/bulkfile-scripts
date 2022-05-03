#!/usr/bin/env python3
import sys
import os
import re

"""
This script accepts a list of gene or transcript symbol/synonyms and converts them into one or more 
FlyBase IDs.

The synonyms file can be download from:

ftp://ftp.flybase.org/releases/current/precomputed_files/synonyms/

Usage:

./symbol_to_id_lookup.py your_symbol_list.txt flybase_synonym_file.tsv > output.tsv

Assumptions:
* Only gene or transcript symbols/synonyms/names.
* Drosophila melanogaster only

Author: Josh Goodman <jogoodma@iu.edu>
"""


def insert_symbol(symbol: str, fbid: str, dict: dict):
    """
    Modifies the dictionary in place by inserting the symbol as the key
    and the fbid as the value.  If the symbol is already present in the
    dictionary the fbid is added to the unique set of FlyBase IDs in the value

    :param symbol:str - A single symbol to insert into the dictionary.
    :param fbid:str -  A single FlyBase ID.
    :param dict:dict - The dictionary reference to modify.
    :return: None
    """
    if symbol and symbol not in dict:
        # If we haven't seen this symbol before initialize the set.
        dict[symbol] = {fbid}
    elif symbol:
        # We have seen this symbol before so we add it to the set.
        dict[symbol].add(fbid)
    return None


def generate_inverted_symbol_dict(sym_file: str):
    """
    Generates an inverted dictionary of all symbols, synonyms, names, etc.
    as keys and a set of FBids as values.

    :param sym_file: str - The FlyBase synonyms file to parse.
    :return: The inverted symbol/synonym dictionary.
    """
    """
     Regex to split name synonyms on commas without spaces.
     Commas without a trailing space indicate a new name.
     Commas with a trailing space indicate a name with a comma in it.
     e.g.
     my gene1,my gene2 -> ['my gene1', 'my gene2']
     my gene1, my gene2 -> ['my gene1, my gene2']
    """
    # Match commas that are not followed by a space.
    comma_ns_re = re.compile(r",(?!\s)")

    # Init the dictionary.
    symbol_dict = {}

    # Open file and loop over lines.
    with open(sym_file, "r") as file:
        for line in file:
            line = line.strip()

            # This script only cares about genes or transcripts ignore the rest.
            if line.startswith("FBgn") or line.startswith("FBtr"):
                # Split out the ID column and all the others.
                fbid, *cols = line.split("\t")
                try:
                    col_len = len(cols)
                    # Dmel only.
                    if cols[0] != "Dmel":
                        continue

                    # Symbol
                    insert_symbol(cols[1], fbid, symbol_dict)

                    # Fullname
                    if col_len >= 3 and cols[2]:
                        insert_symbol(cols[2], fbid, symbol_dict)
                    # Fullname synonyms
                    if col_len >= 4 and cols[3]:
                        [
                            insert_symbol(syn, fbid, symbol_dict)
                            for syn in comma_ns_re.split(cols[3])
                        ]
                    # Symbol synonyms
                    if col_len >= 5 and cols[4]:
                        [
                            insert_symbol(syn, fbid, symbol_dict)
                            for syn in comma_ns_re.split(cols[4])
                        ]
                except IndexError:
                    print(f"Formatting problem found in line:\n{line}", file=sys.stderr)
                    continue
    return symbol_dict


if __name__ == "__main__":
    try:
        # Read in arguments.
        symbols_to_check, fb_synonym = sys.argv[1:3]
        # Generate the inverted dictionary.
        inverted_symbol_dict = generate_inverted_symbol_dict(fb_synonym)

        # Open their symbol file and loop over it.
        with open(symbols_to_check, "r") as file:
            for symbol in file:
                symbol = symbol.strip()
                try:
                    # Fetch the ID set for the symbol in their list.
                    ids = inverted_symbol_dict[symbol]
                    # Print out results.
                    print(f"{symbol}\t{','.join(ids)}")
                except KeyError:
                    # Symbol doesn't exist in dictionary.
                    print(f"{symbol}")
    except ValueError:
        print(
            f"Usage: {os.path.basename(__file__)} your_symbols.txt fb_synonym.tsv",
            file=sys.stderr,
        )
