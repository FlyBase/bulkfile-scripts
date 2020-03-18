#!/usr/bin/env python3
import sys
import os

"""
This script updates a list of FlyBase FBgn IDs to their current IDs
using the FBgn <=> Annotation ID (fbgn_annotation_ID_*.tsv) file.
 
See 

https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#FBgn_.3C.3D.3E_Annotation_ID_.28fbgn_annotation_ID_.2A.tsv.29

Usage:

./fbgn_updater.py your_fbgn_list.txt fbgn_annotation_ID.tsv > output.tsv

Assumptions:
* Accepts only FlyBase gene IDs (FBgn).

Author: Josh Goodman <jogoodma@iu.edu>
"""



def insert_fbid(primary: str = None, secondary: str = None, fbid_dict: dict = {}):
    """
    Modifies the dictionary in place by inserting the secondary FlyBase ID as the key
    and the primary FlyBase ID as the value.  If the secondary ID is already present in the
    dictionary the primary ID is added to the unique set of FlyBase IDs in the value.

    :param primary:str - A single FlyBase ID to insert into the dictionary.
    :param secondary:list -  A list of secondary IDs.
    :param fbid_dict:dict - The dictionary reference to modify.
    :return: None
    """
    if secondary and secondary not in fbid_dict:
        # If we haven't seen this fbid before initialize the set.
        fbid_dict[secondary] = {primary}
    elif secondary:
        # We have seen this fbid before so we add it to the set.
        fbid_dict[secondary].add(primary)
    return None


def generate_inverted_fbid_dict(fbid_file: str):
    """
    Generates an inverted dictionary of all secondary IDs as keys 
    and a set of primary FBids as values.

    :param fbid_file: str - The FlyBase FBgn <=> Annotation ID (fbgn_annotation_ID_*.tsv) file to parse.
    :return: The inverted FlyBase id dictionary and a set of primary FlyBase IDs
    """
    # Init the dictionary and set of current FlyBase ids.
    fbid_dict = {}
    current_fbids = set()

    # Open file and loop over lines.
    with open(fbid_file, "r") as file:
        for line in file:
            line = line.strip()

            # This script only cares about genes for now.
            if not line.startswith('#') and 'FBgn' in line:
                # Split out the ID column and all the others.
                symbol, species, primary_fbid, secondary_fbid_col, *rest = line.split('\t')
                current_fbids.add(primary_fbid)
                secondary_fbid_list = secondary_fbid_col.split(',')
                [insert_fbid(primary_fbid, fbid, fbid_dict) for fbid in secondary_fbid_list]

    return fbid_dict, current_fbids

def main(user_ids, fbid_dict, current_ids):
    # Open the user ID file and loop over it.
    with open(user_ids, 'r') as file:
        for fbid in file:
            fbid = fbid.strip()
            try:
                # If the ID is current, print a line with no additional IDs.
                if fbid in current_ids: 
                    print(f"{fbid}")
                # If the ID is not current, print all possible current IDs.
                else:
                    # Fetch the ID set for the ID in their list.
                    ids = '\t'.join(fbid_dict[fbid])
                    # Print out results.
                    print(f"{fbid}\t{ids}")
            except KeyError:
                # Handle cases when their ID doesn't exist.
                print(f"{fbid}\tNone")
    return None;


if __name__ == '__main__':
    try:
        # Read in arguments.
        user_ids, fbgn_db_acc = sys.argv[1:3]
        # Generate the inverted dictionary.
        inverted_fbid_dict, current_fbids = generate_inverted_fbid_dict(fbgn_db_acc)
        print(f'# FlyBase ID Updater\n# Input = {os.path.abspath(user_ids)}\n# ID Reference = {os.path.abspath(fbgn_db_acc)}')
        print('# Submitted_ID\tUpdated_ID(s)')
        main(user_ids, inverted_fbid_dict, current_fbids)

    except ValueError:
        print(f"Usage: {os.path.basename(__file__)} your_FBgn_ids.txt fbgn_annotation_ID.tsv", file=sys.stderr)

