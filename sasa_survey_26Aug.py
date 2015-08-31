#!/usr/bin/env python
"""
a script to make a SASA Vs. hydrophobicity plot.
"""


def download_pdbs():
    """
    :return: downloads all PDBs for the rostlab database. only the ones actually available...
             print the names of those it failed
    """
    from Bio.PDB import PDBParser, PDBIO, PDBList
    from TMpredict_WinGrade import parse_rostlab_db
    rost_db = parse_rostlab_db()
    pdbl = PDBList()
    failed = []
    for k, v in rost_db.items():
        print k, v
        try:
            pdbl.retrieve_pdb_file(v['pdb'], pdir='PDB')
        except:
            failed.append(v['pdb'])
    print failed

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode')
    args = vars(parser.parse_args())
    if args['mode'] == 'download':
        # use to download the pdbs into ./PDB/
        download_pdbs()

if __name__ == '__main__':
    main()