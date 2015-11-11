#!/usr/bin/env python2.7
# coding=utf-8
"""
a script to make a SASA Vs. hydrophobicity plot.
"""
three_2_one = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
               'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
               'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y', 'UNK': 'U'}


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


def parse_rsa(rsa_file):
    with open(rsa_file, 'r') as fin:
        cont = fin.read().split('\n')
    result = {}
    for l in cont:
        s = l.split()
        if s == []:
            continue
        if s[0] != 'RES':
            continue

        if len(s) < 14:
            t = s[2]
            s[2] = s[2][0]
            s.insert(3, t[1:])
            if s[2].lower() not in result.keys():
                result[s[2].lower()] = {}
        if s[2].lower() not in result.keys():
            result[s[2].lower()] = {}
        result[s[2].lower()][int(s[3])] = {'type': three_2_one[s[1]], 'all_atoms_ABS': float(s[4]), 'all_atoms_REL': float(s[5]),
                             'total_side_ABS': float(s[6]), 'total_side_REL': float(s[7]), 'main_chain_ABS': float(s[8]),
                             'main_chain_REL': float(s[9]), 'non_polar_ABS': float(s[10]), 'non_polar_REL': float(s[11]),
                             'all_polar_ABS': float(s[12]), 'all_polar_REL': float(s[13]), 'chain': s[2]}
        if s[0] == 'END':
            break
    return result


def analyse():
    import pickle
    import os
    import sys
    import operator
    import random
    import matplotlib.pyplot as plt
    from TMpredict_WinGrade import parse_rostlab_db
    from WinGrade import topo_string_to_WGP
    from topo_strings_comparer import prd_parser, spc_parser
    total_sasa_dict = parse_standard_data()
    rost_db = parse_rostlab_db()
    neighbors = 0
    single = 0
    without_neighbours = 0
    accessible_vs_ddg = []
    for k, v in rost_db.items():
        if v['name'] in ['p01730', 'p19054', 'e1c9k9', 'q9qug3', 'p07471']:
            continue
        is_single = os.path.isfile('single_chains/%s_1.pdb' % v['pdb'])
        is_neighbor = os.path.isfile('with_neighbours/%s_%s_with_neighbors.pdb' % (v['pdb'], v['chain'].upper()))
        is_no_neighbour = os.path.isfile('without_neighbours/%s_%s_without_neighbours.pdb' % (v['pdb'], v['chain'].upper()))
        neighbors += 1 if is_neighbor else 0
        without_neighbours += 1 if is_no_neighbour else 0
        single += 1 if is_single else 0

        # before cutting out the spring:
        # prediction = prd_parser('/home/labs/fleishman/elazara/benchmark_paper_new/Mean/Plain', v['name']+'.prd')
        # after cutting out spring, with MSA:
        # prediction = prd_parser('/home/labs/fleishman/elazara/length_21/w_0_with_MSA/', v['name']+'.prd')
        # after cutting out spring, without MSA:
        prediction = prd_parser('/home/labs/fleishman/elazara/length_21/', v['name']+'.prd')
        wgp_pred = topo_string_to_WGP(prediction['best_path_ts'], v['seq'])

        spoc = spc_parser(v['name'])['spoctopus']
        signal = [0, spoc.count('s') + spoc.count('S')]

        if is_single:
            naccess = parse_rsa('single_chains/%s_1.rsa' % v['pdb'])
        else:
            naccess = parse_rsa('with_neighbours/%s_%s_with_neighbors.rsa' % (v['pdb'], v['chain'].upper()))

        wgp_pdbtm = topo_string_to_WGP(v['pdbtm'], v['seq'])
        rost_aln, naccess_aln, score, beg, end = \
            pair_wise_aln_from_seqs(v['seq'], ''.join([a['type'] for a in naccess[v['chain']].values()]))
        for w in wgp_pdbtm.path:
            if w.begin <= signal[1]:
                print 'signes', w.begin, w.end, signal
                continue
            naccess_win = nacces_for_win(naccess, w, naccess_aln, rost_aln, total_sasa_dict, v['chain'])
            predicted = observed_found_in_prediction(w, wgp_pred)
            accessible_vs_ddg.append({'access': naccess_win, 'predicted': predicted, 'grade': w.grade})
    with open('pickled.obj', 'wb') as pkl:
        pickle.dump(accessible_vs_ddg, pkl)
    print 'its pickled'
        # if is_single:
        #     naccess = parse_rsa('single_chains/%s_1.rsa' % v['pdb'])
        #     wgp = topo_string_to_WGP(v['pdbtm'], v['seq'])
        #     rost_aln, naccess_aln, score, beg, end = \
        #         pair_wise_aln_from_seqs(v['seq'], ''.join([a['type'] for a in naccess[v['chain']].values()]))
        #     for w in wgp.path:
        #         naccess_win = nacces_for_win(naccess, w, naccess_aln, rost_aln, total_sasa_dict, v['chain'])
        #         accessible_vs_ddg.append([naccess_win, w.grade])

        # if is_neighbor:
        #     naccess_with = parse_rsa('with_neighbours/%s_%s_with_neighbors.pdb' % (v['pdb'], v['chain'].upper()))
        #     naccess_without = parse_rsa('without_neighbours/%s_%s_without_neighbours.pdb' % (v['pdb'], v['chain'].upper()))
        #     wgp = topo_string_to_WGP(v['pdbtm'], v['seq'])
        #     rost_aln_with, naccess_aln_with, score_with, beg_with, end_with = \
        #         pair_wise_aln_from_seqs(v['seq'], ''.join([a['type'] for a in naccess_with.values()]))
        #     rost_aln_without, naccess_aln_without, score_without, beg_without, end_without = \
        #         pair_wise_aln_from_seqs(v['seq'], ''.join([a['type'] for a in naccess_without.values()]))
        #     for w in wgp.path:
        #         naccess_win_with = nacces_for_win(naccess_with, w, naccess_aln_with, rost_aln_with, total_sasa_dict,
        #                                           v['chain'])
        #         naccess_win_without = nacces_for_win(naccess_without, w, naccess_aln_without, rost_aln_without,
        #                                              total_sasa_dict, v['chain'])
    # print [a[0] for a in accessible_vs_ddg if a[0] is not None]
    # plt.hist([a[0] for a in accessible_vs_ddg if a[0] is not None])
    # ranges = [{'ranges': [a, a+20], 'vals': []} for a in range(0, 60, 20)]


def unpickle_and_plot():
    import pickle
    import matplotlib.pyplot as plt
    import matplotlib
    import random

    matplotlib.rc_context(fname="/home/labs/fleishman/jonathaw/.matplotlib/publishable_matplotlibrc")
    with open('pickled.obj', 'rb') as fin:
        accessible_vs_ddg = pickle.load(fin)

    # Create a figure instance
    fig = plt.figure(1, figsize=(6, 6))


    # Create an axes instance
    ax = fig.add_subplot(111)



    ranges = [{'ranges': [0, 20], 'vals': [], 'pred': []}, {'ranges': [20, 200], 'vals': [], 'pred': []}]
    # ranges.append({'ranges': [60, 100], 'vals': []})
    print 'found a total of %i points' % len(accessible_vs_ddg)
    # plt.subplot(2, 1, 1)
    for a in accessible_vs_ddg:
        for v in ranges:
            if v['ranges'][0] <= a['access'] < v['ranges'][1]:
                v['vals'].append(a['grade'])
                v['pred'].append(a['predicted'])
    # for v in ranges:
    #     print v



    for k in ranges:
        print k['ranges']
        print 'True', [a for i, a in enumerate(k['vals']) if k['pred'][i]]
        print 'False', [a for i, a in enumerate(k['vals']) if not k['pred'][i]]

    bp = ax.boxplot([a['vals'] for a in ranges], positions=[0, 0.5], patch_artist=True, labels=['<20%', '>20%'])
    for box in bp['boxes']:
        # change outline color
        # box.set( color='#7570b3', linewidth=2)
        box.set(color='black', linewidth=2)
        # change fill color
        # box.set( facecolor = '#1b9e77', alpha=0.8 )
        box.set(facecolor=[0,0,0,0])

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        # whisker.set(color='#7570b3', linewidth=2)
        whisker.set(color='black', linewidth=2)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        # cap.set(color='#7570b3', linewidth=2)
        cap.set(color='black', linewidth=3)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        # median.set(color='#b2df8a', linewidth=2)
        median.set(color='red', linewidth=2)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='_', color='#e7298a', markersize=10)

    for i, rng in zip([0, 0.5], ranges):
        x = [i+random.uniform(-0.1, 0.1) for j in range(len(rng['vals']))]
        plt.scatter(x, rng['vals'], alpha=0.5,
                    c=['b' if a else 'r' for a in rng['pred']])

    plt.ylim([-15, 15])
    plt.xlabel('%SASA')
    plt.ylabel(r'$\Delta$G') # [kcal/mol]')
    # plt.subplot(2, 1, 2)
    # scatt = [a for a in accessible_vs_ddg if a['access'] < 120]
    # plt.scatter([a['access'] for a in scatt], [a['grade'] for a in scatt],
    #             c=['b' if a['predicted'] else 'r' for a in scatt])
    plt.show()
    # print 'neighbours', neighbors
    # print 'without_neighbours', without_neighbours
    # print 'single', single


def observed_found_in_prediction(obse_win, predicted_wgp):
    for w in predicted_wgp.path:
        if wingarde_overlap_M(w, obse_win):
            return True
    return False


def wingarde_overlap_M(w1, w2, M=10):
    return abs(w1.begin-w2.begin) <= M and abs(w1.end-w2.end) <= M


def nacces_for_win(naccess, win, naccess_aln, rost_aln, tot_sasa, chain):
    """
    :param naccess: naccess results dictionary for ebntry
    :param win: a WinGrade
    :param naccess_aln: naccess sequnce alignment
    :param rost_aln: rost sequence alignment
    :param tot_sasa: dictionary of total sasa for every residue type
    :return: percent sasa exposed
    """
    typer = 'non_polar_ABS'
    aln_loc = [win.begin+rost_aln[:win.begin].count('-'), win.end+rost_aln[:win.end].count('-')]
    if aln_loc[1] >= max(naccess[chain].keys()):
        return None
    total_sasa = 0.0
    accessible_sasa = 0.
    length = 0.
    for i in range(aln_loc[0], aln_loc[1]+1):
        if i+1 not in naccess[chain].keys():
            continue
        if rost_aln[i] != naccess[chain][i+1]['type']:
            continue
        total_sasa += tot_sasa[rost_aln[i]][typer]
        accessible_sasa += naccess[chain][i+1][typer]
        length += 1.
    if total_sasa == 0:
        return 0
    return 100*accessible_sasa/total_sasa


def pair_wise_aln_from_seqs(seq1, seq2, matrix=None, gap_open=-10, gap_extend=-0.5):
    """
    :param seq1: a fasta sequence
    :param seq2: a fasta sequence
    :param matrix: a substitution matrix
    :param gap_open: gap_open penalty
    :param gap_extend: gap_extend penalty
    :return: aln for seq1, aln for seq2, score, beginning and end for the best alignment
    """
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    if matrix is None:
        matrix = matlist.blosum62
    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)[0]
    return alns


def parse_standard_data():
    with open('/apps/RH6U4/naccess/naccess2.1.1/standard.data', 'r') as fin:
        cont = fin.read().split('\n')
    result = {}
    for l in cont:
        s = l.split()
        if not s:
            continue
        if s[0] != 'ATOM':
            continue
        result[three_2_one[s[3]]] = {'all_atoms_ABS': float(s[4]), 'all_atoms_REL': float(s[5]),
                             'total_side_ABS': float(s[6]), 'total_side_REL': float(s[7]), 'main_chain_ABS': float(s[8]),
                             'main_chain_REL': float(s[9]), 'non_polar_ABS': float(s[10]), 'non_polar_REL': float(s[11]),
                             'all_polar_ABS': float(s[12]), 'all_polar_REL': float(s[13])}
    return result


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode')
    args = vars(parser.parse_args())
    if args['mode'] == 'download':
        # use to download the pdbs into ./PDB/
        download_pdbs()

    elif args['mode'] == 'test':
        a = parse_rsa('with_neighbours/3oe0_A_with_neighbors.rsa')
        for chian, rsa in a.items():
            for k, v in rsa.items():
                print chian, k, v['all_atoms_ABS']

    elif args['mode'] == 'analyse':
        analyse()

    elif args['mode'] == 'unpickle':
        unpickle_and_plot()

    else:
        print 'no mode was found'


if __name__ == '__main__':
    main()