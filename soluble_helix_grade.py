#!/usr/bin/env python2.7
"""
a script bundle to score soluble dataset helical stretches
"""
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde

from WinGrade import grade_segment

import TMConstraint
import ProcessEntry as pe



polyval = pe.MakeHydrophobicityGrade()
work_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/soluble_helix_grading/'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', default='run_job')
    parser.add_argument('-name')
    args = vars(parser.parse_args())

    args['min_length'] = 21
    args['inc_max'] = 10
    args['with_msa'] = False
    args['with_csts'] = False
    args['w'] = 0
    args['z_0'] = 0

    if args['mode'] == 'run_job':
        run_job(args)

    elif args['mode'] == 'collect':
        dGs, ddGs = collect_win_analysis(args)
        for d, dd in zip(dGs, ddGs):
            print d, dd

    elif args['mode'] == 'draw':
        draw_phase(args)

    elif args['mode'] == 'min_dG_hists':
        min_dG_hists(args)

    elif args['mode'] == 'test':
        parse_all_seqs()

    else:
        print 'no mode found'


def min_dG_hists(args):
    """
    :param args:
    :return:
    """
    sol_dGs, sol_ddGs = collect_win_analysis_minimums(args)
    tm_dGs, tm_ddGs = collect_tm_mins(args)

    print 'soluble dGs'
    print '\n'.join([str(a) for a in sol_dGs])
    print 'TM dGs'
    print '\n'.join([str(a) for a in tm_dGs])

    plt.subplot(211)
    plt.hist(sol_dGs, normed=1, bins=100, alpha=0.4, color='r')
    plt.hist(tm_dGs, normed=1, bins=100, alpha=0.4, color='k')

    plt.subplot(212)
    plt.scatter(sol_dGs, sol_ddGs, alpha=0.4, color='r')
    plt.scatter(tm_dGs, tm_ddGs, alpha=0.4, color='k')
    plt.show()


def collect_tm_mins(args):
    from WinGrade import flip_win_grade
    from positive_inside_analysis import parse_prd
    tms_path = '/home/labs/fleishman/elazara/benchmark_paper_new/TOPCONS_dataset/Transmembrane/run_folder/'
    prd_files = [a for a in os.listdir(tms_path) if '.prd' in a and '_msa' not in a and '.txt' not in a]
    full_seqs_dict = parse_all_seqs()

    dGs, ddGs = [], []
    for prd_file in prd_files:
        wgp = parse_prd(tms_path+prd_file)[0]
        if wgp is None:
            continue
        min_w = [w_ for w_ in wgp.path if w_.grade == min([w.grade for w in wgp.path])][0]
        dGs.append(min_w.grade)
        flipped_min_w = flip_win_grade(min_w, full_seqs_dict[prd_file.split('.')[0]])
        ddGs.append(min_w.grade - flipped_min_w.grade)

    return dGs, ddGs


def parse_all_seqs():
    """
    :return: {name: seq}
    """
    all_seqs_path = '/home/labs/fleishman/elazara/benchmark_paper_new/TOPCONS_dataset/Transmembrane/run_folder/All.fasta'
    with open(all_seqs_path, 'r') as fin:
        cont = fin.read().split('>')
    results = {}
    for p in cont:
        s = p.split()
        if len(s) == 2:
            results[s[0]] = s[1]
    return results

def draw_phase(args):
    dGs, ddGs = collect_win_analysis(args)
    values = np.vstack([dGs, ddGs])
    density = gaussian_kde(values)
    z = density([dGs, ddGs])
    plt.scatter(dGs, ddGs, z)
    plt.show()


def collect_win_analysis(args):
    """
    :param args: run argumerns
    :return: dGs and ddGs lists for all soluble proteins.
    """
    wins_files = [a for a in os.listdir(work_path) if '.wns' in a]
    dGs, ddGs = [], []

    for wns_file in wins_files:
        with open(work_path+wns_file, 'r') as fin:
            cont = fin.read().split('\n')
        for l in cont:
            s = l.split()
            if len(s) == 2:
                dGs.append(float(s[0]))
                ddGs.append(float(s[1]))
    return dGs, ddGs


def collect_win_analysis_minimums(args):
    """
    :param args: run argumerns
    :return: min dG for each soluble protein
    """
    wins_files = [a for a in os.listdir(work_path) if '.wns' in a]
    dGs = []
    ddGs = []

    for wns_file in wins_files:
        with open(work_path+wns_file, 'r') as fin:
            cont = fin.read().split('\n')
        res_dG = []
        res_ddG = []
        for l in cont:
            s = l.split()
            if len(s) == 2:
                res_dG.append(float(s[0]))
                res_ddG.append(float(s[1]))
        dGs.append(min(res_dG))
        ddGs.append(res_ddG[np.argmin(res_dG)])
    return dGs, ddGs


def run_job(args):
    seq = read_seq(work_path+args['name']+'.fasta')
    print args['name'], seq
    entry_cst = TMConstraint.TMConstraint(args['name'])
    topo_entry = pe.TopoEntry(args['name'], seq=seq, end_of_SP=0,
                              psipred=pe.parse_psipred(work_path+args['name']+'.ss2'), hydro_polyval=polyval,
                              param_list=args, csts=entry_cst, path_msa=None, original_seq=seq)
    wins = pe.win_grade_generator(topo_entry, mode='all', tm_pos_mode='selective')
    print 'wins', wins
    if wins == []:
        sys.exit()
    else:
        with open(work_path+args['name']+'.wns', 'w+') as fout:
            for w in wins:
                fout.write('%.1f %.1f\n' % (w.grade, w.grade-grade_segment(w.seq[::-1], polyval)))
    # print topo_entry


def read_seq(file_path):
    with open(file_path, 'r') as fin:
            cont = fin.read().split('\n')
    seq = ''
    for i, l in enumerate(cont):
        if len(l) < 1:
            continue
        if l[0] == '>':
           seq = cont[i+1]
    return seq


if __name__ == '__main__':
    main()
