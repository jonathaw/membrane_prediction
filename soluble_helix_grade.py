#!/usr/bin/env python2.7
"""
a script bundle to score soluble dataset helical stretches
"""
import os
import sys
import argparse

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
        collect_win_analysis(args)

    else:
        print 'no mode found'


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
