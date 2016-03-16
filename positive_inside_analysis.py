#!/usr/bin/env python2.7
# coding=utf-8
import os
import argparse
from scipy import stats
from TMpredict_WinGrade import parse_rostlab_db
from WinGrade import parse_WGP, grade_segment, WinGradePath, WinGrade, flip_win_grade
from TMConstraint import TMConstraint, pred2cst
from sasa_survey_26Aug import parse_rsa, nacces_for_win, pair_wise_aln_from_seqs, parse_standard_data
from topo_strings_comparer import spc_parser
from TMpredict_WinGrade import topo_VH_parser
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FuncFormatter
from ProcessEntry import MakeHydrophobicityGrade
import pickle
import pandas as pd
from scipy.stats import gaussian_kde, ttest_1samp, linregress
import numpy as np
from collections import OrderedDict
import sys
import shutil
import pandas as pd

from soluble_helix_grade import collect_win_analysis


matplotlib.rc_context(fname="/home/labs/fleishman/jonathaw/.matplotlib/publishable_matplotlibrc")
# matplotlib.rc_context(fname="/Volumes/labs/fleishman/jonathaw/.matplotlib/publishable_matplotlibrc")
polyval = MakeHydrophobicityGrade()


def to_percent(y, position):
    if y == 0:
        return ''
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(int(100 * round(y, 2)))

    # The percent symbol needs escaping in latex

    if matplotlib.rcParams['text.usetex'] is True:
        return s + r'$\%$'
    else:
        return s + '%'


def analyse_sasa_and_positive(write=False, with_sasa=True, vh=False):
    # for the Rost DB:
    # path_to_nocsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_12Nov/with_msa_no_csts/'
    # path_to_nocsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_12Nov/minus_2/with_msa_no_csts/'
    # for Rost DB with all pos and 0 fidelity
    # path_to_nocsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_12Nov/all_pos_no_fidelity_7Dec/with_msa_no_csts/'
    # path_to_nocsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_29.2/with_msa_no_csts/'
    path_to_nocsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_10.3/with_msa_no_csts/' # msa rost 10.3
    # for the entire VH DB:
    # path_to_nocsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_12Nov/vh_positives/topcons_no_csts/'
    # path_to_nocsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_14.12/with_msa_no_csts/'
    # path_to_nocsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_14.12/with_msa_no_csts/'
    prd_files = [a for a in os.listdir(path_to_nocsts)
                 if a[-4:] == '.prd' and '_msa' not in a]
    rost_db = parse_rostlab_db()
    results = {}
    three_n_over = 0
    total, with_pos = 0., 0.
    for prd_file in prd_files:
        if prd_file.lower() == 'secg.prd':
            continue
        if vh:
            vh_db = topo_VH_parser(prd_file[:-4])
        else:
            rost_prd = rost_db[prd_file[:-4]]
        wgp = parse_prd(path_to_nocsts+prd_file)[0]

        if wgp is None:
            print 'skipping', prd_file, wgp
            continue
        if wgp.total_grade == 0.:
            print 'no path FOUND, F#$% it', prd_file
            continue
        if wgp.win_num >= 3:
            total += 1 if wgp.win_num >= 6 else 0
            three_n_over += 1
            pos_in = [w for w in wgp.path[1:-1] if w.grade-w.inner_tail_score >= -3.0]
            if pos_in:
                with_pos += 1
                # for doing the max win
                # chosen_w = [w for w in pos_in if w.grade == max([a.grade for a in pos_in])][0]
                results[prd_file[:-4]] = {}
                for i, chosen_w in enumerate(pos_in):
                    if with_sasa and not vh:
                        sasas = win_sasa([chosen_w], rost_prd['pdb'], rost_prd['chain'], prd_file[:-4], rost_prd['seq'])
                    else:
                        sasas = [None]
                    if not sasas and not vh:
                        continue
                    sasas = sasas[0]
                    results[prd_file[:-4]][i] = {'wgp': wgp, 'pos_in': pos_in, 'chosen_w': chosen_w, 'sasa': sasas}
                    if write:
                        tm_poses = [[w.begin, w.end, None] for w in wgp.path if w != chosen_w]
                        cst = TMConstraint(name=prd_file.split('.')[0], mode='only', tm_pos=tm_poses, tm_pos_fidelity=0)
                        print cst
                        with open('remove_w_csts/%s_%i.cst' % (prd_file[:-4], i), 'w+') as fout:
                                fout.write(str(cst))
                    # print prd_file
                    # if write:
                    #     try:
                    #         print "took %f with sasa %r at %i-%i" % (chosen_w.grade-chosen_w.inner_tail_score, sasas,
                    #                                                  chosen_w.begin, chosen_w.end)
                    #         cst = parse_prd_csts(prd_file)
                    #         for tm_pos in cst.tm_pos:
                    #             if tm_pos[0]-cst.tm_pos_fidelity <= chosen_w.begin and chosen_w.end <= tm_pos[1]+cst.tm_pos_fidelity:
                    #                 cst.tm_pos.remove(tm_pos)
                    #                 cst.non_tm_pos = [[chosen_w.begin, chosen_w.end]]
                    #                 cst.tm_pos_fidelity = 0
                            # with open('remove_w_csts/%s_%i.cst' % (prd_file[:-4], i), 'w+') as fout:
                            #     fout.write(str(cst))
                        # except:
                        #     pass
    if write and not vh:
        print "found %i with positive in, out of %i with 3 or over" % (len(results.keys()), three_n_over)
        for k, v in results.items():
            print 'AA', k, v
            print 'for %s there are %i positive insides\t %s\t %r' % (k, len(v['pos_in']), rost_db[k]['pdb'], v['pos_in'])
            print 'BB'
    print 'found a total of %i with 3 or more helices' % total
    print 'found %i with a positive helix' % with_pos
    print 'that is %f percent' % (100*with_pos/total)
    return results


def win_sasa(wins, pdb, chain, name, seq):
    total_sasa_dict = parse_standard_data()
    path_d = '/home/labs/fleishman/jonathaw/membrane_prediction_DBs/sasa_survey_26Aug/all_pdbs/'
    is_single = os.path.isfile('%ssingle_chains/%s_1.pdb' % (path_d, pdb))
    is_neighbor = os.path.isfile('%swith_neighbours/%s_%s_with_neighbors.pdb' % (path_d, pdb, chain.upper()))
    is_no_neighbour = os.path.isfile('%swithout_neighbours/%s_%s_without_neighbours.pdb' % (path_d, pdb, chain.upper()))

    spoc = spc_parser(name)['topcons']
    signal = [0, spoc.count('s') + spoc.count('S')]

    if is_single:
        naccess = parse_rsa('%ssingle_chains/%s_1.rsa' % (path_d, pdb))
    else:
        naccess = parse_rsa('%swith_neighbours/%s_%s_with_neighbors.rsa' % (path_d, pdb, chain.upper()))

    rost_aln, naccess_aln, score, beg, end = \
        pair_wise_aln_from_seqs(seq, ''.join([a['type'] for a in naccess[chain].values()]))
    naccess_wins = []
    for w in wins:
        if w.begin <= signal[1]:
            continue
        naccess_wins.append(nacces_for_win(naccess, w, naccess_aln, rost_aln, total_sasa_dict, chain))
    return naccess_wins


def analyse_results():
    notm_prds = [a for a in os.listdir('./') if a[-9:] == '_notm.prd']
    for notm in notm_prds:
        wgp_notm = parse_prd(notm)
        wgp_withtm = parse_prd('prds_with_positive/'+''.join(notm.split('_notm')))
        print 'found %s this with %f, without %f, minus %f' % (notm, wgp_withtm.total_grade, wgp_notm.total_grade,
                                                               wgp_withtm.total_grade-wgp_notm.total_grade)


def parse_prd(file_name):
    text = open(file_name, 'r').read()
    one_liner = text.split('best_path ')[1].split('sec_path')[0].replace('\n', ':')
    if '{ [ }~> total_grade' in one_liner:
        return None, None
    best_wgp = parse_WGP(one_liner)

    one_liner = text.split('sec_path ')[1].split('ddG paths')[0].replace('\n', ':')
    if '{ [ }~> total_grade' in one_liner:
        return None, None
    sec_wgp = parse_WGP(one_liner)
    return best_wgp, sec_wgp


def parse_prd_csts(file_name):
    tm_pos = []
    non_tm_pos = []
    for l in open(file_name, 'r').read().split('\n'):
        s = l.split()
        if not s:
            continue
        if s[0] == 'name':
            name = s[1]
        if s[0] == 'tm_num':
            if s[1] != 'None':
                tm_num = s[1]
            else:
                tm_num = None
        if s[0] == 'tm_pos':
            if s[1] == 'None':
                tm_pos = None
            else:
                tm_pos.append([int(s[1]), int(s[2]), s[3] if s[3] != 'None' else None])
        if s[0] == 'non_tm_pos':
            if s[1] == 'None':
                non_tm_pos = None
            else:
                non_tm_pos.append([int(s[1]), int(s[2])])
        if s[0] == 'tm_pos_fidelity':
            tm_pos_fidelity = int(s[1])
        if s[0] == 'c_term':
            if s[1] != 'None':
                c_term = s[1]
            else:
                c_term = None
        if s[0] == 'n_term':
            if s[1] != 'None':
                n_term = s[1]
            else:
                n_term = None
        if s[0] == 'mode':
            mode = s[1]
    return TMConstraint(name, mode, tm_num, tm_pos, tm_pos_fidelity, c_term, n_term, non_tm_pos)


def from_csts_to_prds():
    """
    :return: after creating the win remove csts. use this to create the folders and jobs to make the with_csts
    prediction
    """
    # work_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_29.2/with_msa_with_csts/'
    work_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_10.3/with_msa_with_csts/' # msa rost 10.3
    csts_files = [a for a in os.listdir(work_path) if '.cst' in a]
    for cst_file in csts_files:
        name = cst_file[:-4]
        os.mkdir(work_path+name)
        shutil.move(cst_file, work_path+name+'/'+name[:-2]+'.cst')
        os.chdir(work_path+name)

        with open('job.%s' % name, 'w+') as fout:
            fout.write('#!/bin/bash\n')
            fout.write('. /usr/share/lsf/conf/profile.lsf\n')
            fout.write('cd %s\n' % (work_path+name))
            fout.write('/apps/RH6U4/python/2.7.6/bin/python '
                       '/home/labs/fleishman/jonathaw/membrane_prediciton/TMpredict_WinGrade.py  '
                       '-name  %s  '
                       '-mode  new  '
                       '-run_type  user_cst  '
                       '-in_path ./ -out_path ./ -db  rost  '
                       '-create_html  False  '
                       '-fidelity  0\n' % name[:-2])

        os.chdir(work_path)


def analyse_3d_results(args):
    """
    analyse the results with all helices (no csts), and with removing the positive (with csts), and that window's SASA
    :return:
    """
    from topo_strings_comparer import parse_rostlab_db, determine_c_term
    rost_db = parse_rostlab_db()
    # for the Rost DB:
    # path_withcsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_12Nov/with_msa_with_csts/'
    # path_withcsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_12Nov/minus_2/with_msa_with_csts/'
    # for the entire VH DB:
    # path_withcsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_12Nov/vh_positives/topcons_with_csts/'
    # for Rost with all positive helices and fidelity=0
    # path_withcsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_12Nov/all_pos_no_fidelity_7Dec/with_msa_with_csts/'
    # path_withcsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_14.12/with_msa_with_csts/'
    # path_withcsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_29.2/with_msa_with_csts/'
    path_withcsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_10.3/with_msa_with_csts/'

    prd_files = [a for a in os.listdir(path_withcsts) if a[-4:] == '.prd' and 'msa' not in a]
    no_csts = analyse_sasa_and_positive(False, with_sasa=args['with_sasa'], vh=args['vh'])
    for k, v in no_csts.items():
        print 'top', k, v.keys()
    min_dists = set()
    all_res = ['R', 'K', 'H', 'D', 'E', 'N', 'Q', 'T', 'S', 'V', 'I', 'L', 'F', 'M', 'P', 'G', 'C', 'Y', 'W', 'A']
    #all_res = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    diff_KR = {k: [] for k in all_res}
    returny, results, wgps = {}, {}, {}
    print 'there are %i prds' % len(prd_files)
    print prd_files
    for prd in prd_files:
        name_full = prd[:-4]
        name = name_full.split('_')[0]
        pos_num = int(name_full.split('_')[1].split('.')[0])
        withcsts_wgp = parse_prd(path_withcsts+prd)[0]
        if withcsts_wgp is None:
            print 'found nothing in the prediction', prd
            continue
        nocsts_wgp = no_csts[name][pos_num]['wgp']

        print 'examining', name
        print nocsts_wgp
        # draw_ddg_fwd_rev(nocsts_wgp, polyval)

        # removing wrongly predicted entries (by c term)
        pdbtm_c = determine_c_term(rost_db[name]['pdbtm'])
        opm_c = determine_c_term(rost_db[name]['opm'])
        topgraph_c = '1' if nocsts_wgp.c_term == 'rev' else '2'
        if topgraph_c != pdbtm_c and topgraph_c != pdbtm_c:
            print 'discarding', name
            print pdbtm_c, opm_c, topgraph_c
            continue

        wgps[name] = nocsts_wgp

        # print name, no_csts[name]['chosen_w'], no_csts[name]['sasa'], withcsts_wgp.total_grade, nocsts_wgp.total_grade
        # print nocsts_wgp
        location = [i+1 for i, w in enumerate(nocsts_wgp.path) if w == no_csts[name][pos_num]['chosen_w']][0]
        perc = 100.*(float(location) / float(nocsts_wgp.win_num))
        min_dist = min([location, nocsts_wgp.win_num-location])
        positives = sum([w.seq.count('K')+w.seq.count('R') for w in nocsts_wgp.path])
        # print location, perc
        # if min_dist < 3:
            # print 'skipping', no_csts[name]['chosen_w']
            # continue
        min_dists.add(min_dist)

        for k in all_res:
            diff_KR[k].append({'RK': arg_lys_assymetry(nocsts_wgp, withcsts_wgp, k),
                               'delta': nocsts_wgp.total_grade-withcsts_wgp.total_grade})

        results[name_full] = {'pos_win': no_csts[name][pos_num]['chosen_w'], 'sasa': no_csts[name][pos_num]['sasa'], 'loc_perc': perc,
                         'withcsts': withcsts_wgp.total_grade, 'nocsts': nocsts_wgp.total_grade, 'min_dist': min_dist,
                         'positives': positives, 'num_pos_win': len([w for w in nocsts_wgp.path if w.grade >= -1]),
                         'KR': arg_lys_assymetry(nocsts_wgp, withcsts_wgp, ['K', 'R']), 'pos_loc': location, 'win_num': nocsts_wgp.win_num}
        returny[name_full] = {'KR': arg_lys_assymetry(nocsts_wgp, withcsts_wgp, ['K', 'R', 'H']),
                         'delta': nocsts_wgp.total_grade-withcsts_wgp.total_grade}
    print '%i entries in results ' % len(results.keys())
    data = [[a['pos_win'].grade-a['pos_win'].inner_tail_score, a['sasa'], a['nocsts']-a['withcsts'], k.upper().split('_')[0],
             a['loc_perc'], a['min_dist'], a['positives'],
             a['num_pos_win'], a['KR'], a['pos_win'].seq if a['pos_win'].direction == 'fwd' else a['pos_win'].seq[::-1],
             a['win_num'], a['pos_loc']] for k, a in results.items()]
    print '%i results in sdata' % len(data)
    #
    print "posWin\t\tSASA\t\tDelta\t\tuniprot\tlocation\tmin_dist\tpositives\tnum_pos_wins\tKRdG\tseq\tpos (win_num)"
    data.sort(key = lambda row: row[2])
    for a in data:
        print "%-2.2f\t%2.2r\t%-2.2f\t%s\t%-2.2f\t%-3i\t%-3i\t%-3i\t%-2.2f\t%-20s\t%i (%i)" % (a[0], a[1], a[2], a[3], a[4]-50, a[5], a[6],
                                                                           a[7], a[8], a[9], a[11], a[10])
    # creates a fig with boxplots around the positive win
    # create_boxplot_array(wgps)

    if False:
        plt.figure()
        x = [a['KR'] for a in results.values()]
        y = [a['nocsts']-a['withcsts'] for a in results.values()]
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

        for md, color in zip(min_dists, ['b', 'r', 'g', 'k', 'm']):
            plt.scatter([a['KR'] for a in results.values() if a['min_dist'] == md],
                        [a['nocsts']-a['withcsts'] for a in results.values() if a['min_dist'] == md],
                        c=color, alpha=0.5, label='min_dist %i' % md)
        plt.plot(x, [a*slope + intercept for a in x], 'r', label='Fitted line')
        plt.text(-5, -25, '%-3.3f*X%-3.3f\nR^2=%-3.3f' % (slope, intercept, r_value))
        plt.show()
        print 'r^2 %f, over %i points' % (r_value, len(x))
        print 'p value %f' % p_value
        print 'std error %f' % std_err
        print 'slope %f, intercept %f' % (slope, intercept)
    elif False:
        for i, RK in enumerate(all_res):
            ax = plt.subplot(4, 5, i)
            x = [a['RK'] for a in diff_KR[RK]]
            y = [a['delta'] for a in diff_KR[RK]]
            plt.xlim([-15, 10])
            plt.xticks(fontsize=4)
            plt.yticks(fontsize=4)
            plt.scatter(x, y, s=2, alpha=0.5)
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
            plt.plot(x, [a*slope + intercept for a in x], 'r')
            txt = '%-3.3f*X%-3.3f\nR^2=%-3.3f' % (slope, intercept, r_value)
            print txt
            plt.text(-14, 10, txt, fontsize=6, verticalalignment='top',
                    horizontalalignment='left')
            ax.set_title('%s' % RK, fontsize=12)
            print 'For ', RK
            print 'r^2 %f, over %i points' % (r_value, len(x))
            print 'p value %f' % p_value
            print 'std error %f' % std_err
            print 'slope %f, intercept %f' % (slope, intercept)
            print '\n'
        plt.legend(loc='upper left')
        plt.savefig('result_KRH.png')
        plt.show()
    else:
        return returny


def create_boxplot_array_max_only(wgps):
    results = OrderedDict([(i, []) for i in range(-20, 21)])
    for name, wgp in wgps.items():
        print name
        pos_wins = [w for w in wgp.path if w.grade >= -1]
        max_pos = max([w.grade for w in pos_wins])
        pos_win = [(i, w) for i, w in enumerate(wgp.path) if w.grade == max_pos][0]
        print wgp
        print pos_win
        for i, w in enumerate(wgp.path):
            print i, w
            results[i-pos_win[0]].append(w.grade)
    print results
    plt.boxplot([v for v in results.values()])
    plt.show()


def create_boxplot_array(wgps, pos_most_pos='pos', hline=0):
    lim = 10
    the_rng = range(-lim, lim+1)

    results = OrderedDict([(i, []) for i in the_rng])
    overall_average = []
    lines = []
    for name, wgp in wgps.items():
        # print name
        if pos_most_pos == 'pos':
            pos_wins = [(i, w) for i, w in enumerate(wgp.path) if w.grade-w.inner_tail_score >= -1]
        else:
            pos_wins = [(i, w) for i, w in enumerate(wgp.path) if w.grade == np.median([t.grade for t in wgp.path])]
        print 'ZZZZZZ %s has %i pos wins' % (name, len(pos_wins))
        print wgp
        for pos_win in pos_wins:
            # if pos_most_pos != 'pos' and pos_win[1].grade > -2:
            #     continue

            # print 'pos_win', pos_win
            if pos_win[0] == 0 or pos_win[0] == len(wgp.path)-1 or pos_win[0] == 1 or pos_win[0] == len(wgp.path)-2:
                # print 'AT EDGE!!!!!!!!!', pos_win[0], len(wgp.path)
                continue

            line_ = []
            for i in range(max(0, pos_win[0]-lim), min(len(wgp.path), pos_win[0]+lim+1)):
                ### disregard wins that are positive or beyond a positive
                for other_pos in [a for a in pos_wins if a != pos_win]:
                    if (i-pos_win[0] >= other_pos[0]-pos_win[0] and i-pos_win[0]>0 and other_pos[0]-pos_win[0]>0) or \
                            (i-pos_win[0] <= other_pos[0]-pos_win[0] and i-pos_win[0]<0 and other_pos[0]-pos_win[0]<0):
                        # print 'diregarding', i, other_pos[0], pos_win[0]
                        continue
                ###
                results[i-pos_win[0]].append(wgp.path[i].grade-grade_segment(wgp.path[i].seq[::-1], polyval))
                line_.append([i-pos_win[0], wgp.path[i].grade-grade_segment(wgp.path[i].seq[::-1], polyval)])
            lines.append(line_)
        [overall_average.append(w.grade - grade_segment(w.seq[::-1], polyval)) for w in wgp.path]
    print results
    print 'AAAAAAAAAAAAAAAAAA'
    # for p, res in results.items():
    #     plt.scatter(np.random.normal(p, 0.05, len(res)), res, s=1)
    # plt.boxplot([v for v in results.values()], positions=the_rng)
    # [plt.text(p, 5, len(results[p])) for p in the_rng]
    # plt.hlines(np.mean(overall_average), xmin=-len(the_rng)-0.5, xmax=len(the_rng)+0.5, color='grey')
    # plt.hlines(0, xmin=-len(the_rng)-0.5, xmax=len(the_rng)+0.5, color='k')
    # plt.hlines(hline, xmin=-len(the_rng)-0.5, xmax=len(the_rng)+0.5, color='g')
    # plt.xlabel('TMH index (positive TMH at 0)')
    # plt.ylabel('ddG true-opposite')
    # plt.xlim([-10.5, 10.5])
    fig = plt.figure(1)
    for i, line_ in enumerate(lines):
        fig.add_subplot(11, 5, i)
        plt.plot([a[0] for a in line_], [a[1] for a in line_])
        print 'AAAAAA', i, line_
        plt.xlabel(str(i))
        plt.xticks([a[0] for a in line_])
    plt.show()


def arg_lys_assymetry(wgp_nocsts, wgp_withcsts, ress):
    """
    :param ress: list or res to count in
    :param wgp_nocsts: wgp with no csts (with the positive helix)
    :param wgp_withcsts: wgp with csts (without the positive helix)
    :return: delta between the energy with no csts only for R and K and with the csts
    """
    nocsts_grade = 0
    for w in wgp_nocsts.path:
        seq = w.seq
        seq_a = ''.join(a if a in ress else 'A' for a in seq)
        nocsts_grade += grade_segment(seq_a, polyval)

    withcsts_grade = 0
    for w in wgp_withcsts.path:
        seq = w.seq
        seq_a = ''.join(a if a in ress else 'A' for a in seq)
        withcsts_grade += grade_segment(seq_a, polyval)

    return nocsts_grade - withcsts_grade


def draw_ddg_fwd_rev(wgp, polyval):
    """
    :param wgp: a WinGradePath
    :return: shows a plot of the dG insertion and ddG^native-opposite(direction) of the WGP
    """
    native_grades, ddg_rev_fwd = [], []
    for w in wgp.path:
        native_grades.append(w.grade)
        opposite_seq = w.seq[::-1]
        opposite_grade = grade_segment(opposite_seq, polyval)
        ddg_rev_fwd.append(opposite_grade) # w.grade-
    plt.scatter(range(1, wgp.win_num+1), native_grades, color='b')
    plt.scatter(range(1, wgp.win_num+1), ddg_rev_fwd, color='r')
    plt.hlines(-1, 1, wgp.win_num, colors='k')
    plt.show()


def extrimities_analysis(args):
    root_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_12Nov/vh_positives/'
    if not os.path.isfile(root_path+'with_pos_grades.obj'):
        print 'creating data'
        path_to_no = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_12Nov/vh_positives/topcons_no_csts/'
        all_prds = [a for a in os.listdir(path_to_no) if '.prd' in a]
        all_else, with_pos = {}, {}
        delta_dict = analyse_3d_results(args)
        for prd in all_prds:
            wgp = parse_prd(path_to_no+prd)
            if prd[:-4] in delta_dict.keys():
                if delta_dict[prd[:-4]]['delta'] < 0:
                    with_pos[prd[:-4]] = {'delta': delta_dict[prd[:-4]]['delta'], 'KR': delta_dict[prd[:-4]]['KR'],
                                          'wgp': wgp}
                    print 'adding ione'
                    continue
            if wgp.win_num >= 3:
                all_else[prd[:-4]] = {'wgp': wgp}

        with_pos_grades, all_else_grades = {'first': [], 'mid': [], 'last': []}, {'first': [], 'mid': [], 'last': []}
        for k, v in all_else.items():
            wgp = v['wgp']
            all_else_grades['first'].append(wgp.path[0].grade)
            all_else_grades['last'].append(wgp.path[-1].grade)
            [all_else_grades['mid'].append(a.grade) for a in wgp.path[1:-1]]
        for k, v in with_pos.items():
            wgp = v['wgp']
            with_pos_grades['first'].append(wgp.path[0].grade)
            with_pos_grades['last'].append(wgp.path[-1].grade)
            [with_pos_grades['mid'].append(a.grade) for a in wgp.path[1:-1]]
        pickle.dump(with_pos_grades, open(root_path+'with_pos_grades.obj', 'wb'))
        pickle.dump(with_pos_grades, open(root_path+'all_else_grades.obj', 'wb'))
    else:
        print 'reading data'
        with_pos_grades = pickle.load(open(root_path+'with_pos_grades.obj', 'rb'))
        all_else_grades = pickle.load(open(root_path+'all_else_grades.obj', 'rb'))

    # plt.hist(all_else_grades['first'], histtype='stepfilled', normed=True, bins=40, alpha=0.5, color='deepskyblue', label='all first')
    # plt.hist(all_else_grades['mid'], histtype='stepfilled', normed=True, bins=40, alpha=0.5, color='blue', label='all mid')
    # plt.hist(all_else_grades['last'], histtype='stepfilled', normed=True, bins=40, alpha=0.5, color='cyan', label='all last')
    # plt.hist(with_pos_grades['first'], histtype='stepfilled', normed=True, bins=40, alpha=0.5, color='orangered', label='pos first')
    # plt.hist(with_pos_grades['mid'], histtype='stepfilled', normed=True, bins=40, alpha=0.5, color='maroon', label='pos mid')
    # plt.hist(with_pos_grades['last'], histtype='stepfilled', normed=True, bins=40, alpha=0.5, color='red', label='pos last')
    # plt.legend()
    # plt.show()

    all_grades = {'all else %s' % k: v for k, v in all_else_grades.items()}
    for k, v in with_pos_grades.items():
        all_grades['with pos %s' % k] = v

    # df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in all_grades.items()]))
    # df.plot(kind='density')

    densities = {}
    xs = np.linspace(-20, 10, 300)
    plt.subplot(2, 1, 1)
    for k, v in all_else_grades.items():
        densities[k] = gaussian_kde(v)
        densities[k].covariance_factor = lambda: .25
        densities[k]._compute_covariance()
        plt.plot(xs, densities[k](xs), label=k)
        plt.legend()
    plt.subplot(2, 1, 2)
    for k, v in with_pos_grades.items():
        densities[k] = gaussian_kde(v)
        densities[k].covariance_factor = lambda: .25
        densities[k]._compute_covariance()
        plt.plot(xs, densities[k](xs), label=k)
        plt.legend()

    plt.xlim([-20, 10])
    plt.title('all else')
    plt.show()


def analyse_positive_same_place():
    pd.set_option('display.width', 1000)
    pd.set_option('display.max_columns', 500)
    # for the Rost DB:
    # path_to_nocsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_12Nov/with_msa_no_csts/'
    # path_to_nocsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_12Nov/minus_2/with_msa_no_csts/'
    # for Rost DB with all pos and 0 fidelity
    # path_to_nocsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_12Nov/all_pos_no_fidelity_7Dec/with_msa_no_csts/'
    # for the entire VH DB:
    # path_to_nocsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_12Nov/vh_positives/topcons_no_csts/'
    path_to_nocsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_14.12/with_msa_no_csts/'
    prd_files = [a for a in os.listdir(path_to_nocsts)
                 if a[-4:] == '.prd' and '_msa' not in a]
    rost_db = parse_rostlab_db()
    results, df = {}, pd.DataFrame(columns=['SASA', 'delta', 'pos_grade', 'KR', 'uniprot',
                                            'loc_(tot)', 'seq'])
    total, with_pos, row_num = 0., 0., 1
    total_over_n, over_n_posw_neg_dalta = 0, 0
    hpolynom = MakeHydrophobicityGrade()
    for prd_file in prd_files:
        rost_prd = rost_db[prd_file[:-4]]
        wgp = parse_prd(path_to_nocsts+prd_file)
        if wgp.total_grade == 0.:
            print 'no path FOUND, F#$% it', prd_file
            continue
        total_over_n += 1 if wgp.win_num > 4 else 0
        this_protein_has_pos_n_neg_delta = False
        if wgp.win_num >= 3:
            total += 1 if wgp.win_num >= 6 else 0
            pos_in = [w for w in wgp.path[1:-1] if w.grade >= -1.0]
            if pos_in:
                with_pos += 1
                results[prd_file[:-4]] = {}
                for i, chosen_w in enumerate(pos_in):
                    sasas = win_sasa([chosen_w], rost_prd['pdb'], rost_prd['chain'], prd_file[:-4], rost_prd['seq'])
                    sasas = sasas[0]
                    no_w_grade, no_w_KR = extract_w_and_calculate_best_path(wgp, chosen_w, hpolynom)
                    # results[prd_file[:-4]][i] = {'wgp': wgp, 'pos_in': pos_in, 'chosen_w': chosen_w, 'sasa': sasas,
                    #                              'KRdelta': no_w_KR, 'no_w_grade': no_w_grade, 'win_num': wgp.win_num,
                    #                              'pos_loc': [i for i, w in enumerate(wgp.path) if w == chosen_w][0]}
                    # print results[prd_file[:-4]][i]
                    df.loc[row_num] = {'SASA': '%2.1f' % sasas if sasas is not None else None, 'delta': '%2.1f' % (wgp.total_grade-no_w_grade),
                                       'pos_grade': '%2.1f' % chosen_w.grade,
                                       'KR': '%2.1f' % no_w_KR, 'uniprot': prd_file[:-4].upper(),
                                       'loc_(tot)': '%i*(%i)' % ([i+1 for i, w in enumerate(wgp.path) if w == chosen_w][0],
                                                                 wgp.win_num), 'seq': chosen_w.seq if chosen_w.direction == 'fwd' else chosen_w.seq[::-1]}
                    if not this_protein_has_pos_n_neg_delta:
                        if wgp.total_grade-no_w_grade < 0:
                            over_n_posw_neg_dalta += 1
                            this_protein_has_pos_n_neg_delta = True
                    row_num += 1
    print df
    print 'total over 4 TMHs %i, all with pos w and negativs ddG %i percentage: %f' % (total_over_n, over_n_posw_neg_dalta,
                                                                                     100*over_n_posw_neg_dalta/total_over_n)


def protter_make():
    path_to_nocsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_14.12/with_msa_no_csts/'
    prd_files = 'o29285.prd'
    wgp = parse_prd(path_to_nocsts+prd_files)
    print 'segments for original:\n', '\n'.join([a.seq if a.direction == 'fwd' else a.seq[::-1] for a in wgp.path])
    extract_w_and_calculate_best_path(wgp, wgp.path[7], MakeHydrophobicityGrade())


def extract_w_and_calculate_best_path(wgp, chosen_win, polyval):
    path_removed = wgp.path[:]
    path_removed.remove(chosen_win)
    all_fwd = [a.seq if a.direction == 'fwd' else a.seq[::-1] for a in path_removed]
    fwd_first = [all_fwd[0]]
    last_is_ = 'fwd'
    for w in all_fwd[1:]:
        if last_is_ == 'fwd':
            fwd_first.append(w[::-1])
            last_is_ = 'rev'
        else:
            fwd_first.append(w)
            last_is_ = 'fwd'
    fwd_first_grade = sum(grade_segment(w, polyval) for w in fwd_first)
    rev_first = [a[::-1] for a in fwd_first]
    rev_first_grade = sum([grade_segment(a, polyval) for a in rev_first])
    if fwd_first_grade < rev_first_grade:
        best_path = fwd_first
        best_grade = fwd_first_grade
        print 'best path without pos win was chosen as first helix FWD'
    else:
        best_path = rev_first
        best_grade = rev_first_grade
        print 'best path without pos win was chosen as first helix REV'

    # for i, k in enumerate(wgp.path[:-1]):
    #     print i, ''.join(a if a in ['K', 'R'] else 'A' for a in k.seq), k.seq
    #     print i, ''.join(a if a in ['K', 'R'] else 'A' for a in best_path[i]), best_path[i]
    KR_grade_best = sum([grade_segment(''.join(a if a in ['K', 'R'] else 'A' for a in w), polyval) for w in best_path])
    KR_grade_original = sum([grade_segment(''.join(a if a in ['K', 'R'] else 'A' for a in w.seq), polyval) for w in wgp.path])
    return best_grade, KR_grade_original-KR_grade_best


def get_ddg_from_prd(prd_file):
    with open(prd_file, 'r') as fin:
        for l in fin.read().split('\n'):
            s = l.split()
            if s[0] == 'ddG':
                return float(s[2].split('|')[1].split('|')[0])


def vh_boxplots():
    vh_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_14.12/VH_topcons_31Jan/'
    prd_files = [a for a in os.listdir(vh_path) if '.prd' in a]
    wgps, wgps_no_pos, no_pos_grade = {}, {}, []
    dG_1, dG_m1, dG_others = [], [], []

    fig = plt.figure()
    p = 1
    no_pos_all, pos_all = [], []
    KRx, KRy = {0: [], 1: [], 2: [], 3:[], 4: []}, {0: [], 1: [], 2: [], 3:[], 4: []}
    pp, mp, mm = [], [], []
    for prd_file in prd_files:
        wgp = parse_prd(vh_path+prd_file)[0]

        if wgp is None:
            continue

        if len(wgp.path) < 5:
            continue

        if get_ddg_from_prd(vh_path+prd_file) < 6:
            continue
        # validation:
        pos_wins = [w for w in wgp.path if w.grade >= -1]
        if len(pos_wins) > 1:
            # print 'found %i pos wins' % len(pos_wins)
            continue

        for w in wgp.path:
            ddG = w.grade-grade_segment(w.seq[::-1], polyval)
            KRx[w.seq.count('K')+w.seq.count('R')].append(w.grade)
            KRy[w.seq.count('K')+w.seq.count('R')].append(ddG)

        if len(pos_wins) == 0:
            [no_pos_grade.append(w.grade-grade_segment(w.seq[::-1], polyval)) for w in wgp.path]
            for w in wgp.path:
                ddG = w.grade-grade_segment(w.seq[::-1], polyval)
                if w.grade > 0 and ddG > 0:
                    pp.append([w.seq, w.grade, ddG])
                if w.grade < 0 and ddG > 0:
                    mp.append([w.seq, w.grade, ddG])
                if w.grade < 0 and ddG < 0:
                    mm.append([w.seq, w.grade, ddG])
                # KRx[w.seq.count('K')+w.seq.count('R')].append(w.grade)
                # KRy[w.seq.count('K')+w.seq.count('R')].append(ddG)

            wgps_no_pos[prd_file[:-4]] = wgp
        wgps[prd_file[:-4]] = wgp

        if pos_wins == []:
            [dG_others.append(w.grade) for w in wgp.path]
            [no_pos_all.append([w.grade, grade_segment(w.seq[::-1], polyval)]) for w in wgp.path]
        # ignore proteins with mTMHs at edges
        elif pos_wins[0] == wgp.path[0] or pos_wins[0] == wgp.path[-1]:
            continue
        else:
            pos_win_i = [i for i, w in enumerate(wgp.path) if w == pos_wins[0]][0]
            # print 'RRR'
            # print str(wgp)
            # print pos_wins, pos_win_i
            # print wgp.path[pos_win_i+1].grade, wgp.path[pos_win_i-1].grade
            [dG_others.append(w.grade) for i, w in enumerate(wgp.path) if (i < pos_win_i-1 or i > pos_win_i+1) and
             (i != 0 and i != len(wgp.path)-1)]
            dG_1.append(wgp.path[pos_win_i+1].grade)
            dG_m1.append(wgp.path[pos_win_i-1].grade)

            [pos_all.append([w.grade, grade_segment(w.seq[::-1], polyval)]) for w in wgp.path]

            fig.add_subplot(11, 5, p)
            p += 1
            wgp_max = wgp.path[-1].end+10
            for w in wgp.path:
                plt.hlines(w.grade, w.begin, w.end, color='r' if w.grade >= -1 else 'k')
                plt.scatter(np.mean([w.begin, w.end]), w.grade-grade_segment(w.seq[::-1], polyval))
                plt.hlines(0, xmin=0, xmax=wgp_max, color='grey')
                plt.xlim([0, wgp_max])
                ddG = w.grade-grade_segment(w.seq[::-1], polyval)
                if w.grade > 0 and ddG > 0:
                    pp.append([w.seq, w.grade, ddG])
                if w.grade < 0 and ddG > 0:
                    mp.append([w.seq, w.grade, ddG])
                if w.grade < 0 and ddG < 0:
                    mm.append([w.seq, w.grade, ddG])
                # KRx[w.seq.count('K')+w.seq.count('R')].append(w.grade)
                # KRy[w.seq.count('K')+w.seq.count('R')].append(ddG)

    print 'found %i non 1/-1 found avg %f' % (len(dG_others), np.mean(dG_others))
    print 'found %i 1 found avg %f' % (len(dG_1), np.mean(dG_1))
    print 't-test diff %f and the p-value is %f' % ttest_1samp(dG_1, np.mean(dG_others))
    print 'found %i -1 found avg %f' % (len(dG_m1), np.mean(dG_m1))
    print 't-test diff %f and the p-value is %f' % ttest_1samp(dG_m1, np.mean(dG_others))
    # plt.boxplot([dG_m1, dG_1, dG_others])
    plt.show()

    # print 'x with pos', [a[0] for a in pos_all]
    # print 'y with pos', [a[1] for a in pos_all]
    # print 'x no pos', [a[0] for a in no_pos_all]
    # print 'y no pos', [a[1] for a in no_pos_all]

    plt.scatter([a[0] for a in no_pos_all], [a[1] for a in no_pos_all], color='r')
    plt.scatter([a[0] for a in pos_all], [a[1] for a in pos_all], color='b')

    all_x = [a[0] for a in no_pos_all] + [a[0] for a in pos_all]
    all_y = [a[1] for a in no_pos_all] + [a[1] for a in pos_all]
    slope, intercept, r_value, p_value, std_err = stats.linregress(all_x, all_y)
    plt.plot(np.arange(-10, 20), [intercept+slope*a for a in np.arange(-10, 20)])
    print 'found %f r value and %f p-value' % (r_value**2, p_value)
    plt.hlines(0, -20, 10)
    plt.vlines(0, -20, 15)
    plt.show()

    # create_boxplot_array(wgps, pos_most_pos='pos', hline=np.mean(no_pos_grade))
    # create_boxplot_array(wgps_no_pos, pos_most_pos='all',hline=np.mean(no_pos_grade))

    # print 'seqs plus plus', pp
    # print 'seqs minus plus', mp
    # print 'seqs minus minus', mm

    for i, k in enumerate(KRx.keys()):
        # print 'moving to k %i' % k
        for x, y in zip(KRx[k], KRy[k]):
            print '%.3f%s%.3f' % (x, str(''.join(['\t']*(i+1))), y)


def topgraphMSA_phase_diagram():
    from WinGrade import topo_string_to_WGP
    from topo_strings_comparer import overlappM
    # msa_path = '/home/labs/fleishman/elazara/length_21/w_0_with_MSA/'
    msa_path = '/home/labs/fleishman/elazara/TopGraph/TopGarph_symmetric/vh/msa_no_fildel/for_jonathan/' # vh msa symmetric
    prd_files = [a for a in os.listdir(msa_path) if '.prd' in a and 'msa' not in a]
    rost_db = parse_rostlab_db()
    vh_names_lower_to_wierd = vh_name_dict()

    result = []
    # true, false = [], []
    for prd_file in prd_files:
        wgp, clue = parse_prd(msa_path+prd_file)
        if clue is None:
            continue
        try:
            spoc = spc_parser(rost_db[prd_file[:-4]]['name'])['topcons']
            wgp_pdbtm = topo_string_to_WGP(rost_db[prd_file[:-4]]['pdbtm'], rost_db[prd_file[:-4]]['seq'])
        except:
            spoc = spc_parser(vh_names_lower_to_wierd[prd_file[:-4]])['topcons']
        signal = [0, spoc.count('s') + spoc.count('S')]

        for w in wgp.path:
            if w.begin <= signal[1]:
                continue
            # predicted = overlappM([w.begin, w.end], [[a.begin, a.end] for a in wgp_pdbtm.path])
            # if predicted:
            #     true.append([w.grade, w.grade-grade_segment(w.seq[::-1], polyval)])
            # else:
            #     print 'fail', w
            #     print wgp_pdbtm
            #     false.append([w.grade, w.grade-grade_segment(w.seq[::-1], polyval)])
            result.append([w.grade, w.grade-grade_segment(w.seq[::-1], polyval)])

    if args['mode'] != 'msa':
        return [a[0] for a in result], [a[1] for a in result]

    soluble_dGs, soluble_ddGs = collect_win_analysis(dict())

    plt.scatter(soluble_dGs, soluble_ddGs, edgecolors='g', alpha=0.2)
    plt.scatter([a[0] for a in result], [a[1] for a in result], color='r', alpha=0.8)

    # plt.scatter([a[0] for a in true], [a[1] for a in true], color='k')
    # plt.scatter([a[0] for a in false], [a[1] for a in false], color='r')
    xmin = min([a[0] for a in result]+soluble_dGs)  # +[a[0] for a in false])
    xmax = max([a[0] for a in result]+soluble_dGs)  # +[a[0] for a in false])
    ymin = min([a[1] for a in result]+soluble_ddGs)  # +[a[1] for a in false])
    ymax = max([a[1] for a in result]+soluble_ddGs)  # +[a[1] for a in false])
    #
    plt.hlines(0, xmin, xmax)
    plt.vlines(0, ymin, ymax)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    # plt.show()

    # for a in result:
    #     print a[0], '\t', a[1]


def draw_random_combined_scatter_hist():
    import random
    data_1_x = [int(1000*random.random()) for i in xrange(100)]
    data_1_y = [int(1000*random.random()) for i in xrange(100)]
    data_2_x = [int(1000*random.random()) for i in xrange(100)]
    data_2_y = [int(1000*random.random()) for i in xrange(100)]

    ax_scat = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)
    ax_hist_dG = plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=1)
    ax_hist_ddG = plt.subplot2grid((3, 3), (1, 2), colspan=1, rowspan=2)

    ax_scat.scatter(data_1_x, data_1_y, alpha=0.2, color='grey')
    ax_scat.scatter(data_2_x, data_2_y, alpha=0.8, color='red')

    ax_hist_dG.hist(data_1_x, color='grey', alpha=0.2)
    ax_hist_dG.hist(data_2_x, color='red', alpha=0.8)

    ax_hist_ddG.hist(data_1_y, color='grey', alpha=0.2, orientation='horizontal')
    ax_hist_ddG.hist(data_2_y, color='red', alpha=0.8, orientation='horizontal')

    plt.show()


def draw_combined_scatter_hist():
    import matplotlib.gridspec as gridspec
    import pickle
    import os
    wp = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/soluble_helix_grading/data.obj'

    if not os.path.isfile(wp):
        print 'gathering data from soluble'
        soluble_dGs, soluble_ddGs = collect_win_analysis(dict())
        print 'gathering TM data'
        tm_dGs, tm_ddGs = topgraphMSA_phase_diagram()
        print 'finished gathering data'
        pickle.dump((soluble_dGs, soluble_ddGs, tm_dGs, tm_ddGs), open(wp, 'wb'))
    else:
        print 'reading pickled data'
        soluble_dGs, soluble_ddGs, tm_dGs, tm_ddGs = pickle.load(open(wp, 'rb'))

    soluble_color = 'cornflowerblue'
    tm_color = 'indianred'

    formatter = FuncFormatter(to_percent)

    xmin, xmax = min(soluble_dGs+tm_dGs)-1, max(soluble_dGs+tm_dGs)+1
    ymin, ymax = min(soluble_ddGs+tm_ddGs)-1, max(soluble_ddGs+tm_ddGs)+1

    # arrang
    gs1 = gridspec.GridSpec(3, 3)
    gs1.update(wspace=0.0, hspace=0.0)
    ax_scat = plt.subplot(gs1[1:, :-1])
    ax_hist_dG = plt.subplot(gs1[0, :-1])
    ax_hist_ddG = plt.subplot(gs1[1:, -1])

    ax_scat.axhline(color='k', linestyle='dotted', markeredgewidth=8)
    ax_scat.axvline(color='k', linestyle='dotted', markeredgewidth=8)
    ax_scat.scatter(soluble_dGs, soluble_ddGs, alpha=0.2, color=soluble_color)
    ax_scat.scatter(tm_dGs, tm_ddGs, alpha=0.4, color=tm_color)
    ax_scat.set_xlim([xmin, xmax])
    ax_scat.set_ylim([ymin, ymax])
    ax_scat.set_xlabel('dG')
    ax_scat.set_ylabel('ddG')

    ax_hist_dG.hist(soluble_dGs, color=soluble_color, alpha=0.2, normed=1, bins=50)
    ax_hist_dG.hist(tm_dGs, color=tm_color, alpha=0.4, normed=1, bins=50)
    ax_hist_dG.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off')
    ax_hist_dG.set_xlim([xmin, xmax])
    ax_hist_dG.yaxis.set_major_formatter(formatter)
    ax_hist_dG.set_ylabel('Frequency')

    ax_hist_ddG.hist(soluble_ddGs, color=soluble_color, alpha=0.2, orientation='horizontal', normed=1, bins=50)
    ax_hist_ddG.hist(tm_ddGs, color=tm_color, alpha=0.4, orientation='horizontal', normed=1, bins=50)
    ax_hist_ddG.tick_params(axis='both', which='both', bottom='off', top='off', labelleft='off', right='off', left='off')
    ax_hist_ddG.set_ylim([ymin, ymax])
    ax_hist_ddG.xaxis.set_major_formatter(formatter)
    ax_hist_ddG.set_xlabel('Frequency')

    plt.show()


def analyse_positive_inside_results(args):
    score_threshold = -1
    # work_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_10.3/with_msa_no_csts/'
    work_path = '/home/labs/fleishman/elazara/TopGraph/TopGarph_symmetric/rost/temp_msa_50/correct_only/correct_c_term_non_conservative/' # rostr msa correct by C term
    prd_files = [a for a in os.listdir(work_path) if '.prd' in a]

    rost_db = parse_rostlab_db()
    df = pd.DataFrame(columns=['uniprot', 'pos_win_grade', 'pos_win_num', 'delta_wgp_flipped', 'SASA', 'pos_total', 'seq'])

    for prd_file in prd_files:
        name = prd_file.split('.')[0]
        full_wgp = parse_prd(work_path+prd_file)[0]

        # skipping if no wgp, or less than 3 wins
        if full_wgp is None:
            continue
        # print full_wgp
        if full_wgp.win_num < 4:
            # print '%s has only %i windows. skipping' % (prd_file, full_wgp.win_num)
            continue

        pos_wins = [w for w in full_wgp.path if w.grade > score_threshold and
                    full_wgp.path[-1] != w != full_wgp.path[0]]

        rost_prd = rost_db[name]

        for i, pos_w in enumerate(pos_wins):
            # print 'positive window', pos_w
            delta_wgp_flipped, flipped_grade, flipped_c_term, flipped_n_term = \
                find_constrained_delta(full_wgp, pos_w, rost_prd['seq'])
            w_sasa = win_sasa([pos_w], rost_prd['pdb'], rost_prd['chain'], name, rost_prd['seq'])[0]
            pos_win_index = [i+1 for i, w in enumerate(full_wgp.path) if w == pos_w][0]
            df = df.append({'uniprot': name,
                            'pos_win_grade': pos_w.grade,
                            'pos_win_num': pos_win_index,
                            'delta_wgp_flipped': delta_wgp_flipped,
                            'flipped_grade': flipped_grade,
                            'SASA': w_sasa,
                            'pos_total': '%i (%i)' % (pos_win_index, full_wgp.win_num),
                            'seq': pos_w.seq if pos_w.direction == 'fwd' else pos_w.seq[::-1],
                            'flipped_c_term': 'in' if flipped_c_term == 'rev' else 'out',
                            'flipped_n_term': 'in' if flipped_n_term == 'rev' else 'out'},
                           ignore_index=True)
    print df


def find_constrained_delta(wgp, w_pos, full_seq, verbose=False):
    """
    :param verbose: whether to print stuff or not
    :param wgp: the WinGradePath
    :param full_seq: full protein sequence
    :param w_pos: a positive window to remove
    :return: the delta between, and the best grade
    """
    pos_win_ind = [i for i, w in enumerate(wgp.path) if w == w_pos][0]

    # crerate wgp where all wins after pos win are flipped
    opt_1_list = []
    for i, w in enumerate(wgp.path):
        if i < pos_win_ind:
            opt_1_list.append(w)
        elif i == pos_win_ind:
            continue
        elif i > pos_win_ind:
            # opt_1_list.append(WinGrade(w.begin, w.end, 'fwd' if w.direction == 'rev' else 'rev', w.seq[::-1], polyval,
            #                            inner_tail=find_inner_tail(full_seq, w, 'fwd' if w.direction == 'rev' else 'rev')
            #                            ))
            opt_1_list.append(flip_win_grade(w, full_seq))
    opt_1_wgp = WinGradePath(opt_1_list)

    # crerate wgp where all wins before pos win are flipped
    opt_2_list = []
    for i, w in enumerate(wgp.path):
        if i > pos_win_ind:
            opt_2_list.append(w)
        elif i == pos_win_ind:
            continue
        elif i < pos_win_ind:
            # opt_2_list.append(WinGrade(w.begin, w.end, 'fwd' if w.direction == 'rev' else 'rev', w.seq[::-1], polyval,
            #                            inner_tail=find_inner_tail(full_seq, w, 'fwd' if w.direction == 'rev' else 'rev')
            #                            ))
            opt_2_list.append(flip_win_grade(w, full_seq))
    opt_2_wgp = WinGradePath(opt_2_list)

    best_wgp = opt_1_wgp if opt_1_wgp.total_grade < opt_2_wgp.total_grade else opt_2_wgp

    if verbose:
        print 'wgp\n', wgp
        print 'wgp1\n', opt_1_wgp
        print 'wgp2\n', opt_2_wgp
        print 'best\n', best_wgp
        print 'seq %s\n' % full_seq

    return wgp.total_grade - best_wgp.total_grade, best_wgp.total_grade, best_wgp.c_term, best_wgp.path[0].direction


def find_inner_tail(full_seq, w, direction):
    if direction == 'fwd':
        return full_seq[max([0, w.begin-5]): w.begin]
    elif direction == 'rev':
        return full_seq[w.end: min([w.end+5, len(full_seq)])]



def vh_name_dict():
    """
    :return: {lower_case: wierd case}
    """
    list_ = [a.split('.')[0] for a in os.listdir('/home/labs/fleishman/jonathaw/membrane_prediction_DBs/spoctopus_SPDB')
            if 3 <= len(a.split('.')[0]) <= 4 and '.spc' in a]
    return {k.lower(): k for k in list_}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode')
    parser.add_argument('-with_sasa', default=True)
    parser.add_argument('-vh', default=False)
    args = vars(parser.parse_args())

    if args['mode'] == 'analyse_sasa_and_positive':
        analyse_sasa_and_positive(True, args['with_sasa'], args['vh'])

    elif args['mode'] == 'analyse_results':
        analyse_results()

    elif args['mode'] == '3d_results':
        analyse_3d_results(args)

    elif args['mode'] == 'extrimities_analysis':
        extrimities_analysis(args)

    elif args['mode'] == 'positive_inside':
        analyse_positive_same_place()

    elif args['mode'] == 'protter':
        protter_make()

    elif args['mode'] == 'vh_boxplot':
        vh_boxplots()

    elif args['mode'] == 'msa':
        topgraphMSA_phase_diagram()

    elif args['mode'] == 'from_csts_to_prds':
        from_csts_to_prds()

    elif args['mode'] == 'combined':
        draw_combined_scatter_hist()

    elif args['mode'] == 'new_pos':
        analyse_positive_inside_results(args)

    elif args['mode'] == 'test':
        draw_combined_scatter_hist()

    else:
        print 'no mode chosen'
