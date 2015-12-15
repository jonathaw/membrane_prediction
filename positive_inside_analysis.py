#!/usr/bin/env python2.7
import os
import argparse
from scipy import stats
from TMpredict_WinGrade import parse_rostlab_db
from WinGrade import parse_WGP
from TMConstraint import TMConstraint
from sasa_survey_26Aug import parse_rsa, nacces_for_win, pair_wise_aln_from_seqs, parse_standard_data
from topo_strings_comparer import spc_parser
from TMpredict_WinGrade import topo_VH_parser
import matplotlib.pyplot as plt
from TMpredict_WinGrade import MakeHydrophobicityGrade
import pickle
import pandas as pd
from scipy.stats import gaussian_kde
import numpy as np
import sys

polyval = MakeHydrophobicityGrade()


def analyse_sasa_and_positive(write=False, with_sasa=True, vh=False):
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
        wgp = parse_prd(path_to_nocsts+prd_file)
        if wgp.total_grade == 0.:
            print 'no path FOUND, F#$% it', prd_file
            continue
        if wgp.win_num >= 3:
            total += 1 if wgp.win_num >= 6 else 0
            three_n_over += 1
            pos_in = [w for w in wgp.path[1:-1] if w.grade >= -1.0]
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
                    # print prd_file
                    if write:
                        try:
                            print "took %f with sasa %r at %i-%i" % (chosen_w.grade, sasas, chosen_w.begin, chosen_w.end)
                            cst = parse_prd_csts(prd_file)
                            for tm_pos in cst.tm_pos:
                                if tm_pos[0]-cst.tm_pos_fidelity <= chosen_w.begin and chosen_w.end <= tm_pos[1]+cst.tm_pos_fidelity:
                                    cst.tm_pos.remove(tm_pos)
                                    cst.non_tm_pos = [[chosen_w.begin, chosen_w.end]]
                                    cst.tm_pos_fidelity = 0
                                    # print str(cst)
                            with open('remove_w_csts/%s_%i.cst' % (prd_file[:-4], i), 'w+') as fout:
                                fout.write(str(cst))
                        except:
                            pass
    if write and not vh:
        print "found %i with positive in, out of %i with 3 or over" % (len(results.keys()), three_n_over)
        for k, v in results.items():
            print 'for %s there are %i positive insides\t %s\t %r' % (k, len(v['pos_in']), rost_db[k]['pdb'], v['pos_in'])
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

    spoc = spc_parser(name)['spoctopus']
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
        return None
    wgp = parse_WGP(one_liner)
    return wgp


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
    path_withcsts = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/positive_inside_14.12/with_msa_with_csts/'
    prd_files = [a for a in os.listdir(path_withcsts) if a[-4:] == '.prd' and 'msa' not in a]
    no_csts = analyse_sasa_and_positive(False, with_sasa=args['with_sasa'], vh=args['vh'])
    for k, v in no_csts.items():
        print k, v.keys()
    min_dists = set()
    all_res = ['R', 'K', 'H', 'D', 'E', 'N', 'Q', 'T', 'S', 'V', 'I', 'L', 'F', 'M', 'P', 'G', 'C', 'Y', 'W', 'A']
    #all_res = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    diff_KR = {k: [] for k in all_res}
    returny, results = {}, {}
    print 'there are %i prds' % len(prd_files)
    print prd_files
    for prd in prd_files:
        name_full = prd[:-4]
        name = name_full.split('_')[0]
        pos_num = int(name_full.split('_')[1].split('.')[0])
        withcsts_wgp = parse_prd(path_withcsts+prd)
        if withcsts_wgp is None:
            print 'found nothing in the prediction', prd
            continue
        nocsts_wgp = no_csts[name][pos_num]['wgp']

        # removing wrongly predicted entries (by c term)
        pdbtm_c = determine_c_term(rost_db[name]['pdbtm'])
        opm_c = determine_c_term(rost_db[name]['opm'])
        topgraph_c = '1' if nocsts_wgp.c_term == 'rev' else '2'
        if topgraph_c != pdbtm_c and topgraph_c != pdbtm_c:
            print 'discarding', name
            print pdbtm_c, opm_c, topgraph_c
            continue

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
    data = [[a['pos_win'].grade, a['sasa'], a['nocsts']-a['withcsts'], k.upper().split('_')[0], a['loc_perc'], a['min_dist'], a['positives'],
             a['num_pos_win'], a['KR'], a['pos_win'].seq if a['pos_win'].direction == 'fwd' else a['pos_win'].seq[::-1],
             a['win_num'], a['pos_loc']] for k, a in results.items()]
    print '%i results in sdata' % len(data)
    #
    print "posWin\t\tSASA\t\tDelta\t\tuniprot\tlocation\tmin_dist\tpositives\tnum_pos_wins\tKRdG\tseq\tpos (win_num)"
    data.sort(key = lambda row: row[2])
    for a in data:
        print "%-2.2f\t%2.2r\t%-2.2f\t%s\t%-2.2f\t%-3i\t%-3i\t%-3i\t%-2.2f\t%-20s\t%i (%i)" % (a[0], a[1], a[2], a[3], a[4]-50, a[5], a[6],
                                                                           a[7], a[8], a[9], a[11], a[10])

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


def grade_segment(seq, polyval):
        membrane_position = np.linspace(-20, 20, endpoint=True, num=len(seq))
        grade = 0
        for i, aa in enumerate(seq):
            if aa in polyval.keys():
                grade += np.polyval(polyval[aa], membrane_position[i])
        return grade


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
    else:
        best_path = rev_first
        best_grade = rev_first_grade

    # for i, k in enumerate(wgp.path[:-1]):
    #     print i, ''.join(a if a in ['K', 'R'] else 'A' for a in k.seq), k.seq
    #     print i, ''.join(a if a in ['K', 'R'] else 'A' for a in best_path[i]), best_path[i]
    KR_grade_best = sum([grade_segment(''.join(a if a in ['K', 'R'] else 'A' for a in w), polyval) for w in best_path])
    KR_grade_original = sum([grade_segment(''.join(a if a in ['K', 'R'] else 'A' for a in w.seq), polyval) for w in wgp.path])
    return best_grade, KR_grade_original-KR_grade_best

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

    else:
        print 'no mode chosen'
