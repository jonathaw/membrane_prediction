#!/usr/bin/env python2.7
import os
import re
import argparse
import matplotlib.pyplot as plt

from positive_inside_analysis import parse_prd
from topo_strings_comparer import spc_parser
from TMpredict_WinGrade import parse_rostlab_db

predictors = ['polyphobius', 'topcons', 'spoctopus', 'philius', 'octopus', 'scampi']


def main_rost():
    prd_files = [a for a in os.listdir('./') if '.prd' in a and '_msa' not in a]
    rost_db = parse_rostlab_db()
    new_old = rost_new_old()
    topgraph_none = []

    follow = 'q8dkp6'

    old_new_totals = {'new': 0, 'old': 0}
    results = {}
    for prd_file in prd_files:
        name = prd_file.split('.')[0].lower()
        best_wgp, sec_wgp = parse_prd(prd_file)

        if best_wgp is None:
            topgraph_none.append(name)
            continue

        topc = spc_parser(name)

        signal_peptide = topc['topcons'].count('S') + topc['topcons'].count('s')

        best_wgp_loc_list = wgp_to_loc_list(best_wgp, signal_peptide)
        sec_wgp_loc_list = wgp_to_loc_list(sec_wgp, signal_peptide)

        old_new_totals[new_old[name]] += 1

        if name == follow:
            print 'at %s found loc list %r' % (name, best_wgp_loc_list)

        best_tgr_qok, best_tgr_ovm = qok_pdbtm_opm(rost_db[name], best_wgp_loc_list, signal_peptide, verbose=name==follow)
        sec_tgr_qok, sec_tgr_ovm = qok_pdbtm_opm(rost_db[name], sec_wgp_loc_list, signal_peptide)

        best_or_sec_qok = best_tgr_qok or sec_tgr_qok
        best_or_sec_ovm = best_tgr_ovm or sec_tgr_ovm

        results[name] = {'old_new': new_old[name],
                         'tm_num': len(pdbtm_opm_loc_list(rost_db[name]['pdbtm'], signal_peptide)),
                         'topgraph': {'qok': best_tgr_qok, 'ovm': best_tgr_ovm},
                         'best_or_sec': {'qok': best_or_sec_qok, 'ovm': best_or_sec_ovm}}

        for predictor in predictors:
            prd_qok, prd_ovm = qok_pdbtm_opm(rost_db[name], ts_loc_list(topc[predictor], signal_peptide), signal_peptide)
            results[name][predictor] = {'qok': prd_qok, 'ovm': prd_ovm}

    # prints resutls sliced by old/new
    print_results_by_old_new(results, predictors, old_new_totals)

    # prints results sliced by 1, 2-4 >4 TMHs
    print_results_by_tm_num(results)

    # print names TopGraph got wrong
    print_names_topgraph_got_wrong(results)

    # prints namse TopGraph got wrong by both best and sec best
    print_names_topgraph_got_wrong_best_and_sec(results)

    # print total percentage correct for TopGraph, TopGraph best or sec, and TOPCONS
    print_total_results(results)


def print_total_results(results):
    print 'Total results'
    total, total_right, total_best_or_sec, total_topc = 0, 0, 0, 0
    for v in results.values():
        total += 1
        total_right += 1 if v['topgraph']['ovm'] else 0
        total_best_or_sec += 1 if v['best_or_sec']['ovm'] else 0
        total_topc += 1 if v['topcons']['ovm'] else 0
    print 'TopGrapgh got %f right' % (100.*total_right/total)
    print 'TopGrapgh best or sec got %f right' % (100.*total_best_or_sec/total)
    print 'TOPCONS got %f right' % (100.*total_topc/total)


def print_names_topgraph_got_right(results):
    print 'these are the entries topgraph got RIGHT by overlapM:'
    for k, v in results.items():
        if v['topgraph']['ovm']:
            print k


def print_names_topgraph_got_wrong(results):
    print 'these are the entries topgraph got WRONG by overlapM:'
    l = []
    for k, v in results.items():
        if not v['topgraph']['ovm']:
            l.append(k)
            print k
    print l


def print_names_topgraph_got_wrong_best_and_sec(results):
    print 'these are the entries TopGraph got WRONG at best and second best paths:'
    l = []
    for k, v in results.items():
        if not v['best_or_sec']['ovm']:
            l.append(k)
            print k
    print l


def print_results_by_old_new(results, predictors, old_new_totals):
    print '########## results by old / new:'
    new_old_sums_qok = {predictor: {'new': 0, 'old': 0} for predictor in predictors+['topgraph']}
    new_old_sums_ovm = {predictor: {'new': 0, 'old': 0} for predictor in predictors+['topgraph']}
    for v in results.values():
        for predictor in predictors+['topgraph']:
            new_old_sums_qok[predictor][v['old_new']] += 1 if v[predictor]['qok'] else 0
            new_old_sums_ovm[predictor][v['old_new']] += 1 if v[predictor]['ovm'] else 0
    new_percentages_qok = {predictor: 100.0*new_old_sums_qok[predictor]['new']/old_new_totals['new'] for predictor in
                       predictors+['topgraph']}
    old_percentages_qok = {predictor: 100.0*new_old_sums_qok[predictor]['old']/old_new_totals['old'] for predictor in
                       predictors+['topgraph']}
    new_percentages_ovm = {predictor: 100.0*new_old_sums_ovm[predictor]['new']/old_new_totals['new'] for predictor in
                       predictors+['topgraph']}
    old_percentages_ovm = {predictor: 100.0*new_old_sums_ovm[predictor]['old']/old_new_totals['old'] for predictor in
                       predictors+['topgraph']}

    print 'new_percentages_qok', new_percentages_qok
    print 'old_percentages_qok', old_percentages_qok
    print 'new_percentages_ovm', new_percentages_ovm
    print 'old_percentages_ovm', old_percentages_ovm
    print '########## end results old / new'


def print_results_by_tm_num(results):
    print '%%%%%%%%%% results by 1 2-4 5:'
    totals = {'1': 0, '2-4': 0, '5': 0}
    sums_qok = {predictor: {'1': 0, '2-4': 0, '5': 0} for predictor in predictors+['topgraph', 'best_or_sec']}
    sums_ovm = {predictor: {'1': 0, '2-4': 0, '5': 0} for predictor in predictors+['topgraph', 'best_or_sec']}

    # either_sums_qok = {predictor: {'1': 0, '2-4': 0, '5': 0} for predictor in predictors+['topgraph', 'best_or_sec']}
    # either_sums_ovm = {predictor: {'1': 0, '2-4': 0, '5': 0} for predictor in predictors+['topgraph', 'best_or_sec']}

    for v in results.values():
        if v['tm_num'] == 1:
            tm_num = '1'
        elif 2 <= v['tm_num'] <= 4:
            tm_num = '2-4'
        elif v['tm_num'] >= 5:
            tm_num = '5'

        # either_sums_qok['topgraph'][tm_num] += 1 if v['best_or_sec']['qok'] else 0
        # either_sums_ovm['topgraph'][tm_num] += 1 if v['best_or_sec']['ovm'] else 0

        totals[tm_num] += 1

        for predictor in predictors+['topgraph', 'best_or_sec']:
            sums_qok[predictor][tm_num] += 1 if v[predictor]['qok'] else 0
            sums_ovm[predictor][tm_num] += 1 if v[predictor]['ovm'] else 0

    tm_num_1 = {predictor: 100.0*sums_ovm[predictor]['1']/totals['1'] for predictor in
                predictors+['topgraph', 'best_or_sec']}
    tm_num_2_4 = {predictor: 100.0*sums_ovm[predictor]['2-4']/totals['2-4'] for predictor in
                  predictors+['topgraph', 'best_or_sec']}
    tm_num_5 = {predictor: 100.0*sums_ovm[predictor]['5']/totals['5'] for predictor in
                predictors+['topgraph', 'best_or_sec']}

    print 'the totals are:'
    for k, v in totals.items():
        print k, v

    print '1 TM results:'
    for k, v in tm_num_1.items():
        print k, v

    print '2-4 TM results:'
    for k, v in tm_num_2_4.items():
        print k, v

    print '5 TM results:'
    for k, v in tm_num_5.items():
        print k, v
    print '%%%%%%%%%% end results by 1 2-4 5'


def pred_loc_list_clean_by_uuus(loc_list, uuus_list, verbose=False):
    """
    :param loc_list: prediction loc list
    :param uuus_list: unknown loc list
    :param verbose: whehter to talk
    :return: a clean loc list. if there are uuus or uuus cut a loc in the prediction to under 10 loc list, it is cut out
    """
    to_cancel = []
    if verbose:
        print 'in cleaner'
    for w in loc_list:
        for uuu in uuus_list:
            if verbose:
                print 'examining', w, uuu
            if uuu[0] <= w[0] <= uuu[1] and uuu[0] <= w[1] <= uuu[1]:
                if verbose:
                    print 'canceling due to uuus', w, uuu
                to_cancel.append(w)
            elif w[0] <= uuu[0] <= w[1] and w[0] <= uuu[1] <= w[1] and uuu[1]-uuu[0] > 1:
                if verbose:
                    print 'canceling due to overlap', w, uuu
                to_cancel.append(w)
            elif uuu[0] <= w[0] <= uuu[1] or uuu[0] <= w[1] <= uuu[1]:
                new_range = [a for a in range(w[0], w[1]+1) if a not in range(uuu[0], uuu[1]+1)]
                new_loc = [new_range[0], new_range[-1]]
                if new_loc[1]-new_loc[0] <= 10:
                    if verbose:
                        print 'canceling due to length', new_loc
                    to_cancel.append(w)
                    continue
                w = new_loc
    return [a for a in loc_list if a not in to_cancel]


def qok_pdbtm_opm(rost_entry, loc_list, signal_peptide, verbose=False):
    """
    :param rost_entry: Rost database entry about query
    :param loc_list: [[begin, end]...] for prediction
    :return: True if either pdbtm or opm agree by Qok, True/False by overlapM
    """
    pdbtm_uuus = pdbtm_opm_uuus_list(rost_entry['pdbtm'], signal_peptide, verbose=verbose)
    opm_uuus = pdbtm_opm_uuus_list(rost_entry['opm'], signal_peptide, verbose=verbose)

    cleaned_loc_list = pred_loc_list_clean_by_uuus(loc_list, pdbtm_uuus+opm_uuus, verbose=verbose)

    pdbtm_loc_list = pdbtm_opm_loc_list(rost_entry['pdbtm'], signal_peptide)
    opm_loc_list = pdbtm_opm_loc_list(rost_entry['opm'], signal_peptide)

    clean_pdbtm_loc_list = pred_loc_list_clean_by_uuus(pdbtm_loc_list, pdbtm_uuus)
    clean_opm_loc_list = pred_loc_list_clean_by_uuus(opm_loc_list, opm_uuus)

    if verbose:
        print 'checker'
        print 'pdbtm loc list %r, uuu list %r' % (pdbtm_loc_list, pdbtm_uuus)
        print 'opm loc list %r, uuu list %r' % (opm_loc_list, opm_uuus)
        print 'clean pred', cleaned_loc_list

    if len(clean_pdbtm_loc_list) != len(cleaned_loc_list) and len(clean_opm_loc_list) != len(cleaned_loc_list):
        if verbose:
            print 'failing on tm num. prediction %i, pdbtm %i, opm %i' % \
                  (len(cleaned_loc_list), len(pdbtm_loc_list), len(opm_loc_list))
        return False, False

    pdbtm_qok = qok_protein(pdbtm_loc_list, cleaned_loc_list)
    opm_qok = qok_protein(opm_loc_list, cleaned_loc_list)

    if verbose:
        print 'pdbtm loc list' % pdbtm_loc_list
        print 'opm loc list' % opm_loc_list

    pdbtm_overlapm = overlapM_loc_lists(pdbtm_loc_list, cleaned_loc_list, verbose=verbose)
    opm_overlapm = overlapM_loc_lists(opm_loc_list, cleaned_loc_list, verbose=verbose)

    if verbose:
        print 'pdbtm result', pdbtm_overlapm
        print 'opm result', opm_overlapm

    return any([pdbtm_qok, opm_qok]), any([pdbtm_overlapm, opm_overlapm])


def overlapM_loc_lists(observed_loc_list, predicted_loc_list, M=10, verbose=False):
    """
    :param observed_loc_list: [[begin, end]...]
    :param predicted_loc_list: [[begin, end]...]
    :param M: overlap threshold
    :return: True/False if overlap for every windw is over/equal M
    >>> overlapM_loc_lists([[10, 20], [30, 40]], [[10, 20], [30, 40]])
    True
    >>> overlapM_loc_lists([[10, 20]], [[10, 20]])
    True
    >>> overlapM_loc_lists([[10, 20]], [[8, 22]])
    True
    >>> overlapM_loc_lists([[10, 20]], [[15, 20]])
    False
    >>> overlapM_loc_lists([[10, 20]], [[10, 200]])
    True
    >>> overlapM_loc_lists([[10, 20]], [[10, 18]])
    False
    """
    if verbose:
        print 'overlapM_loc_lists verbose!!!'
    for obs_w, prd_w in zip(observed_loc_list, predicted_loc_list):
        if verbose:
            print 'looking at observed win %r and predicted win %r' % (obs_w, prd_w)
        obs_w_range = range(obs_w[0], obs_w[1]+1)
        prd_w_range = range(prd_w[0], prd_w[1]+1)
        overlap = len([a for a in prd_w_range if a in obs_w_range])
        if overlap < M:
            if verbose:
                print 'failed at observed win %r and predicted win %r' % (obs_w, prd_w)
            return False
    return True


def wgp_to_loc_list(wgp, signal_peptide):
    """
    :param wgp: wgp
    :param signal_peptide: int, signal peptide size
    :return: [[begin, end]...] of wgp
    >>> from WinGrade import WinGrade, WinGradePath
    >>> w1 = WinGrade(0, 10, 'fwd', 'A', grade=1, length_element=1, charges=1)
    >>> w2 = WinGrade(20, 30, 'rev', 'B', grade=2, length_element=2, charges=2)
    >>> w3 = WinGrade(30, 40, 'fwd', 'C', grade=3, length_element=3, charges=3)
    >>> wgp_to_loc_list(WinGradePath([w1, w2, w3]), 0)
    [[0, 10], [20, 30], [30, 40]]
    >>> wgp_to_loc_list(WinGradePath([w1, w2, w3]), 15)
    [[20, 30], [30, 40]]
    """
    result = []
    for w in wgp.path:
        if w.begin >= signal_peptide:
            result.append([w.begin, w.end])
    return result


def qok_protein(observed_list, predicted_list):
    """
    :param observed_ts: the observed topo string as TMH locations list
    :param predicted_ts: the predicted topo string as TMH locations list
    :return: True/Flase if stands by Qok. both ends within 5 residues from observe,
            and all and only observed TMHs are predicted
    >>> qok_protein([[10, 20], [30, 40]], [[10, 20], [30, 40]])
    True
    >>> qok_protein([[10, 20], [30, 40]], [[10, 20], [30, 50]])
    False
    >>> qok_protein([[10, 20], [30, 40]], [[10, 26], [30, 40]])
    False
    """
    for obs_w, prd_w in zip(observed_list, predicted_list):
        if not qok_helix(obs_w, prd_w):
            return False
    return True


def qok_helix(observed_w, predicted_w):
    """
    :param observed_w: observed window in [begin, end] format
    :param predicted_w: predicted window in [begin, end] format
    :return: True or False by Qok
    >>> qok_helix([10, 20], [10, 20])
    True
    >>> qok_helix([10, 20], [5, 15])
    True
    >>> qok_helix([10, 20], [25, 35])
    False
    >>> qok_helix([10, 20], [14, 35])
    False
    """
    return observed_w[0]-5 <= predicted_w[0] <= observed_w[0]+5 and \
           observed_w[1]-5 <= predicted_w[1] <= observed_w[1]+5


def ts_loc_list(ts, signal_peptide):
    """
    :param ts: a topo string
    :return: [[begin, end]...]
    >>> ts = 'iiiiimmmmmooooommmmmiiiii'
    >>> ts_loc_list(ts, 0)
    [[5, 9], [15, 19]]
    >>> ts_loc_list(ts, 11)
    [[15, 19]]
    """
    hhh = re.compile('[Mm]*')
    hhh_list = [[a.span()[0], a.span()[1]-1] for a in re.finditer(hhh, ts) if a.span()[1]-a.span()[0] > 1
                and a.span()[0] >= signal_peptide]
    return hhh_list


def pdbtm_opm_uuus_list(ts, signal_peptide, verbose=False):
    """
    :param ts: a topo string from either PDBTM or OPM
    :param signal_peptide: int of signal peptide length
    :return: [[begin, end]... ]
    >>> ts = 'iiiiihhhhhoooooHHHHHiiiii'
    >>> pdbtm_opm_loc_list(ts, 0)
    [[5, 9], [15, 19]]
    >>> pdbtm_opm_loc_list(ts, 11)
    [[15, 19]]
    >>> pdbtm_opm_uuus_list('uuuaaabbb', 0)
    [[0, 5]]
    """
    uuu = re.compile(r'([uULl])(\1*)')
    uuu_list = [[a.span()[0], a.span()[1]-1] for a in re.finditer(uuu, ts) if a.span()[1]-a.span()[0] > 1]
    if verbose:
        print 'uuu finder found', uuu_list
    return uuu_list


def pdbtm_opm_loc_list(ts, signal_peptide):
    """
    :param ts: a topo string from either PDBTM or OPM
    :param signal_peptide: int of signal peptide length
    :return: [[begin, end]... ]
    >>> ts = 'iiiiihhhhhoooooHHHHHiiiii'
    >>> pdbtm_opm_loc_list(ts, 0)
    [[5, 9], [15, 19]]
    >>> pdbtm_opm_loc_list(ts, 11)
    [[15, 19]]
    """
    hhh = re.compile(r'([hH])(\1*)')
    hhh_list = [[a.span()[0], a.span()[1]-1] for a in re.finditer(hhh, ts) if a.span()[1]-a.span()[0] > 1
                and a.span()[0] >= signal_peptide]
    return hhh_list


def rost_new_old():
    """
    :return: uniprot(lower case): new/old in the Rost database
    """
    result = {}
    with open('/home/labs/fleishman/jonathaw/membrane_prediction_DBs/rost_old_names.txt') as fin:
        for a in fin.read().split('\n'):
            result[a.lower().rstrip()] = 'old'
    with open('/home/labs/fleishman/jonathaw/membrane_prediction_DBs/rost_new_names.txt') as fin:
        for a in fin.read().split('\n'):
            result[a.lower().rstrip()] = 'new'
    return result


def observed_lengths(args):
    """
    :param args: run arguments
    :return: draws histograms of lengths of either PDBTM or OPM win lengths and inter helix loops
    """
    rost_db = parse_rostlab_db()
    lengths = []
    between_helices_lengths = []

    for k, v in rost_db.items():
        topc = spc_parser(k)
        signal_peptide = topc['topcons'].count('S') + topc['topcons'].count('s')
        obs_loc_list = pdbtm_opm_loc_list(v[args['data_base']], signal_peptide)

        for i, w in enumerate(obs_loc_list):
            lengths.append(w[1]-w[0])

            if i+1 in range(0, len(obs_loc_list)):
                between_helices_lengths.append(obs_loc_list[i+1][0] - w[1])

    plt.hist(lengths, 30, normed=1, facecolor='green', alpha=0.75)
    plt.hist(between_helices_lengths, 100, normed=1, facecolor='blue', alpha=0.5)
    plt.xlabel('Window lengths in %s dataset' % args['data_base'])
    plt.ylabel('Frequency')
    plt.xlim([0, 100])
    plt.grid(True)
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', default='rost')
    parser.add_argument('-data_base', default='pdbtm')

    args = vars(parser.parse_args())

    if args['mode'] == 'rost':
        main_rost()

    elif args['mode'] == 'observed_lengths':
        observed_lengths(args)

    else:
        print 'no mode found'