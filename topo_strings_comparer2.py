#!/usr/bin/env python2.7
import os
import re
import argparse
import matplotlib.pyplot as plt

from positive_inside_analysis import parse_prd
from topo_strings_comparer import spc_parser
from TMpredict_WinGrade import parse_rostlab_db


def main_rost():
    predictors = ['polyphobius', 'topcons', 'spoctopus', 'philius', 'octopus', 'scampi']
    prd_files = [a for a in os.listdir('./') if '.prd' in a and '_msa' not in a]
    rost_db = parse_rostlab_db()
    new_old = rost_new_old()
    topgraph_none = []

    old_new_totals = {'new': 0, 'old': 0}
    results = {}
    for prd_file in prd_files:
        name = prd_file.split('.')[0].lower()
        wgp = parse_prd(prd_file)[0]
        # print name
        if wgp is None:
            topgraph_none.append(name)
            continue
        topc = spc_parser(name)

        signal_peptide = topc['topcons'].count('S') + topc['topcons'].count('s')

        wgp_loc_list = wgp_to_loc_list(wgp, signal_peptide)
        # print wgp, topc, rost_db[name]
        old_new_totals[new_old[name]] += 1
        tgr_qok, tgr_ovm = qok_pdbtm_opm(rost_db[name], wgp_loc_list, signal_peptide)
        results[name] = {'old_new': new_old[name],
                         'tm_num': len(pdbtm_opm_loc_list(rost_db[name]['pdbtm'], signal_peptide)),
                         'topgraph': {'qok': tgr_qok, 'ovm': tgr_ovm}}
        for predictor in predictors:
            prd_qok, prd_ovm = qok_pdbtm_opm(rost_db[name], ts_loc_list(topc[predictor], signal_peptide), signal_peptide)
            results[name][predictor] = {'qok': prd_qok, 'ovm': prd_ovm}

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


def qok_pdbtm_opm(rost_entry, loc_list, signal_peptide):
    """
    :param rost_entry: Rost database entry about query
    :param loc_list: [[begin, end]...] for prediction
    :return: True if either pdbtm or opm agree by Qok, True/False by overlapM
    """
    pdbtm_loc_list = pdbtm_opm_loc_list(rost_entry['pdbtm'], signal_peptide)
    opm_loc_list = pdbtm_opm_loc_list(rost_entry['opm'], signal_peptide)

    if len(pdbtm_loc_list) != len(loc_list) and len(opm_loc_list) != len(loc_list):
        return False, False

    pdbtm_qok = qok_protein(pdbtm_loc_list, loc_list)
    opm_qok = qok_protein(opm_loc_list, loc_list)

    pdbtm_overlapm = overlapM_loc_lists(pdbtm_loc_list, loc_list)
    opm_overlapm = overlapM_loc_lists(opm_loc_list, loc_list)

    return any([pdbtm_qok, opm_qok]), any([pdbtm_overlapm, opm_overlapm])


def overlapM_loc_lists(observed_loc_list, predicted_loc_list, M=10):
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
    for obs_w, prd_w in zip(observed_loc_list, predicted_loc_list):
        obs_w_range = range(obs_w[0], obs_w[1]+1)
        prd_w_range = range(prd_w[0], prd_w[1]+1)
        overlap = len([a for a in prd_w_range if a in obs_w_range])
        if overlap < M:
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
    hhh = re.compile('[hH]*')
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', default='rost')

    args = vars(parser.parse_args())

    if args['mode'] == 'rost':
        main_rost()