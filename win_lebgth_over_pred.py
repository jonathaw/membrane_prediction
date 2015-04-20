"""
this script is for going through prediction results, determining false and true TMHs, and testing their hp moment.
trying to determine a threshold of hp moment for which to throw away potential windows
"""


def main():
    import os
    import re
    from WinGrade import hp_moment
    from TMpredict_WinGrade import MakeHydrophobicityGrade
    import matplotlib.pyplot as plt
    import numpy as np
    from random import uniform

    polyval = MakeHydrophobicityGrade()
    prd_list = [x for x in os.listdir(os.getcwd()) if re.match('.*\.prd', x)]
    correct_moments_over = []
    incorrect_moments_over = []
    correct_moments_qok = []
    incorrect_moments_qok = []
    obs_moments_pdbtm = []
    obs_moments_opm = []
    max_mom = 0
    min_mom = 1000
    lengths = {'co_qok': [], 'in_qok': [], 'co_over': [], 'in_over': [], 'pdbtm': [], 'opm': []}
    c1 = np.arange(-10, 10, 1)
    c2 = np.arange(-1, 1, 0.1)
    c3 = np.arange(-0.1, 0.1, 0.01)
    for one in c1:
        for two in c2:
            for three in c3:
                results = {(one, two, three): {'pass': 0, 'fail': 0}}

    for prd_file in prd_list:
        prd = prd_parser(prd_file)
        pred_helices = split_to_helices(prd['seq'], prd['pred'])
        pdbtm_helices = split_to_helices(prd['seq'], prd['obs_pdbtm'])
        opm_helices = split_to_helices(prd['seq'], prd['obs_opm'])
        # print prd['seq']
        # print prd['obs_pdbtm']
        # print prd['obs_opm']
        # print prd['pred']
        for pred_h in pred_helices:
            qok = 'not_found'
            over = 'not_found'
            for obs_h in pdbtm_helices:
                if tuples_overlap_m(pred_h[0], obs_h[0]):
                    correct_moments_over.append(hp_moment(pred_h[1], polyval=polyval, poly_param={'c0': 1}))
                    lengths['co_over'].append(pred_h[0][1]-pred_h[0][0])
                    over = 'found'
                if tuples_Qok(pred_h[0], obs_h[0]):
                    correct_moments_qok.append(hp_moment(pred_h[1], polyval=polyval, poly_param={'c0': 1}))
                    lengths['co_qok'].append(pred_h[0][1]-pred_h[0][0])
                    qok = 'found'
            if over == 'not_found':
                lengths['in_over'].append(pred_h[0][1]-pred_h[0][0])
                incorrect_moments_over.append(hp_moment(pred_h[1], polyval=polyval, poly_param={'c0': 1}))
            if qok == 'not_found':
                lengths['in_qok'].append(pred_h[0][1]-pred_h[0][0])
                incorrect_moments_qok.append(hp_moment(pred_h[1], polyval=polyval, poly_param={'c0': 1}))

        for obs_h in pdbtm_helices:
            lengths['pdbtm'].append(pred_h[0][1]-pred_h[0][0])
            obs_moments_pdbtm.append(hp_moment(obs_h[1], polyval, {'c0': 1}))
        for obs_h in opm_helices:
            lengths['opm'].append(pred_h[0][1]-pred_h[0][0])
            obs_moments_opm.append(hp_moment(obs_h[1], polyval, {'c0': 1}))

    for k, v in lengths.items():
        print k, max(v), min(v), sum(v)/len(v), np.percentile(v, 25), np.percentile(v, 50), np.percentile(v, 75)
    positions = np.arange(6)
    plt.figure()
    plt.boxplot(lengths.values(), labels=lengths.keys(), showmeans=True, meanline=True, positions=positions)
    plt.ylim((16, 27))
    for list_, i in zip(lengths.values(), positions):
        x = [i+uniform(-0.1, 0.1) for a in list_]
        plt.scatter(x, list_, marker='.', s=8, c='grey')
    # plt.scatter([1.1]*len(correct_moments_over), [y for y in correct_moments_over])
    # plt.scatter([0.1]*len(incorrect_moments_over), [y for y in incorrect_moments_over])
    # plt.scatter([1]*len(correct_moments_qok), [y for y in correct_moments_qok], marker='+', color='r')
    # plt.scatter([0]*len(incorrect_moments_qok), [y for y in incorrect_moments_qok], marker='+', color='r')
    plt.show()


def tuples_overlap_m(tup1, tup2, m=10):
    return len([a for a in range(tup1[0], tup1[1]+1) if a in range(tup2[0], tup2[1]+1)]) >= m


def tuples_Qok(tup1, tup2):
    return len([a for a in range(tup1[0], tup1[1]+1) if a not in range(tup2[0], tup2[1]+1)]) <= 5


def split_to_helices(seq, ts):
    """
    :param seq: AA sequence
    :param ts: topo-string (rostlab format) from either database
    :param psipred: psipred results as list of helix propensity
    :return: tuples of seq and psipred results list for that seq, corresponding to helices in the ts.
    """
    import re
    result = []
    hhh = re.compile('[hH]*')
    helices = [(a.start(), a.end()) for a in hhh.finditer(ts) if a.end()-a.start() > 1]
    for helix in helices:
        result.append(((helix[0], helix[1]), seq[helix[0]: helix[1]]))
    return result


def prd_parser(file_name):
    result = {}
    with open(file_name, 'r') as f:
        cont = f.read().split('\n')
    for line in cont:
        split = line.split()
        if split == []: continue
        if split[0] == 'Results': tech = split[3]
        if split[0] == 'sequence': result['seq'] = split[1]
        if split[0] == 'obs': result['obs_'+tech] = split[2]
        if split[0] == 'pre': result['pred'] = split[2]
        if split[0] == 'PDB': result['pdb'] = split[2]
        if split[0] == 'PDB': result['chain'] = split[3]
    return result


if __name__ == '__main__':
    main()