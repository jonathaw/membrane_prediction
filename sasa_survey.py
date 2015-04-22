from Bio.PDB import *


def read_rostdb_entries():
    f = open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/datasetEntryList.tsv', 'r')
    cont = f.read().split('\n')
    result = {}
    for line in cont:
        split = line.split()
        if split[1].split(':')[0] not in result.keys():
            result[split[1].split(':')[0]] = [split[1].split(':')[1]]
        result[split[1].split(':')[0]].append(split[1].split(':')[1])
    return result


def download_and_get_chains():
    from Bio.PDB import PDBParser, PDBIO
    failed = []
    pdbs_dict = read_rostdb_entries()
    io = PDBIO()
    pdbl = PDBList()
    for pdb_e, chains in pdbs_dict.items():
        for chain_e in chains:
            try:
                pdbl.retrieve_pdb_file(pdb_e, pdir='./')
                pdb = PDBParser().get_structure(pdb_e, 'pdb'+pdb_e.lower()+'.ent')
                for chain in pdb.get_chains():
                    if chain.get_id() == chain_e:
                        io.set_structure(chain)
                        io.save(pdb.get_id() + '_' + chain.get_id() + '.pdb')
            except:
                failed.append((pdb_e, chain_e))
    print "failures:", failed


def naccess_parser(name):
    with open('/home/labs/fleishman/jonathaw/membrane_prediction_DBs/sasa_survey/'+name+'.rsa', 'r') as f:
    # with open('/home/labs/fleishman/jonathaw/membrane_prediction_DBs/sasa_survey/single_bio_unit/pdb'+name+'.rsa', 'r') as f:
        cont = f.read().split('\n')
    result = []
    three_to_one = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}
    for line in cont:
        split = line.split()
        if len(split) == 0: continue
        if split[0] == 'RES':
            if split[1] not in three_to_one.keys(): continue
            if len(split[2]) > 1:
                result.append({three_to_one[split[1]]: float(split[4])})
            result.append({three_to_one[split[1]]: float(split[5])})
    return result


def helix_rsa(seq, naccess):
    result = []
    naccess_seq = ''.join(a.keys()[0] for a in naccess)
    loc = (naccess_seq.find(seq), naccess_seq.find(seq)+len(seq))
    for i in range(loc[0], loc[1]):
        result.append(naccess[i].values()[0])
    return sum(result) / len(result)


def observed_TMHs_hp_burial():
    import os
    import re
    import matplotlib.pyplot as plt
    import numpy as np
    from WinGrade import grade_segment
    from TMpredict_WinGrade import MakeHydrophobicityGrade
    from hp_moment_over_pred import prd_parser, split_to_helices


    polyval = MakeHydrophobicityGrade()
    prd_list = [x for x in os.listdir(os.getcwd()) if re.match('.*\.prd', x)]
    maxi = 0
    mini = 100000
    results = {'pdbtm': [], 'opm': [], 'correct': [], 'incorrect': [], 'pred': []}

    for prd_file in prd_list:
        prd_pars = prd_parser(prd_file)
        prd_nacces = naccess_parser(prd_pars['pdb'] +'_' +prd_pars['chain'].upper())
        print 'aaa', prd_nacces, prd_pars
        pdbtm_helices = split_to_helices(prd_pars['seq'], prd_pars['obs_pdbtm'])
        pred_helices = split_to_helices(prd_pars['seq'], prd_pars['pred'])
        for helix in pdbtm_helices:
            rsa = helix_rsa(helix[1], prd_nacces)
            hp = grade_segment(helix[1], polyval)
            results['pdbtm'].append((hp, rsa))
            if rsa > maxi:
                maxi = rsa
                max_rsa = helix
                max_pdb = prd_pars['pdb']
                max_chian = prd_pars['chain']
                max_grade = hp
            if rsa < mini:
                mini = rsa
                mini_rsa = helix
                mini_pdb = prd_pars['pdb']
                mini_chian = prd_pars['chain']
                mini_grade = hp
        for helix in pred_helices:
            rsa = helix_rsa(helix[1], prd_nacces)
            hp = grade_segment(helix[1], polyval)
            results['pred'].append((hp, rsa))
            if pred_over_m_obs(helix, pdbtm_helices):
                results['correct'].append((hp, rsa))
            else:
                results['incorrect'].append((hp, rsa))
    print 'max', maxi, max_rsa, max_pdb, max_chian, max_grade
    print 'min', mini, mini_rsa, mini_pdb, mini_chian, mini_grade
    print 'total number of data points in pdbtm', len(results['pdbtm'])
    print 'total number of data points in pred', len(results['pred'])

    axis_labels = ['% Burial', 'dG']
    labels = ['<5', '5-25', '25-45', '>45']
    ranges = [(0, 5), (5, 25), (25, 45), (45, np.inf)]
    labels = ranges2labels(ranges)
    # pdbtm_data_box = [[a[0] for a in pdbtm_data if a[1]<5.], [a[0] for a in pdbtm_data if 25.>a[1]>5.],
    #         [a[0] for a in pdbtm_data if a[1]>25]]
    # pred_data_box = [[a[0] for a in pred_data if a[1]<5.], [a[0] for a in pred_data if 25.>a[1]>5.],
    #         [a[0] for a in pred_data if a[1]>25]]

    pdbtm_data_box = break2bons(results['pdbtm'], ranges)
    pred_data_box = break2bons(results['pred'], ranges)
    correct_data_box = break2bons(results['correct'], ranges)
    incorrect_data_box = break2bons(results['incorrect'], ranges)


    boxplot_with_scatter(pdbtm_data_box, labels, 221, 'pdbtm', axis_labels)
    boxplot_with_scatter(pred_data_box, labels, 222, 'prediction', axis_labels)
    boxplot_with_scatter(correct_data_box, labels, 223, 'correct', axis_labels)
    boxplot_with_scatter(incorrect_data_box, labels, 224, 'incorrect', axis_labels)

    plt.show()


def pred_over_m_obs(pred, obs_s, m=10):
    for obs in obs_s:
        if tuples_overlap_m(pred[0], obs[0], m):
            return True
    return False


def tuples_overlap_m(tup1, tup2, m=10):
    return len([a for a in range(tup1[0], tup1[1]+1) if a in range(tup2[0], tup2[1]+1)]) >= m


def break2bons(data, ranges):
    result = []
    for r in ranges:
        result.append([a[0] for a in data if r[0] < a[1] < r[1]])
    return result


def ranges2labels(ranges):
    return [str(r[0])+'-'+str(r[1]) for r in ranges]


def boxplot_with_scatter(data, labels, fig_num, title='default', axis_labels=('x axis', 'y axis')):
    import numpy as np
    import matplotlib.pyplot as plt
    from random import uniform
    positions = np.arange(len(data))
    plt.subplot(fig_num)
    plt.boxplot(data, labels=labels, positions=positions)
    min_y = min([i for o in data for i in o])
    for list_, i in zip(data, positions):
        x = [i+uniform(-0.1, 0.1) for a in list_]
        plt.scatter(x, list_, marker='.', s=8, c='grey')
        plt.text(i, min_y-3, len(list_))
    plt.title(title)
    plt.xlabel(axis_labels[0])
    plt.ylabel(axis_labels[1])


if __name__ == '__main__':
    # download_and_get_chains()
    observed_TMHs_hp_burial()