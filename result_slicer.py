# ## a script to run through TMpredict results from the Rostlab DB, and present it's Qok by
### New/Old, Procaryotic/Eukaryotic, and length.
import matplotlib.pyplot as plt
import numpy as np


def main(prd_path):
    import os
    import re

    global have_conf
    table_data = rostlab_table_parser()
    prd_file_list = [x for x in os.listdir(prd_path) if re.match('.*\.prd', x)]
    prd_results_dict = {name.split('.')[0]: pred_reader(prd_path, name) for name in prd_file_list}
    if all(a['conf'] == False for a in prd_results_dict.values()):
        have_conf = False
    else:
        have_conf = True
    plot_euk_prok(table_data, prd_results_dict)
    plot_new_old(table_data, prd_results_dict)
    plot_tm_num(table_data, prd_results_dict)
    plt.figlegend((rects_all_, rects_over10, rects_topo), ('Qok', 'over10', 'topo'), 'lower right')
    plt.suptitle(os.getcwd(), fontsize=14, fontweight='bold')
    plt.show()


def plot_new_old(table_data, prd_results_dict):
    global rects_all_, rects_over10, rects_topo
    new = {'tot_conf': 0, 'Qok_conf': 0, 'tot_all': 0, 'Qok_all': 0, 'over10': 0, 'topo': 0}
    old = {'tot_conf': 0, 'Qok_conf': 0, 'tot_all': 0, 'Qok_all': 0, 'over10': 0, 'topo': 0}
    known_uni = prd_results_dict.keys()
    for uni, table_entry in table_data.items():
        if uni not in known_uni:
            continue
        if table_entry['old_new'] == 'new':
            if prd_results_dict[uni]['conf']:
                new['tot_conf'] += 1
                new['Qok_conf'] += 1 if prd_results_dict[uni]['Qok'] else 0
            new['over10'] += 1 if prd_results_dict[uni]['over10'] else 0
            new['topo'] += 1 if prd_results_dict[uni]['topo'] else 0
            new['tot_all'] += 1
            new['Qok_all'] += 1 if prd_results_dict[uni]['Qok'] else 0
        if table_entry['old_new'] == 'old':
            if prd_results_dict[uni]['conf']:
                old['tot_conf'] += 1
                old['Qok_conf'] += 1 if prd_results_dict[uni]['Qok'] else 0
            old['over10'] += 1 if prd_results_dict[uni]['over10'] else 0
            old['topo'] += 1 if prd_results_dict[uni]['topo'] else 0
            old['tot_all'] += 1
            old['Qok_all'] += 1 if prd_results_dict[uni]['Qok'] else 0
    all_results = {'Qok_conf': new['Qok_conf'] + old['Qok_conf'], 'tot_conf': new['tot_conf'] + old['tot_conf'],
                   'Qok_all': new['Qok_all'] + old['Qok_all'], 'tot_all': new['tot_all'] + old['tot_all'],
                   'over10': new['over10']+old['over10'], 'topo': new['topo']+old['topo']}
    if have_conf:
        conf = (100 * all_results['Qok_conf'] / all_results['tot_conf'], 100 * new['Qok_conf'] / new['tot_conf'],
                100 * old['Qok_conf'] / old['tot_conf'])
    all_ = (100 * all_results['Qok_all'] / all_results['tot_all'], 100 * new['Qok_all'] / new['tot_all'],
            100 * old['Qok_all'] / old['tot_all'])
    all_over10 = (100 * all_results['over10'] / all_results['tot_all'], 100 * new['over10'] / new['tot_all'],
                  100 * old['over10'] / old['tot_all'])
    all_topo = (100 * all_results['topo'] / all_results['tot_all'], 100 * new['topo'] / new['tot_all'],
                  100 * old['topo'] / old['tot_all'])
    N = 3
    width = 0.35
    ind = np.arange(N)
    plt.subplot(221)

    if have_conf:
        rects_conf = plt.bar(ind, conf, width/2, color='r')
    rects_all_ = plt.bar(ind + width/2, all_, width/2, color='b')
    rects_over10 = plt.bar(ind, all_over10, width/2, color='k')
    rects_topo = plt.bar(ind-width/2, all_topo, width/2, color='g')

    plt.ylim((0, max(all_over10)+5))
    plt.ylabel('Qok')
    plt.title('Qok New Vs. Old')
    plt.xticks(ind + width, ['Full', 'New', 'Old'])

    if have_conf:
        autolabel(rects_conf)
    autolabel(rects_all_)
    autolabel(rects_over10)
    autolabel(rects_topo)


def plot_euk_prok(table_data, prd_results_dict):
    euk_ = {'tot_conf': 0, 'Qok_conf': 0, 'tot_all': 0, 'Qok_all': 0, 'over10': 0, 'topo': 0}
    prok = {'tot_conf': 0, 'Qok_conf': 0, 'tot_all': 0, 'Qok_all': 0, 'over10': 0, 'topo': 0}
    known_uni = prd_results_dict.keys()
    for uni, table_entry in table_data.items():
        if uni not in known_uni:
            continue
        if table_entry['euk_prok'] == 'euk':
            if prd_results_dict[uni]['conf']:
                euk_['tot_conf'] += 1
                euk_['Qok_conf'] += 1 if prd_results_dict[uni]['Qok'] else 0
            euk_['tot_all'] += 1
            euk_['Qok_all'] += 1 if prd_results_dict[uni]['Qok'] else 0
            euk_['over10'] += 1 if prd_results_dict[uni]['over10'] else 0
            euk_['topo'] += 1 if prd_results_dict[uni]['topo'] else 0
        if table_entry['euk_prok'] == 'prok':
            if prd_results_dict[uni]['conf']:
                prok['tot_conf'] += 1
                prok['Qok_conf'] += 1 if prd_results_dict[uni]['Qok'] else 0
            prok['tot_all'] += 1
            prok['Qok_all'] += 1 if prd_results_dict[uni]['Qok'] else 0
            prok['over10'] += 1 if prd_results_dict[uni]['over10'] else 0
            prok['topo'] += 1 if prd_results_dict[uni]['topo'] else 0
    if have_conf:
        conf = (100 * euk_['Qok_conf'] / euk_['tot_conf'], 100 * prok['Qok_conf'] / prok['tot_conf'])
    all_ = (100 * euk_['Qok_all'] / euk_['tot_all'], 100 * prok['Qok_all'] / prok['tot_all'])
    over10 = (100 * euk_['over10'] / euk_['tot_all'], 100 * prok['over10'] / prok['tot_all'])
    topo = (100 * euk_['topo'] / euk_['tot_all'], 100 * prok['topo'] / prok['tot_all'])

    N = 2
    width = 0.35
    ind = np.arange(N)
    plt.subplot(222)

    if have_conf:
        rects_conf = plt.bar(ind, conf, width/2, color='r')
    rects_all_ = plt.bar(ind + width/2, all_, width/2, color='b')
    rects_over10 = plt.bar(ind, over10, width/2, color='k')
    rects_topo = plt.bar(ind-width/2, topo, width/2, color='g')

    plt.ylabel('percent')
    plt.title('Eukaryotes Vs. Prokaryotes')
    plt.xticks(ind + width, ['Eukaryotes', 'Prokaryotes'])
    plt.ylim((0, max(over10) + 5.))
    if have_conf:
        autolabel(rects_conf)

    autolabel(rects_all_)
    autolabel(rects_over10)
    autolabel(rects_topo)


def plot_tm_num(table_data, prd_results_dict):
    tmh_1 = {'tot_conf': 0, 'Qok_conf': 0, 'tot_all': 0, 'Qok_all': 0, 'over10': 0, 'topo': 0}
    tmh_2_5 = {'tot_conf': 0, 'Qok_conf': 0, 'tot_all': 0, 'Qok_all': 0, 'over10': 0, 'topo': 0}
    tmh_5 = {'tot_conf': 0, 'Qok_conf': 0, 'tot_all': 0, 'Qok_all': 0, 'over10': 0, 'topo': 0}
    known_uni = prd_results_dict.keys()
    for uni, table_entry in table_data.items():
        if uni not in known_uni:
            continue
        if table_entry['tm_pdbtm'] == 1:
            if prd_results_dict[uni]['conf']:
                tmh_1['tot_conf'] += 1
                tmh_1['Qok_conf'] += 1 if prd_results_dict[uni]['Qok'] else 0
            tmh_1['tot_all'] += 1
            tmh_1['Qok_all'] += 1 if prd_results_dict[uni]['Qok'] else 0
            tmh_1['over10'] += 1 if prd_results_dict[uni]['over10'] else 0
            tmh_1['topo'] += 1 if prd_results_dict[uni]['topo'] else 0
        if 1 < table_entry['tm_pdbtm'] <= 5:
            if prd_results_dict[uni]['conf']:
                tmh_2_5['tot_conf'] += 1
                tmh_2_5['Qok_conf'] += 1 if prd_results_dict[uni]['Qok'] else 0
            tmh_2_5['tot_all'] += 1
            tmh_2_5['Qok_all'] += 1 if prd_results_dict[uni]['Qok'] else 0
            tmh_2_5['over10'] += 1 if prd_results_dict[uni]['over10'] else 0
            tmh_2_5['topo'] += 1 if prd_results_dict[uni]['topo'] else 0
        if table_entry['tm_pdbtm'] > 5:
            if prd_results_dict[uni]['conf']:
                tmh_5['tot_conf'] += 1
                tmh_5['Qok_conf'] += 1 if prd_results_dict[uni]['Qok'] else 0
            tmh_5['tot_all'] += 1
            tmh_5['Qok_all'] += 1 if prd_results_dict[uni]['Qok'] else 0
            tmh_5['over10'] += 1 if prd_results_dict[uni]['over10'] else 0
            tmh_5['topo'] += 1 if prd_results_dict[uni]['topo'] else 0
    if have_conf:
        conf = (100 * tmh_1['Qok_conf'] / tmh_1['tot_conf'], 100 * tmh_2_5['Qok_conf'] / tmh_2_5['tot_conf'],
                100 * tmh_5['Qok_conf'] / tmh_5['tot_conf'])
    all_ = (100 * tmh_1['Qok_all'] / tmh_1['tot_all'], 100 * tmh_2_5['Qok_all'] / tmh_2_5['tot_all'],
            100 * tmh_5['Qok_all'] / tmh_5['tot_all'])
    all_over10 = (100 * tmh_1['over10'] / tmh_1['tot_all'], 100 * tmh_2_5['over10'] / tmh_2_5['tot_all'],
            100 * tmh_5['over10'] / tmh_5['tot_all'])
    all_topo = (100 * tmh_1['topo'] / tmh_1['tot_all'], 100 * tmh_2_5['topo'] / tmh_2_5['tot_all'],
            100 * tmh_5['topo'] / tmh_5['tot_all'])

    N = 3
    width = 0.35
    ind = np.arange(N)
    plt.subplot(223)
    if have_conf:
        rects_conf = plt.bar(ind, conf, width/2, color='r')
    rects_all_ = plt.bar(ind + width/2, all_, width/2, color='b')
    rects_over10 = plt.bar(ind, all_over10, width/2, color='k')
    rects_topo = plt.bar(ind-width/2, all_topo, width/2, color='g')

    plt.ylabel('%')
    plt.title('#TMH')
    plt.xticks(ind + width, ['1', '2-5', '5'])
    plt.ylim((0, max(all_over10) + 5.))

    if have_conf:
        autolabel(rects_conf)
    autolabel(rects_all_)
    autolabel(rects_over10)
    autolabel(rects_topo)


def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        plt.text(rect.get_x() + rect.get_width() / 2., 1.05 * height, '%d' % int(height),
                 ha='center', va='bottom')


def rostlab_table_parser():
    results = {}
    with open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/datasetEntryList.tsv') as f:
        cont = f.read().split('\n')
    for line in cont:
        split = line.lower().split()
        results[split[0]] = {'uniprot': split[0], 'pdb': split[1].split(':')[0], 'chain': split[1].split(':')[1],
                             'old_new': split[2], 'euk_prok': split[3], 'tm_opm': int(split[4]),
                             'tm_pdbtm': int(split[5])}
    return results


def pred_reader(prd_path, file_name):
    result = {}
    with open(prd_path + '/' + file_name, 'r') as f:
        cont = f.read().split('\n')
    for line in cont:
        split = line.split()
        if len(split) > 1:
            if split[0] == 'uniprot':
                result['uniprot'] = split[2]
            elif split[0] == 'Results':
                tech = split[3]
            elif split[1] == 'best':
                result['best_score'] = float(split[3])
            elif split[1] == 'second':
                result['second_score'] = float(split[4])
            elif split[0] == 'Qok':
                result['Qok_' + tech] = True if split[3] == 'True' else False
            elif split[0] == 'overlap10':
                result['over10'] = True if split[3] == 'True' else False
            elif split[0] == 'hp_threshold':
                result['hp_threshold'] = float(split[1])
            elif split[0] == 'min_length':
                result['min_length'] = int(split[1])
            elif split[0] == 'psi_helix':
                result['psi_helix'] = float(split[1])
            elif split[0] == 'psi_res_nume':
                result['psi_res_num'] = float(split[1])
            elif split[0] == 'Topo':
                result[tech + 'topo'] = True if split[2] == 'True' else False
    result['Qok'] = result['Qok_pdbtm'] or result['Qok_opm']
    result['topo'] = result['pdbtmtopo'] or result['opmtopo']
    result['best_minus_second'] = result['best_score'] - result['second_score']
    result['conf'] = True if abs(result['best_minus_second']) >= 2.5 else False
    return result


if __name__ == '__main__':
    import os

    main(os.getcwd())