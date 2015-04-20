"""
a script to create prd_sum_ files for ROC directories, under mode summerize.
under mode compare compares ROC files
"""
import matplotlib.pyplot as plt
import numpy as np


def prd_summarize(dir_name):
    import os
    import re
    import sys
    from result_slicer import pred_reader
    prd_file_list = [x for x in os.listdir(dir_name) if re.match('.*\.prd', x)]
    if prd_file_list == []:
        sys.exit()
    results = {'Qok': 0, 'over10': 0, 'topo': 0, 'proteins': 0}
    for prd_file in prd_file_list:
        prd = pred_reader(os.getcwd()+'/'+dir_name, prd_file)
        results['Qok'] += 1 if prd['Qok'] else 0
        results['over10'] += 1 if prd['over10'] else 0
        results['topo'] += 1 if prd['topo'] else 0
        results['proteins'] += 1
    qok_percent = percent(results['Qok'], results['proteins'])
    over10_percent = percent(results['over10'], results['proteins'])
    topo_percent = percent(results['topo'], results['proteins'])

    if args['polynom']:
        name_split = dir_name.split('/')[-1].split('_')
        name = '_'.join(name_split[1:])
        c0 = float(name_split[1])
        c1 = float(name_split[2])
        c2 = float(name_split[3])
        c3 = float(name_split[4])
    else:
        name = str(prd['hp_threshold'])+'_'+str(prd['min_length'])+'_'+str(prd['psi_helix'])+'_'+str(prd['psi_res_num'])
    o = open(os.getcwd()+'/prd_sum_'+name+'.roc', 'wr+')
    o.write('hp_threshold %f\n' % prd['hp_threshold'])
    o.write('min_length %i\n' % prd['min_length'])
    o.write('psi_helix %f\n' % prd['psi_helix'])
    o.write('psi_res_nume %i\n\n' % prd['psi_res_num'])
    if args['polynom']:
        o.write('c0 %f\n' % c0)
        o.write('c1 %f\n' % c1)
        o.write('c2 %f\n' % c2)
        o.write('c3 %f\n\n' % c3)
    o.write('# proteins %i\n' % results['proteins'])
    o.write('qok_percent %f\n' % qok_percent)
    o.write('over10_percent %f\n' % over10_percent)
    o.write('topo_percent %f\n' % topo_percent)
    o.close()


def percent(up, down):
    if down == 0:
        return None
    else:
        return float(100 * float(up) / float(down))


def compare_ROC():
    import os
    import re
    roc_results = []
    best_by = {'qok': (0., {}), 'over10': (0., {}), 'topo': (0., {})}
    roc_file_list = [x for x in os.listdir(os.getcwd()) if re.match('prd.*\.roc', x)]
    for roc_file in roc_file_list:
        print roc_file
        roc_res = roc_parser(roc_file)
        print roc_res
        if roc_res == {}: continue
        if roc_res['proteins'] < 180: continue
        roc_results.append(roc_res)
        if roc_res['qok_percent'] > best_by['qok'][0]:
            best_by['qok'] = (roc_res['qok_percent'], roc_res)
        if roc_res['over10_percent'] > best_by['over10'][0]:
            best_by['over10'] = (roc_res['over10_percent'], roc_res)
        if roc_res['topo_percent'] > best_by['topo'][0]:
            best_by['topo'] = (roc_res['topo_percent'], roc_res)
    for k, v in best_by.items():
        print 'best by %s: %f -> %r' % (k, v[0], v[1])
    i = 1
    for some in ['c0', 'c1', 'c2', 'c3']:
        plot_by_something(roc_results, best_by, some, i)
        i += 1
    # plot_by_hp_threshold(roc_results, best_by)
    # plot_by_min_length(roc_results, best_by)
    # plot_by_psi_helix(roc_results, best_by)
    # plot_by_psi_res_num(roc_results, best_by)
    plt.figlegend((qok, over10, topo), ('Qok', 'over10', 'topo'), 'upper left')
    plt.show()


def plot_by_something(roc_results, best_by, some, plot_num):
    global qok, over10, topo
    plt.subplot(220+plot_num)
    qok = plt.scatter([x[some] for x in roc_results], [x['qok_percent'] for x in roc_results], marker='o')
    over10 = plt.scatter([x[some] for x in roc_results], [x['over10_percent'] for x in roc_results], marker='+')
    topo = plt.scatter([x[some] for x in roc_results], [x['topo_percent'] for x in roc_results], marker='*')
    for perc, mark in zip(['qok', 'over10', 'topo'], ['o', '+', '*']):
        plt.scatter(best_by[perc][1][some], best_by[perc][0],marker=mark, s=100)
    plt.xlabel(some)
    plt.ylabel('percent')
    plt.title(some+' Vs. percent')


def plot_by_psi_res_num(roc_results, best_by):
    global qok, over10, topo
    plt.subplot(224)
    qok = plt.scatter([x['psi_res_nume'] for x in roc_results], [x['qok_percent'] for x in roc_results], marker='o')
    over10 = plt.scatter([x['psi_res_nume'] for x in roc_results], [x['over10_percent'] for x in roc_results], marker='+')
    topo = plt.scatter([x['psi_res_nume'] for x in roc_results], [x['topo_percent'] for x in roc_results], marker='*')
    for perc, mark in zip(['qok', 'over10', 'topo'], ['o', '+', '*']):
        plt.scatter(best_by[perc][1]['psi_res_nume'], best_by[perc][0],marker=mark, s=100)
    plt.xlabel('psi_res_num')
    plt.ylabel('percent')
    plt.title('psi_res_num Vs. percent')
    
    
def plot_by_psi_helix(roc_results, best_by):
    plt.subplot(223)
    plt.scatter([x['psi_helix'] for x in roc_results], [x['qok_percent'] for x in roc_results], marker='o')
    plt.scatter([x['psi_helix'] for x in roc_results], [x['over10_percent'] for x in roc_results], marker='+')
    plt.scatter([x['psi_helix'] for x in roc_results], [x['topo_percent'] for x in roc_results], marker='*')
    for perc, mark in zip(['qok', 'over10', 'topo'], ['o', '+', '*']):
        plt.scatter(best_by[perc][1]['psi_helix'], best_by[perc][0],marker=mark, s=100)
    plt.xlabel('psi_helix')
    plt.ylabel('percent')
    plt.title('psi_helix Vs. percent')
    

def plot_by_min_length(roc_results, best_by):
    plt.subplot(222)
    plt.scatter([x['min_length'] for x in roc_results], [x['qok_percent'] for x in roc_results], marker='o')
    plt.scatter([x['min_length'] for x in roc_results], [x['over10_percent'] for x in roc_results], marker='+')
    plt.scatter([x['min_length'] for x in roc_results], [x['topo_percent'] for x in roc_results], marker='*')
    for perc, mark in zip(['qok', 'over10', 'topo'], ['o', '+', '*']):
        plt.scatter(best_by[perc][1]['min_length'], best_by[perc][0],marker=mark, s=100)
    plt.xlabel('min_length')
    plt.ylabel('percent')
    plt.title('min_length Vs. percent')
    
    
def plot_by_hp_threshold(roc_results, best_by):
    plt.subplot(221)
    plt.scatter([x['hp_threshold'] for x in roc_results], [x['qok_percent'] for x in roc_results], marker='o')
    plt.scatter([x['hp_threshold'] for x in roc_results], [x['over10_percent'] for x in roc_results], marker='+')
    plt.scatter([x['hp_threshold'] for x in roc_results], [x['topo_percent'] for x in roc_results], marker='*')
    for perc, mark in zip(['qok', 'over10', 'topo'], ['o', '+', '*']):
        plt.scatter(best_by[perc][1]['hp_threshold'], best_by[perc][0],marker=mark, s=100)
    plt.xlabel('hp threshold')
    plt.ylabel('percent')
    plt.title('hp threshold Vs. percent')


def roc_parser(file_name):
    result = {}
    with open(file_name, 'r') as f:
        cont = f.read().split('\n')
    for line in cont:
        split = line.split()
        if len(split) == 0:
            continue
        if split[0] == '#':
            result['proteins'] = int(split[2])
        else:
            result[split[0]] = float(split[1])

    return result

if __name__ == '__main__':
    import argparse
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode')
    parser.add_argument('-name')
    parser.add_argument('-polynom', default=True)
    args = vars(parser.parse_args())
    if args['mode'] == 'summerize':
        prd_summarize(args['name'])
    elif args['mode'] == 'compare':
        compare_ROC()