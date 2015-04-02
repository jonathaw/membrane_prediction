import os
from operator import itemgetter
import matplotlib.pyplot as plt

file_list = [a for a in os.listdir(os.getcwd()) if a[-4:] == '.roc']

results = []
for roc in file_list:
    f = open(roc, 'r')
    cont = f.read().split('\n')
    f.close()
    res = {}
    for line in cont:
        split = line.split()
        if len(split) < 2:
            continue
        if split[0] == 'hp_threshold':
           res['hp_threshold'] = float(split[1])
        if split[0] == 'min_length':
           res['min_length'] = int(split[1])
        if split[0] == 'psi_helix':
           res['psi_helix'] = float(split[1])
        if split[0] == 'psi_res_nume':
           res['psi_res_num'] = int(split[1])
        if split[0] == 'Q_ok':
           res['Q_ok'] = split[1]
        if split[1] == 'correct':
           res['topo_percent'] = float(split[5])
    results.append({k: v for k, v in res.items()})
    # break

sort_results = sorted(results, key=itemgetter('Q_ok'))
print sort_results[-1]
print sort_results[0]
f, p = plt.subplots(2, 2)
p[0, 0].scatter([a['hp_threshold'] for a in sort_results], [a['Q_ok'] for a in sort_results], color='b', label='Qok')
p[0, 0].scatter([a['hp_threshold'] for a in sort_results], [a['topo_percent'] for a in sort_results], marker='+', color='r', label='Topo')
p[0, 0].set_title('Q_ok Vs. hp_threshold')
p[0, 1].scatter([a['min_length'] for a in sort_results], [a['Q_ok'] for a in sort_results], color='b')
p[0, 1].scatter([a['min_length'] for a in sort_results], [a['topo_percent'] for a in sort_results], marker='+', color='r')
p[0, 1].set_title('Q_ok Vs. min_length')
p[1, 0].scatter([a['psi_helix'] for a in sort_results], [a['Q_ok'] for a in sort_results],color='b')
p[1, 0].scatter([a['psi_helix'] for a in sort_results], [a['topo_percent'] for a in sort_results], marker='+',color='r')
p[1, 0].set_title('Q_ok Vs. psi_helix')
p[1, 1].scatter([a['psi_res_num'] for a in sort_results], [a['Q_ok'] for a in sort_results],color='b')
p[1, 1].scatter([a['psi_res_num'] for a in sort_results], [a['topo_percent'] for a in sort_results], marker='+',color='r')
p[1, 1].set_title('Q_ok Vs. psi_res_num')

plt.show()