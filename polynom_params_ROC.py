import os
import sys
from operator import itemgetter
import matplotlib.pyplot as plt

file_list = [a for a in os.listdir(os.getcwd()) if a[:4] == 'ROC_']
results = []
for roc in file_list:
    f = open(roc, 'r')
    cont = f.read().split('\n')
    f.close()
    name_split = roc.split('_')
    res = {'c0': name_split[2]}
    for line in cont:
        split = line.split()
        if len(split) < 2:
            continue
        print 'split', split



    # results.append({k: v for k, v in res.items()})
    results.append(res)
    # break
sys.exit()
sort_pdbtm = sorted(results, key=itemgetter('Qok_pdbtm'))
sort_opm = sorted(results, key=itemgetter('Qok_opm'))
print 'best by pdbtm', sort_pdbtm[-1]
print 'worst by pdbtm', sort_pdbtm[0]
print 'best by opm', sort_opm[-1]
print 'worst by opm', sort_opm[0]
best_pdbtm = sort_pdbtm[-1]
best_opm = sort_opm[-1]

print best_pdbtm

f, p = plt.subplots(2, 2)
poses = [(0, 0), (0, 1), (1, 0), (1, 1)]
twins = {}
plts = []
jitters = [0.08, 0.05, 0.01, 0.1]
for differ, pos, delta in zip(['hp_threshold', 'min_length', 'psi_helix', 'psi_res_num'], poses, jitters):
    # delta = (max([a[differ] for a in sort_pdbtm])-min([a[differ] for a in sort_pdbtm]))/100
    # delta = 0.01
    plts.append(p[pos[0], pos[1]].scatter([a[differ]+delta for a in sort_pdbtm], [a['Qok_pdbtm'] for a in sort_pdbtm], color='b', marker='.', label='Qok_pdbtm '+differ))
    plts.append(p[pos[0], pos[1]].scatter([a[differ]-delta for a in sort_opm], [a['Qok_opm'] for a in sort_opm], color='r', marker='.',label='Qok_opm '+differ))
    # p[pos[0], pos[1]].set_ylim([10, 30])

    twins[(pos[0], pos[1])] = p[pos[0], pos[1]].twinx()
    plts.append(twins[(pos[0], pos[1])].scatter([a[differ]+delta for a in sort_pdbtm], [a['topo_percent_pdbtm'] for a in sort_pdbtm], marker='+', color='b', label='Topo '+differ))
    plts.append(twins[(pos[0], pos[1])].scatter([a[differ]-delta for a in sort_opm], [a['topo_percent_opm'] for a in sort_opm], marker='+', color='r', label='Topo' +differ))
    # twins[(pos[0], pos[1])].set_ylim([20, 50])

    plts.append(p[pos[0], pos[1]].scatter(best_pdbtm[differ], best_pdbtm['Qok_pdbtm'], color='k', marker='o', label='best_pdbtm '+differ, s=32))
    plts.append(p[pos[0], pos[1]].scatter(best_opm[differ], best_opm['Qok_opm'], color='deeppink', marker='o', label='best_opm '+differ, s=32))
    p[pos[0], pos[1]].set_title(differ)
    p[pos[0], pos[1]].set_xlabel(differ)
    p[pos[0], pos[1]].set_ylabel('Qok')
    twins[(pos[0], pos[1])].set_ylabel('Topo Correct')
plt.legend((a for a in plts),('Qok pdbtm','Qok opm','Topo pdbtm','Topo opm','best pdbtm', 'best opm'),loc = 'right', bbox_to_anchor = (0,-0.1,1,1),
            bbox_transform = plt.gcf().transFigure)
plt.show()
print 'AAAA'
