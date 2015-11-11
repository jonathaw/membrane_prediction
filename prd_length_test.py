#!/usr/bin/env python2.7
import os
from TMpredict_WinGrade import parse_rostlab_db
rost_db = parse_rostlab_db()

in_path = '/home/labs/fleishman/elazara/length_21/w_0_with_MSA/correct_prds'
prds = [a for a in os.listdir(in_path) if a[-4:] == '.prd' and '_msa' not in a]
result = {a: [] for a in range(10)}
for p in prds:
    lengths = []
    name = p.split('.')[0]

    cont = open(p, 'r').read()

    best_path_txt = cont.split('{')[1].split('}')[0]
    for l in best_path_txt.split('\n'):
        l = l.replace('[', '').replace(']', '')
        s = l.split()
        if not s:
            continue
        length = int(s[2])-int(s[0])
        lengths.append(length)
    # if 21 in lengths and 30 in lengths:
    #     print 'found', name, lengths
    if len(lengths) == 6:
        print 'I got ot 6 ', name, rost_db[name]['pdb']
    result[max(lengths)-min(lengths)].append(name)
print result

