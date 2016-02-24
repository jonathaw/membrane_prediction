#!/usr/bin/env python2.7
"""
a script to analyse GO profiles of transmembrane proteins
"""
import urllib2
import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
from bioservices import QuickGO

from TMpredict_WinGrade import parse_rostlab_db
from positive_inside_analysis import parse_prd

quickGO = QuickGO(verbose=False)
terms_of_interest = ['receptor', 'transport', 'symport', 'antiport', 'channel', 'anchor']
topgraph_msa_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/GO_topgraph_MSA/'
work_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/GO_topgraph_MSA/analysis/'
vh_files_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/GO_VH_topgraph_MSA/'
vh_work_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/GO_VH_topgraph_MSA/analysis/'


def main_rost():
    rost_db = parse_rostlab_db()
    results = {}
    if not os.path.isfile(work_path+'analysis.obj'):
        print 'creating data'
        for name in rost_db.keys():
            # name = 'P20963'
            # print name
            best_wgp, sec_wgp = parse_prd(topgraph_msa_path+name.lower()+'.prd')
            if best_wgp is None:
                continue
            go_terms = get_GO_for_uniprot(name)
            # print 'for %s found these terms %s' % (name, ' '.join(go_terms))
            total_grade = best_wgp.total_grade
            avg_grade = np.mean([w.grade for w in best_wgp.path])
            ddg_best_sec = best_wgp.total_grade - sec_wgp.total_grade
            results[name] = {'go_terms': list(go_terms), 'total_grade': total_grade, 'avg_grade': avg_grade,
                             'ddg_best_sec': ddg_best_sec}
            if len(go_terms) > 1:
                print '%s has %i go terms [%s]' % (name, len(go_terms), ' '.join(go_terms))
            # break
        with open(work_path+'analysis.obj', 'wb') as fout:
            pickle.dump(results, fout)
    else:
        with open(work_path+'analysis.obj', 'rb') as fin:
            print 'reading data'
            results = pickle.load(fin)
    for param in terms_of_interest:
        print 'found %i proteins at %s' % (len([a for a in results.values() if param in a['go_terms']]), param)
    for i, param in enumerate(['total_grade', 'avg_grade', 'ddg_best_sec']):
        plt.subplot(1, 3, i+1)
        plt.hist([v[param] for v in results.values() if 'receptor' in v['go_terms']], alpha=0.5, label='receptor')
        plt.hist([v[param] for v in results.values() if 'transport' in v['go_terms']], alpha=0.5,
                 label='transport')
        plt.hist([v[param] for v in results.values() if 'symport' in v['go_terms']], alpha=0.5,
                 label='symport')
        plt.hist([v[param] for v in results.values() if 'antiport' in v['go_terms']], alpha=0.5,
                 label='antiport')
        plt.hist([v[param] for v in results.values() if 'channel' in v['go_terms']], alpha=0.5, label='channel')
        # plt.hist([v[param] for v in results.values() if 'anchor' in v['go_terms']], alpha=0.5, label='anchor')
        plt.legend(loc='upper right')
    plt.show()


def main_vh():
    colibri_to_uniprot = parse_colibri2uniprot()
    prd_files = [a for a in os.listdir(vh_files_path) if '.prd' in a]
    results = {}
    if not os.path.isfile(vh_files_path+'vh_analysis.obj'):
        print 'creating data'
        for prd_file in prd_files:
            # name = 'P20963'
            # print name
            name = prd_file.split('.')[0]
            best_wgp, sec_wgp = parse_prd(vh_files_path+name+'.prd')
            if best_wgp is None:
                continue
            go_terms = get_GO_for_uniprot(colibri_to_uniprot[name.lower()])
            # print 'for %s found these terms %s' % (name, ' '.join(go_terms))
            total_grade = best_wgp.total_grade
            avg_grade = np.mean([w.grade for w in best_wgp.path])
            ddg_best_sec = best_wgp.total_grade - sec_wgp.total_grade
            results[name] = {'go_terms': list(go_terms), 'total_grade': total_grade, 'avg_grade': avg_grade,
                             'ddg_best_sec': ddg_best_sec}
            if len(go_terms) > 1:
                print '%s has %i go terms [%s]' % (name, len(go_terms), ' '.join(go_terms))
            # break
        with open(work_path+'vh_analysis.obj', 'wb') as fout:
            pickle.dump(results, fout)
    else:
        with open(vh_files_path+'vh_analysis.obj', 'rb') as fin:
            print 'reading data'
            results = pickle.load(fin)
    for param in terms_of_interest:
        print 'found %i proteins at %s' % (len([a for a in results.values() if param in a['go_terms']]), param)
    for i, param in enumerate(['total_grade', 'avg_grade', 'ddg_best_sec']):
        plt.subplot(1, 3, i+1)
        plt.hist([v[param] for v in results.values() if 'receptor' in v['go_terms']], alpha=0.5, label='receptor',
                 normed=True)
        plt.hist([v[param] for v in results.values() if 'transport' in v['go_terms']], alpha=0.5,
                 label='transport', normed=True)
        # plt.hist([v[param] for v in results.values() if 'symport' in v['go_terms']], alpha=0.5,
        #          label='symport', normed=True)
        # plt.hist([v[param] for v in results.values() if 'antiport' in v['go_terms']], alpha=0.5,
        #          label='antiport', normed=True)
        plt.hist([v[param] for v in results.values() if 'channel' in v['go_terms']], alpha=0.5, label='channel',
                 normed=True)
        # plt.hist([v[param] for v in results.values() if 'anchor' in v['go_terms']], alpha=0.5, label='anchor')
        plt.legend(loc='upper right')
    plt.show()


def get_GO_for_uniprot(uniprot):
    url = 'http://www.uniprot.org/uniprot/%s' % uniprot.upper()
    try:
        data = urllib2.urlopen(url)
    except:
        return set()
    go_names = []
    for line in data: # files are iterable
        if 'GO:' in line:
            GO_id = 'GO:'+line.split('GO:')[1].split('"')[0]
            go_names += [a for a in quickGO.Term(GO_id, frmt='obo').split('\n') if 'name' in a]
    terms = set()
    found = False
    for go_name in go_names:
        for term in terms_of_interest:
            if term in go_name:
                print uniprot, 'known NAME', term
                found = True
                terms.add(term)
        if not found:
            print '%s unknown go name %s' % (uniprot, go_name)
    data.close()

    new_terms = set()
    for trm in terms:
        if trm in ['transport', 'antiport', 'symport']:
            new_terms.add('transport')
        else:
            new_terms.add(trm)
    return new_terms


def parse_colibri2uniprot():
    with open('/home/labs/fleishman/jonathaw/membrane_prediction_DBs/VH_colibri2uniprot.txt', 'r') as fin:
        cont = fin.read().split('\n')
    results = {}
    for l in cont:
        if len(l) == 0:
            continue
        if l[0] != 'b':
            continue
        s = l.split()
        num_i = [i for i, a in enumerate(s) if a.isdigit()][0]
        for i in s[num_i+1:]:
            results[i.lower().split(';')[0]] = s[3]
    return results



if __name__ == '__main__':
    main_vh()
