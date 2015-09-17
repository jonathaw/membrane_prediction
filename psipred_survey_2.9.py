#!/usr/bin/env python
"""
a script to analyse the psipred results compared to whether the AA is membranal or not in the Rost database
"""


def main():
    # check_all_aa_points_as_boxplot()
    from TMpredict_WinGrade import parse_rostlab_db
    import matplotlib.pyplot as plt
    import numpy as np
    rostlab_dict = parse_rostlab_db()
    passed = []
    didnt_pass = []
    could_pass = []
    for k, v in rostlab_dict.items():
        # print k, v
        psipred = parse_psipred(k, v['seq'])
        tms = ts2hp_seq(v['seq'], v['pdbtm'])
        # print v['seq']
        # print v['pdbtm']
        # print ''.join([str(a) for a in range(10)]*50)
        # print tms
        for tm in tms:
            # print 'testing', tm
            if is_not_helical(v['seq'], tm, psipred):
                # print tm, 'didnt pass'
                didnt_pass.append(tm)
                if try_to_pass(tm, v['seq'], psipred):
                    could_pass.append(tm)
                else:
                    print v['name'], tm, np.mean([psipred[a]['e'] for a in range(tm[0], tm[1]+1)]), np.mean([psipred[a]['c'] for a in range(tm[0], tm[1]+1)]), np.mean([psipred[a]['h'] for a in range(tm[0], tm[1]+1)])
                    for i in range(tm[0], tm[1]+1):
                        print 'c: %f h: %f e: %f' % (psipred[i]['c'], psipred[i]['h'], psipred[i]['e'])
            else:
                # print tm, 'has passed'
                passed.append(tm)

        # break
    print 'failed %i times' % len(didnt_pass)
    print 'could have passed %i' % len(could_pass)
    print 'succeeded %i times' % len(passed)


def try_to_pass(tm, seq, psi):
    if not is_not_helical(seq, [tm[0]+7, tm[1]], psi):
        return True
    if not is_not_helical(seq, [tm[0], tm[1]-7], psi):
        return True


# def is_not_helical(seq, pos, psi, verbose=False):
#     """
#     >>> is_not_helical('AAA', [0, 10], 'a')
#     """
#     win_size = 5
#     for i in range(pos[0], pos[1]-win_size+2):
#         print range(i, i+win_size)
#         if all(psi[i]['h'] <= 0.4 for i in range(i, i+win_size)):
#             return True
#     else:
#         return False

def is_not_helical(seq, pos, psi, verbose=False):
    win_size = 6
    for i in range(pos[0], pos[1]-win_size+2):
        if all(psi[j]['e'] >= 0.5 for j in range(i, i+win_size))or \
                all(psi[j]['c'] >= 0.5 for j in range(i, i+win_size)) or \
                all(psi[j]['h'] <= 0.1 for j in range(i, i+win_size)):
            return True
    cs = []
    es = []
    hs = []
    for i in range(pos[0], pos[0]+3):
        cs.append(psi[i]['c'] >= 0.5)
        es.append(psi[i]['e'] >= 0.5)
        hs.append(psi[i]['h'] <= 0.1)
    if all(cs) or all(es) or all(hs):
        return True
    cs = []
    es = []
    hs = []
    for i in range(pos[1]-3, pos[1]):
        cs.append(psi[i]['c'] >= 0.5)
        es.append(psi[i]['e'] >= 0.5)
        hs.append(psi[i]['h'] <= 0.1)
    if all(cs) or all(es) or all(hs):
        return True

# def is_not_helical(seq, pos, psi, verbose=False):
#     win_size = 6
#     for i in range(pos[0], pos[1]-win_size+2):
#         if all(psi[j]['e'] >= 0.5 for i in range(i, i+win_size))or all(psi[i]['c'] >= 0.5 for i in range(i, i+win_size)) :
#             return True
#     cs = []
#     es = []
#     for i in range(pos[0], pos[0]+3):
#         cs.append(psi[i]['c'] > 0.5)
#         es.append(psi[i]['e'] > 0.5)
#     if all(cs) or all(es):
#         return True
#     cs = []
#     es = []
#     for i in range(pos[1]-3, pos[1]):
#         cs.append(psi[i]['c'] > 0.5)
#         es.append(psi[i]['e'] > 0.5)
#     if all(cs) or all(es):
#         return True
#     return all(psi[i]['c'] > 0.5 for i in range(pos[0], pos[0]+3)) or \
#             all(psi[i]['c'] > 0.5 for i in range(pos[1]-3, pos[1])) or \
    #         all(psi[i]['e'] > 0.5 for i in range(pos[0], pos[0]+3)) or \
    #         all(psi[i]['e'] > 0.5 for i in range(pos[1]-3, pos[1]))

# def is_not_helical(seq, pos, psi, verbose=False):
#     """
#     :param pos:start and end positions of the corresponding window
#     :param psi: a dict of c/e/h propensities
#     :return:True if the average e or c propensities are above certain thresholds determined in
#     psipred_vs_mm_nomm.py
#     """
#     import numpy as np
#     from WinGrade import sequential_coiled
#     assert all([psi[i]['aa'] == seq[i] for i in range(pos[0], pos[1])])
#     if sequential_coiled(pos, psi, verbose):
#         if verbose:
#             print 'failing sequential'
#         return True
#     if verbose:
#         print 'didnt fail sequental'
#         print 'found avgs: e: %f, c: %f, h: %f' % \
#               (np.mean([psi[a]['e'] for a in range(pos[0], pos[1])]), np.mean([psi[a]['c'] for a in range(pos[0], pos[1])]), np.mean([psi[a]['h'] for a in range(pos[0], pos[1])]))
#     for a in range(pos[0], pos[1]):
#         print 'in %i found e: %f c: %f h: %f' % (a, psi[a]['e'], psi[a]['c'], psi[a]['h'])
    # return True if (np.mean([psi[a]['e'] for a in range(pos[0], pos[1]+1)]) >= 0.3 or
    #                 np.mean([psi[a]['c'] for a in range(pos[0], pos[1]+1)]) >= 0.48) else False#  or
    #                 np.mean([psi[a]['h'] for a in range(pos[0], pos[1]+1)]) <= 0.3) else False


def ts2hp_seq(seq, ts):
    '''
    :param seq: AA seq
    :param ts: topo string of prediction/observation
    :return: list of AA sequences of the TMs in ts
    '''
    import re
    hp_seqs = []
    hhh = re.compile('[hH]*')
    for it in hhh.finditer(ts):
        if it.end() - it.start() > 1:
            hp_seqs.append([it.start(), it.end()-1])
    return hp_seqs


def parse_psipred(name, seq):
    results = {}
    with open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/psipred/'+name+'.ss2', 'r') as f:
        cont = f.read().split('\n')
    for line in cont:
        split = line.split()
        if len(split) != 6:
            continue
        assert split[1] == seq[int(split[0])-1], 'residues dont match %i, %s %s' % (int(split[0])-1, split[1], seq[int(split[0])-1])
        results[int(split[0])-1] = {'aa': split[1], 'c': float(split[3]), 'h': float(split[4]), 'e': float(split[5])}
    return results


if __name__ == '__main__':
    main()