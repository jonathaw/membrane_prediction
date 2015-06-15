"""
a script to analyse the psipred results compared to whether the AA is membranal or not in the Rost database
"""


def main():
    # check_all_aa_points_as_boxplot()
    from TMpredict_WinGrade import parse_rostlab_db
    import matplotlib.pyplot as plt
    import numpy as np
    IS_BETA_CUTOFF = 0.3
    BETA_NUM_CUTOFF = 5
    missed_h = 0
    rostlab_dict = parse_rostlab_db()
    for k, v in rostlab_dict.items():
        psipred = parse_psipred(k)
        # print v['seq']
        # print v['pdbtm']
        tms = ts2hp_seq(v['seq'], v['pdbtm'])
        for tm in tms:
            are_beta = 0
            for i in range(tm[0], tm[1]+1):
                if psipred[i]['e'] >= IS_BETA_CUTOFF:
                    are_beta += 1
            if are_beta >= BETA_NUM_CUTOFF:
                missed_h += 1
                print "MISSED ME!!!!", k, tm, [psipred[a]['e'] for a in range(tm[0], tm[1]+1)], are_beta
    print "total misses (helices)", missed_h


def check_beta_average():
    '''
    main function here. tests every TM in the Rost data base for its average sheet propensity, also every non-TM.
    outpuits the number of windows that will be discarded for each.
    :return:
    '''
    from TMpredict_WinGrade import parse_rostlab_db
    import matplotlib.pyplot as plt
    import numpy as np
    global IS_BETA_CUTOFF, IS_COIL_CUTOFF, IS_HELIX_CUTOFF
    IS_BETA_CUTOFF = 0.3
    IS_COIL_CUTOFF = 0.48
    IS_HELIX_CUTOFF = 0.3

    tm_missed_h = 0
    non_ym_missed_h = 0
    how_many_non_tms = 0
    tot_passed_of_tmh = 0
    tot_NOT_passed_of_tmh = 0
    tot_passed_of_NON_tmh = 0
    tot_NOT_passed_of_NON_tmh = 0
    rostlab_dict = parse_rostlab_db()
    for k, v in rostlab_dict.items():
        psipred = parse_psipred(k)
        tms = ts2hp_seq(v['seq'], v['pdbtm'])
        for tm in tms:
            # avg = sum([psipred[a]['e'] for a in range(tm[0], tm[1]+1)]) / len(range(tm[0], tm[1]+1))
            # avg = psipred_avg(range(tm[0], tm[1]+1), psipred, 'e')
            # avg_c = psipred_avg(range(tm[0], tm[1]+1), psipred, 'c')
            # avg_h = psipred_avg(range(tm[0], tm[1]+1), psipred, 'h')
            # med = psipred_median(range(tm[0], tm[1]+1), psipred)
            # avgs.append(avg)
        # if avg >= IS_BETA_CUTOFF or avg_c >= IS_COIL_CUTOFF or avg_h <= IS_HELIX_CUTOFF:
        #     print avg,avg_c,avg_h
            if not pass_thresholds(psipred, tm[0], tm[1]):
                tm_missed_h += 1

        ## check how many non_tms will be canceled thanks to threshold
        non_tms = ts2non_tms(v['seq'], v['pdbtm'])
        for non_tm in non_tms:
            rng = range(non_tm[0], non_tm[1]+1)
            for i in rng:
                if i+20 in rng:
                    # avg = psipred_avg(range(i, i+20), psipred, 'e')
                    # avg_c = psipred_avg(range(i, i+20), psipred, 'c')
                    # avg_h = psipred_avg(range(i, i+20), psipred, 'h')
                    # med = psipred_median(range(i, i+20), psipred)
                    how_many_non_tms += 1
                    # if (avg >= IS_BETA_CUTOFF or avg_c >= IS_COIL_CUTOFF or avg_h >= IS_HELIX_CUTOFF):
                    if not pass_thresholds(psipred, i, i+20):
                        non_ym_missed_h += 1

        for i in range(len(v['seq'])-20):
            if do_range_overlap_ranges(range(i+1, i+21), tms):
                if pass_thresholds(psipred, i+1, i+21):
                    tot_passed_of_tmh += 1
                else:
                    tot_NOT_passed_of_tmh += 1
            else:
                if pass_thresholds(psipred, i+1, i+21):
                    tot_passed_of_NON_tmh += 1
                else:
                    tot_NOT_passed_of_NON_tmh += 1
        # break

    print 'TM totalt misses (helices)', tm_missed_h
    print 'NON TM total misses (helices)', non_ym_missed_h
    print "overall %i non tms examined" % how_many_non_tms
    print "\nTotal helices passed and are TMHs %i, Total helices not pass and are TMHs %i" % (tot_passed_of_tmh, tot_NOT_passed_of_tmh)
    print "Total helices passed and are not TMHs %i, Total helices not pass and not TMHs %i" % (tot_passed_of_NON_tmh, tot_NOT_passed_of_NON_tmh)
    # plt.hist(avgs)
    # plt.show()


def do_ranges_overlap(rng1, rng2):
    for i in rng1:
        if i in rng2:
            return True
    else:
        return False


def do_range_overlap_ranges(rng, rngs):
    for rng2 in rngs:
        if do_ranges_overlap(rng, rng2):
            return True
    return False


def pass_thresholds(psi, start, end):
    beta = psipred_avg(range(start+1, end), psi, 'e')
    coil = psipred_avg(range(start+1, end), psi, 'c')
    hlix = psipred_avg(range(start+1, end), psi, 'h')
    return hlix >= IS_HELIX_CUTOFF and beta <= IS_BETA_CUTOFF and coil <= IS_COIL_CUTOFF



def psipred_avg(rng, psipred, t):
    from numpy import mean
    return mean([psipred[a][t] for a in rng])


def psipred_median(rng, psipred, t):
    from numpy import median
    return median([psipred[a][t] for a in rng])


def ts2non_tms(seq, ts):
    '''
    :param seq: AA seq
    :param ts: topo string of prediction/observation
    :return: list of AA sequences of the TMs in ts
    '''
    import re
    hp_seqs = []
    hhh = re.compile('[120lUu]*')
    for it in hhh.finditer(ts):
        if it.end() - it.start() > 1:
            hp_seqs.append([it.start()+1, it.end()-1])
    return hp_seqs


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
            hp_seqs.append([it.start(), it.end()])
    return hp_seqs


def check_all_aa_points_as_boxplot():
    from TMpredict_WinGrade import parse_rostlab_db
    import matplotlib.pyplot as plt
    import numpy as np
    rostlab_dict = parse_rostlab_db()
    membranal = []
    non_membranal = []
    for k, v in rostlab_dict.items():
        psipred = parse_psipred(k)
        for i, aa in enumerate(v['pdbtm']):
            if aa.lower() == 'h':
                membranal.append(psipred[i+1])
            else:
                non_membranal.append(psipred[i+1])
    positions = np.arange(6)
    plt.subplot(111)
    labels = ['MM Coil', 'MM sheet', 'MM Helix', 'no Coil', 'no sheet', 'no Helix']
    mm_coil = [a['c'] for a in membranal]
    mm_sheet = [a['e'] for a in membranal]
    mm_helix = [a['h'] for a in membranal]
    no_coil = [a['c'] for a in non_membranal]
    no_sheet = [a['e'] for a in non_membranal]
    no_helix = [a['h'] for a in non_membranal]
    plt.boxplot([mm_coil, mm_sheet, mm_helix, no_coil, no_sheet, no_helix], labels=labels, positions=positions)
    plt.ylim([-0.5, 1.5])
    plt.show()

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


def parse_psipred(name):
    ss2_file = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/psipred/'+name+'.ss2'
    results = {}
    with open(ss2_file, 'r') as f:
        cont = f.read().split('\n')
    for line in cont:
        split = line.split()
        if len(split) != 6: continue
        results[int(split[0])] = {'aa': split[1], 'c': float(split[3]), 'h': float(split[4]), 'e': float(split[5])}
    return results


if __name__ == '__main__':
    # main()
    check_beta_average()