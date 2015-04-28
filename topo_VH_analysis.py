"""
a script to summerize topo prediction results over the Von-Heijne and phobius data.
"""


def main():
    import os
    import re
    import matplotlib.pyplot as plt
    import numpy as np
    from sasa_survey import boxplot_with_scatter
    energy_gap = -1.5
    prd_list = [x for x in os.listdir(os.getcwd()) if re.match('.*\.prd', x)]
    results = {'entries': 0, 'pred_correct': 0, 'phobius_correct': 0, 'pred_correct_gap': 0, 'pred_incorrect_gap': 0,
               'scampi_correct': 0}
    pred_correct, pred_incorrect = [], []
    for prd in prd_list:
        prd_pars = topo_prd_parse(prd)
        # print prd_pars
        if prd_pars['name'] == 'AraJ': continue
        results['entries'] += 1
        results['pred_correct'] += 1 if prd_pars['pred_correct'] == 'True' else 0

        topcons = topcons_topo_parser('/home/labs/fleishman/jonathaw/membrane_topcons/topo_VH_topcons/all_results/'+
                         prd_pars['name'].lower()+'.prd')

        if prd_pars['pred_correct'] == 'True':
            pred_correct.append(prd_pars['best_sec_best_delta'])
            if prd_pars['best_sec_best_delta'] <= energy_gap:
                results['pred_correct_gap'] += 1
                results['phobius_correct'] += 1 if prd_pars['phobius_correct'] == 'True' else 0
                if (prd_pars['obs_c_term'] == 'out' and topcons['top'][-1] == '2') or (prd_pars['obs_c_term'] == 'in'
                                                                                       and topcons['top'][-1] == '1'):
                    results['scampi_correct'] += 1

        else:
            pred_incorrect.append(prd_pars['best_sec_best_delta'])
            if prd_pars['best_sec_best_delta'] <= energy_gap:
                results['pred_incorrect_gap'] += 1




    print 'energy gap', energy_gap
    print 'num correct', len(pred_correct)
    print 'num incorrect', len(pred_incorrect)
    plt.figure()
    # plt.subplot(121)
    print [results['pred_correct'], results['entries']-results['pred_correct']]
    # boxplot_with_scatter([pred_correct, pred_incorrect], ['correct', 'incorrect'], 121)

    print 'result gap correct', results['pred_correct_gap']
    print 'result gap incorrect', results['pred_incorrect_gap']
    print 'result phobius correct', results['phobius_correct']
    print 'result scampi correct', results['scampi_correct']
    plt.subplot(111)
    the_range = [(a, a+1) for a in np.arange(-6, 0, 1)]
    correct_bins = [float(len(b)) for b in break2bins(pred_correct, the_range)]
    incorrect_bins = [float(len(b)) for b in break2bins(pred_incorrect, the_range)]
    correct_prec = [100*float(correct_bins[i])/float(incorrect_bins[i]+correct_bins[i]) for i in range(len(correct_bins))]
    Ns = [correct_bins[i]+incorrect_bins[i] for i in range(len(correct_bins))]
    bars = plt.bar([a for a in np.arange(len(correct_prec))], correct_prec)
    # plt.xlabel(ranges2labels(the_range))
    plt.xlabel('dG ranges')
    plt.xticks([a+0.5 for a in np.arange(len(the_range))], ranges2labels(the_range), rotation='vertical')
    plt.ylabel('C terminus correct prediction percent')
    plt.ylim((0, 101))
    for i, rect in enumerate(bars):
        height = rect.get_height()
        print height, Ns[i], correct_bins[i], incorrect_bins[i]
        plt.text(rect.get_x() + rect.get_width() / 2., 1.005 * height, '%d' % int(Ns[i]),
                 ha='center', va='bottom')
    # plt.subplot(122)
    # positions = np.arange(2)
    # plt.bar(positions, results)
    plt.show()


def ranges2labels(ranges):
    return [str(r[0])+' , '+str(r[1]) for r in ranges]


def break2bins(data, ranges):
    result = []
    for r in ranges:
        result.append([a for a in data if r[0] < a < r[1]])
    return result


def topcons_topo_parser(name):
    result = {}
    with open(name, 'r') as f:
        cont = f.read().split('\n')
    for line in cont:
        split = line.split()
        if len(split) <= 1: continue
        result[split[0]] = split[1]
    return result

def topo_prd_parse(name):
    with open(name, 'r') as f:
        cont = f.read().split('\n')
    result = {}
    for line in cont:
        split = line.split()
        if split == []: continue
        try:
            result[split[0]] = float(split[1])
        except:
            result[split[0]] = split[1]
    return result


if __name__ == '__main__':
    main()