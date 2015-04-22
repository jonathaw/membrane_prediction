"""
a script to summerize topo prediction results over the Von-Heijne and phobius data.
"""


def main():
    import os
    import re
    import matplotlib.pyplot as plt
    from sasa_survey import boxplot_with_scatter
    prd_list = [x for x in os.listdir(os.getcwd()) if re.match('.*\.prd', x)]
    results = {'entries': 0, 'pred_correct': 0, 'phobius_correct': 0, 'pred_correct_gap': 0, 'pred_incorrect_gap': 0}
    pred_correct, pred_incorrect = [], []
    for prd in prd_list:
        prd_pars = topo_prd_parse(prd)
        results['entries'] += 1
        results['pred_correct'] += 1 if prd_pars['pred_correct'] == 'True' else 0
        results['phobius_correct'] += 1 if prd_pars['phobius_correct'] == 'True' else 0
        if prd_pars['pred_correct'] == 'True':
            pred_correct.append(prd_pars['best_sec_best_delta'])
            if prd_pars['best_sec_best_delta'] <= -3:
                results['pred_correct_gap'] += 1
        else:
            pred_incorrect.append(prd_pars['best_sec_best_delta'])
            if prd_pars['best_sec_best_delta'] <= -3:
                results['pred_incorrect_gap'] += 1
    print 'num correct', len(pred_correct)
    print 'num incorrect', len(pred_incorrect)
    # plt.subplot(121)
    # plt.boxplot([pred_correct, pred_incorrect])
    # plt.ylim((-20,5))
# def boxplot_with_scatter(data, labels, fig_num, title='default', axis_labels=('x axis', 'y axis')):
    print [results['pred_correct'], results['entries']-results['pred_correct']]
    boxplot_with_scatter([pred_correct, pred_incorrect], ['correct', 'incorrect'], 111)

    print 'result gap correct', results['pred_correct_gap']
    print 'result gap incorrect', results['pred_incorrect_gap']
    # plt.subplot(122)
    # positions = np.arange(2)
    # plt.bar(positions, results)
    plt.show()


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