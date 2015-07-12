"""
A script to analyse TopoGraph prediciton runs. useful both for ROC and single folders
"""


def compare_ROC(path):
    '''
    :param path: a path to a directory to analyse
    :return: runs over all ROC_ folders in path and analyses their predcitions if they have > 180 results (.prd files).
    prints the name of the folder with the best results
    '''
    import re
    roc_dirs = [x for x in os.listdir(path) if re.match('ROC_.*', x)]
    best_grade = 0
    for dir in roc_dirs:
        results, totals = prd_directory(dir)
        # print totals
        # print results
        # print sum(totals.values())
        # print sum(results.values())
        if sum(totals.values()) >= 170 and sum(results.values()) > best_grade:
            print 'inside', totals, results, dir
            best_grade = sum(results.values())
            best_name = dir
    print best_grade
    print best_name


def prd_directory(dir_path):
    """
    :param dir_path: path to directory to analyse
    :return: if in ROC mode returns prediction results. if in single mode, shows a graph of the results
    """
    import re, os
    from TMpredict_WinGrade import parse_rostlab_db
    from topcons_result_parser import topcons2rostlab_ts_format
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    M = 10
    file_list = [x for x in os.listdir(dir_path) if re.match('.*\.prd', x)]
    if len(file_list) < args['num_prd']: return {'tm_1': 0, 'tm_2_5': 0, 'tm_5': 0}, {'tm_1': 0, 'tm_2_5': 0, 'tm_5': 0}
    rostlab_db_dict = parse_rostlab_db()
    predictors = ['polyphobius', 'topcons', 'spoctopus', 'philius', 'octopus', 'scampi', 'pred_ts']
    results = {a: {'tm_1': 0, 'tm_2_5': 0, 'tm_5': 0} for a in predictors}
    totals = {'tm_1': 0, 'tm_2_5': 0, 'tm_5': 0}

    errors = {'over': 0, 'miss': 0, 'exact': 0, 'total': 0}

    for file_name in file_list:
        pred = prd_parser(dir_path, file_name)
        # if pred['name'] != 'p0c7b7': continue
        obse = rostlab_db_dict[pred['name']]
        topc = spc_parser(pred['name'])
        predictors = {k: topcons2rostlab_ts_format(v) for k, v in topc.items() if k not in ['name', 'seq']}
        predictors['pred_ts'] = pred['pred_ts']
        first_passage = True
        for predictor in predictors:
            comp_pdbtm = comparer(obse['pdbtm'], predictors[predictor], M)
            comp_opm = comparer(obse['opm'], predictors[predictor], M)

            if predictor == 'pred_ts' and not comp_pdbtm['overlapM_ok']:
                # print 'AAAAAAAHHHHHHH !!!!!! :('
                # print 'obse_tm_num', comp_pdbtm['obse_tm_num']
                # print 'pred_tm_num', comp_pdbtm['pred_tm_num']
                # print 'ok', comp_pdbtm['overlapM_ok_helices']
                if comp_pdbtm['obse_tm_num'] > comp_pdbtm['pred_tm_num']:
                    # print 'MISS', pred['name'], comp_pdbtm['obse_tm_num']
                    errors['miss'] += 1
                elif comp_pdbtm['obse_tm_num'] < comp_pdbtm['pred_tm_num']:
                    # print 'OVER', pred['name'], comp_pdbtm['obse_tm_num']
                    errors['over'] += 1
                else:
                    errors['exact'] += 1
                errors['total'] += 1

            if comp_pdbtm['obse_tm_num'] == 0 or comp_opm['obse_tm_num'] == 0: continue
            if comp_pdbtm['obse_tm_num'] != comp_opm['obse_tm_num']:
                if comp_pdbtm['overlapM_ok']:
                    results[predictor][tm_num2range(comp_pdbtm['obse_tm_num'])] += 1
                    if first_passage: totals[tm_num2range(comp_pdbtm['obse_tm_num'])] += 1
                elif comp_opm['overlapM_ok']:
                    results[predictor][tm_num2range(comp_opm['obse_tm_num'])] += 1
                    if first_passage: totals[tm_num2range(comp_opm['obse_tm_num'])] += 1
            else:
                results[predictor][tm_num2range(comp_opm['obse_tm_num'])] += 1 \
                    if (comp_pdbtm['overlapM_ok'] or comp_opm['overlapM_ok']) else 0
                if first_passage: totals[tm_num2range(comp_pdbtm['obse_tm_num'])] += 1
            first_passage = False
    if args['mode'] == 'ROC':
        data = {k: v for k, v in results['pred_ts'].items()}
        return data, totals
    else:
        print 'results', results
        print 'totals', totals
        print 'errors', errors

        plt.figure()
        data = {}
        for predictor, results_d in results.items():
            data[predictor] = {k: 100*float(v)/float(totals[k]) for k, v in results_d.items()}
        print 'pps', results['polyphobius']
        font = {'family': 'normal', 'size': 22}
        matplotlib.rc('font', **font)
        print data
        print 'range', np.arange(0, 1./3., 1./(7.*3.)), len(np.arange(0, 1./3., 1./(7.*3.)))
        ind = np.arange(3)
        width = 1./3. * (1./7.)
        incs = np.arange(0, 1./3., 1./(7.*3.))
        colors = ['red', 'blue', 'green', 'black', 'orange', 'pink', 'grey']
        print ind
        plots = {}
        for predictor, details, inc, col in zip(data.keys(), data.values(), incs, colors):
            print predictor, details, inc
            plots[predictor] = plt.bar(ind + inc, details.values(), width, color=col)
        plt.ylim((0, 105))
        plt.xlim((-0.15, 3.7))
        plt.xticks(np.arange(3)+0.15, ['1', '2-5', '5<'])
        plt.xlabel('Number of TMH')
        plt.ylabel('Overlap 10 Accuracy (%)')
        plt.title('TMH prediction comparison')
        names = [k for k in plots.keys()]
        names[0] = 'TopoGraph'
        plt.legend(plots.values(), names, 'upper right')
        plt.show()


def devide_with_zero(num1, num2):
    if num2 == 0: return 0
    else: return float(num1)/float(num2)


def tm_num2range(num):
    if num == 1:
        return 'tm_1'
    elif 2 <= num <= 5:
        return 'tm_2_5'
    elif num > 5:
        return 'tm_5'


def comparer(obse, pred, M):
    """
    :param obse: observed topo string
    :param pred: predicted topo string
    :param M: overlap M.
    :return: dictionary with observe/predicted #TM, and if it passes overlap M.
    """
    import re
    assert len(obse) == len(pred), 'observed and predicted strings lengths do not match'
    # print obse
    # print pred
    allowed = ['1', '2', 'h', 'l']
    obse_cln, pred_cln = '', ''
    for i, c in enumerate(obse):
        obse_cln += c if (pred[i].lower() in allowed and c.lower() in allowed) else 'u'
    for i, c in enumerate(pred):
        pred_cln += c if (obse[i].lower() in allowed and c.lower() in allowed) else 'u'
    obse_cln = obse_cln.replace('l', 'h')
    obse_cln = obse_cln.replace('L', 'h')
    # print obse_cln
    # print pred_cln
    result = {}

    hhh = re.compile('[hH]*')
    obse_list = [(a.start(), a.end()) for a in hhh.finditer(obse_cln) if a.end()-a.start() > 1]
    pred_list = [(a.start(), a.end()) for a in hhh.finditer(pred_cln) if a.end()-a.start() > 1]

    # print obse_list[0][0], obse_list[0][1]
    # print pred_list
    # print obse[obse_list[0][0]:obse_list[0][1]]
    # print pred[pred_list[0][0]:pred_list[0][1]]
    result['obse_tm_num'] = len(obse_list)
    result['pred_tm_num'] = len(pred_list)

    result['overlapM_ok_helices'] = 0
    for pred_helix in pred_list:
        result['overlapM_ok_helices'] += 1 if overlappM(pred_helix, obse_list, M) else 0
    # all helices must be predicted correctly, no extra or "missing" helices
    result['overlapM_ok'] = True if result['overlapM_ok_helices'] == result['obse_tm_num'] == result['pred_tm_num'] \
        else False

    return result


def overlappM(pred_h, obse_set, M=10):
    for obse_h in obse_set:
        result = len([a for a in range(pred_h[0], pred_h[1]) if a in range(obse_h[0], obse_h[1])])
        if result >= M: return True
    return False


def spc_parser(name):
    """
    :param name: protein name
    :return: dictionary with data from topcons regarding name
    """
    with open('/home/labs/fleishman/jonathaw/membrane_prediction_DBs/spoctopus_SPDB/'+name+'.spc', 'r') as f:
        cont = f.read().split('\n')
    result = {}
    for line in cont:
        split = line.split()
        if split == []: continue
        result[split[0]] = split[1]
    return result


def prd_parser(file_path, file_name):
    with open(file_path+'/'+file_name, 'r') as f:
        cont = f.read().split('\n')
    result = {}
    for line in cont:
        split = line.split()
        if split == [] or len(split) < 2 or split[1] == 'None': continue
        result[split[0]] = split[1]
    result['name'] = result['name'].translate(None, "'")
    return result


def compare_just_one():
    """
    mode: one
    path: path to .prd
    name: protein id
    :return: compares a single prediction to it's database input.
    """
    from TMpredict_WinGrade import parse_rostlab_db
    from topcons_result_parser import topcons2rostlab_ts_format
    M = 10
    rostlab_db_dict = parse_rostlab_db()
    pred = prd_parser(args['path'], args['name'])
    obse = rostlab_db_dict[pred['name']]
    topc = spc_parser(pred['name'])
    predictors = {k: topcons2rostlab_ts_format(v) for k, v in topc.items() if k not in ['name', 'seq']}
    predictors_results = {k: None for k in topc.keys() if k not in ['name', 'seq']}
    for predictor in predictors:
        comp_pdbtm = comparer(obse['pdbtm'], predictors[predictor], M)
        comp_opm = comparer(obse['opm'], predictors[predictor], M)
        predictors_results[predictor] = comp_pdbtm['overlapM_ok'] or comp_opm['overlapM_ok']
    comp_pdbtm = comparer(obse['pdbtm'], pred['pred_ts'], M)
    comp_opm = comparer(obse['opm'], pred['pred_ts'], M)
    if comp_opm['overlapM_ok'] and comp_pdbtm['overlapM_ok']: print 'TopoGraph is correct by both'
    elif comp_opm['overlapM_ok']: print 'TopoGraph is correct ONLY by OPM'
    elif comp_pdbtm['overlapM_ok']: print 'TopoGraph is correct ONLY by PDBTM'
    print predictors_results

if __name__ == '__main__':
    import argparse
    import os
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', default='single', type=str)
    parser.add_argument('-path', default=os.getcwd(), type=str)
    parser.add_argument('-name', type=str)
    parser.add_argument('-num_prd', default=180, type=int)
    args = vars(parser.parse_args())

    if args['mode'] == 'single':
        prd_directory(args['path'])
    elif args['mode'] == 'ROC':
        compare_ROC(args['path'])
    elif args['mode'] == 'one':
        compare_just_one()
