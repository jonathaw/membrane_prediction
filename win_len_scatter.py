#!/usr/bin/env python
def main():
    import os
    import re
    from topo_strings_comparer import prd_parser, parse_rostlab_db
    import matplotlib.pyplot as plt
    prd_files = [x for x in os.listdir(os.getcwd()) if re.match('.*prd', x)]
    rost_db = parse_rostlab_db()

    result = []

    for prd_file in prd_files:
        prd_name = prd_file.split('.')[0]
        prd = prd_parser(os.getcwd(), prd_file)
        hhh = re.compile('[h]*')
        HHH = re.compile('[H]*')
        obse_list = [(a.start(), a.end()) for a in hhh.finditer(rost_db[prd_name]['pdbtm']) if a.end()-a.start() > 1
                     and a.end()-a.start() >= 10]
        obse_list.extend([(a.start(), a.end()) for a in HHH.finditer(rost_db[prd_name]['pdbtm']) if a.end()-a.start() > 1
                          and a.end()-a.start() >= 10])
        pred_list = [(a.start(), a.end()) for a in HHH.finditer(prd['pred_ts']) if a.end()-a.start() > 1
                     and a.end()-a.start() >= 10]
        for p_win in pred_list:
            res = win_overlappM(p_win, obse_list)
            if res:
                result.append((p_win[1]-p_win[0]+1, res[1]-res[0]+1))
    plt.scatter([a[0] for a in result], [a[1] for a in result])
    plt.xlim([0, 40])
    plt.ylim([0, 40])
    # plt.hist([[a[0] for a in result], [a[1] for a in result]])
    plt.show()


def win_overlappM(pred_h, obse_set, M=10):
    for obse_h in obse_set:
        result = len([a for a in range(pred_h[0], pred_h[1]) if a in range(obse_h[0], obse_h[1])])
        if result >= M:
            return obse_h
    return False


if __name__ == '__main__':
    main()