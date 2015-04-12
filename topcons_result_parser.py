"""
A script to read topcons output files, get the SCAMPI prediction, and write .ord files for it, similar
to TMpredict.
"""


def main():
    import os
    import re
    from TMpredict_WinGrade import result_comparer, results_writer, parse_rostlab_db, result_comparer_10overlap
    topcons_path = '/home/labs/fleishman/jonathaw/membrane_topcons/all_results/'
    rostlab_db_dict = parse_rostlab_db()
    file_list = [x for x in os.listdir(topcons_path)
                 if re.match('.*\.txt', x)]
    for file_i in file_list:
        entry = topcons_parser(topcons_path, file_i)
        topo_string = entry['rost_format_scampi']
        pred_tm = len(re.findall('[1u20LU]h', entry['rost_format_scampi']))

        rost_results = rostlab_db_dict[entry['name'].lower()]

        entry_results = {}
        entry_results['pdbtm'] = result_comparer(rost_results['pdbtm'], topo_string, pred_tm)
        entry_results['opm'] = result_comparer(rost_results['opm'], topo_string, pred_tm)
        overlap10 = {'pdbtm': result_comparer_10overlap(rost_results['pdbtm'], topo_string),
                     'opm': result_comparer_10overlap(rost_results['pdbtm'], topo_string)}
        results_writer(entry_results, rost_results, topo_string, None, None, pred_tm, topcons_path, overlap10)


def topcons2rostlab_ts_format(ts):
    """
    :param ts: a topcons format topo string
    :return: a rostlab format topo string
    """
    ts_list = list(ts)
    topcons2rostlab_dict = {'i': '1', 'o': '2', 'M': 'h', 'm': 'h'}
    for i, val in enumerate(ts_list):
        ts_list[i] = topcons2rostlab_dict[val]
    return ''.join(ts_list)


def topcons_parser(path, file_name):
    """
    :param path: path to topcons output files
    :param file_name: file name
    :return: dictionary with name and SCAMPI topo string
    """
    result = {}
    with open(path+file_name, 'r') as f:
        cont = f.read().split('\n')
    for i, line in enumerate(cont):
        split = line.split()
        if len(split) <= 1:
            continue
        if split[1] == 'name:':
            result['name'] = split[2]
        if split[0] == 'SCAMPI':
            result['scampi'] = cont[i+1]
    result['rost_format_scampi'] = topcons2rostlab_ts_format(result['scampi'])
    return result

if __name__ == '__main__':
    main()