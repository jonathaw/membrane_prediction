def read_entries(name=False, num1=False, num2=False):
    topdb_db = open('data_sets/top_ap_TOPDB.txt', 'r')
    cont = topdb_db.read().split('>')
    topdb_dict = [{'name': i.split('\n')[0], 'seq': i.split('\n')[1], 'topo': i.split('\n')[2],
                   'ss2': 'data_sets/psipred/'+i.split('\n')[0]+'.ss2'} for i in cont if len(i) > 0]
    if name:
        return [i for i in topdb_dict if i['name'] == name]
    elif num2 and num1:
        return topdb_dict[num1:num2]
    else:
        return topdb_dict[:num2]


def topo_define(topdb_dict):
    for i in topdb_dict['topo']:
        if i == 'I':
            return 'in'
        elif i == 'O':
            return 'out'


def topo_compare(query, topdb):
    results = {'tm_overlap': 0, 'tm_topdb': 0, 'tm_query': 0, 'query_num_tm': 0, 'topdb_num_tm': 0,
               'query_n_term': query[0], 'topdb_n_term': topdb[0]}
    print query
    print topdb
    last_q = query[0]
    last_t = topdb[0]
    for q, t in zip(query[1:], topdb[1:]):
        if q == 'M' and t == 'M':
            results['tm_overlap'] += 1
        elif q == 'M' and t != 'M':
            results['tm_query'] += 1
        elif q != 'M' and t == 'M':
            results['tm_topdb'] += 1
        if q == 'M' and last_q != 'M':
            results['query_num_tm'] += 1
        if t == 'M' and last_t != 'M':
            results['topdb_num_tm'] += 1
        last_q = q
        last_t = t
    print results


if __name__ == '__main__':
    read_entries()