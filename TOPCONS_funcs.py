#!/usr/bin/env python2.7
"""
TOPCONS  related functions
"""


def parse_topcons_output(args):
    with open(args['in_path']+args['topcons_file'], 'r') as fin:
        cont = fin.read().split('\n')
    result = {}
    for i, l in enumerate(cont):
        s = l.split()
        try:
            if s[1] == 'name:':
                l_name = l
                name = ''.join(s[2:])
        except:
            pass
        try:
            if s[0] in ['TOPCONS', 'OCTOPUS', 'Philius', 'PolyPhobius', 'SCAMPI', 'SPOCTOPUS']:
                result[s[0].lower()] = cont[i+1]
        except:
            pass

    # print name
    # print result
    # if args['special']:
    #     name_list = parse_assaf()
    #     for n in name_list:
    #         if n in l_name.lower():
    #             name = n
    #             break
    if args['special']:
        name = l_name.split('GN=')[1].split()[0]
    with open(args['out_path']+name+'.topc', 'wr+') as fout:
        fout.write('name %s\n' % name)
        for k, v in result.items():
            fout.write('%s %s\n' % (k, v))

def parse_assaf():
    with open('/home/labs/fleishman/elazara/benchmark_paper_new/eSOL/Soluble/list.list', 'r') as cin:
        return cin.read().lower().split('\n')


if __name__ == '__main__':
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', type=str, help='job name')
    parser.add_argument('-topcons_file', type=str, help='name of TOPCONS result file')
    parser.add_argument('-in_path', type=str, default=os.getcwd()+'/', help='in path')
    parser.add_argument('-out_path', type=str, default=os.getcwd()+'/', help='out path')
    parser.add_argument('-special', type=bool)
    args = vars(parser.parse_args())
    parse_topcons_output(args)