#!/usr/bin/env python2.7
import os
import matplotlib.pyplot as plt
import numpy as np


def parse_psipred_timeit(file_name):
    with open(file_name, 'r') as fin:
        cont = fin.read().split('\n')
    for l in cont:
        s = l.split()
        if len(s) == 0:
            continue
        if s[0] == 'real':
            min = float(s[1].split('m')[0])
            sec = float(s[1].split('m')[1].split('s')[0])
            return min*60.+sec


def parse_fasta(file_name):
    with open(file_name, 'r') as fin:
        cont = fin.read().split('\n')
    return ''.join(cont[1:])


def main():
    result = {}
    with open('all_names.txt', 'r') as fin:
        for l in fin.read().split('\n'):
            if l == '':
                continue
            result[l] = {}
            result[l]['seq'] = parse_fasta('%s/%s.fasta' % (l, l))
            result[l]['psipred_time'] = parse_psipred_timeit('%s/%s.err' % (l, l))
    plt.figure()
    # plt.scatter([len(a['seq']) for a in result.values()], [a['psipred_time'] for a in result.values()])
    plt.boxplot([a['psipred_time'] for a in result.values()])
    plt.savefig('plot.png')
    plt.show()
    print 'the median is', np.median([a['psipred_time'] for a in result.values()])
    print 'the mean is', np.mean([a['psipred_time'] for a in result.values()])


if __name__ == '__main__':
    main()