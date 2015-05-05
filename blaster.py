'''
A script to make blast files (make_blast), make MSA db files (multiple fastas) and MSA files.
'''


def ncbiXML_parser(file_name):
    import re
    with open(file_name, 'r') as f:
        cont = f.read().split()
    results = {}
    result = {}
    for line in cont:
        if len(line) == 0: continue
        split = re.split('<?>?', line)
        if len(split) < 2: continue
        # if split[1] == 'Hsp_qseq':
        #     query = split[2].translate(None, '-')
        if split[1] == 'Hit_id':
            result['hit_id'] = split[2]
        if split[1] == 'Hsp_hseq':
            result['hit_seq'] = split[2].translate(None, '-')
            results[result['hit_id']] = result
            result = {}
    return results#, query


def make_MSA():
    import sys
    from TMpredict_WinGrade import parse_rostlab_db
    name = sys.argv[1].split('/')[-1].split('.')[0]
    path_bl = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/blasts/'
    output_bl = path_bl + name + '_blast.xml'
    # seq_dict, query = ncbiXML_parser(output_bl)
    seq_dict = ncbiXML_parser(output_bl)
    query = {k: v for k, v in parse_rostlab_db().items() if k == name}.values()[0]
    print query
    with open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/MSA_fastas/'
                      +name+'_MSA_fastas.fa', 'wr+') as o:
        o.writelines('>%s\n' % name)
        o.writelines('%s\n' % query['seq'])
        for k, v in seq_dict.items():
            o.writelines('>%s\n' % k)
            o.writelines('%s\n' % v['hit_seq'])


def make_blasts():
    import os
    import sys
    from Bio.Blast.Applications import NcbiblastpCommandline
    name = sys.argv[1].split('/')[-1].split('.')[0]
    print 'blasting ', name
    path_fa = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/fastas/'
    path_bl = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/blasts/'
    input_fa = path_fa + name + '.fasta'
    output_bl = path_bl + name + '_blast.xml'
    comline = NcbiblastpCommandline(query=input_fa, db="/shareDB/nr/nr", evalue=0.001, outfmt=5, out=output_bl)
    print comline
    os.system(str(comline))


if __name__ == '__main__':
    # make_blasts()
    make_MSA()
    pass