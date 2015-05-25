'''
A script to make blast files (make_blast), make MSA db files (multiple fastas) and MSA files.
'''


def ncbiXML_parser(file_name):
    """
    :param file_name: a blast .xml file name with path
    :return: a dictionary of dictionaries. target-name:{hit_id: name, hit_seq: seq (no gaps)}
    """
    import re
    with open(file_name, 'r') as f:
        cont = f.read().split('<Hit>')
    results = {}
    result = {}
    for para in cont:
        if len(para) == 0: continue
        for line in para.split('\n'):
            split = re.split('<?>?', line)
            if len(split) < 2: continue
            if split[1] == 'Hit_id':
                result['hit_id'] = split[2]
            if split[1] == 'Hsp_hseq' and 'hit_id' in result.keys():
                result['hit_seq'] = split[2].translate(None, '-')
                results[result['hit_id']] = result
                result = {}
    return results


def blast2fasta():
    '''
    :return: takes one blast .xml result from rost_msa_prep/blast and makes a multiple fasta file of the same sequences
    in the same folder
    '''
    from TMpredict_WinGrade import parse_rostlab_db
    name = args['name']
    path_bl = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/rost_msa_prep/blast/'
    output_bl = path_bl + name + '_blast.xml'
    seq_dict = ncbiXML_parser(output_bl)
    print name
    query = {k: v for k, v in parse_rostlab_db().items() if k == name.split('_')[0]}.values()[0]
    # print query
    # print seq_dict
    with open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/rost_msa_prep/blast/'
                      +name+'_blast.fa', 'wr+') as o:
        o.writelines('>%s\n' % name)
        o.writelines('%s\n' % query['seq'])
        for k, v in seq_dict.items():

            o.writelines('>%s\n' % k)
            o.writelines('%s\n' % v['hit_seq'])
    with open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/rost_msa_prep/names.txt',
              'a') as o:
        o.write(name+'\n')


def make_blasts_local():
    '''
    not working very well
    :return: run blastp locally
    '''
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


def make_blast_remote():
    '''
    not working very well
    :return: run blastp remotely
    '''
    from Bio.Blast import NCBIWWW
    # from Bio import SeqIO
    # from Bio import Seq
    query = fasta_parser()
    print 'blasting', args['name']
    result = NCBIWWW.qblast('blastp', 'nr', query['seq'], auto_format=5, hitlist_size=1000, expect=0.0001)
    print 'writing to', '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/blasts/'\
                        +args['name']+'_remote_blast.xml'
    with open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/blasts/'+args['name']+
                      '_remote_blast.xml', 'wr+') as o:
        o.write(result.read())
    result.close()


def fasta_parser():
    '''
    :return: {name: name, seq: seq} dictionary of fastas
    '''
    with open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/fastas/'
                      +args['name']+'.fasta', 'r') as f:
        cont = f.read().split('\n')
    return {'name': cont[0][1:], 'seq': cont[1]}


def msa_cleaner():
    '''
    :return:goes over a MSA, and prints a "clean" multiple fasta file. clean means all targets have very little gaps in
    segments where the query has a "block" (long segment with no gaps). used to reduce sie of MSA for TMpredict. insures
    targets are from same family (more or less) of the query.
    '''
    from MSA_for_TMpredict import read_fasta_msa
    name = args['name']
    msa = read_fasta_msa('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/rost_msa_prep/msa_unclean/'
                         +args['name']+'_blast_cluster_msa.fa')

    query = [a for a in msa if a.name.lower() == args['name'].lower()][0]
    query_blocks = find_seq_blocks(query.seq)
    passed_cond = []
    for target in msa:
        if blocks_overlap(query_blocks, target.seq):
            passed_cond.append(target)
    print 'processing %s, passed %i, out of %i' % (name, len(passed_cond), len(msa))
    with open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/rost_msa_prep/clean/'
                      +name+'_blast_cluster_msa_cln.fa', 'wr+') as o:
        o.writelines('>%s\n' % name)
        o.writelines('%s\n' % query.seq.translate(None, '-'))
        for target in passed_cond:
            if target.name == '': continue
            o.writelines('>%s\n' % target.name)
            o.writelines('%s\n' % target.seq.translate(None, '-'))


def blocks_overlap(query_blocks, target_seq):
    '''
    :param query_blocks: list of tuples describing blocks in the query sequence (stretches with no gaps)
    :param target_seq: sequence of blast target to be tested
    :return: True if all blocks in the query correspond to less than 10% gaps in the target sequence
    '''
    for bl in query_blocks:
        M = (bl[1]-bl[0]) / 10 + 1
        if not len([a for a in target_seq[bl[0]:bl[1]] if a == '-']) <= M:
            return False
    return True


def find_seq_blocks(seq):
    '''
    :param seq: a sequence
    :return: list of tuples describing blocks, such that they are longer than 14, and have no gaps
    '''
    import re
    gaps = re.compile('[ACDEFGHIKLMNPQRSTVWY]+')
    return [(a.start(), a.end()-1) for a in gaps.finditer(seq) if a.end()-a.start() > 14]


def add_query2clusters():
    """
    :return: goes over a cluster result file and adds the query fasta if it is missing
    """
    path_cluster = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/rost_msa_prep/cluster/'
    query = fasta_parser()
    with open(path_cluster+args['name']+'_blast_cluster.fa', 'r') as f:
        cont = f.read().split('\n')
    if '>'+args['name'] not in cont:
        with open(path_cluster+args['name']+'_blast_cluster.fa', 'a') as o:
            o.write('>%s\n%s\n' % (query['name'], query['seq']))


if __name__ == '__main__':
    import argparse
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', type=str)
    parser.add_argument('-mode', type=str)
    args = vars(parser.parse_args())
    args['name'] = args['name'].lower().split('/')[-1].split('_')[0]
    if args['mode'] == 'blast_remote':
        make_blast_remote()
    elif args['mode'] == 'blast_local':
        make_blasts_local()
    elif args['mode'] == 'blast2fasta':
        blast2fasta()
    elif args['mode'] == 'clean_msa':
        msa_cleaner()
    elif args['mode'] == 'add_query':
        add_query2clusters()