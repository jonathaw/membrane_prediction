#!/usr/bin/env python2.7
import urllib2
from os import system
import shutil
from collections import OrderedDict
import threading
from Bio import Entrez


def make_blast_db(args):
    """
    :param args: run arguments
    :return:
    """
    print 'making blast DB'
    system("makeblastdb -in %s_blast.fasta -dbtype 'prot' -out %s_blastdb" % (args['name'], args['name']))


def run_blastp_ncbi(args):
    from Bio.Blast import NCBIWWW
    from Bio import SeqIO, Seq

    seq = Seq.Seq(args['seq'])
    print 'running blastp'
    result = NCBIWWW.qblast('blastp', 'nr', seq, expect=args['blast_evalue'])
    print 'finished running BLASTP'
    with open(args['name'] + args['blast_suffix'], 'wr+') as fout:
        fout.write(result.read())


def run_psipred(args):
    print 'running psipred'
    if args['psipred_db'] == 'normal':
        system('/apps/RH6U4/psipred/3.5/bin/runpsipred-1 %s.fasta' % args['name'])
    elif args['psipred_db'] == 'nr':
        system('/apps/RH6U4/psipred/3.5/bin/runpsipred-nr %s.fasta' % args['name'])


def concurrent_map(func, data):
    """
    picked up from http://code.activestate.com/recipes/577360-a-multithreaded-concurrent-version-of-map/
    Similar to the bultin function map(). But spawn a thread for each argument
    and apply `func` concurrently.

    Note: unlike map(), we cannot take an iterable argument. `data` should be an
    indexable sequence.
    """

    N = len(data)
    result = [None] * N

    # wrapper to dispose the result in the right slot
    def task_wrapper(i):
        result[i] = func(data[i])

    threads = [threading.Thread(target=task_wrapper, args=(i,)) for i in xrange(N)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    return result


def fetch(url):
    return urllib2.urlopen(url).read()


def gather_psipred_seqs(args):
    import os
    print 'gathering sequences from psipred .blast'
    blast_file = [a for a in os.listdir('./') if a[-6:] == '.blast'][0]
    with open(blast_file, 'r') as fin:
        cont = fin.read().split('\n')
    with open('%s_blast.fasta' % args['name'], 'wr+') as fout:
        names = set()
        url_dict = {}
        for l in cont:
            if args['psipred_db'] == 'normal':
                if l[:9] == 'UniRef90_':
                    s = l.split()
                    name = s[0].split('_')[1]
                    try:
                        evalue = float(s[-1])
                    except:
                        evalue = float(str(1)+s[-1])
                    if evalue <= args['blast_evalue'] and name[:2] != 'UP' and name not in names:
                        # score = int(s[-2])
                        url_dict[name] = 'http://www.uniprot.org/uniprot/%s.fasta' % name
                    names.add(name)
            elif args['psipred_db'] == 'nr':
                if l.count('|') >= 2:
                    s = l.split()
                    name = s[0].split('|')[1]
                    if name not in names:
                        try:
                            evalue = float(s[-1])
                        except:
                            evalue = float(str(1)+s[-1])
                        if evalue <= args['blast_evalue'] and name not in names:
                            url_dict[name] = 'http://www.uniprot.org/uniprot/?sort=score&desc=&compress=no&query=%s&fil=&format=fasta' % name
                            # try:
                            #     seq = ''.join(html.read().split('>')[1].split('\n')[1:])
                            # except:
                            #     continue
                    names.add(name)
        gathered = 0
        for res in concurrent_map(fetch, list(url_dict.values())):
            s = res.split('\n')
            if len(s) > 1:
                seq = ''.join(res.split('>')[1].split('\n')[1:])
                name = s[0].split('|')[1]
                if 0.8*len(args['seq']) <= len(seq) <= 1.2*len(args['seq']):
                    fout.write('>%s\n%s\n' % (name, seq))
                    gathered += 1
    print 'gathered %i sequences' % gathered


def retrive_psipred_efetch(args):
    import os
    print 'gathering sequences from psipred .blast'
    blast_file = [a for a in os.listdir('./') if a[-6:] == '.blast'][0]
    with open(blast_file, 'r') as fin:
        cont = fin.read().split('\n')
    with open('%s_blast.fasta' % args['name'], 'wr+') as fout:
        names = set()
        wanted_names = []
        for l in cont:
            if args['psipred_db'] == 'normal':
                if l[:9] == 'UniRef90_':
                    s = l.split()
                    name = s[0].split('_')[1]
                    try:
                        evalue = float(s[-1])
                    except:
                        evalue = float(str(1)+s[-1])
                    if evalue <= args['blast_evalue'] and name[:2] != 'UP' and name not in names:
                        wanted_names.append(name)
                    names.add(name)
            elif args['psipred_db'] == 'nr':
                if l.count('|') >= 2:
                    s = l.split()
                    name = s[0].split('|')[1]
                    if name not in names:
                        try:
                            evalue = float(s[-1])
                        except:
                            evalue = float(str(1)+s[-1])
                        if evalue <= args['blast_evalue'] and name not in names:
                            wanted_names.append(name)
                    names.add(name)
    Entrez.email = 'jonathan.weinstein2012@gmail.com'
    handle = Entrez.efetch(db='protein', id=','.join(wanted_names), rettype='fasta', retmode='text')
    with open('%s_blast.fasta' % args['name'], 'wr+') as fout:
        for entry in handle.read().split('>'):
            if entry != '':
                s = entry.split('\n')
                # name = '|'.join(s[0].split('|')[0:2])+'|'
                name = s[0].split('|')[3]
                seq = ''.join(s[1:])
                print name, seq
                # if 0.8*len(args['seq']) <= len(seq) <= 1.2*len(args['seq']):
                fout.write('>%s\n%s\n' % (name, seq))
    handle.close()


def run_blastp(args):
    from os import system
    print 'running BLASTP'
    print 'blastp -query %s -db %s -evalue %f -max_target_seqs %i -outfmt 5 -out %s' % (args['name'] +
                                                                                        args['seq_file_suffix'],
                                                                                        args['blastp_db'],
                                                                                        args['blast_evalue'],
                                                                                        args['blast_max_target_seqs'],
                                                                                        args['name'] +
                                                                                        args['blast_suffix'])
    system('blastp -query %s -db %s -evalue %f -max_target_seqs %i -outfmt 5 -out %s' %
           (args['path']+args['name'] + args['seq_file_suffix'], args['blastp_db'], args['blast_evalue'],
            args['blast_max_target_seqs'], args['path']+args['name'] + args['blast_suffix']))


def retrive_pfam_matches(args):
    """
    :param args: run arguments. uniprot is essential
    :return: {accession: [start, end]} from the PFAM database
    """
    import re
    import urllib2

    html = urllib2.urlopen('http://pfam.xfam.org/protein/%s?output=xml' % args['uniprot'])

    triangles_out = re.compile('>.*<')
    accession_re = re.compile('accession="[a-zA-Z0-9]*"')
    start_re = re.compile('start="[0-9]*"')
    end_re = re.compile('end="[0-9]*"')

    html_split = html.read().split('\n')
    matches = {}
    for i, l in enumerate(html_split):
        s = l.split()
        if len(s) == 0:
            continue
        if s[0] == '<sequence':
            html_seq = triangles_out.findall(l)[0][1:-1]
            assert args['seq'] == html_seq, 'sequences do not match'
        if s[0] == '<match':
            segment = '\n'.join(html_split[i - 1:i + 3])
            accession = accession_re.findall(segment)[0].split('=')[1].replace('"', '')
            start = start_re.findall(segment)[0].split('=')[1].replace('"', '')
            end = end_re.findall(segment)[0].split('=')[1].replace('"', '')
            matches[accession] = {'segment': [int(start), int(end)]}
    return matches


def retrieve_pfam_MSA_add_query(args, match, tupl):
    """
    :param match: accession at Pfam
    :return:None. downloads an MSA file with no -
    """
    # import prody
    # res = prody.fetchPfamMSA(match, gaps='none', format='fasta', alignment=args['pfam_db'])
    # print 'prody returned', res
    html = urllib2.urlopen(
        'http://pfam.xfam.org/family/%s/alignment/%s/format?format=fasta&alnType=%s&order=t&case=l&gaps=dash'
        % (match, args['pfam_db'], args['pfam_db']))
    with open('%s_%s.fasta' % (match, args['pfam_db']), 'wr+') as fout:
        sp = html.read().upper().split('>')
        head = '>'.join(sp[0:min(500, len(sp))])
        fout.write(head)
    write_seq_to_file({'name': match, 'seq': args['seq'][tupl[0]:tupl[1] + 1], 'seq_file_suffix': '.fasta'})
    system('muscle -profile -in1 %s -in2 %s -out %s' % ('%s_%s.fasta' % (match, args['pfam_db']),
                                                        match + args['seq_file_suffix'], 'temp.fasta'))
    shutil.move('temp.fasta', '%s_%s.fasta' % (match, args['pfam_db']))


def remove_bad_seqs(args, tup):
    fastas = read_multi_fastas(args['name'] + '_' + args['pfam_db'] + args['seq_file_suffix'])
    PF_seq = fastas[args['name']]['seq']
    query_len = tup[1] - tup[0]
    removed, kept = 0, 0
    with open(args['name'] + '_temp.fasta', 'wr+') as fout:
        for k, v in fastas.items():
            if query_len * 0.9 <= len(v['seq'].replace('-', '')) <= query_len * 1.1 and \
                            calc_aln_identity(PF_seq, v['seq']) >= 0.35:
                fout.write('>%s\n%s\n' % (k, v['seq']))
                kept += 1
            else:
                removed += 1
    shutil.move(args['name'] + '_temp.fasta', '%s_%s.fasta' % (args['name'], args['pfam_db']))
    print 'removed %i sequences, kept %i sequences' % (removed, kept)


def calc_aln_identity(seq1, seq2):
    counter = 0
    for i, aa in enumerate(seq1):
        counter += 1 if aa == seq2[i] else 0
    return float(counter) / len(seq1)


def concatenate_msas(args, matches):
    """
    :param args: run arguments
    :param matches: {accession: {segment: [start, end], msa: {name: {name: ,seq: }}
    :return:creates a MSA with the query at the top.
    """
    GAPS = '------------------------------'
    matches_by_start = OrderedDict({v['segment'][0]: {'accession': k, 'msa': v['msa'], 'start': v['segment'][0],
                                                      'end': v['segment'][1]} for k, v in matches.items()})
    length_to_now, last_end = 0, 0
    seqs = {args['name']: {'name': args['name'], 'seq': args['seq'][0]}}
    sorted_keys = sorted(matches_by_start.keys())
    # if sorted_keys[0] > 0:
    # seqs[args['name']] = {'name': args['name'], 'seq': args['seq'][0:sorted_keys[0]]+GAPS}
    #     length_to_now = sorted_keys[0]+len(GAPS)
    for start in sorted_keys:
        print matches_by_start[start]['start'], matches_by_start[start]['end'], matches_by_start[start]['accession']
        length_added = 0
        if start > last_end:
            seqs[args['name']]['seq'] += args['seq'][last_end + 1:start] + GAPS
            length_to_now = len(seqs[args['name']]['seq'])
            length_added = len(args['seq'][last_end + 1:start] + GAPS)
        for val in matches_by_start[start]['msa'].values():
            nm = val['name'].split('/')[0]
            if nm not in seqs.keys():
                if nm == matches_by_start[start]['accession']:
                    seqs[args['name']]['seq'] += val['seq'] + GAPS
                else:
                    seqs[nm] = {'name': nm, 'seq': ''.join(['-'] * length_to_now) + val['seq'] + GAPS}
            else:
                seqs[nm]['seq'] += ''.join(['-'] * length_added) + val['seq'] + GAPS
        last_end = matches_by_start[start]['end']
        length_to_now = len(seqs[args['name']]['seq'])

    ### adding the last segment, if there's one for the query, and the appropriate gaps
    if last_end < len(args['seq']):
        seqs[args['name']]['seq'] += args['seq'][last_end + 1:]
        end_len = len(args['seq'][last_end:])
        for v in seqs.values():
            if v['name'] != args['name']:
                seqs[v['name']]['seq'] += ''.join(['-'] * end_len)
    with open('%s.msa' % args['name'], 'wr+') as fout:
        fout.write('>%s\n%s\n' % (args['name'], seqs[args['name']]['seq']))
        for v in seqs.values():
            if v['name'] != args['name']:
                fout.write('>%s\n%s\n' % (v['name'], v['seq']))


def write_seq_to_file(args):
    print 'writing sequence to file'
    with open(args['path']+args['name'] + args['seq_file_suffix'], 'wr+') as fout:
        fout.write('>%s\n%s' % (args['name'], args['seq']))


def parse_blast_xml(args):
    """
    :param args: run arguments
    :return: makes a file with query's name with fastas from the blast. includes the query as the first entry
    """
    from Bio.Blast.NCBIXML import read

    print 'parsing BLASTP reuslts'
    xml_handle = open(args['path']+args['name'] + args['blast_suffix'], 'r')
    xml = read(xml_handle)
    fout = open(args['path']+args['name'] + args['fastas_suffix'], 'wr+')
    fout.write(">%s\n" % args['name'])
    fout.write("%s\n" % args['seq'])
    for aln in xml.alignments:
        for hsp in aln.hsps:
            # coverage = float(hsp.align_length) / float(len(args['seq'])) #  Adi's coverage.
            identity = float(hsp.identities) / float(hsp.align_length)
            evalue = float(hsp.expect)
            percent_query_len = (float(hsp.query_end) - float(hsp.query_start)) / float(len(args['seq']))
            seq = str(hsp.sbjct).translate(None, "-")
            if percent_query_len > args['blast_percent_query_len'] and identity > args['blast_identity'] \
                    and evalue < args['blast_evalue'] and len(args['seq']) * 0.8 <= len(seq) <= len(args['seq']) * 1.2:
                # fout.write(">%s_%f\n" % (aln.accession, identity))
                fout.write(">%s_%f\n" % (aln.hit_def, identity))
                fout.write("%s\n" % seq)
    xml_handle.close()
    fout.close()


def cd_hit_cluster(args):
    from os import system

    print 'RUNNING CDHIT %s' % args['path']+args['name'] + args['fastas_suffix']
    system('cd-hit -i %s -o %s -c %f ' % (args['path']+args['name'] + args['fastas_suffix'], args['path']+args['name']
                                          + args['cdhit_suffix'], args['cdhit_threshold']))


def read_multi_fastas(file):
    print 'reading fastas'
    with open(file, 'r') as f:
        cont = f.read().split('>')
    result = {}
    for entry in cont:
        split_entry = entry.split('\n')
        if len(split_entry) < 2:
            continue
        name = '_'.join(split_entry[0].rstrip().split())
        seq = ''.join(a.rstrip() for a in split_entry[1:])
        result[name] = {'name': name, 'seq': seq}
    return result


def insure_query_in_cd_hit(args, segment=None):
    print 'insuring the query is in the cdhit results'
    cdhit_fastas = read_multi_fastas(args['path']+args['name'] + args['cdhit_suffix'])
    with open(args['path']+args['name'] + args['cdhit_suffix'], 'w') as fout:
        fout.write(">%s\n" % args['name'])
        if segment is None:
            fout.write("%s\n" % args['seq'])
        else:
            fout.write("%s\n" % args['seq'][segment[0] - 1:segment[1]])
        for hit in cdhit_fastas.values():
            if hit['name'] == args['name']:
                continue
            else:
                fout.write(">%s\n" % hit['name'])
                fout.write("%s\n" % hit['seq'])


def run_muscle(args):
    import os

    print 'RUNNIN MUSCLE %s' % args['path']+args['name'] + args['cdhit_suffix']
    os.system('muscle -in ' + args['path']+args['name'] + args['cdhit_suffix'] + ' -out ' + args['path']+args['name'] +
              args['msa_suffix'])


if __name__ == '__main__':
    import argparse
    import os
    import sys
    import timeit

    parser = argparse.ArgumentParser()
    parser.add_argument('-name', type=str, help='a name for the entry')
    parser.add_argument('-uniprot', type=str, default=None, help='if available a uniprot ID can save a lot of time')
    parser.add_argument('-use_pfam', type=str, default=False,
                        help='retrieve sequences from pfam if a uniprot ID is given')
    parser.add_argument('-seq', type=str, default=None, help='entry AA sequence')
    parser.add_argument('-mode', type=str, help='which mode to run')
    parser.add_argument('-path', type=str, default=os.getcwd())
    parser.add_argument('-blast_evalue', type=float, default=0.0001, help='what evalue threshold to set for blastp')
    parser.add_argument('-blast_max_target_seqs', type=int, default=1500, help='maximal number of returned blastp seqs')
    parser.add_argument('-blast_suffix', type=str, default='_blast.xml', help='what suffix to make blast results under')
    parser.add_argument('-seq_file_suffix', type=str, default='.fasta', help='what suffix to make the lone fasta file')
    parser.add_argument('-fastas_suffix', type=str, default='_blast.fasta', help='suffix for blast results fastas')
    parser.add_argument('-blast_percent_query_len', type=float, default=0.9, help='query length in blast')
    parser.add_argument('-blast_identity', type=float, default=0.4, help='blast identity threshold')
    parser.add_argument('-cdhit_threshold', type=float, default=0.97, help='threshold for cdhit')
    parser.add_argument('-cdhit_suffix', type=str, default='.cdhit', help='suffix for cdhit results')
    parser.add_argument('-msa_suffix', type=str, default='.msa', help='suffix fro MSA file')
    parser.add_argument('-blastp_db', type=str, default='/shareDB/nr/nr', help='the database for BLASTP')
    parser.add_argument('-pfam_db', type=str, default='rp15', help='the Pfam database to use')
    parser.add_argument('-psipred_db', type=str, default='nr', help='what BLASTP DB to use for psipred')

    args = vars(parser.parse_args())

    if args['path'][-1] != '/':
        args['path'] += '/'

    if args['uniprot'] is not None and args['use_pfam']:
        if args['seq'] is None:
            args['seq'] = read_multi_fastas(args['name'] + args['seq_file_suffix'])[args['name']]['seq']
        matches = retrive_pfam_matches(args)
        args_copy = args.copy()
        print 'matches', matches
        for acc, val in matches.items():
            print 'PROCESSING MATCH %s' % acc
            tup = val['segment']
            args_copy['name'] = acc
            args_copy['fastas_suffix'] = '_%s.fasta' % args['pfam_db']
            retrieve_pfam_MSA_add_query(args, acc, tup)
            remove_bad_seqs(args_copy, tup)
            # cd_hit_cluster(args_copy)
            # insure_query_in_cd_hit(args_copy, tup)
            # run_muscle(args_copy)
            if len(matches.keys()) == 1:
                fastas = read_multi_fastas('%s_%s.fasta' % (acc, args['pfam_db']))
                with open(args['name'] + args['msa_suffix'], 'wr+') as fout:
                    fout.write('>%s\n%s\n' % (args['name'], fastas[acc]['seq']))
                    for v in fastas.values():
                        if v['name'] != acc:
                            fout.write('>%s\n%s\n' % (v['name'], v['seq']))
            else:
                matches[acc]['msa'] = read_multi_fastas('%s_%s.fasta' % (acc, args['pfam_db']))
        if len(matches.keys()) > 1:
            concatenate_msas(args, matches)
    if args['seq'] is None:
        args['seq'] = read_multi_fastas(args['name'] + args['seq_file_suffix'])[args['name']]['seq']

    if args['mode'] == 'seq2msa':
        write_seq_to_file(args)
        run_blastp(args)
        tic = timeit.default_timer()
        parse_blast_xml(args)
        cd_hit_cluster(args)
        insure_query_in_cd_hit(args)
        run_muscle(args)
        print 'process took %f seconds' % (timeit.default_timer() - tic)

    elif args['mode'] == 'seq_file2msa':
        tic_blast = timeit.default_timer()
        run_blastp_ncbi(args)
        tic = timeit.default_timer()
        parse_blast_xml(args)
        cd_hit_cluster(args)
        insure_query_in_cd_hit(args)
        run_muscle(args)
        print 'blast took %f seconds' % (tic - tic_blast)
        print 'process took %f seconds' % (timeit.default_timer() - tic)
        print 'all took %f seconds' % (timeit.default_timer() - tic_blast)

    elif args['mode'] == 'xml2msa':
        tic = timeit.default_timer()
        parse_blast_xml(args)
        cd_hit_cluster(args)
        insure_query_in_cd_hit(args)
        run_muscle(args)
        print 'process took %f seconds' % (timeit.default_timer() - tic)

    # elif args['mode'] == 'seq_file2msa_psipred':
    #     tic = timeit.default_timer()
    #     run_psipred(args)
    #     toc_psipred = timeit.default_timer()
    #     gather_psipred_seqs(args)
    #     cd_hit_cluster(args)
    #     insure_query_in_cd_hit(args)
    #     run_muscle(args)
    #     with open('%s.log' % args['name'], 'wr+') as fout:
    #         fout.write('process took %f seconds\n' % (timeit.default_timer()-tic))
    #         fout.write('psipred took %f seconds\n' % (toc_psipred-tic))
    #         fout.write('processing after psipred took %f seconds\n' % (timeit.default_timer()-toc_psipred))

    elif args['mode'] == 'psipred2msa':
        tic = timeit.default_timer()
        gather_psipred_seqs(args)
        cd_hit_cluster(args)
        insure_query_in_cd_hit(args)
        run_muscle(args)
        with open('%s.log' % args['name'], 'wr+') as fout:
            fout.write('process took %f seconds\n' % (timeit.default_timer()-tic))

    elif args['mode'] == 'psipred2msa_efetch':
        tic = timeit.default_timer()
        retrive_psipred_efetch(args)
        cd_hit_cluster(args)
        insure_query_in_cd_hit(args)
        run_muscle(args)
        with open('%s.log' % args['name'], 'wr+') as fout:
            fout.write('process took %f seconds\n' % (timeit.default_timer()-tic))

    elif args['mode'] == 'seq_file2msa_psipred':
        tic = timeit.default_timer()
        run_psipred(args)
        toc_psipred = timeit.default_timer()
        retrive_psipred_efetch(args)
        cd_hit_cluster(args)
        insure_query_in_cd_hit(args)
        run_muscle(args)
        with open('%s.log' % args['name'], 'wr+') as fout:
            fout.write('process took %f seconds\n' % (timeit.default_timer()-tic))
            fout.write('psipred took %f seconds\n' % (toc_psipred-tic))
            fout.write('processing after psipred took %f seconds\n' % (timeit.default_timer()-toc_psipred))

    elif args['mode'] == 'psipred2msa_blast':
        tic = timeit.default_timer()
        retrive_psipred_efetch(args)
        make_blast_db(args)
        args['blastp_db'] = args['name']+'_blastdb'
        run_blastp(args)
        parse_blast_xml(args)
        cd_hit_cluster(args)
        insure_query_in_cd_hit(args)
        run_muscle(args)
        [os.remove(a) for a in os.listdir('./') if a.split('.')[1] in ['blast', 'phr', 'pin', 'psq']and a[:3] != 'psi']
        with open('%s.log' % args['name'], 'wr+') as fout:
            fout.write('process took %f seconds\n' % (timeit.default_timer()-tic))

    elif args['mode'] == 'test':
        retrive_psipred_efetch(args)

    else:
        print 'no mode given'
