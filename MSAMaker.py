#!/usr/bin/env python
def run_blastp(args):
    from os import system
    print 'running BLASTP'
    system('blastp -query %s -db /shareDB/nr/nr -evalue %f -max_target_seqs %i -outfmt 5 -out %s' %
           (args['name']+args['seq_file_suffix'], args['blast_evalue'], args['blast_max_target_seqs'],
            args['name']+args['blast_suffix']))


def write_seq_to_file(args):
    print 'writing sequence to file'
    with open(args['name']+args['seq_file_suffix'], 'wr+') as fout:
        fout.write('>%s\n%s' % (args['name'], args['seq']))


def parse_blast_xml(args):
    """
    :param args: run arguments
    :return: makes a file with query's name with fastas from the blast. includes the query as the first entry
    """
    from Bio.Blast.NCBIXML import read
    print 'parsing BLASTP reuslts'
    xml_handle = open(args['name']+args['blast_suffix'], 'r')
    xml = read(xml_handle)
    fout = open(args['name']+args['fastas_suffix'], 'wr+')
    fout.write(">%s\n" % args['name'])
    fout.write("%s\n" % args['seq'])
    for aln in xml.alignments:
        for hsp in aln.hsps:
            coverage = float(hsp.align_length) / float(len(args['seq'])) #  Adi's coverage.
            identity = float(hsp.identities) / float(hsp.align_length)
            evalue = float(hsp.expect)
            percent_query_len = (float(hsp.query_end)-float(hsp.query_start)) / float(len(args['seq']))
            seq = str(hsp.sbjct).translate(None, "-")
            if percent_query_len > args['blast_percent_query_len'] and identity > args['blast_identity'] \
                    and evalue < args['blast_evalue'] and len(args['seq'])*0.9 <= len(seq) <= len(args['seq'])*1.1:
                fout.write(">%s_%f\n" % (aln.accession, identity))
                fout.write("%s\n" % seq)
    xml_handle.close()
    fout.close()


def cd_hit_cluster(args):
    from os import system
    print 'running cdhit'
    system('cd-hit -i %s -o %s -c %f ' % (args['name']+args['fastas_suffix'], args['name']+args['cdhit_suffix'],
                                          args['cdhit_threshold']))


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


def insure_query_in_cd_hit(args):
    print 'insuring the query is in the cdhit results'
    cdhit_fastas = read_multi_fastas(args['name']+args['cdhit_suffix'])
    with open(args['name']+args['cdhit_suffix'], 'w') as fout:
        fout.write(">%s\n" % args['name'])
        fout.write("%s\n" % args['seq'])
        for hit in cdhit_fastas.values():
            if hit['name'] == args['name']:
                continue
            else:
                fout.write(">%s\n" % hit['name'])
                fout.write("%s\n" % hit['seq'])


def run_muscle(args):
    import os
    print 'running muscle'
    os.system('muscle -in ' + args['name']+args['cdhit_suffix'] + ' -out ' + args['name'] + args['msa_suffix'])


if __name__ == '__main__':
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', type=str, help='a name for the entry')
    parser.add_argument('-seq', type=str, help='entry AA sequence')
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

    args = vars(parser.parse_args())

    if args['mode'] == 'seq2msa':
        write_seq_to_file(args)
        run_blastp(args)
        parse_blast_xml(args)
        cd_hit_cluster(args)
        insure_query_in_cd_hit(args)
        run_muscle(args)

    elif args['mode'] == 'seq_file2msa':
        run_blastp(args)
        parse_blast_xml(args)
        cd_hit_cluster(args)
        insure_query_in_cd_hit(args)
        run_muscle(args)

    elif args['mode'] == 'xml2msa':
        parse_blast_xml(args)
        cd_hit_cluster(args)
        insure_query_in_cd_hit(args)
        run_muscle(args)
