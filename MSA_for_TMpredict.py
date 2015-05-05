class single_fasta():
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.gaps = [i for i, aa in enumerate(seq) if aa == '-']

    def __str__(self):
        return 'name: %s\nseq: %s\n' % (self.name, self.seq)

    def __repr__(self):
        return 'name: %s\nseq: %s\n' % (self.name, self.seq)

    def nogap2withgap(self, pos):
        gaps = 0
        for i, aa in enumerate(self.seq):
            gaps += 1 if aa == '-' else 0
            if i-gaps == pos-1:
                return i

    def withgap2nogap(self, pos):
        gaps = 0
        for i, aa in enumerate(self.seq):
            gaps += 1 if aa == '-' else 0
            if i == pos:
                return i-gaps


def target_has_gaps_in_query_stretch(query, target, start_wg, end_wg):
    query_gaps_in_stretch = [i for i in query.gaps if start_wg <= i <= end_wg]
    target_gaps_in_stretch = [i for i in target.gaps if start_wg <= i <= end_wg]
    return True if query_gaps_in_stretch == target_gaps_in_stretch else False



class TMpredict_MSA():
    def __init__(self, name):
        path_msa = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/MSA_fastas/'
        msa_file_name = name+'_MSA_fastas_muscle.fa'
        self.stack = read_fasta_msa(path_msa+msa_file_name)
        self.query = [a for a in self.stack if a.name == name][0]
        # print '\n\nquery:', self.query
        # print self.query.nogap2withgap(9)
        # print self.query.withgap2nogap(4)
        self.retrieve_seqs(0, 74)

    def retrieve_seqs(self, start, end):
        start_wg = self.query.nogap2withgap(start)
        end_wg = self.query.nogap2withgap(end)
        for target in self.stack:
            print 'q seq', self.query.seq[start_wg:end_wg]
            print 't seq', target.seq[start_wg:end_wg]
            print target_has_gaps_in_query_stretch(self.query, target, start_wg, end_wg)


def read_fasta_msa(file_name):
    with open(file_name) as f:
        cont = f.read().split('>')
    result = []
    for para in cont:
        split = para.split('\n')
        if split == []: continue
        result.append(single_fasta(split[0], ''.join(split[1:]).rstrip()))
    return result


if __name__ == '__main__':
    import sys
    # read_fasta_msa(sys.argv[1])
    TMpredict_MSA(sys.argv[1])