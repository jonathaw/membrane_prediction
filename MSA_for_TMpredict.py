class single_fasta():
    """
    a class for describing an AA sequence
    """
    def __init__(self, name, seq):
        """
        :param name: sequence name
        :param seq: sequence
        :return: a class instance
        """
        self.name = name
        self.seq = seq
        self.gaps = [i for i, aa in enumerate(seq) if aa == '-']

    def __str__(self):
        return 'name: %s\nseq: %s\n' % (self.name, self.seq)

    def __repr__(self):
        return 'name: %s\nseq: %s\n' % (self.name, self.seq)

    def nogap2withgap(self, pos):
        """
        :param pos: a position coordinate in the original, no gaps, sequence
        :return: the position corresponding to the same pos in the seq with gaps
        """
        gaps = 0
        for i, aa in enumerate(self.seq):
            gaps += 1 if aa == '-' else 0
            if i-gaps == pos-1:
                return i

    def withgap2nogap(self, pos):
        """
        :param pos: a position coordinate in the aligned seq with gaps
        :return: the position corresponding to the same pos in the seq without gaps
        """
        gaps = 0
        for i, aa in enumerate(self.seq):
            gaps += 1 if aa == '-' else 0
            if i == pos:
                return i-gaps


def target_has_gaps_in_query_stretch(query, target, start_wg, end_wg):
    """
    :param query: query sequence
    :param target: target sequence
    :param start_wg: start position for segemnt, with gaps
    :param end_wg: end position for segment, with gaps
    :return: True if gaps pattern in target matches that of the query in the given segment, else False
    """
    query_gaps_in_stretch = [i for i in query.gaps if start_wg <= i <= end_wg]
    target_gaps_in_stretch = [i for i in target.gaps if start_wg <= i <= end_wg]
    xs = target.seq[start_wg:end_wg].count('X') == 0
    return True if query_gaps_in_stretch == target_gaps_in_stretch and xs else False


class TMpredict_MSA():
    '''
    A class for handling MSA input in fasta format. reads in the sequences as single_fasta objects
    '''
    def __init__(self, name, polyval, poly_param):
        # path_msa = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/rost_msa_prep/msa_ready/'
        # msa_file_name = name+'_ready.fa'
        path_msa = '/home/labs/fleishman/jonathaw/membrane_prediction_DBs/blasts_adi/'
        msa_file_name = name+'.fasta_msa.aln'
        self.polyval = polyval
        self.poly_param = poly_param
        self.stack = read_fasta_msa(path_msa+msa_file_name)
        self.query = [a for a in self.stack if a.name == name][0]

    def retrieve_seqs(self, start, end, direction):
        from WinGrade import WinGrade
        start_wg = self.query.nogap2withgap(start)
        end_wg = self.query.nogap2withgap(end)
        if direction == 'fwd':
            best_win = WinGrade(start, end, direction, gap_remover(self.query.seq[start_wg:end_wg]), self.polyval,
                                self. poly_param)
        if direction == 'rev':
            best_win = WinGrade(start, end, direction, gap_remover(self.query.seq[start_wg:end_wg])[::-1],
                                self.polyval, self. poly_param)
        query_grade = best_win.grade
        name = 'query'
        for target in self.stack:
            # print 'q seq', self.query.seq[start_wg:end_wg]
            # print 't seq', target.seq[start_wg:end_wg]
            if target_has_gaps_in_query_stretch(self.query, target, start_wg, end_wg):
                # print 'processing', target.seq[start_wg:end_wg]
                if direction == 'fwd':
                    temp_win_grade = WinGrade(start, end, direction, gap_remover(target.seq[start_wg:end_wg]),
                                              self.polyval, self.poly_param)
                elif direction == 'rev':
                    temp_win_grade = WinGrade(start, end, direction, gap_remover(target.seq[start_wg:end_wg])[::-1],
                                              self.polyval, self.poly_param)
                # if target.name == self.query.name:
                    # print 'found the query!!!', temp_win_grade
                if temp_win_grade.grade < best_win.grade and abs(temp_win_grade.grade-query_grade) < 7:
                    best_win = temp_win_grade
                    name = target.name
        if direction == 'fwd':
            return WinGrade(start, end, direction, gap_remover(self.query.seq[start_wg:end_wg]), self.polyval,
                            self.poly_param, name, gap_remover(best_win.seq))
        elif direction == 'rev':
            return WinGrade(start, end, direction, gap_remover(self.query.seq[start_wg:end_wg])[::-1], self.polyval,
                            self.poly_param, name, gap_remover(best_win.seq)[::-1])


def gap_remover(seq):
    return ''.join([aa for aa in seq if aa != '-'])


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
    from TMpredict_WinGrade import MakeHydrophobicityGrade

    # read_fasta_msa(sys.argv[1])
    TMpredict_MSA(sys.argv[1], MakeHydrophobicityGrade(), {'c0': 3., 'c1': 6.8, 'c2': 0.7, 'c3': -0.03})