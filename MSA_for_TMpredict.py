class single_fasta:
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
        if pos == 0:
            return 0
        gaps = 0
        for i, aa in enumerate(self.seq):
            gaps += 1 if aa == '-' else 0
            if i - gaps == pos - 1:
                return i
        if i + 1 == pos:
            return i - gaps - 1

        return len(self.seq) # if not found, returns the end of the seq

    def withgap2nogap(self, pos):
        """
        :param pos: a position coordinate in the aligned seq with gaps
        :return: the position corresponding to the same pos in the seq without gaps
        """
        if pos == 0:
            return 0
        gaps = 0
        for i, aa in enumerate(self.seq):
            gaps += 1 if aa == '-' else 0
            if i == pos:
                return i - gaps
        if i + 1 == pos:
            return i - gaps - 1


def target_has_gaps_in_query_stretch(query, target, start_wg, end_wg):
    """
    :param query: query sequence
    :param target: target sequence
    :param start_wg: start position for segemnt, with gaps
    :param end_wg: end position for segment, with gaps
    :return: True if gaps pattern in target matches that of the query in the given segment, else False
    >>> q = single_fasta('query', 'ABC--DE-G')
    >>> t = single_fasta('targe', 'AB--FDEG-')
    >>> target_has_gaps_in_query_stretch(q, t, 0, 8)
    False
    >>> q = single_fasta('query', 'ABC--DE-FG')
    >>> t = single_fasta('targe', 'DFG--YU-LK')
    >>> target_has_gaps_in_query_stretch(q, t, 0, 8)
    True
    """
    query_gaps_in_stretch = [i for i in query.gaps if start_wg <= i <= end_wg]
    target_gaps_in_stretch = [i for i in target.gaps if start_wg <= i <= end_wg]
    xs = target.seq[start_wg:end_wg].count('X') == 0
    return True if query_gaps_in_stretch == target_gaps_in_stretch and xs else False


class TMpredict_MSA:
    '''
    A class for handling MSA input in fasta format. reads in the sequences as single_fasta objects
    '''
    def __init__(self, name, polyval, poly_param, PERCENTILE, path_msa):
        # path_msa = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/rost_msa_prep/msa_ready/'
        # msa_file_name = name+'_ready.fa'
        #path_msa = '/home/labs/fleishman/jonathaw/membrane_prediction_DBs/blasts_adi/'
        path_msa = path_msa
        msa_file_name = name + '.msa'
        self.polyval = polyval
        self.poly_param = poly_param
        self.stack = read_fasta_msa(path_msa + msa_file_name)
        self.query = [a for a in self.stack if a.name == name][0]
        self.percentile = PERCENTILE

    def retrieve_seqs_old(self, start, end, direction):
        from WinGrade import WinGrade

        start_wg = self.query.nogap2withgap(start)
        end_wg = self.query.nogap2withgap(end)
        if direction == 'fwd':
            best_win = WinGrade(start, end, direction, gap_remover(self.query.seq[start_wg:end_wg]), self.polyval,
                                self.poly_param)
        if direction == 'rev':
            best_win = WinGrade(start, end, direction, gap_remover(self.query.seq[start_wg:end_wg])[::-1],
                                self.polyval, self.poly_param)
        query_grade = best_win.grade
        name = 'query'
        for target in self.stack:
            if target_has_gaps_in_query_stretch(self.query, target, start_wg, end_wg):
                if direction == 'fwd':
                    temp_win_grade = WinGrade(start, end, direction, gap_remover(target.seq[start_wg:end_wg]),
                                              self.polyval, self.poly_param)
                elif direction == 'rev':
                    temp_win_grade = WinGrade(start, end, direction, gap_remover(target.seq[start_wg:end_wg])[::-1],
                                              self.polyval, self.poly_param)
                    # if target.name == self.query.name:
                    # print 'found the query!!!', temp_win_grade
                if temp_win_grade.grade < best_win.grade and abs(temp_win_grade.grade - query_grade) < 0:
                    best_win = temp_win_grade
                    name = target.name
        if direction == 'fwd':
            return WinGrade(start, end, direction, gap_remover(self.query.seq[start_wg:end_wg]), self.polyval,
                            self.poly_param, name, gap_remover(best_win.seq))
        elif direction == 'rev':
            return WinGrade(start, end, direction, gap_remover(self.query.seq[start_wg:end_wg])[::-1], self.polyval,
                            self.poly_param, name, gap_remover(best_win.seq)[::-1])


def retrieve_seqs(msa, start, end, direction, inner_tail_length=5, outer_tail_length=5):
    from WinGrade import WinGrade

    start_wg = msa.query.nogap2withgap(start)
    end_wg = msa.query.nogap2withgap(end)

    # calculate inner tail KR addition
    if direction == 'fwd':
        inner_tail_start_wg = max([0, msa.query.nogap2withgap(start-inner_tail_length)])
        inner_tail_end_wg = msa.query.nogap2withgap(start)

        outer_tail_start_wg = msa.query.nogap2withgap(start)+1
        outer_tail_end_wg = min([msa.query.nogap2withgap(len(msa.query.seq)),
                                 msa.query.nogap2withgap(start)+outer_tail_length])+1
        # print 'AAAAAAAAAAAAAAA', inner_tail_start_wg, inner_tail_end_wg, start
    else:
        inner_tail_start_wg = msa.query.nogap2withgap(end)+1
        inner_tail_end_wg = min([msa.query.nogap2withgap(len(msa.query.seq)),
                                 msa.query.nogap2withgap(end+inner_tail_length)])+1

        outer_tail_start_wg = max(0, msa.query.nogap2withgap(start)-outer_tail_length)
        outer_tail_end_wg = msa.query.nogap2withgap(start)
        # print 'BBBBBBBBBBBBBBB', inner_tail_start_wg, inner_tail_end_wg

    if direction == 'fwd':
        query_grade = WinGrade(start, end, direction, gap_remover(msa.query.seq[start_wg:end_wg]), msa.polyval,
                            msa.poly_param,
                               inner_tail=gap_remover(msa.query.seq[inner_tail_start_wg:inner_tail_end_wg]),
                               outer_tail=gap_remover(msa.query.seq[outer_tail_start_wg:outer_tail_end_wg])).grade
    elif direction == 'rev':
        query_grade = WinGrade(start, end, direction, gap_remover(msa.query.seq[start_wg:end_wg])[::-1],
                            msa.polyval, msa.poly_param,
                               inner_tail=gap_remover(msa.query.seq[inner_tail_start_wg:inner_tail_end_wg]),
                               outer_tail=gap_remover(msa.query.seq[outer_tail_start_wg:outer_tail_end_wg])).grade
    grade_stack = {}

    for target in msa.stack:
        if target_has_gaps_in_query_stretch(msa.query, target, start_wg, end_wg):
            if direction == 'fwd':
                temp_win_grade = WinGrade(start, end, direction, gap_remover(target.seq[start_wg:end_wg]),
                                          msa.polyval, msa.poly_param,
                                          inner_tail=gap_remover(target.seq[inner_tail_start_wg:inner_tail_end_wg]),
                                          outer_tail=gap_remover(msa.query.seq[outer_tail_start_wg:outer_tail_end_wg]))
                # print 'A', inner_tail_start_wg, inner_tail_end_wg, gap_remover(target.seq[inner_tail_start_wg:inner_tail_end_wg])
            elif direction == 'rev':
                temp_win_grade = WinGrade(start, end, direction, gap_remover(target.seq[start_wg:end_wg])[::-1],
                                          msa.polyval, msa.poly_param,
                                          inner_tail=gap_remover(target.seq[inner_tail_start_wg:inner_tail_end_wg]),
                                          outer_tail=gap_remover(msa.query.seq[outer_tail_start_wg:outer_tail_end_wg]))
                # print 'B', inner_tail_start_wg, inner_tail_end_wg, gap_remover(target.seq[inner_tail_start_wg:inner_tail_end_wg])

            if abs(temp_win_grade.grade-query_grade) <= msa.poly_param['msa_threshold'] and temp_win_grade.seq != '':
                grade_stack[temp_win_grade.grade] = {'win': temp_win_grade, 'name': target.name}

    # med_grade = median_low(grade_stack.keys())
    med_grade = precentile_for_wins(grade_stack.keys(), msa.percentile)

    # med_grade = min(grade_stack.keys())
    med_win = grade_stack[med_grade]['win']
    med_name = grade_stack[med_grade]['name']

    if direction == 'fwd':
        return WinGrade(start, end, direction, gap_remover(msa.query.seq[start_wg:end_wg]), msa.polyval,
                        msa.poly_param, med_name, gap_remover(med_win.seq),
                        inner_tail=gap_remover(msa.query.seq[inner_tail_start_wg:inner_tail_end_wg]),
                        outer_tail=gap_remover(msa.query.seq[outer_tail_start_wg:outer_tail_end_wg]))

    elif direction == 'rev':
        return WinGrade(start, end, direction, gap_remover(msa.query.seq[start_wg:end_wg])[::-1], msa.polyval,
                        msa.poly_param, med_name, gap_remover(med_win.seq),
                        inner_tail=gap_remover(msa.query.seq[inner_tail_start_wg:inner_tail_end_wg]),
                        outer_tail=gap_remover(msa.query.seq[outer_tail_start_wg:outer_tail_end_wg]))


def median_low(list_):
    if len(list_) % 2 != 0:
        from numpy import median
        return float(median(list_))
    else:
        list_.sort()
        return list_[int(len(list_)/2 - 0.5)]


def precentile_for_wins(list_, PERCENTILE):
    import numpy as np
    a = np.array(list_)

    return np.percentile(a, PERCENTILE, interpolation='lower')


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
    # from TMpredict_WinGrade import MakeHydrophobicityGrade
    pass
    # read_fasta_msa(sys.argv[1])
    # TMpredict_MSA(sys.argv[1], MakeHydrophobicityGrade(), {'c0': 3., 'c1': 6.8, 'c2': 0.7, 'c3': -0.03})