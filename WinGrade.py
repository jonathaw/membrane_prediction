class WinGrade():
    '''
    A class to parameterise a window in hydrophobicity manners.
    '''
    def __init__(self, begin, end, direction, seq, polyval, poly_param):
        """
        :param begin: seq position (from 0) where the window begins
        :param end: seq position (from 0) where the window ends
        :param grade: the hydrophobicity grade calculated using the energy function
        :param direction: either fwd or rev
        :param seq: the sequence of the segment
        """
        self.begin = begin
        self.end = end
        self.seq = seq
        self.length = end-begin
        # self.poly_param = {'c0': 0.27, 'c1': 9.29, 'c2': -0.645, 'c3': 0.00822}

        self.poly_param = poly_param
        self.length_element = self.length_polynom()
        self.hp_moment = hp_moment(seq, polyval, poly_param)
        self.hp_sum = self.grade_segment(polyval)
        self.grade = self.hp_sum + self.hp_moment + self.length_element

        self.direction = direction
        self.span = range(self.begin, self.end+1)
        # self.grade_norm = self.grade / (self.end-self.begin-19)

    def __str__(self):
        return '%-4i to %-4i in %3s => %10f %-35s' % (self.begin, self.end, self.direction, self.grade,
                                                      self.seq)

    def __repr__(self):
        return '%-4i to %-4i in %3s => %10f %-35s' % (self.begin, self.end, self.direction, self.grade, self.seq)

    def grade_grade_colliding(self, other):
        return True if len(set(self.span) & set(other.span)) != 0 else False

    def set_grade_colliding(self, grade_set):
        for grade in grade_set:
            if self.grade_grade_colliding(grade):
                return True
        return False

    def print_wingrade(self):
        print '%i to %i => %f' % (self.begin, self.end, self.grade)

    def grade_segment(self, polyval):
        import numpy as np
        membrane_position = np.linspace(-15, 15, endpoint=True, num=self.length)
        grade = 0
        for i, aa in enumerate(self.seq):
            grade += np.polyval(polyval[aa], membrane_position[i])
        return grade    #+self.hp_moment(polyval)+self.length_polynom()

    def length_polynom(self):
        """
        :return: length polynomial as in hp_moment article
        """
        poly_param = self.poly_param
        l = self.length
        return poly_param['c1'] + poly_param['c2']*l + poly_param['c3']*l**2

    def pos_in_wingrade(self, pos):
        return True if pos >= self.begin and pos <= self.end else False

    def set_grade(self, val):
        self.grade = self.grade + val


def hp_moment(seq, polyval, poly_param):
    """
    :param polyval: the hydrophobicity scal polynom values
    :return: the hydrophobicity moment of the helix, calculated as in
    von Heijne, G. (2007) Molecular code for transmembrane-helix recognition by the Sec61 translocon.
    doi:10.1038/nature06387
    """
    import numpy as np
    deg100_in_rad = (100.*2.*np.pi) / 360.
    # membrane_position = np.linspace(-15, 15, endpoint=True, num=len(seq))
    sin_elem = 0.
    cos_elem = 0.
    for i, aa in enumerate(seq):
        if aa == 'G' or aa == 'P' or np.polyval(polyval[aa], 0) < 0.:
            continue
        sin_elem += np.polyval(polyval[aa], 0) * np.sin(deg100_in_rad*i)
        cos_elem += np.polyval(polyval[aa], 0) * np.cos(deg100_in_rad*i)
    return poly_param['c0']*np.sqrt(sin_elem**2 + cos_elem**2)


def grade_segment(seq, polyval):
    import numpy as np
    membrane_position = np.linspace(-15, 15, endpoint=True, num=len(seq))
    grade = 0
    for i, aa in enumerate(seq):
        grade += np.polyval(polyval[aa], membrane_position[i])
    return grade