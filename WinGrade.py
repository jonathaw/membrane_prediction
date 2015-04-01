class WinGrade():
    '''
    A class to parameterise a window in hydrophobicity manners.
    '''
    def __init__(self, begin, end, direction, seq, polyval):
        '''
        :param begin: seq position (from 0) where the window begins
        :param end: seq position (from 0) where the window ends
        :param grade: the hydrophobicity grade calculated using the energy function
        :param direction: either fwd or rev
        :param seq: the sequence of the segment
        '''
        self.begin = begin
        self.end = end
        self.seq = seq
        self.length = end-begin
        self.grade = self.grade_segment(polyval)
        self.direction = direction
        self.span = range(self.begin, self.end+1)
        # self.grade_norm = self.grade / (self.end-self.begin-19)

    def __str__(self):
        return '%-4i to %-4i in %3s => %10f %10f %-35s' % (self.begin, self.end, self.direction, self.grade,
                                                           self.grade_norm, self.seq)

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
        return grade

    def pos_in_wingrade(self, pos):
        return True if pos >= self.begin and pos <= self.end else False

    def set_grade(self, val):
        self.grade = self.grade + val