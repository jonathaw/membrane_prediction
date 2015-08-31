#!/usr/bin/env python
class WinGrade():
    '''
    A class to parameterise a window in hydrophobicity manners.
    '''
    def __init__(self, begin, end, direction, seq, polyval=None, poly_param=None, msa_name=None, msa_seq=None,
                 grade=None, length_element=None, charges=None):
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
        self.length = len(seq)
        if poly_param is not None:
            self.poly_param = poly_param
        else:
            self.poly_param = {'w': 0, 'z_0': 0}
        if polyval is None:
            from TMpredict_WinGrade import MakeHydrophobicityGrade
            polyval = MakeHydrophobicityGrade()
        # self.length_element = self.length_polynom()
        if length_element is None:
            self.length_element = membrane_deformation(self.length, self.poly_param['w'], self.poly_param['z_0'])
        else:
            self.length_element = length_element
        # self.hp_moment = hp_moment(seq, polyval, poly_param)
        # self.hp_sum = self.grade_segment(polyval)
          # self.hp_sum + self.hp_moment + self.length_element
        if charges is None:
            self.charges = count_charges(self.seq)
        else:
            self.charges = charges
        # if self.charges >= 2:
        #     print 'found %i charges' % self.charges
            # raise Exception()
        if grade is None:
            self.grade = self.grade_segment(polyval) + self.length_element
        else:
            self.grade = grade + self.length_element
        self.direction = direction
        self.span = range(self.begin, self.end+1)
        # self.grade_norm = self.grade / (self.end-self.begin-19)

        self.msa_name = msa_name
        if msa_name != None:
            self.msa_seq = msa_seq
            self.msa_grade = grade_segment(msa_seq, polyval) + length_polynom(msa_seq, poly_param) #  + hp_moment(msa_seq, polyval, poly_param) \
                             #  + length_polynom(msa_seq, poly_param)
        if self.seq == 'SOURCE' or self.seq == 'SINK' or self.seq == 'SOURCEPATH':
            self.grade = 0.0

    def __repr__(self):
        if self.msa_name == None:
            return '%-4i to %-4i in %3s => %10f %-35s %i %f' % (self.begin, self.end, self.direction, self.grade, self.seq, self.charges, self.length_element)
        else:
            return '%-4i to %-4i in %3s => %10f %-35s %s %s %f' % (self.begin, self.end, self.direction, self.grade,
                                                               self.seq, self.msa_name, self.msa_seq, self.msa_grade)

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
        membrane_position = np.linspace(-20, 20, endpoint=True, num=self.length)
        grade = 0
        for i, aa in enumerate(self.seq):
            if aa in polyval.keys():
                grade += np.polyval(polyval[aa], membrane_position[i])
        return grade    #+self.hp_moment(polyval)+self.length_polynom()

    def length_polynom(self):
        """
        :return: length polynomial as in hp_moment article
        """
        poly_param = self.poly_param
        #print 'length', poly_param
        l = self.length
        return poly_param['c1'] + poly_param['c2']*l + poly_param['c3']*(l**2)

    def pos_in_wingrade(self, pos):
        return True if pos >= self.begin and pos <= self.end else False

    def set_grade(self, val):
        self.grade = self.grade + val

    def same_as_other(self, other):
        """
        :param other:
        :return:
        >>> w1 = WinGrade(0, 10, 'fwd', 'A', grade=1, length_element=1, charges=1)
        >>> w2 = WinGrade(20, 30, 'rev', 'B', grade=2, length_element=2, charges=2)
        >>> w1.same_as_other(w1)
        True
        >>> w1.same_as_other(w2)
        False
        """
        return self.__dict__ == other.__dict__

    def within_segment(self, segment, flaps=0):
        return self.begin >= segment[0]-flaps and self.end <= segment[1]+flaps

    def within_segments(self, segments, flaps=0):
        for seg in segments:
            if self.within_segment(seg, flaps):
                return True
        return False
# a=WinGradePath([WinGrade(1,2,'',''), WinGrade(3,4,'',''), WinGrade(100,110,'',''), WinGrade(50,60,'','')])
# a.add_win(WinGrade(30, 40,'fwd', 'AAAAAAAAAAAAAAAA'))
class WinGradePath():
    def __init__(self, win_list):
        import operator
        self.path = sorted(win_list, key=operator.attrgetter('begin'))
        self.total_grade = self.grade_path()
        self.win_num = self.path_length()
        if win_list != []:
            self.c_term = win_list[-1].direction
        else:
            self.c_term = None

    def __repr__(self):
        msg = '{ ['
        for i, win in enumerate(self.path):
            if win.seq == '':
                continue
            msg += str(win)
            if i < len(self.path)-1:
                msg += ' : '
            else:
                msg += ']'
        msg += ' }~> total_grade %10f win_num %2i c_term %s' % (self.total_grade, self.win_num, self.c_term)
        return msg

    def add_win(self, win):
        import operator
        wins = self.path[:] + [win]
        self.path = sorted(wins, key=operator.attrgetter('begin'))
        self.total_grade = self.grade_path()
        self.win_num = self.path_length()
        self.c_term = self.path[-1].direction

    def grade_path(self):
        return sum(a.grade for a in self.path if a.seq != '')

    def path_length(self):
        return len([a for a in self.path if a.seq != ''])

    def same_as_other(self, other):
        for k1, v1 in self.__dict__.items():
            for k2, v2 in other.__dict__.items():
                if k1 == k2:
                    if k1 == 'path':
                        for w1 in v1:
                            booli = False
                            for w2 in v2:
                                if w1.same_as_other(w2):
                                    booli = True
                                    break
                            if not booli:
                                return False

                    else:
                        if v1 != v2:
                            return False
        return True

    def first(self):
        return self.path[0]

    def last(self):
        return self.path[-1]

    def overlap_me_first(self, other):
        """
        :param other: another path
        :return: True iff the last win of self is the same as the first of other
        """
        return self.last().same_as_other(other.first())

    def overlap_me_last(self, other):
        """
        :param other: another path
        :return: True iff the first win of self is the same as the last of other
        """
        return self.first().same_as_other(other.last())


def count_charges(seq):
    import numpy as np
    membrane_position = np.linspace(-20, 20, endpoint=True, num=len(seq))
    charges = 0
    for i, aa in enumerate(seq):
        if -12.5 <= membrane_position[i] <= 12.5:
            if aa in ['E', 'D', 'K', 'R', 'N', 'Q', 'H']:
                charges += 1
    return charges


def membrane_deformation(l, w, z_0):
    '''
    :param l: helix length in AAs
    :param w: spring constant
    :param z_0: neutral membrane width
     taken from Ben-Tal - Kessel http://bental.tau.ac.il/MCPep/overview.html
    :return:membrane deformation penalty, dGdef
    >>> membrane_deformation(40, 0.001, 20)
    1.6
    '''
    return w * ((1.5*l - z_0) ** 2)


def hp_moment(seq, polyval, poly_param):
    """
    :param polyval: the hydrophobicity scal polynom values
    :return: the hydrophobicity moment of the helix, calculated as in
    von Heijne, G. (2007) Molecular code for transmembrane-helix recognition by the Sec61 translocon.
    doi:10.1038/nature06387
    """
    import numpy as np
    #print 'hp_moment', poly_param
    deg100_in_rad = (100.*2.*np.pi) / 360.
    # membrane_position = np.linspace(-15, 15, endpoint=True, num=len(seq))
    sin_elem = 0.
    cos_elem = 0.
    for i, aa in enumerate(seq):
        if aa not in polyval.keys():
            continue
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


def length_polynom(seq, poly_param):
    """
    :return: length polynomial as in hp_moment article
    """
    # poly_param = self.poly_param
    #print 'length', poly_param
    l = len(seq)
    return poly_param['c1'] + poly_param['c2']*l + poly_param['c3']*(l**2)


def flatten_path_list(path_list):
    result = []
    for lst in path_list:
        for win in lst.path:
            if win not in result:
                result.append(win)
    return result


def sequential_coiled(pos, psi, verbose=False):
    """
    :param pos: tuple of start and end of win
    :param psi: dictionary with psipred resutls
    :return: True iff a stretch of 3 residues at either end of the win is coiled > 0.5
    >>> pos = (0, 10)
    >>> psi = {0: {'c': 0.6}, 1: {'c': 0.6}, 2: {'c': 0.6}, 3: {'c': 0.6}, 4: {'c': 0.6}, 5: {'c': 0.6}, 6: {'c': 0.6}, 7: {'c': 0.6}, 8: {'c': 0.6}, 9: {'c': 0.6}, 10: {'c': 0.6}}
    >>> sequential_coiled(pos, psi)
    True
    >>> psi = {0: {'c': 0.4}, 1: {'c': 0.4}, 2: {'c': 0.4}, 3: {'c': 0.6}, 4: {'c': 0.6}, 5: {'c': 0.6}, 6: {'c': 0.6}, 7: {'c': 0.6}, 8: {'c': 0.4}, 9: {'c': 0.4}, 10: {'c': 0.4}}
    >>> sequential_coiled(pos, psi)
    False
    >>> psi = {0: {'c': 0.4}, 1: {'c': 0.4}, 2: {'c': 0.4}, 3: {'c': 0.9}, 4: {'c': 0.9}, 5: {'c': 0.9}, 6: {'c': 0.9}, 7: {'c': 0.9}, 8: {'c': 0.4}, 9: {'c': 0.4}, 10: {'c': 0.4}}
    >>> sequential_coiled(pos, psi)
    False
    >>> pos = (120, 141)
    >>> psi = {120: {'c': 0.078}, 121: {'c': 0.111}, 122: {'c': 0.087}, 123: {'c': 0.077}, 124: {'c': 0.067}, 125: {'c': 0.074}, 126: {'c': 0.042}, 127: {'c': 0.02}, 128: {'c': 0.029}, 129: {'c': 0.041}, 130: {'c': 0.068}, 131: {'c': 0.234}, 132: {'c': 0.5}, 133: {'c': 0.778}, 134: {'c': 0.848}, 135: {'c': 0.819}, 136: {'c': 0.792}, 137: {'c': 0.806}, 138: {'c': 0.615}, 139: {'c': 0.458}, 140: {'c': 0.419}, 141: {'c': 0.403}}
    >>> sequential_coiled(pos, psi)
    True
    """
    if verbose:
        print 'sequential:'
        print all(psi[i]['c'] > 0.4 for i in range(pos[0], pos[0]+3)), [psi[i]['c'] for i in range(pos[0], pos[0]+3)]
        print all(psi[i]['c'] > 0.4 for i in range(pos[1]-3, pos[1])), [psi[i]['c'] for i in range(pos[1]-3, pos[1])]
        print all(psi[i]['e'] > 0.4 for i in range(pos[1]-3, pos[1])), [psi[i]['e'] for i in range(pos[1]-3, pos[1])]
        print all(psi[i]['e'] > 0.4 for i in range(pos[0], pos[0]+3)), [psi[i]['e'] for i in range(pos[0], pos[0]+3)]
    for i in range(pos[0], pos[1]-10+1):
        if all(psi[j]['h'] < 0.4 for j in range(i, i+10+1)):
            return True
    return all(psi[i]['c'] > 0.4 for i in range(pos[0], pos[0]+3)) or \
           all(psi[i]['c'] > 0.4 for i in range(pos[1]-3, pos[1])) or \
           all(psi[i]['e'] > 0.4 for i in range(pos[0], pos[0]+3)) or \
           all(psi[i]['e'] > 0.4 for i in range(pos[1]-3, pos[1])) or \
           all(psi[i]['h'] < 0.4 for i in range(pos[0], pos[0]+3)) or \
           all(psi[i]['h'] < 0.4 for i in range(pos[1]-3, pos[1]))



def parse_WGP(text):
    """
    :param text:
    :return:
    >>> test = "{ [0   to 10   in fwd => 1 A            1 1 : 20   to 30   in rev =>  2 B              2 2 : 30  to 40  in fwd => 3 C             3 3] }~> total_grade 6 win_num  3 c_term fwd"
    >>> w1 = WinGrade(0, 10, 'fwd', 'A', grade=1, length_element=1, charges=1)
    >>> w2 = WinGrade(20, 30, 'rev', 'B', grade=2, length_element=2, charges=2)
    >>> w3 = WinGrade(30, 40, 'fwd', 'C', grade=3, length_element=3, charges=3)
    >>> parse_WGP(test).same_as_other(WinGradePath([w1, w2, w3]))
    True
    """
    import re
    wins = re.search('\[.*\]', text).group()
    spl = wins[1:-1].split(':')
    win_list = []
    for w in spl:
        ws = w.split()
        temp_win = WinGrade(int(ws[0]), int(ws[2]), ws[4], ws[7], grade=float(ws[6]), length_element=float(ws[9]), charges=int(ws[8]))
        win_list.append(temp_win)
    return WinGradePath(win_list)