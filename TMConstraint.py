class TMConstraint():
    """
    a class to describe constraints on a protein's topology
    """
    def __init__(self, name, tm_num=None, tm_pos=None, tm_pos_fidelity=None, c_term=None, n_term=None):
        self.name = name
        self.tm_num = tm_num
        self.tm_pos = tm_pos
        self.c_term = c_term
        self.n_term = n_term
        self.tm_pos_fidelity = tm_pos_fidelity

    def __repr__(self):
        msg = ''
        msg += 'name %s\n' % self.name
        msg += 'tm_num %r\n' % self.tm_num
        for t in self.tm_pos:
            msg += 'tm_pos %i %i %s\n' % (t[0], t[1], t[2])
        msg += 'tm_pos_fidelity %r\n' % self.tm_pos_fidelity
        msg += 'c_term %r\n' % self.c_term
        msg += 'n_term %r\n' % self.n_term
        return msg

    def test_tm_num(self, win_list):
        """
        :param win_list: a list of predicted windows
        :return: True if tm nums agree
        """
        return len(win_list) == self.tm_num

    def test_c_term(self, win_list):
        '''
        :param win_list: a list of predicted windows
        :return: True if c terminus agrees with constraints
        '''
        return (win_list[-1].direction == 'fwd' and self.c_term == 'out') or \
               (win_list[-1].direction == 'rev' and self.c_term == 'in')

    def test_n_term(self, win_list):
        """
        :param win_list: a list of predicted windows
        :return: True if c terminus agrees with constraints
        """
        return (win_list[0].direction == 'fwd' and self.n_term == 'in') or \
               (win_list[0].direction == 'rev' and self.n_term == 'out')

    def test_tm_pos(self, win_list):
        """
        :param win_list: a list of predicted windows
        :return: True iff all positions in cst have a win in win_list closer by both ends than fidelity
        """
        for cst_pos in self.tm_pos:
            if not win_list_near_pos_by_fidel(cst_pos, win_list, self.tm_pos_fidelity):
                return False
        return True

    def test_manager(self, win_list, verbose=False):
        """
        :param win_list: a list of predicted windows
        :return: True iff win_list passes all non None constraints
        """
        if self.tm_num is not None:
            if not self.test_tm_num(win_list):
                if verbose:
                    print 'failed tm_num'
                return False
        if self.c_term is not None:
            if not self.test_c_term(win_list):
                if verbose:
                    print 'failed c_term'
                return False
        if self.n_term is not None:
            if not self.test_n_term(win_list):
                if verbose:
                    print 'failed n_term'
                return False
        if self.tm_pos is not None:
            if not self.test_tm_pos(win_list):
                if verbose:
                    print 'failed tm_pos'
                return False
        return True

    def test_manager_segment(self, win_list, segment, verbose=False):
        """
        :param win_list: a list of win grades
        :param segment: a segment on the sequence to be tested for non-satisfied constraints
        :param verbose: wheher to print or not
        :return: True iff all poses that are within segment are satisfied by win_list
        """
        if self.tm_pos is None:
            return True
        for pos in self.tm_pos:
            print 'checking pos', pos
            print 'with win_list', win_list
            if pos[0] >= segment[0] and pos[1] <= segment[1]:
                print 'is within segment', segment
                if not win_list_near_pos_by_fidel(pos, win_list, self.tm_pos_fidelity):
                    if verbose:
                        print 'failed at pos', pos
                    return False
        return True

    def unsat_cst(self, win_list):
        """
        :param win_list: a list of win grades
        :return: a list of WinGrades that are the most negative and satisfy tm_poses that cannot be satisfied by a
                    negative win in win_list
        """
        if self.tm_pos is None:
            return []
        unsat = self.tm_pos[:]
        for win in win_list:
            if win.grade < 0:
                for pos in unsat:
                    if win_near_pos_by_fidel(pos, win, self.tm_pos_fidelity):
                        unsat.remove(pos)
        print "unsat_cst found ", unsat
        sat_wins = []
        for pos in unsat:
            sat_wins.append(self.satisfying_win(pos, win_list))
        return sat_wins

    def satisfying_win(self, pos, win_list):
        """
        :param pos: constraint pos
        :param win_list: a list of WinGrades
        :return: the most negative fwd and rev windows that satisfy pos
        """
        best_win_fwd, best_win_rev = win_list[0], win_list[0]
        for win in win_list:
            if win_near_pos_by_fidel(pos, win, self.tm_pos_fidelity) and win.direction == 'fwd' and \
                            win.grade < best_win_fwd.grade:
                best_win_fwd = win
            elif win_near_pos_by_fidel(pos, win, self.tm_pos_fidelity) and win.direction == 'rev' and \
                            win.grade < best_win_rev.grade:
                best_win_rev = win
        if pos[2] == 'fwd':
            return [best_win_fwd]
        elif pos[2] == 'rev':
            return [best_win_rev]
        else:
            return [best_win_fwd, best_win_rev]


def win_list_near_pos_by_fidel(pos, win_list, fidel):
    """
    :param pos: a cst position
    :param win_list: a list of predicted windows
    :param fidel: fidelity
    :return: True if any win in win_list agrees with win_near_pos_by_fidel
    """
    for win in win_list:
        if win.direction == pos[2] or pos[2] is None:
            if win_near_pos_by_fidel(pos, win, fidel):
                return True
    return False


def win_near_pos_by_fidel(pos, win, fidel):
    """
    :param pos: a cst position
    :param win: a WinGrade
    :param fidel: fidelity
    :return: True if both WinGrade ends within pos+fidel
    """
    return pos[0]-fidel <= win.begin and pos[1]+fidel >= win.end


def parse_cst(name, in_path):
    """
    :param name: protein name
    :param in_path: path to .cst
    :return: a TMConstraint object with only position constraints
    """
    with open(in_path+'/'+name+'.cst', 'r') as f:
        cont = f.read().split('\n')
    result = {'tm_pos': []}
    for ln in cont:
        sp = ln.split()
        if len(sp) <= 0:
            continue
        if sp[0] == 'tm_pos':
            if sp[1] != 'None':
                result['tm_pos'].append((int(sp[1]), int(sp[2]), sp[3] if sp[3] != 'None' else None))
            else:
                result['tm_pos'] = None
            continue
        try:
            result[sp[0]] = int(sp[1])
        except:
            if sp[1] == 'None':
                result[sp[0]] = None
            else:
                result[sp[0]] = sp[1]
    return TMConstraint(result['name'], result['tm_num'], result['tm_pos'], result['tm_pos_fidelity'], result['c_term'],
                        result['n_term'])


def pred2cst(name, path, ts):
    """
    :param name: protein name
    :param path: path for writing .cst
    :param ts: topo string
    :return: writes a file with TMconstraint instance, and returns the instance
    """
    import re
    hhh = re.compile('[hH]*')
    tms = [(a.start(), a.end(), None) for a in hhh.finditer(ts) if a.end()-a.start() > 1 ]
    msg = str(TMConstraint(name, tm_pos=tms, tm_pos_fidelity=5))
    with open(path + '/' + name + '.cst', 'wr+') as fout:
        fout.write(msg)
    return TMConstraint(name, tm_pos=tms, tm_pos_fidelity=5)


if __name__ == '__main__':
    import argparse
    import os
    from topo_strings_comparer import prd_parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', type=str)
    parser.add_argument('-mode', default='pred2cst', type=str)
    parser.add_argument('-path', default=os.getcwd())
    args = vars(parser.parse_args())
    if args['mode'] == 'pred2cst':
        prd = prd_parser(args['path'], args['name'].lower() + '.prd')
        print pred2cst(args['name'].lower(), args['path'], prd['pred_ts'])
    elif args['mode'] == 'cst2TMC':
        tmc = parse_cst(args['name'].lower(), args['path'])
        print tmc