#!/usr/bin/env python2.7
class TMConstraint():
    """
    a class to describe constraints on a protein's topology
    """
    def __init__(self, name, mode=None, tm_num=None, tm_pos=None, tm_pos_fidelity=0, c_term=None, n_term=None, non_tm_pos=None):
        self.name = name
        self.tm_num = tm_num
        if tm_pos is None:
            self.tm_pos = None
        else:
            self.tm_pos = sorted(tm_pos)
        self.c_term = c_term
        self.n_term = n_term
        self.tm_pos_fidelity = tm_pos_fidelity
        self.non_tm_pos = non_tm_pos
        self.mode = mode

    def __repr__(self):
        msg = ''
        msg += 'name %s\n' % self.name
        msg += 'tm_num %r\n' % self.tm_num
        if self.tm_pos is not None:
            for t in self.tm_pos:
                msg += 'tm_pos %i %i %s\n' % (t[0], t[1], t[2])
        else:
            msg += 'tm_pos None\n'
        msg += 'tm_pos_fidelity %r\n' % self.tm_pos_fidelity
        if self.non_tm_pos is not None:
            for t in self.non_tm_pos:
                msg += 'non_tm_pos %i %i\n' % (t[0], t[1])
        else:
            msg += 'non_tm_pos None\n'
        msg += 'c_term %r\n' % self.c_term
        msg += 'n_term %r\n' % self.n_term
        msg += 'mode %s\n' % self.mode
        return msg

    def test_tm_num(self, win_list):
        """
        :param win_list: a list of predicted windows
        :return: True if tm nums agree
        """
        return len([a for a in win_list if a.seq != '']) == self.tm_num

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

    def test_tm_pos(self, win_list, verbose=False):
        """
        :param win_list: a list of predicted windows
        :return: True iff all positions in cst have a win in win_list closer by both ends than fidelity
        """
        for cst_pos in self.tm_pos:
            if not win_list_near_pos_by_fidel(cst_pos, win_list, self.tm_pos_fidelity):
                if verbose:
                    print 'tm_pos fail at', cst_pos
                return False
        return True

    def test_non_tm_pos(self, win_list):
        """
        :param win_list:
        :return:
        """
        for win in win_list:
            if not win_not_in_segments(win, self.non_tm_pos):
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
            if not self.test_tm_pos(win_list, verbose):
                if verbose:
                    print 'failed tm_pos'
                return False
        if self.non_tm_pos is not None:
            if not self.test_non_tm_pos(win_list):
                if verbose:
                    print 'failed non_tm_pos'
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

    def pos_cst_segment(self, pos):
        """
        :param pos: a pos from tm_pos
        :return: semgent beginning and end (-/+ fidelity)
        """
        return (pos[0]-self.tm_pos_fidelity, pos[1]+self.tm_pos_fidelity)


def win_not_in_segments(win, segments):
    """
    :param win: a WinGrade instance
    :param segments: a list of segments to test
    :return: True iff the win does not overlap with ANY of the segments
    """
    for seg in segments:
        if (seg[0] <= win.begin <= seg[1]) or (seg[0] <= win.begin <= seg[1]):
            return False
    return True


def seg_within_tm_pos(seg, tm_pos, fidelity):
    """
    :param seg: a segment tuple, for a potential WinGrade (start, end)
    :param tm_pos: a pos list of type [(begin, end, direction)]
    :param fidelity:
    :return:True if there's no tm_pos, or if the segment is within any of the tm_poses
    >>> seg = (10, 20)
    >>> tm_pos = [(0, 5), (30, 40)]
    >>> fidelity = 5
    >>> seg_within_tm_pos(seg, tm_pos, fidelity)
    False
    >>> tm_pos = [(0, 5), (10, 20),(30, 40)]
    >>> seg_within_tm_pos(seg, tm_pos, fidelity)
    True
    >>> tm_pos = [(0, 5), (11, 19),(30, 40)]
    >>> seg_within_tm_pos(seg, tm_pos, fidelity)
    True
    >>> tm_pos = [(0, 5), (21, 29),(30, 40)]
    >>> seg_within_tm_pos(seg, tm_pos, fidelity)
    False
    """
    if tm_pos is None:
        return False
    for pos in tm_pos:
        if pos[0]-fidelity <= seg[0] <= pos[1]+fidelity and pos[0]-fidelity <= seg[1] <= pos[1]+fidelity:
            return True
    return False


def win_in_tm_pos(pos, tm_pos, fidelity):
    """
    :param pos: a pos for a potential window
    :param tm_pos: a list of tm_pos constraints
    :param fidelity: fidelity
    :return: False if none, or a tm_pos constraint that contains the pos
    >>> pos = [20, 30]
    >>> tm_pos = [[0, 10], [15, 25], [90, 100]]
    >>> fidelity = 5
    >>> win_in_tm_pos(pos, tm_pos, fidelity)
    [15, 25]
    >>> pos = [100, 110]
    >>> win_in_tm_pos(pos, tm_pos, fidelity)
    False
    >>> tm_pos = None
    >>> win_in_tm_pos(pos, tm_pos, fidelity)
    False
    """
    if tm_pos is None:
        return False, None
    for i, tm in enumerate(sorted(tm_pos)):
        if tm[0]-fidelity <= pos[0] and tm[1]+fidelity >= pos[1]:
            return tm, i
    return False, None


def pos_in_non_tm(pos, non_tm_pos):
    """
    :param pos: a tm pos [begin, end]
    :param non_tm_pos: list of non tm poses [[begin, end],...]
    :return: True iff pos is within any of the non_tm_pos
    >>> pos_in_non_tm([10, 20], [[0, 5]])
    False
    >>> pos_in_non_tm([10, 20], [[0, 15]])
    True
    >>> pos_in_non_tm([10, 20], [[9, 21]])
    True
    >>> pos_in_non_tm([100, 120], [[9, 21]])
    False
    >>> pos_in_non_tm([10, 20], [[12, 15]])
    True
    """
    if non_tm_pos is None or non_tm_pos is []:
        return False
    pos_range = range(pos[0], pos[1]+1)
    for non in non_tm_pos:
        non_range = range(non[0], non[1]+1)
        if any([a in non_range for a in pos_range]):
            return True
    return False


def all_satisfying_wins(pos, win_list, fidelity):
    """
    :param pos: constraint pos
    :param win_list: a list of WinGrades
    :return: all wins in win_list that satisfy pos
    """
    result = []
    for win in win_list:
        if win.direction == pos[2] or pos[2] is None:
            if win_near_pos_by_fidel(pos, win, fidelity):
                result.append(win)
    return result


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
    if in_path[-1] != '/':
        in_path += '/'
    with open(in_path+name+'.cst', 'r') as f:
        cont = f.read().split('\n')
    result = {'tm_pos': [], 'non_tm_pos': []}
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
        if sp[0] == 'non_tm_pos':
            if sp[1] != 'None':
                result['non_tm_pos'].append((int(sp[1]), int(sp[2])))
            else:
                result['non_tm_pos'] = None
            continue

        try:
            result[sp[0]] = int(sp[1])
        except:
            if sp[1] == 'None':
                result[sp[0]] = None
            else:
                result[sp[0]] = sp[1]
    return TMConstraint(result['name'], result['mode'], result['tm_num'], result['tm_pos'], result['tm_pos_fidelity'],
                        result['c_term'], result['n_term'], non_tm_pos=result['non_tm_pos'])


def pred2cst(name, path, ts, cst_mode, tm_pos_fidelity, write=True):
    """
    :param name: protein name
    :param path: path for writing .cst
    :param ts: topo string
    :return: writes a file with TMconstraint instance, and returns the instance
    """
    import re
    hhh = re.compile('[hH]*')
    tms = [(a.start(), a.end(), None) for a in hhh.finditer(ts) if a.end()-a.start() > 1]
    # msg = str(TMConstraint(name, tm_pos=tms, tm_pos_fidelity=tm_pos_fidelity, mode=cst_mode))
    tmc = TMConstraint(name, tm_pos=tms, tm_pos_fidelity=tm_pos_fidelity, mode=cst_mode)
    if write:
        # with open(path + name + '.cst', 'wr+') as fout:
        #     fout.write(msg)
        write_cst(path, name, tmc)
    return tmc
    # return TMConstraint(name, tm_pos=tms, tm_pos_fidelity=tm_pos_fidelity, mode=cst_mode)


def write_cst(path, name, tmc):
    with open(path + name + '.cst', 'wr+') as fout:
            fout.write(str(tmc))


def rost2cst_tm_num(args):
    import re
    from topo_strings_comparer import parse_rostlab_db, spc_parser
    topc = spc_parser(args['name'].lower())
    signle_peptide = topc['spoctopus'].count('S') + topc['spoctopus'].count('s')
    rostdb = parse_rostlab_db()[args['name'].lower()]
    pdbtm = rostdb['pdbtm']
    pdbtm_cln = 's'*signle_peptide + pdbtm[signle_peptide:]
    opm = rostdb['opm']
    opm_cln = 's'*signle_peptide + opm[signle_peptide:]
    hhh = re.compile('[hH]*')
    pdbtm_num = len([(a.start(), a.end(), None) for a in hhh.finditer(pdbtm_cln) if a.end()-a.start() > 1])
    opm_num =len([(a.start(), a.end(), None) for a in hhh.finditer(opm_cln) if a.end()-a.start() > 1])
    tmc = TMConstraint(args['name'].lower(), tm_num=max(pdbtm_num, opm_num))
    with open(args['path']+args['name'].lower()+'.cst', 'wr+') as fout:
        fout.write(str(tmc))
    return tmc


def rost2cst_tm_pos(args):
    import re
    from topo_strings_comparer import parse_rostlab_db, spc_parser
    topc = spc_parser(args['name'].lower())
    signle_peptide = topc['spoctopus'].count('S') + topc['spoctopus'].count('s')
    rostdb = parse_rostlab_db()[args['name'].lower()]
    pdbtm = rostdb['pdbtm']
    pdbtm_cln = 's'*signle_peptide + pdbtm[signle_peptide:]
    # opm = rostdb['opm']
    # opm_cln = 's'*signle_peptide + opm[signle_peptide:]
    hhh = re.compile('[hH]*')
    pdbtm = [[a.start(), a.end(), None] for a in hhh.finditer(pdbtm_cln) if a.end()-a.start() > 1]
    # opm = [(a.start(), a.end(), None) for a in hhh.finditer(opm_cln) if a.end()-a.start() > 1]
    for h in pdbtm:
        if h[1] < len(pdbtm_cln):
            if pdbtm_cln[h[1]] == '1':
                h[2] = 'rev'
                continue
            elif pdbtm_cln[h[1]] == '2':
                h[2] = 'fwd'
                continue
        if pdbtm_cln[h[0]-1] == '1':
            h[2] = 'fwd'
        else:
            h[2] = 'rev'

    # for i, h in enumerate(pdbtm[:-1]):
    #     if h[2] == pdbtm[i+1][2]:
    #         print 'HELP!!!!!!!!!', args['name']
    #         import sys
    #         sys.exit()

    tmc = TMConstraint(args['name'].lower(), tm_pos=pdbtm, mode='only', tm_pos_fidelity=7)
    with open(args['path']+args['name'].lower()+'.cst', 'wr+') as fout:
        fout.write(str(tmc))
    return tmc


def topcons2cst(args):
    import re
    from math import ceil, floor
    from topo_strings_comparer import spc_parser
    try:
        topc = spc_parser(args['name'].lower(), path=args['path'])
    except:
        topc = spc_parser(args['name'], path=args['path'])
    single_peptide = topc['spoctopus'].count('S') + topc['spoctopus'].count('s')
    topcons = topc['topcons']
    topcons_cln = 's'*single_peptide + topcons[single_peptide:]
    hhh = re.compile('[Mm]*')
    topc_hhh = [(a.start(), a.end()-1, None) for a in hhh.finditer(topcons_cln) if a.end()-a.start() > 1]
    final_topc_hhh = []
    for seg in topc_hhh:
        diff = 20-(seg[1]-seg[0])
        if diff > 0:
            if seg[0] <= single_peptide:
                continue
            s = [1, 1, 1]
            s[0] = max(0, seg[0]-ceil(diff/2.))
            s[1] = min(len(topcons), seg[1]+ceil(diff/2.))
            s[2] = seg[2]
            final_topc_hhh.append(tuple(s))
        else:
            final_topc_hhh.append(seg)
    tmc = TMConstraint(args['name'].lower(), tm_num=None, tm_pos=final_topc_hhh, tm_pos_fidelity=args['tm_pos_fidelity']
                       , mode='only')
    with open(args['path']+args['name'].lower()+'.cst', 'wr+') as fout:
        fout.write(str(tmc))
    return tmc


if __name__ == '__main__':
    import argparse
    import os
    from topo_strings_comparer import prd_parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', type=str)
    parser.add_argument('-mode', default='pred2cst', type=str)
    parser.add_argument('-path', default=os.getcwd()+'/')
    parser.add_argument('-cst_mode', default=None)
    parser.add_argument('-tm_pos_fidelity', type=int, default=3)
    args = vars(parser.parse_args())
    if args['mode'] == 'pred2cst':
        prd = prd_parser(args['path'], args['name'].lower() + '.prd')
        print pred2cst(args['name'].lower(), args['path'], prd['pred_ts'], args['cst_mode'], args['tm_pos_fidelity'])

    elif args['mode'] == 'cst2TMC':
        tmc = parse_cst(args['name'].lower(), args['path'])
        print tmc

    elif args['mode'] == 'rost2cst':
        rost2cst(args)
        print rost2cst(args)

    elif args['mode'] == 'topcons2cst':
        tmc = topcons2cst(args)
        print tmc

    elif args['mode'] == 'rost2cst_tm_pos':
        rost2cst_tm_pos(args)

    else:
        print 'no mode recived'