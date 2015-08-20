#!/usr/bin/env python
class HphobicityScore():
    def __init__(self, name, seq, ss2_file, hydro_polyval, param_list, csts):
        '''
        :param name: protein's name
        :param uniprot: the entrie's uniprot code
        :param seq: protein's sequence
        :param hydro_polyval: polynom valued dictionary
        :return: a stack of WinGrade instances, and their utilities
        '''
        global HP_THRESHOLD, INC_MAX, MIN_LENGTH, PSI_HELIX, PSI_RES_NUM, known_tm_num, poly_param, output_path
        HP_THRESHOLD = param_list['hp_threshold']
        INC_MAX = 10  # 8
        MIN_LENGTH = param_list['min_length']
        PSI_HELIX = param_list['psi_helix']
        PSI_RES_NUM = param_list['psi_res_num']
        known_tm_num = param_list['known_tm_num']
        poly_param = {k: v for k, v in param_list.items() if k in ['c0', 'c1', 'c2', 'c3', 'z_0', 'w']}
        print 'poly', poly_param
        output_path = param_list['result_path']
        self.csts = csts
        self.percentile = param_list['msa_percentile']
        self.name = name
        self.seq = seq
        self.ss2_file = ss2_file
        self.seq_length = len(seq)
        # self.psipred = self.PsiReaderHelix()
        self.psipred = self.parse_psipred()
        self.polyval = hydro_polyval
        self.with_msa = param_list['with_msa']
        # self.topo = self.topo_greedy_chooser()
        # self.n_term_orient = self.topo[0].direction
        self.after_SP = 0 if self.seq[0] != 'u' else max([i for i, aa in enumerate(self.seq) if aa == 'u']) + 1
        self.WinGrades = self.win_grade_generator(self.after_SP, self.seq_length, 'both')
        # self.topo_string = self.make_topo_string()
        # print 'before sort', self.WinGrades
        # self.sorted_grade = self.sort_WinGrades()
        # print 'after sort', self.sorted_grade
        # self.minimas = self.local_minima_finder(direction='both')
        # self.fwd_minimas = self.local_minima_finder(direction='fwd')
        # self.rev_minimas = self.local_minima_finder(direction='rev')
        # self.topo = self.topo_determine_new()
        # self.topo_minimas = self.topo_determine()
        # self.sorted_grade_norm = self.sort_WinGrades_norm()
        # self.minimas_norm = self.local_minima_finder_norm()
        # self.topo = self.topo_brute()
        self.topo_best, self.topo_best_val, self.topo_sec_best, self.topo_sec_best_val = self.topo_graph
        try:
            self.best_c_term = 'out' if self.topo_best[-1].direction == 'fwd' else 'in'
        except:
            self.best_c_term = None
        try:
            self.sec_best_c_term = 'out' if self.topo_sec_best[-1].direction == 'fwd' else 'in'
        except:
            self.sec_best_c_term = None

    def __str__(self):
        """
        :return: A message with all Fwd/Rev minimas, and the selected topology
        """
        return 'Selectod topology:\n' + '\n'.join([str(i) for i in self.best_topo])

    # def win_grade_generator(self, pos1, pos2, fwd_or_rev):
    #     '''
    #     :return:grades all segments of self.seq, and aggregates them as WinGrades
    #     '''
    #     from WinGrade import WinGrade
    #     PSI_HP_THRESHOLD = -8.
    #     psi = self.psipred
    #     grades = []
    #     for i in range(pos1, pos2):
    #         for inc in range(min(INC_MAX, self.seq_length - 20 - i)):
    #             is_not_helical = self.is_not_helical((i, i+20+inc), psi)
    #             fwd_temp = WinGrade(i, i+20+inc, 'fwd', self.seq[i:i+20+inc], self.polyval)
    #             rev_temp = WinGrade(i, i+20+inc, 'rev', self.seq[i:i+20+inc][::-1], self.polyval)
    #             if (fwd_or_rev == 'both' or fwd_or_rev == 'fwd') and (not is_not_helical or fwd_temp.grade < PSI_HP_THRESHOLD):
    #                 grades.append(WinGrade(i, i+20+inc, 'fwd', self.seq[i:i+20+inc], self.polyval))
    #             if (fwd_or_rev == 'both' or fwd_or_rev == 'rev') and (not is_not_helical or rev_temp.grade < PSI_HP_THRESHOLD):
    #                 grades.append(WinGrade(i, i+20+inc, 'rev', self.seq[i:i+20+inc][::-1], self.polyval))
    #     return grades

    def win_grade_generator(self, pos1, pos2, fwd_or_rev):
        '''
        :return:grades all segments of self.seq, and aggregates them as WinGrades
        '''
        from WinGrade import WinGrade, count_charges
        import sys
        if self.with_msa:
            from MSA_for_TMpredict import TMpredict_MSA
            msa_obj = TMpredict_MSA(self.name, self.polyval, poly_param, self.percentile)
        psi = self.psipred
        grades = []
        # print "psi", psi
        # setup toolbar
        bar_width = len(range(pos1, pos2))
        sys.stdout.write("[%s]" % (" " * bar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (bar_width+1))  # return to start of line, after '['
        for i in range(pos1, pos2+1):
            for inc in range(min(INC_MAX, self.seq_length - MIN_LENGTH - i)):
                if not self.is_not_helical((i, i+MIN_LENGTH+inc), psi) and \
                                count_charges(self.seq[i:i+MIN_LENGTH+inc]) < 3:
                    if self.with_msa:
                        grades.append(msa_obj.retrieve_seqs(i, i+MIN_LENGTH+inc, 'fwd'))
                        grades.append(msa_obj.retrieve_seqs(i, i+MIN_LENGTH+inc, 'rev'))
                    else:
                        # print 'creating win %s' % self.seq[i:i+MIN_LENGTH+inc]
                        grades.append(WinGrade(i, i+MIN_LENGTH+inc, 'fwd', self.seq[i:i+MIN_LENGTH+inc], self.polyval,
                                            poly_param))
                        grades.append(WinGrade(i, i+MIN_LENGTH+inc, 'rev', self.seq[i:i+MIN_LENGTH+inc][::-1], self.polyval,
                                            poly_param))
            # writes to the progress bar
            sys.stdout.write("-")
            sys.stdout.flush()
        sys.stdout.write("]\n")
        return grades

    def win_grade_generator_old(self, pos1, pos2, fwd_or_rev):
        '''
        :return:grades all segments of self.seq, and aggregates them as WinGrades
        '''
        from WinGrade import WinGrade
        import sys
        if self.with_msa:
            from MSA_for_TMpredict import TMpredict_MSA
            msa_obj = TMpredict_MSA(self.name, self.polyval, poly_param)
        PSI_HP_THRESHOLD = -8.
        psi = self.psipred
        grades = []

        # setup toolbar
        bar_width = len(range(pos1, pos2))
        sys.stdout.write("[%s]" % (" " * bar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (bar_width+1))  # return to start of line, after '['

        for i in range(pos1, pos2):
            for inc in range(min(INC_MAX, self.seq_length - MIN_LENGTH - i)):
                # print i, inc
                is_not_helical = self.is_not_helical((i, i+MIN_LENGTH+inc), psi)
                if self.with_msa:
                    fwd_temp = msa_obj.retrieve_seqs(i, i+MIN_LENGTH+inc, 'fwd')
                    rev_temp = msa_obj.retrieve_seqs(i, i+MIN_LENGTH+inc, 'rev')
                else:
                    fwd_temp = WinGrade(i, i+MIN_LENGTH+inc, 'fwd', self.seq[i:i+MIN_LENGTH+inc], self.polyval,
                                        poly_param)
                    rev_temp = WinGrade(i, i+MIN_LENGTH+inc, 'rev', self.seq[i:i+MIN_LENGTH+inc][::-1], self.polyval,
                                        poly_param)
                if (fwd_or_rev == 'both' or fwd_or_rev == 'fwd') and \
                        (not is_not_helical or fwd_temp.grade < PSI_HP_THRESHOLD) and fwd_temp.hp_moment <= 6:
                    grades.append(fwd_temp)
                if (fwd_or_rev == 'both' or fwd_or_rev == 'rev') and \
                        (not is_not_helical or rev_temp.grade < PSI_HP_THRESHOLD) and rev_temp <= 6:
                    grades.append(rev_temp)
            # writes to the progress bar
            sys.stdout.write("-")
            sys.stdout.flush()
        sys.stdout.write("\n")
        return grades    

    def print_HphobicityScore(self):
        '''
        :return:a utility funciton to print things from the class
        '''
        print 'printing ', self.name
        for i in self.WinGrades:
            i.print_WinGrade()

    def sort_WinGrades(self):
        '''
        :return:a list of the WinGrades, sorted by grade
        '''
        import operator
        dict_win_grades = {i.grade: i for i in self.WinGrades}
        sort_dict = sorted(dict_win_grades.items(), key=operator.itemgetter(0))
        return [{i[0]: i[1]} for i in sort_dict]

    def local_minima_finder(self, direction):
        """
        :direction: the direction requested to m=find minimas. rev/fwd/both
        :return:returns a list of the minimal grade non-clashing WinGrades
        """
        chosen_grades = []
        i = 0
        while len(chosen_grades) == 0:
            if self.sorted_grade[i].values()[0].direction == direction or direction == 'both':
                chosen_grades = [self.sorted_grade.pop(i).values()[0]]
            i += 1
        for grade in self.sort_WinGrades():
            if direction == 'both' or (direction == grade.values()[0].direction):
                if not grade.values()[0].set_grade_colliding(chosen_grades) \
                        and grade.values()[0].grade <= HP_THRESHOLD:
                    chosen_grades.append(grade.values()[0])
        return chosen_grades

    def sort_WinGrades_norm(self):
        import operator
        dict_win_grades = {i.grade_norm: i for i in self.WinGrades}
        sort_dict = sorted(dict_win_grades.items(), key=operator.itemgetter(0))
        return [{i[0]: i[1]} for i in sort_dict]

    def local_minima_finder_norm(self):
        chosen_grades = [self.sort_WinGrades_norm().pop(0).values()[0]]
        for grade in self.sort_WinGrades_norm():
            if not grade.values()[0].set_grade_colliding(chosen_grades):
                chosen_grades.append(grade.values()[0])
        return chosen_grades

    def plot_win_grades(self):
        import matplotlib.pyplot as plt
        import matplotlib.lines as mlines
        plt.figure()
        for win_grade in self.WinGrades:
            plt.plot((win_grade.begin, win_grade.end), (win_grade.grade, win_grade.grade),
                     'k--' if win_grade.direction == 'fwd' else 'grey', )
        print 'topo_best'
        for minima in self.topo_best:
            print minima
            plt.plot((minima.begin, minima.end), (minima.grade, minima.grade),
                     'b--' if minima.direction == 'fwd' else 'r--', lw=4)
        print 'topo_sec_best'
        for minima in self.topo_sec_best:
            print minima
            # if minima.grade > 0: continue
            # plt.plot((minima.begin, minima.end), (minima.grade, minima.grade), 'b--', lw=4)
            plt.plot((minima.begin, minima.end), (minima.grade, minima.grade),
                     'c--' if minima.direction == 'fwd' else 'm--', lw=4)
        # for minima in self.rev_minimas:
        #     if minima.grade > 0: continue
        #     plt.plot((minima.begin, minima.end), (minima.grade, minima.grade), 'r--')
        black_line = mlines.Line2D([], [], 'k--', marker='', lw=2, label='Fwd grade')
        grey_line = mlines.Line2D([], [], color='grey', marker='', lw=2, label='Rev grade')
        blue_line = mlines.Line2D([], [], color='blue', marker='', lw=2, label='Fwd minima')
        red_line = mlines.Line2D([], [], color='red', marker='', lw=2, label='Rev minima')
        green_line = mlines.Line2D([], [], color='green', marker='', lw=2, label='Fwd topo minima')
        purple_line = mlines.Line2D([], [], color='purple', marker='', lw=2, label='Rev topo minima')
        # plt.legend(handles=[black_line, grey_line, blue_line, red_line, green_line, purple_line], ncol=2)
        plt.xlabel('Sequence Position')
        plt.ylabel('Energy')
        plt.title('Win Grades Energy Plot for %s' % self.name)
        # plt.show()
        print 'saving fig to', output_path + self.name + '.pdf'
        plt.savefig(output_path + '/' + self.name + '.pdf')

    def topo_determine(self):
        """
        :return:a sorted list of minimas of flipping direction which has the lowest
        global energy
        """
        fwd_sorted = sorted(self.fwd_minimas, key=lambda x: x.begin)
        rev_sorted = sorted(self.rev_minimas, key=lambda x: x.begin)
        fwd_start = []
        rev_start = []
        for i in range(max([len(fwd_sorted), len(rev_sorted)])):
            fwd_start.append(fwd_sorted[i] if i % 2 == 0 and i < len(fwd_sorted) else None)
            fwd_start.append(rev_sorted[i] if i % 2 != 0 and i < len(rev_sorted) else None)
            rev_start.append(rev_sorted[i] if i % 2 == 0 and i < len(rev_sorted) else None)
            rev_start.append(fwd_sorted[i] if i % 2 != 0 and i < len(fwd_sorted) else None)
        fwd_start = [a for a in fwd_start if a is not None]
        rev_start = [a for a in rev_start if a is not None]
        fwd_score = sum([a.grade for a in fwd_start])
        rev_score = sum([a.grade for a in rev_start])
        print fwd_start
        print rev_start
        return fwd_start if fwd_score < rev_score and len(fwd_start) >= len(rev_start)\
            else rev_start

    def topo_determine_new(self):
            """
            :return:a sorted list of minimas of flipping direction which has the lowest
            global energy
            """
            fwd_sorted = sorted(self.fwd_minimas, key=lambda x: x.begin)
            rev_sorted = sorted(self.rev_minimas, key=lambda x: x.begin)
            fwd_start = []
            rev_start = []

            # for i in range(max([len(fwd_sorted), len(rev_sorted)])):
            #     fwd_start.append(fwd_sorted[i] if i % 2 == 0 and i < len(fwd_sorted) else None)
            #     fwd_start.append(rev_sorted[i] if i % 2 != 0 and i < len(rev_sorted) else None)
            #     rev_start.append(rev_sorted[i] if i % 2 == 0 and i < len(rev_sorted) else None)
            #     rev_start.append(fwd_sorted[i] if i % 2 != 0 and i < len(fwd_sorted) else None)
            fwd_start = [a for a in fwd_start if a is not None]
            rev_start = [a for a in rev_start if a is not None]
            fwd_score = sum([a.grade for a in fwd_start])
            rev_score = sum([a.grade for a in rev_start])
            print fwd_start
            print rev_start
            return fwd_start if fwd_score < rev_score and len(fwd_start) >= len(rev_start)\
                else rev_start

    def is_not_helical_old(self, pos, psi):
        """
        :param pos:start and end positions of the corresponding window
        :param psi: a list of alpha-helical propensity grades for the entire sequence
        :return:True if there are more than x AAs that are less helicla than y
        """
        non_helical = 0
        for i in range(pos[0], pos[1]):
            if psi[i] <= PSI_HELIX:
                non_helical += 1
        return False if non_helical < PSI_RES_NUM else True

    def is_not_helical(self, pos, psi):
        """
        :param pos:start and end positions of the corresponding window
        :param psi: a dict of c/e/h propensities
        :return:True if the average e or c propensities are above certain thresholds determined in
        psipred_vs_mm_nomm.py
        """
        import numpy as np
        assert all([psi[i]['aa'] == self.seq[i] for i in range(pos[0], pos[1])])
        return True if (np.mean([psi[a]['e'] for a in range(pos[0], pos[1])]) >= 0.3 or
                        np.mean([psi[a]['c'] for a in range(pos[0], pos[1])]) >= 0.48 or
                        np.mean([psi[a]['h'] for a in range(pos[0], pos[1])]) <= 0.3) else False

    def PsiReaderHelix(self):
        """
        :return:reads the sequence's psipred results. returns them as a string
        """
        ss2_file = open(self.ss2_file)
        result = []
        for line in ss2_file:
            split = line.split()
            if len(split) > 5:
                result.append(float(split[4]))
        ss2_file.close()
        return result

    def parse_psipred(self):
        results = {}
        with open(self.ss2_file, 'r') as f:
            cont = f.read().split('\n')
        for line in cont:
            split = line.split()
            if len(split) != 6: continue
            results[int(split[0])-1] = {'aa': split[1], 'c': float(split[3]), 'h': float(split[4]), 'e': float(split[5])}
        return results

    def PsiReaderHelix_old(self):
        """
        :return:reads the sequence's psipred results. returns them as a string
        """
        import re
        ss2_file = open(self.ss2_file)
        line_re = re.compile('^\s*([0-9]*)\s*([A-Z]*)\s*([A-Z]*)\s*(0\.[0-9]*)\s*(0\.[0-9]*)\s*(0\.[0-9]*)')
        result = []
        for line in ss2_file:
            if line_re.search(line):
                result.append(float(line_re.search(line).group(5)))
        ss2_file.close()
        return result

    def make_topo_string(self):
        in_or_out = 'out' if self.n_term_orient == 'rev' else 'in'
        tm_count = 0
        result = 'x'
        for i in range(len(self.seq)):
            if tm_count == len(self.topo):
                result += 'I' if in_or_out == 'in' else 'O'
            elif not self.topo[tm_count].pos_in_wingrade(i) and result[-1] != 'M': # not TM, and not emmidiately after TM
                result += 'I' if in_or_out == 'in' else 'O'
            elif self.topo[tm_count].pos_in_wingrade(i) and result[-1] != 'M': # first in TM
                result += 'M'
                in_or_out = 'in' if in_or_out == 'out' else 'out'
            elif self.topo[tm_count].pos_in_wingrade(i) and result[-1] == 'M': # in TM, not first
                result += 'M'
            elif not self.topo[tm_count].pos_in_wingrade(i) and result[-1] == 'M': # not in TM, emmidiatley after TM
                tm_count += 1
                result += 'I' if in_or_out == 'in' else 'O'
        return result[1:]

    def topo_greedy(self, direction):
        """
        :param direction:the direction for hte first (N') TM
        :return: a list of WinGrades starting in the desired direction, describing the minimas
        """
        big_win_len = 40
        result = []
        fwd_or_rev = direction
        i = 0
        last_end = -1
        while i < range(len(self.seq)-20):
            temp_win = self.win_grade_generator(i, i+20, fwd_or_rev)
            if len(temp_win) == 0:
                i += 1
                continue
            if temp_win[0].grade > -4.0:
                i += 1
            else:
                if i <= 10:
                    t = 0
                elif 10 < i < last_end + 1:
                    t = last_end+1
                else:
                    t = i
                win_grade_set = self.win_grade_generator(t, t+big_win_len, fwd_or_rev)
                min_grd = win_grade_set[0]
                for grd in win_grade_set:
                    if grd.grade < min_grd.grade:
                        min_grd = grd
                result.append(min_grd)
                fwd_or_rev = 'rev' if fwd_or_rev == 'fwd' else 'fwd'
                i = win_grade_set[-1].end + 1
                last_end = win_grade_set[-1].end
            if i+big_win_len > self.seq_length:
                break
        return result

    def topo_greedy_chooser(self):
        fwd = self.topo_greedy('fwd')
        rev = self.topo_greedy('rev')
        fwd_tot = sum([i.grade for i in fwd])
        rev_tot = sum([i.grade for i in rev])
        print 'fwd', fwd_tot, len(fwd)
        print 'rev', rev_tot, len(rev)
        if len(fwd) != len(rev):
            return fwd if len(fwd) > len(rev) else rev
        else:
            return fwd if fwd_tot < rev_tot else rev

    def topo_brute(self):
        from random import shuffle
        print "in brute"
        chosen_set = []
        chosen_grade = 1000
        for i in range(1000):
            all_wins = self.WinGrades[:]
            shuffle(all_wins)
            temp_set = [all_wins.pop()]
            while len(all_wins) != 0:
                temp_win = all_wins.pop()
                if not temp_win.set_grade_colliding(temp_set): temp_set.append(temp_win)
            if sum(a.grade for a in temp_set) < chosen_grade:
                chosen_set = temp_set[:]
                chosen_grade = sum(a.grade for a in temp_set)
        print chosen_set
        print chosen_grade
        return chosen_set

    def find_graph_path(self, pred, path, source_node):
        # find all wins in the minimal energy path from source to last win
        # print 'source_node', source_node
        while path[-1].seq != source_node.seq:
            for k, v in pred.items():
                if v is None:
                    continue
                if k.seq == path[-1].seq:
                    path.append(v)
        path.pop(-1)    # get rid of source_node
        return path[::-1]   # revert the path to begin->end

    def find_path_frompaths(self, pred, last_path, source_node_path):
        comp_path = [last_path]
        # print 'in find fro path'
        # print 'pred', pred
        # print 'last_path', last_path
        # print 'source node path', source_node_path
        # print 'last_path.first().seq', last_path.first().seq
        # print 'source_node_path.first().seq', source_node_path.first().seq
        # print 'comp_path[-1].first().seq', comp_path[-1].first().seq
        # print 'source_node_path.first().seq', source_node_path.first().seq
        while comp_path[-1].first().seq != source_node_path.first().seq:
            for k, v in pred.items():
                if v is None:
                    continue
                if k.first().seq == comp_path[-1].first().seq:
                    # print 'app', v
                    comp_path.append(v)
            # print 'noow its', comp_path, len(comp_path)
        # print comp_path[::-1], type(comp_path)
        print 'AAAA'
        print 'BBB', comp_path
        # comp_path.pop(-1)
        comp_path.remove(source_node_path)
        return [item for sublist in comp_path[::-1] for item in sublist.path][1:]

    @property
    def topo_graph(self):
        """
        :return: a topology
        """
        import networkx as nx
        from WinGrade import WinGrade, WinGradePath
        import operator
        # print self.WinGrades
        if self.csts.mode != 'only':
            win_list = [a for a in self.WinGrades if a.grade < HP_THRESHOLD and a.charges < 3]  # make copy of negative wingrades list
        else:
            win_list = [a for a in self.WinGrades if a.grade < HP_THRESHOLD and a.charges < 3
                        and a.within_segments(self.csts.tm_pos, self.csts.tm_pos_fidelity)]
        print 'constructing graph'
        G = nx.DiGraph()
        source_node = WinGrade(0, 0, 'fwd', '', self.polyval, poly_param)   # define source win
        if self.with_msa:
            if self.csts.tm_pos is None:
                [G.add_edge(source_node, a, weight=a.msa_grade) for a in win_list]  # add all win to source edges with msa
            else:
                [G.add_edge(source_node, a, weight=a.msa_grade) for a in win_list if a.end <= self.csts.tm_pos[0][1]]  # add all win to source edges with msa
        else:
            if self.csts.tm_pos is None:
                [G.add_edge(source_node, a, weight=a.grade) for a in win_list]  # add all win to source edges
            else:
                [G.add_edge(source_node, a, weight=a.grade) for a in win_list if a.end <= self.csts.tm_pos[0][1]]  # add all win to source edges
        for win1 in win_list:   # add all win to win edges where condition applies
            cst_after = self.cst_after(win1)
            for win2 in win_list:
                if cst_after is None or win2.begin < cst_after[1]: #  make sure no edges skip a cst
                    # condition: non-overlapping, 2 is after 1, opposite directions
                    if not win1.grade_grade_colliding(win2) and win2.begin > win1.end and win1.direction != win2.direction \
                            and win2.begin-win1.end >= 2:
                        if not self.with_msa:
                            G.add_edge(win1, win2, weight=win2.grade)
                        elif self.with_msa:
                            G.add_edge(win1, win2, weight=win2.msa_grade)
        if self.csts.tm_num is None:
            print "About to Bellman-Ford"
            pred, dist = nx.bellman_ford(G, source_node)
            print "Finished Bellman-Fording"
            sorted_dist = sorted(dist.items(), key=operator.itemgetter(1))
            for k, v in sorted_dist:
                print k, v
            temp_direction = 'A'
            best_path, sec_best_path = [], []
            best_score, sec_best_score = None, None
            for last_path, total_grade in sorted_dist:
                temp_path = self.find_graph_path(pred, [last_path], source_node)
                if temp_path == [] or temp_path == [source_node]:
                    continue
                print 'tp', temp_path, total_grade
                if temp_direction != 'A' and self.csts.test_manager(temp_path) \
                        and temp_path[-1].direction is not temp_direction:
                    # print 'looking at sec', temp_path[-1].direction, temp_direction
                    sec_best_path = temp_path[:]
                    if sec_best_path[0].seq == 'SOURCE':
                        sec_best_path.pop(0)
                    sec_best_score = sum([a.grade for a in sec_best_path])
                    break
                elif self.csts.test_manager(temp_path) and temp_path[-1].direction != temp_direction:
                    # print 'looking at BEST', temp_path
                    best_path = temp_path[:]
                    best_score = sum([a.grade for a in best_path])
                    temp_direction = best_path[-1].direction
                    # print 'now temp_direction is:', temp_direction
            # print 'found best'
            print best_path
            print best_score
            # print 'found sec best',
            # print sec_best_path[-1].direction
            # print sec_best_score
        else:
            """
            a method to find minimal path with only tm_num nodes. a BFS algorithm that stops every path building once
            it gets to tm_num nodes. only saves the best solutions for c_term fwd & rev, and only if they pass all
            constraints
            """
            from collections import deque
            from copy import deepcopy
            import timeit
            tic = timeit.default_timer()
            print 'looking for %i wins' % self.csts.tm_num
            best_fwd = WinGradePath([])
            best_rev = WinGradePath([])
            Q = deque([WinGradePath([source_node])])
            while Q:  # continue as long as there si something in Q
                path = Q.popleft()  # get the FIFO out
                # print len(Q), path.win_num, self.csts.tm_num, path
                if path.win_num == self.csts.tm_num:  # the path has enough nodes
                    if self.csts.test_manager(path.path, True):
                        if path.c_term == 'fwd' and path.total_grade < best_fwd.total_grade:
                            best_fwd = deepcopy(path)
                            print 'found besf fwd', best_fwd
                        elif path.c_term == 'rev' and path.total_grade < best_rev.total_grade:
                            best_rev = deepcopy(path)
                            print 'found besf rev', best_rev
                    continue
                for neighbor in G.successors_iter(path.last()):  # go over all successors of last win of path
                    temp_path = deepcopy(path)
                    temp_path.add_win(neighbor)
                    Q.append(temp_path)
            print 'best_fwd', best_fwd
            print 'best_rev', best_rev
            if best_fwd.total_grade < best_rev.total_grade:
                best_path = best_fwd.path
                best_score = best_fwd.total_grade
                sec_best_path = best_rev.path
                sec_best_score = best_rev.total_grade
            else:
                best_path = best_rev.path
                best_score = best_rev.total_grade
                sec_best_path = best_fwd.path
                sec_best_score = best_fwd.total_grade
            print 'elapsed time for BFS', timeit.default_timer() - tic
        return best_path, best_score, sec_best_path, sec_best_score

    def cst_between(self, win1, win2):
        """
        :param win1: a WinGrade
        :type win1: WinGrade
        :param win2: a WinGrade
        :type win2: WinGrade
        :return:
        """
        print 'cst between', self.csts.tm_pos
        print win1
        print win2
        if self.csts.tm_pos is None:
            print 'None so False'
            return False
        for cst in self.csts.tm_pos:
            if cst[0]-self.csts.tm_pos_fidelity > win1.begin:
                if cst[1]+self.csts.tm_pos_fidelity < win2.begin:
                    print 'found', cst, ' so True'
                    return True
            else:
                print 'found', cst, ' so False'
                return False
        print 'non csts found'
        return False

    def cst_after(self, win):
        """
        :param win: a WinGrade instance
        :return: a (beginning, end, direction) tuple from tm_pos that is after win
        """
        if self.csts.tm_pos is None:
            return None
        for cst in sorted(self.csts.tm_pos):
            if cst[0]-self.csts.tm_pos_fidelity > win.begin:
                return cst
        return None

    # @property
    def topo_graph_old(self):
        """
        topology determination using graph theory. wins are nodes
        :return: list of WinGrades describing best topologies
        """
        import networkx as nx
        from WinGrade import WinGrade
        import operator
        # print self.WinGrades
        win_list = [a for a in self.WinGrades if a.grade < HP_THRESHOLD and a.charges < 3]  # make copy of negative wingrades list
        G = nx.DiGraph()
        source_node = WinGrade(0, 0, 'fwd', '', self.polyval, poly_param)   # define source win
        if self.with_msa:
            [G.add_edge(source_node, a, weight=a.msa_grade) for a in win_list]  # add all win to source edges with msa
        else:
            [G.add_edge(source_node, a, weight=a.grade) for a in win_list]  # add all win to source edges
        for win1 in win_list:   # add all win to win edges where condition applies
            for win2 in win_list:
                # condition: non-overlapping, 2 is after 1, opposite directions
                if not win1.grade_grade_colliding(win2) and win2.begin > win1.end and win1.direction != win2.direction \
                        and win2.begin-win1.end >= 2:
                    if not self.with_msa:
                        G.add_edge(win1, win2, weight=win2.grade)
                    elif self.with_msa:
                        G.add_edge(win1, win2, weight=win2.msa_grade)
        # use bellman-ford algorithm to find minimum paths to all nodes
        # print win_list
        print "About to Bellman-Ford"
        pred, dist = nx.bellman_ford(G, source_node)
        print "Finished Bellman-Fording"
        # force #TM in topology to be known_tm_num
        sorted_dist = sorted(dist.items(), key=operator.itemgetter(1))
        best_path_val = [sorted_dist[0][0]]
        min_val = sorted_dist[0][1]
        if best_path_val == [source_node]:
            print 'topology not found'
            return best_path_val, min_val, best_path_val, min_val
        best_path = self.find_graph_path(pred, best_path_val, source_node)
        copy_sorted_dist = sorted_dist[:]
        # print 'pred', pred
        # print 'dist', dist
        # print 'sorted_dict', sorted_dist
        # print 'best_path_val', best_path_val
        # print 'min_val', min_val
        # print 'best_path', best_path
        while known_tm_num != len(best_path) and known_tm_num != -100:
            copy_sorted_dist = copy_sorted_dist[1:]
            best_path_val = [copy_sorted_dist[0][0]]
            min_val = copy_sorted_dist[0][1]
            best_path = self.find_graph_path(pred, best_path_val, source_node)
        # if no best path is found, the minimal energy win grade is chosen as a single window
        if len(best_path) < known_tm_num != -100:
            best_path_val = [sorted_dist[0][0]]
            min_val = sorted_dist[0][1]
            best_path = self.find_graph_path(pred, best_path_val, source_node)
        # force secondary best topo to have known_tm_num if known, and opposite topology
        sec_best_path_val = [sorted_dist[0][0]]
        sec_min_val = sorted_dist[0][1]
        sec_best_path = self.find_graph_path(pred, sec_best_path_val, source_node)
        copy_sorted_dist = sorted_dist[:]
        while (known_tm_num != len(sec_best_path) or known_tm_num == -100) and \
                        sec_best_path[-1].direction == best_path[-1].direction:
            copy_sorted_dist = copy_sorted_dist[1:]
            sec_best_path_val = [copy_sorted_dist[0][0]]
            sec_min_val = copy_sorted_dist[0][1]
            sec_best_path = self.find_graph_path(pred, sec_best_path_val, source_node)
            # a condition to solve cases where the sec gets to an empty solution
            if sec_best_path == []:
                sec_best_path_val = [sorted_dist[0][0]]
                sec_best_path = self.find_graph_path(pred, sec_best_path_val, source_node)
                break
        # if no best path is found, the minimal energy win grade is chosen as a single window
        if len(sec_best_path) < known_tm_num and known_tm_num != -100:
            sec_best_path_val = [sorted_dist[0][0]]
            sec_min_val = sorted_dist[0][1]
            sec_best_path = self.find_graph_path(pred, sec_best_path_val, source_node)

        return best_path, min_val, sec_best_path, sec_min_val

    @property
    def topo_graph_rel_old(self):
        """
        topology determination using graph theory. wins are nodes
        :return: list of WinGrades describing best topologies
        """
        import networkx as nx
        from WinGrade import WinGrade, WinGradePath
        import operator
        # print self.WinGrades

        win_list = [a for a in self.WinGrades if a.grade < HP_THRESHOLD and a.charges < 3]  # make copy of negative wingrades list
        G = nx.DiGraph()
        print win_list
        source_node = WinGrade(0, 0, 'fwd', 'SOURCE', self.polyval, poly_param)   # define source win
        sink_node = WinGrade(0, 0, 'fwd', 'SINK', self.polyval, poly_param)   # define sink win
        print 'SOURCE NODE', source_node
        for a in win_list:
            print a
            G.add_edge(source_node, a, weight=a.grade)
        if self.with_msa:
            [G.add_edge(source_node, a, weight=a.msa_grade) for a in win_list]  # add all win to source edges with msa
        else:
            [G.add_edge(source_node, a, weight=a.grade) for a in win_list]  # add all win to source edges
        [G.add_edge(a, sink_node, weight=0) for a in win_list]
        for win1 in win_list:   # add all win to win edges where condition applies
            for win2 in win_list:
                # condition: non-overlapping, 2 is after 1, opposite directions
                if not win1.grade_grade_colliding(win2) and win2.begin > win1.end and win1.direction != win2.direction \
                        and win2.begin-win1.end >= 2:
                    if not self.with_msa:
                        G.add_edge(win1, win2, weight=win2.grade)
                    elif self.with_msa:
                        G.add_edge(win1, win2, weight=win2.msa_grade)

        unsat_csts = self.csts.unsat_cst(win_list)
        print 'UNSAT CSTS', unsat_csts
        path_segs = [[source_node]] + unsat_csts + [[sink_node]]
        best_paths = []
        for seg in range(len(path_segs)-1):
            print 'at seg', seg
            for first in path_segs[seg]:
                for second in path_segs[seg+1]:
                    print "looking at segment", first, second
                    best_paths.append(self.shortest_path(G, first, second))
        print 'best paths', best_paths

        G_segs = nx.DiGraph()
        # connect a "source node list" to all first nodes
        source_path = WinGradePath([WinGrade(0, 0, 'fwd', 'SOURCEPATH', self.polyval, poly_param)])
        print 'source_path type', type(source_path)
        [G_segs.add_edge(source_path, a, weight=a.total_grade) for a in best_paths]
        for path1 in best_paths:
            print 'path1', path1
            for path2 in best_paths:
                print 'path1', path1
                print 'path2', path2
                if path1.last().same_as_other(path2.first()):
                    print 'connecting'
                    G_segs.add_edge(path1, path2, weight=path2.total_grade-path1.first().grade)
        print "MADE THE SECOND GRAPH", G_segs
        print 'source_path type', type(source_path)
        # use bellman-ford algorithm to find minimum paths to all nodes
        # print win_list
        print "About to Bellman-Ford"
        pred, dist = nx.bellman_ford(G_segs, source_path)
        print "Finished Bellman-Fording"
        # force #TM in topology to be known_tm_num
        sorted_dist = sorted(dist.items(), key=operator.itemgetter(1))
        best_path_last = [sorted_dist[0][0]]
        min_val = sorted_dist[0][1]

        if best_path_last == [source_node]:
            print 'topology not found'
            return best_path_last, min_val, best_path_last, min_val

        temp_dir = 'A'

        print 'SORTED DIST', sorted_dist

        for last_path, val in sorted_dist:
            if last_path is None:
                continue
            print 'last_path', last_path, val
            # print self.find_graph_path(pred, [last_win], source_node)
            temp_path = self.find_path_frompaths(pred, last_path, source_path)
            if temp_path == []:
                continue
            print 'ABBA', temp_path
            continue
            import sys
            sys.exit()
            if temp_dir != 'A' and self.csts.test_manager(temp_path) and temp_path[-1].direction != temp_dir:
                sec_best_path = temp_path
                sec_best_score = val
                break
            elif self.csts.test_manager(temp_path) and temp_path[-1].direction != temp_dir:
                # print 'I PASSED !!!!'
                best_path = temp_path
                best_score = val
                temp_dir = temp_path[-1].direction
        return best_path, best_score, sec_best_path, sec_best_score

    @property
    def topo_graph_12_8(self):
        """
        topology determination using graph theory. wins are nodes
        :return: list of WinGrades describing best topologies
        """
        import networkx as nx
        from WinGrade import WinGrade, WinGradePath
        import operator
        # print self.WinGrades

        win_list = [a for a in self.WinGrades if a.grade < HP_THRESHOLD and a.charges < 3]  # make copy of negative wingrades list
        G = nx.DiGraph()
        print win_list
        source_node = WinGrade(0, 0, 'fwd', 'SOURCE', self.polyval, poly_param)   # define source win
        sink_node = WinGrade(0, 0, 'fwd', 'SINK', self.polyval, poly_param)   # define sink win
        print 'SOURCE NODE', source_node
        for a in win_list:
            print a
            G.add_edge(source_node, a, weight=a.grade)
        if self.with_msa:
            [G.add_edge(source_node, a, weight=a.msa_grade) for a in win_list]  # add all win to source edges with msa
        else:
            [G.add_edge(source_node, a, weight=a.grade) for a in win_list]  # add all win to source edges
        [G.add_edge(a, sink_node, weight=0) for a in win_list]
        for win1 in win_list:   # add all win to win edges where condition applies
            for win2 in win_list:
                # condition: non-overlapping, 2 is after 1, opposite directions
                if not win1.grade_grade_colliding(win2) and win2.begin > win1.end and win1.direction != win2.direction \
                        and win2.begin-win1.end >= 2:
                    if not self.with_msa:
                        G.add_edge(win1, win2, weight=win2.grade)
                    elif self.with_msa:
                        G.add_edge(win1, win2, weight=win2.msa_grade)

        unsat_csts = self.csts.unsat_cst(win_list)
        print 'UNSAT CSTS', unsat_csts

        import sys
        # sys.exit()

        # use bellman-ford algorithm to find minimum paths to all nodes
        # print win_list
        print "About to Bellman-Ford"
        pred, dist = nx.bellman_ford(G, source_node)
        print "Finished Bellman-Fording"
        # force #TM in topology to be known_tm_num
        sorted_dist = sorted(dist.items(), key=operator.itemgetter(1))


        path_segs = [[source_node]] + unsat_csts + [[sink_node]]
        print 'path_segs', path_segs
        best_paths = []
        for seg in range(len(path_segs)-1):
            print 'at seg', seg
            for first in path_segs[seg]:
                for second in path_segs[seg+1]:
                    print "looking at segment", first, second
                    best_paths.append(self.path_source_sink(pred, sorted_dist, first, second))
                    # best_paths.append(self.shortest_path(G, first, second))

        sys.exit()
        '''
        THE IDEA IS THIS: do Bellman-Ford once. go through all segments (source->cst->cst->...->sink),
        use path_source_sink to find the minimal path between them. join these paths to a secondary graph and BF again.
        go over these results, and look for ones that pass csts.
        I got to creating path_source_sink, have an infinite loop in find_graph_path. fix this.

        THIS CONCEPT IS WRONG BECAUSE BF WILL NEVE ATTACHTO TTHE TARGET, THE SINK OF THE CONSTRAINT

        I have to force find_graph_path to add the sink if necessary (the cst win)
        '''
        best_path_last = [sorted_dist[0][0]]
        min_val = sorted_dist[0][1]

        if best_path_last == [source_node]:
            print 'topology not found'
            return best_path_last, min_val, best_path_last, min_val

        temp_dir = 'A'

        print 'SORTED DIST', sorted_dist

        for last_path, val in sorted_dist:
            if last_path is None:
                continue
            print 'last_path', last_path, val
            # print self.find_graph_path(pred, [last_win], source_node)
            temp_path = self.find_path_frompaths(pred, last_path, source_path)
            if temp_path == []:
                continue
            print 'ABBA', temp_path
            continue
            import sys
            sys.exit()
            if temp_dir != 'A' and self.csts.test_manager(temp_path) and temp_path[-1].direction != temp_dir:
                sec_best_path = temp_path
                sec_best_score = val
                break
            elif self.csts.test_manager(temp_path) and temp_path[-1].direction != temp_dir:
                # print 'I PASSED !!!!'
                best_path = temp_path
                best_score = val
                temp_dir = temp_path[-1].direction
        return best_path, best_score, sec_best_path, sec_best_score

    def path_source_sink(self, pred, sorted_dist, source, sink):
        """
        :param pred: pred dict from bellman ford
        :param sorted_dist: sorted dist dict from bellman ford
        :param source: a source
        :param sink: a sink
        :return: the most negative path between source and sink
        """
        for win, val in sorted_dist:
            if win.same_as_other(sink):
                return self.find_graph_path(pred, [win], source)

    def shortest_path(self, g, source, sink):
        import networkx as nx
        from numpy import inf as inf
        from WinGrade import WinGradePath
        all_paths = nx.all_shortest_paths(g, source, sink, 'grade')
        best_path, best_path_grade = ([], inf)
        for path in all_paths:
            # print 'path', len(path), sum([a.grade for a in path])
            path_grade = sum([a.grade for a in path])
            if path_grade < best_path_grade:
                print 'paht better energy', path
                if self.csts.test_manager_segment(path, (source.begin, sink.end), True):
                    best_path_grade = path_grade
                    best_path = path[:]
        print 'returning', best_path
        if best_path == []:
            print 'NOTHING PASSED CONSTRAINTS', source, sink
        # if best_path[-1].seq == 'SINK':
        #     print 'last is SINK'
        #     best_path.pop(-1)
        return WinGradePath(best_path)

    def shortest_path_old(self, g, source, sink):
        import networkx as nx
        from numpy import inf as inf
        from WinGrade import WinGradePath
        all_paths = nx.all_shortest_paths(g, source, sink, 'grade')
        best_path, best_path_grade = ([], inf)
        for path in all_paths:
            # print 'path', len(path), sum([a.grade for a in path])
            path_grade = sum([a.grade for a in path])
            if path_grade < best_path_grade:
                print 'paht better energy', path
                if self.csts.test_manager_segment(path, (source.begin, sink.end), True):
                    best_path_grade = path_grade
                    best_path = path[:]
        print 'returning', best_path
        if best_path == []:
            print 'NOTHING PASSED CONSTRAINTS', source, sink
        # if best_path[-1].seq == 'SINK':
        #     print 'last is SINK'
        #     best_path.pop(-1)
        return WinGradePath(best_path)

    # @property
    def topo_graph_by_segs(self):
        """
        :return: best and second best paths and their scores
        segments the sequence into [source - cst_1 - cst_2 - ... - cst_i - sink] such that every cst_i is a constraint
        of the sequence. will create a graph from all WinGrades within every segment, and use all available paths for
        that segment to create a subsequent graph of paths. Bellman-Ford is used agin to find the minimum of the
        resulting graph
        """
        import networkx as nx
        from WinGrade import WinGrade, WinGradePath
        from WinGrade import flatten_path_list
        import operator
        # print self.WinGrades

        win_list = [a for a in self.WinGrades if a.grade < HP_THRESHOLD and a.charges < 3]  # make copy of negative wingrades list
        # G = nx.DiGraph()
        # print win_list
        # source_node = WinGrade(0, 0, 'fwd', 'SOURCE', self.polyval, poly_param)   # define source win
        # sink_node = WinGrade(0, 0, 'fwd', 'SINK', self.polyval, poly_param)   # define sink win
        # print 'SOURCE NODE', source_node
        # [G.add_edge(source_node, a, weight=a.grade) for a in win_list]
        # if self.with_msa:
        #     [G.add_edge(source_node, a, weight=a.msa_grade) for a in win_list]  # add all win to source edges with msa
        # else:
        #     [G.add_edge(source_node, a, weight=a.grade) for a in win_list]  # add all win to source edges
        # [G.add_edge(a, sink_node, weight=0) for a in win_list]
        # for win1 in win_list:   # add all win to win edges where condition applies
        #     for win2 in win_list:
        #         # condition: non-overlapping, 2 is after 1, opposite directions
        #         if not win1.grade_grade_colliding(win2) and win2.begin > win1.end and win1.direction != win2.direction \
        #                 and win2.begin-win1.end >= 2:
        #             if not self.with_msa:
        #                 G.add_edge(win1, win2, weight=win2.grade)
        #             elif self.with_msa:
        #                 G.add_edge(win1, win2, weight=win2.msa_grade)
        # unsat_csts = self.csts.unsat_cst(win_list)
        # print 'UNSAT CSTS', unsat_csts

        print 'CSTS', self.csts.tm_pos

        source_cst = (-2, -1, None)
        sink_cst = (self.seq_length+1, self.seq_length+2, None)
        if self.csts.tm_pos is not None:
            path_segs = [source_cst] + self.csts.tm_pos + [sink_cst]
        else:
            path_segs = [source_cst, sink_cst]
        print 'path_segs', path_segs
        best_paths = []
        for seg in range(len(path_segs)-1):
            print 'at seg', seg, path_segs[seg], path_segs[seg+1]
            best_paths.extend(self.paths_source2sink(path_segs[seg], path_segs[seg+1], win_list))

        # print 'best_paths', best_paths
        print 'Creating graph of Paths TATATATA@@@@@'
        print 'there are %i potential nodes' % len(best_paths)

        source_path_node = WinGradePath([WinGrade(0, 0, 'fwd', 'SOURCE', self.polyval, poly_param)])
        G = nx.DiGraph()
        [G.add_edge(source_path_node, a, weight=a.total_grade) for a in best_paths]
        for path1 in best_paths:
            for path2 in best_paths:
                if path1.overlap_me_first(path2):
                    G.add_edge(path1, path2, weight=path2.total_grade-path1.first().grade)

        print 'MADE THE MOTHA_FUCKING GRAPH', G

        pred, dist = nx.bellman_ford(G, source_path_node)
        sorted_dist = sorted(dist.items(), key=operator.itemgetter(1))
        print 'finished BF, and sorted dist is:'
        temp_direction = 'A'
        best_path, sec_best_path = [], []
        best_score, sec_best_score = None, None

        for last_path, total_grade in sorted_dist:
            temp_path_list = self.path_list_from_pred(pred, last_path, source_path_node)
            temp_path_list_flat = flatten_path_list(temp_path_list[:])
            # print 'temp_path_list', temp_path_list_flat
            if temp_direction != 'A' and self.csts.test_manager(temp_path_list_flat) \
                    and temp_path_list_flat[-1].direction is not temp_direction:
                print 'looking at sec', temp_path_list_flat[-1].direction, temp_direction
                sec_best_path = temp_path_list_flat[:]
                if sec_best_path[0].seq == 'SOURCE':
                    sec_best_path.pop(0)
                sec_best_score = sum([a.grade for a in sec_best_path])
                break
            elif self.csts.test_manager(temp_path_list_flat) and temp_path_list_flat[-1].direction != temp_direction:
                print 'looking at BEST'
                best_path = temp_path_list_flat[:]
                if best_path[0].seq == 'SOURCE':
                    best_path.pop(0)
                best_score = sum([a.grade for a in best_path])
                temp_direction = best_path[-1].direction
                print 'now temp_direction is:', temp_direction
        print 'found best'
        print best_path[-1].direction
        print best_score
        print 'fiound sec best',
        print sec_best_path[-1].direction
        print sec_best_score
        return best_path, best_score, sec_best_path, sec_best_score

    def path_list_from_pred(self, pred, last_path, source):
        """
        :param pred: Bellman-Ford pred dictionary (win: predecessor)
        :param last_path: last win in Bellman-Ford result
        :param source: the source node, where to stop
        :return: the path, made of WinGRadePaths, that ends at last_path
        """
        # find all wins in the minimal energy path from source to last win
        # print 'source_node', source_node
        path = [last_path]
        while not path[-1].same_as_other(source):
            for k, v in pred.items():
                if v is None:
                    continue
                if k.same_as_other(path[-1]):
                    path.append(v)
        path.pop(-1)    # get rid of source_node
        return path[::-1]   # revert the path to begin->end

    def paths_source2sink(self, source, sink, wins):
        """
        :param source: a source segment (start, end, direction) contstraint (or actual source/sink)
        :param sink: a source segment (start, end, direction) contstraint (or actual source/sink)
        :param wins: the availabel list of WinGrades
        :return: a list of all Bellman-Ford minimal paths between all possible WinGrades in source and sink

        MAYBE: make return only best 100, so that the graphs are much smaller...

        """
        from WinGrade import WinGrade, WinGradePath
        from TMConstraint import all_satisfying_wins
        import networkx as nx
        seg = (source[0]-self.csts.tm_pos_fidelity, sink[1]+self.csts.tm_pos_fidelity) #  the segment for the graph
        seg_wins = [a for a in wins if a.within_segment(seg)]
        G = nx.DiGraph()
        source_node = WinGrade(0, 0, 'fwd', 'SOURCE', self.polyval, poly_param)   # define source win

        if self.with_msa:
            [G.add_edge(source_node, a, weight=a.msa_grade) for a in seg_wins]
        else:
            [G.add_edge(source_node, a, weight=a.grade) for a in seg_wins]
        for win1 in seg_wins:
            for win2 in seg_wins:
                if not win1.grade_grade_colliding(win2) and \
                                win2.begin > win1.end \
                        and win1.direction != win2.direction \
                        and win2.begin-win1.end >= 2:
                    if not self.with_msa:
                        G.add_edge(win1, win2, weight=win2.grade)
                    elif self.with_msa:
                        G.add_edge(win1, win2, weight=win2.msa_grade)
        if source == (-2, -1, None):
            sources = [source_node]
        else:
            sources = all_satisfying_wins(source, wins, self.csts.tm_pos_fidelity)
        all_paths = []
        for src in sources:
            pred, dist = nx.bellman_ford(G, src)
            for k, v in pred.items():
                if k.within_segment((sink[0], sink[1]), self.csts.tm_pos_fidelity) or \
                                sink == (self.seq_length+1, self.seq_length+2, None):
                    temp_path = self.find_graph_path(pred, [k], src)
                    if source == (-2, -1, None) or \
                            temp_path[0].within_segment((source[0], source[1], None), self.csts.tm_pos_fidelity):
                        all_paths.append(WinGradePath(temp_path))

        return all_paths