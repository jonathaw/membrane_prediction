class HphobicityScore():
    def __init__(self, name, seq, ss2_file, hydro_polyval, param_list):
        '''
        :param name: protein's name
        :param uniprot: the entrie's uniprot code
        :param seq: protein's sequence
        :param hydro_polyval: polynom valued dictionary
        :return: a stack of WinGrade instances, and their utilities
        '''
        global HP_THRESHOLD, INC_MAX, MIN_LENGTH, PSI_HELIX, PSI_RES_NUM, known_tm_num, poly_param, output_path
        HP_THRESHOLD = param_list['hp_threshold']
        INC_MAX = 11  # 8
        MIN_LENGTH = param_list['min_length']
        PSI_HELIX = param_list['psi_helix']
        PSI_RES_NUM = param_list['psi_res_num']
        known_tm_num = param_list['known_tm_num']
        poly_param = {k: v for k, v in param_list.items() if k in ['c0', 'c1', 'c2', 'c3', 'z_0', 'w']}
        print 'poly', poly_param
        output_path = param_list['result_path']
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
        self.best_c_term = 'out' if self.topo_best[-1].direction == 'fwd' else 'in'
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
        from WinGrade import WinGrade
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
                if not self.is_not_helical((i, i+MIN_LENGTH+inc), psi):
                    if self.with_msa:
                        grades.append(msa_obj.retrieve_seqs(i, i+MIN_LENGTH+inc, 'fwd'))
                        grades.append(msa_obj.retrieve_seqs(i, i+MIN_LENGTH+inc, 'rev'))
                    else:
                        try:
                            # print 'creating win %s' % self.seq[i:i+MIN_LENGTH+inc]
                            grades.append(WinGrade(i, i+MIN_LENGTH+inc, 'fwd', self.seq[i:i+MIN_LENGTH+inc], self.polyval,
                                                poly_param))
                            grades.append(WinGrade(i, i+MIN_LENGTH+inc, 'rev', self.seq[i:i+MIN_LENGTH+inc][::-1], self.polyval,
                                                poly_param))
                        except:
                            continue
            # writes to the progress bar
            sys.stdout.write("-")
            sys.stdout.flush()
        sys.stdout.write("\n")
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
        while path[-1].seq != source_node.seq:
            for k, v in pred.items():
                if v is None:
                    continue
                if k.seq == path[-1].seq:
                    path.append(v)
        path.pop(-1)    # get rid of source_node
        return path[::-1]   # revert the path to begin->end

    @property
    def topo_graph_draft(self):
        """
        topology determination using graph theory. wins are nodes
        :return: list of WinGrades describing best topologies
        """
        import networkx as nx
        from WinGrade import WinGrade
        import operator
        win_list = [a for a in self.WinGrades if a.grade < HP_THRESHOLD]  # make copy of negative wingrades list
        G = nx.DiGraph()
        source_node = WinGrade(0, 0, 'fwd', '', self.polyval, poly_param)   # define source win
        [G.add_edge(source_node, a, weight=a.grade) for a in win_list]  # add all win to source edges
        for win1 in win_list:   # add all win to win edges where condition applies
            for win2 in win_list:
                # condition: non-overlapping, 2 is after 1, opposite directions
                if not win1.grade_grade_colliding(win2) and win2.begin > win1.end and win1.direction != win2.direction:
                    G.add_edge(win1, win2, weight=win2.grade)
        # use bellman-ford algorithm to find minimum paths to all nodes
        # print win_list
        pred, dist = nx.bellman_ford(G, source_node)
        # force #TM in topology to be known_tm_num
        sorted_dist = sorted(dist.items(), key=operator.itemgetter(1))
        # print known_tm_num, sorted_dist
        booli = True
        while booli:
            # min_val = min(dist.values())
            # best_path_val = [k for k, v in dist.items() if v == min_val]
            best_path_val = [sorted_dist[0][0]]
            min_val = sorted_dist[0][1]
            best_path = self.find_graph_path(pred, best_path_val, source_node)
            if len(best_path) == known_tm_num or known_tm_num == -100:
                booli = False
            else:
                print 'not long enough', len(best_path), known_tm_num
                # dist = {k: v for k, v in dist.items() if v != min_val}
                sorted_dist.pop()
        if len(best_path) == 1 and known_tm_num != 1 and known_tm_num != -100:
            print 'Could not find %i TMHs, outputing best option' % known_tm_num
            sorted_dist = sorted(dist.items(), key=operator.itemgetter(1))
            best_path_val = [sorted_dist[0][0]]
            min_val = sorted_dist[0][1]
            best_path = self.find_graph_path(pred, best_path_val, source_node)

        # if no best path is found, the minimal energy win grade is chosen as a single window
        if best_path == []:
            best_path_val = min([v for k, v in dist.items() if k.seq != ''])
            best_path = [k for k, v in dist.items() if v == best_path_val]

        # force secondary best topology to have known_tm_num
        for k, v in sorted_dist:
            if k.direction != best_path[-1].direction:
                sec_best_path = [k]
                sec_best_val = v
                if len(sec_best_path) == known_tm_num or known_tm_num == -100:
                    break
        if len(sec_best_path) == 1 and known_tm_num != 1 and known_tm_num != -100:
            print 'Could not find %i TMHs for second best, outputing best option' % known_tm_num
            for k, v in sorted(dist.items(), key=operator.itemgetter(1)):
                if k.direction != best_path[-1].direction:
                    sec_best_path = [k]
                    sec_best_val = v
                    break

        # sec_best_path, sec_best_val = [(k, v) for k, v in sorted(dist.items(), key=operator.itemgetter(1)) if k.direction != best_path[-1].direction]
        sec_best_path = self.find_graph_path(pred, sec_best_path, source_node)
        # print 'topo finished'
        return best_path, min_val, sec_best_path, sec_best_val

    @property
    def topo_graph(self):
        """
        topology determination using graph theory. wins are nodes
        :return: list of WinGrades describing best topologies
        """
        import networkx as nx
        from WinGrade import WinGrade
        import operator
        # print self.WinGrades
        win_list = [a for a in self.WinGrades if a.grade < HP_THRESHOLD]  # make copy of negative wingrades list
        G = nx.DiGraph()
        source_node = WinGrade(0, 0, 'fwd', '', self.polyval, poly_param)   # define source win
        if self.with_msa:
            [G.add_edge(source_node, a, weight=a.msa_grade) for a in win_list]  # add all win to source edges with msa
        else:
            [G.add_edge(source_node, a, weight=a.grade) for a in win_list]  # add all win to source edges
        for win1 in win_list:   # add all win to win edges where condition applies
            for win2 in win_list:
                # condition: non-overlapping, 2 is after 1, opposite directions
                if not win1.grade_grade_colliding(win2) and win2.begin > win1.end and win1.direction != win2.direction:
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