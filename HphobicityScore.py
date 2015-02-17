class HphobicityScore():
    def __init__(self, name, seq, hydro_polyval):
        '''
        :param name: protein's name
        :param seq: protein's sequence
        :param hydro_polyval: polynom valued dictionary
        :return: a stack of WinGrade instances, and their utilities
        '''
        global HP_THRESHOLD
        HP_THRESHOLD = 0.
        self.name = name
        self.seq = seq
        self.seq_length = len(seq)
        self.polyval = hydro_polyval
        self.WinGrades = self.win_grade_generator()
        self.sorted_grade = self.sort_WinGrades()
        self.minimas = self.local_minima_finder(direction='both')
        self.topo_minimas = self.topo_determine()
        self.n_term_orient = self.topo_minimas[0].direction
        # self.sorted_grade_norm = self.sort_WinGrades_norm()
        # self.minimas_norm = self.local_minima_finder_norm()

    def win_grade_generator(self):
        '''
        :return:grades all segments of self.seq, and aggregates them as WinGrades
        '''
        from WinGrade import WinGrade
        grades = []
        for i in range(self.seq_length):
            for inc in range(16):
                if i+20+inc > self.seq_length:
                    continue
                grades.append(WinGrade(i, i+20+inc, 'fwd', self.seq[i:i+20+inc], self.polyval))
                grades.append(WinGrade(i, i+20+inc, 'rev', self.seq[i:i+20+inc][::-1], self.polyval))
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
                     color='black' if win_grade.direction == 'fwd' else 'grey')
        for minima in self.minimas:
            plt.plot((minima.begin, minima.end), (minima.grade, minima.grade),
                     color='blue' if minima.direction == 'fwd' else 'red', lw=2)
        print 'found the following minimas:'
        for minima in self.topo_minimas:
            print minima
            plt.plot((minima.begin, minima.end), (minima.grade, minima.grade),
                     color='green' if minima.direction == 'fwd' else 'purple', lw=4)
        black_line = mlines.Line2D([], [], color='black', marker='', lw=2, label='Fwd grade')
        grey_line = mlines.Line2D([], [], color='grey', marker='', lw=2, label='Rev grade')
        blue_line = mlines.Line2D([], [], color='blue', marker='', lw=2, label='Fwd minima')
        red_line = mlines.Line2D([], [], color='red', marker='', lw=2, label='Rev minima')
        green_line = mlines.Line2D([], [], color='green', marker='', lw=2, label='Fwd topo minima')
        purple_line = mlines.Line2D([], [], color='purple', marker='', lw=2, label='Rev topo minima')
        plt.legend(handles=[black_line, grey_line, blue_line, red_line, green_line, purple_line], ncol=2)
        plt.xlabel('Sequence Position')
        plt.ylabel('Energy')
        plt.title('Win Grades Energy Plot for %s' % self.name)
        plt.show()

    def plot_energy_landscape(self):
        # import matplotlib.pyplot as plt
        # from mpl_toolkits.mplot3d import Axes3D
        # from scipy.interpolate import griddata
        # import numpy as np
    #     data = [(pos, 20+inc, grade) for pos, windows in enumerate(self.fwd_grades) for inc, grade in enumerate(windows)]
    #     print data
    #     data = [(pos, 20+inc, grade) for pos, windows in enumerate(self.rev_grades) for inc, grade in enumerate(windows)]
    #     print data
        print [(i.keys()[0], 20+int(inc), grade) for i in self.fwd_grades for inc, grade in enumerate(i.values()[0])]
        print '\n\n\n'
        print [(i.keys()[0], 20+int(inc), grade) for i in self.rev_grades for inc, grade in enumerate(i.values()[0])]
    #     x, y, z = zip(*data)
    #     z = map(float, z)
    #     grid_x, grid_y = np.mgrid[min(x):max(x):100j, min(y):max(y):100j]
    #     grid_z = griddata((x, y), z, (grid_x, grid_y), method='cubic')
    #
    #     fig = plt.figure()
    #     ax = fig.gca(projection='3d')
    #     ax.plot_surface(grid_x, grid_y, grid_z, cmap=plt.cm.Spectral)
    #     plt.show()
        pass

    def topo_determine(self):
        """
        :return:a sorted list of minimas of flipping direction which has the lowest
        global energy
        """
        fwd_minimas = self.local_minima_finder(direction='fwd')
        rev_minimas = self.local_minima_finder(direction='rev')
        fwd_sorted = sorted(fwd_minimas, key=lambda x: x.begin)
        rev_sorted = sorted(rev_minimas, key=lambda x: x.begin)
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
        return fwd_start if fwd_score < rev_score and len(fwd_start) >= len(rev_start)\
            else rev_start