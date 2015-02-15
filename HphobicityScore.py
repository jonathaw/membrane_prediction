class HphobicityScore():
    def __init__(self, name, seq, hydro_polyval):
        '''
        :param name: protein's name
        :param seq: protein's sequence
        :param hydro_polyval: polynom valued dictionary
        :return: a stack of WinGrade instances, and their utilities
        '''
        self.name = name
        self.seq = seq
        self.seq_length = len(seq)
        self.polyval = hydro_polyval
        self.WinGrades = self.win_grade_generator()
        self.sorted_grade = self.sort_WinGrades()
        self.minimas = self.local_minima_finder()
        self.sorted_grade_norm = self.sort_WinGrades_norm()
        self.minimas_norm = self.local_minima_finder_norm()

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

    def local_minima_finder(self):
        '''
        :return:returns a list of the minimal grade non-clashing WinGrades
        '''
        chosen_grades = [self.sort_WinGrades().pop(0).values()[0]]
        for grade in self.sort_WinGrades():
            if not grade.values()[0].set_grade_colliding(chosen_grades):
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
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        for minima in self.local_minima_finder():
            print minima.begin, minima.end, minima.grade
            ax1.plot((minima.begin, minima.end), (minima.grade, minima.grade),
                     color='red' if minima.direction == 'fwd' else 'blue')
        plt.show()

    def plot_energy_landscape(self):
    #     import matplotlib.pyplot as plt
    #     from mpl_toolkits.mplot3d import Axes3D
    #     from scipy.interpolate import griddata
    #     import numpy as np
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


