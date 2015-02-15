def main():
    global hydrophobicity_polyval
    hydrophobicity_polyval = MakeHydrophobicityGrade()
    # temp = HphobicityScore('temp', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    # temp = HphobicityScore('temp', 'YSYRFVWWAISTAAMLYILY')
    temp = HphobicityScore('1E12', 'MSITSVPGVVDAGVLGAQSAAAVRENALLSSSLWVNVALAGIAILVFVYMGRTIRPGRPRLIWGATLMIPLVSISSYLGLLSGLTVGMIEMPAGHALAGEMVRSQWGRYLTWALSTPMILLALGLLADVDLGSLFTVIAADIGMCVTGLAAAMTTSALLFRWAFYAISCAFFVVVLSALVTDWAASASSAGTAEIFDTLRVLTVVLWLGYPIVWAVGVEGLALVQSVGVTSWAYSVLDVFAKYVFAFILLRWVANNERTVAVAGQTLGTMSSDD')
    # temp = HphobicityScore('1IWG', 'MPNFFIDRPIFAWVIAIIIMLAGGLAILKLPVAQYPTIAPPAVTISASYPGADAKTVQDTVTQVIEQNMNGIDNLMYMSSNSDSTGTVQITLTFESGTDADIAQVQVQNKLQLAMPLLPQEVQQQGVSVEKSSSSFLMVVGVINTDGTMTQEDISDYVAANMKDAISRTSGVGDVQLFGSQYAMRIWMNPNELNKFQLTPVDVITAIKAQNAQVAAGQLGGTPPVKGQQLNASIIAQTRLTSTEEFGKILLKVNQDGSRVLLRDVAKIELGGENYDIIAEFNGQPASGLGIKLATGANALDTAAAIRAELAKMEPFFPSGLKIVYPYDTTPFVKISIHEVVKTLVEAIILVFLVMYLFLQNFRATLIPTIAVPVVLLGTFAVLAAFGFSINTLTMFGMVLAIGLLVDDAIVVVENVERVMAEEGLPPKEATRKSMGQIQGALVGIAMVLSAVFVPMAFFGGSTGAIYRQFSITIVSAMALSVLVALILTPALCATMLKPIAKGDHGEGKKGFFGWFNRMFEKSTHHYTDSVGGILRSTGRYLVLYLIIVVGMAYLFVRLPSSFLPDEDQGVFMTMVQLPAGATQERTQKVLNEVTHYYLTKEKNNVESVFAVNGFGFAGRGQNTGIAFVSLKDWADRPGEENKVEAITMRATRAFSQIKDAMVFAFNLPAIVELGTATGFDFELIDQAGLGHEKLTQARNQLLAEAAKHPDMLTSVRPNGLEDTPQFKIDIDQEKAQALGVSINDINTTLGAAWGGSYVNDFIDRGRVKKVYVMSEAKYRMLPDDIGDWYVRAADGQMVPFSAFSSSRWEYGSPRLERYNGLPSMEILGQAAPGKSTGEAMELMEQLASKLPTGVGYDWTGMSYQERLSGNQAPSLYAISLIVVFLCLAALYESWSIPFSVMLVVPLGVIGALLAATFRGLTNDVYFQVGLLTTIGLSAKNAILIVEFAKDLMDKEGKGLIEATLDAVRMRLRPILMTSLAFILGVMPLVISTGAGSGAQNAVGTGVMGGMVTATVLAIFFVPVFFVVVRRRFSRKNEDIEHSHTVDHH')
    # temp = HphobicityScore('1BRX', 'EAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSAGDGAAATS')
    # temp.print_HphobicityScore()
    # temp.local_minima_finder()
    temp.plot_energy_landscape()
    pass


class HphobicityScore():
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.fwd_grades = self.grade_seq(self.seq)
        self.rev_grades = self.grade_seq(self.seq[::-1])[::-1]
        self.rev_grades = [{len(self.seq)-x.keys()[0]: x.values()[0]} for x in self.rev_grades]

        # self.fwd_min = [min(x) for x in self.fwd_grades]
        # self.rev_min = [min(x) for x in self.rev_grades]

    def grade_segment(self, segment):
        import numpy as np
        membrane_position = np.linspace(-15, 15, endpoint=True, num=len(segment))
        grade = 0
        for i, aa in enumerate(segment):
            grade += np.polyval(hydrophobicity_polyval[aa], membrane_position[i])
        return grade

    def grade_window(self, window):
        new_grade = []
        for inc in range(0, min([17, len(window)-19])):
            new_grade.append(self.grade_segment(window[:20+inc]))
        return new_grade

    def grade_seq(self, seq):
        window_grades = []
        for i in range(0, len(seq)-19):
            window_grades.append({i: self.grade_window(seq[i:i+35])})
        return window_grades

    def print_HphobicityScore(self):
        print 'printing ', self.name
        # print self.fwd_grades
        # print self.rev_grades
        # for i, window in enumerate(self.fwd_grades):
        #     print window, self.fwd_min[i]
            # print i, window
        for i in self.rev_grades:
            print i, type(i)
            for j, grd in enumerate(i.values()[0]):
                print j, grd
        # for i, window in enumerate(self.rev_grades):
        #     print i, window

    def is_clash(self, min_point, point_list):
        # supposed to check if a point (start, length) is within the range of a list of points of the same structure
        for point in point_list:
            if range(point['start'], point['start']+point['length']) in \
                    range(min_point['start'], min_point['start']+min_point['length']):
                print 'aaa'
                pass

    def local_minima_finder(self):
        from collections import OrderedDict
        data = {grade: {'start': pos, 'length': 20+inc} for pos, windows in enumerate(self.fwd_grades)
                for inc, grade in enumerate(windows)}
        data_sorted = OrderedDict(sorted(data.items(), key=lambda x: x[0]))
        chosen_minimas = []
        # this loop will test each point in the sorted dictionary, and will keep the mionimums that are not clashing
        for grade, val in data_sorted.items():
            print grade, val['start'], val['length']
            print self.is_clash(val, chosen_minimas)

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


def MakeHydrophobicityGrade():
    # hydrophobicity_grade = open('Poly_Values.txt', 'r')
    hydrophobicity_polyval = {}
    hydrophobicity_grade = open('poly_value_11.2.txt', 'r')
    for line in hydrophobicity_grade:
        split = line.split()
        hydrophobicity_polyval[split[0]] = [float(n) for n in split[1:6]]
    hydrophobicity_grade.close()
    return hydrophobicity_polyval


if __name__ == '__main__':
    main()