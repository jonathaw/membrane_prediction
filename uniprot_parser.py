import re
import numpy as np
import matplotlib.pyplot as plt
import csv
import math
import subprocess
# import os

# MEMEBRANE_RANGE = range(-15, 15, 30 / 21)
MEMBRANE_FULL = range(0, 35)
MEMBRANE_THIRD = [0, 1, 4, 8, 11, 15, 18, 19, 22, 26, 29, 33, 36]
MEMBRANE_TWO_THIRDS = [0, 1, 2, 4, 5, 8, 9, 11, 12, 13, 15, 16, 18, 19, 20, 22, 23, 26, 27, 29, 30, 31, 33, 34, 36]
hydrophobicity_polyval = {}
window_grades = {}
uniprot_entris = 0
main_dict = {}
SMOOTH_SIZE = 1
LOOP_BYPASS = 3
MIN_WIN = 20
SECONDARY_MINIMA_THRESHOLD = 7
PRIMARY_MINIMA_THRESHOLD = 4


def MakeHydrophobicityGrade():
    global hydrophobicity_polyval
    hydrophobicity_grade = open('Poly_Values.txt', 'r')
    for line in hydrophobicity_grade:
        split = line.split(' ')
        hydrophobicity_polyval[split[0]] = [float(n) for n in split[1:6]]
    hydrophobicity_grade.close()


def grade_seq(seq, mode):
    sumer = 0
    seq_inverse = seq[::-1]
    sumer_inverse = 0
    not_span = 0
    not_span_inverse = 0
    membrane_position = np.linspace(-15, 15, endpoint=True, num=len(seq))
    if mode == 'full':
        membrane_span = MEMBRANE_FULL
    elif mode == 'third':
        membrane_span = MEMBRANE_THIRD
    elif mode == 'two_thirds':
        membrane_span = MEMBRANE_TWO_THIRDS
    for i, aa in enumerate(seq):
        if i in membrane_span:
            sumer += np.polyval(hydrophobicity_polyval[aa], membrane_position[i])
        else:
            not_span += np.polyval(hydrophobicity_polyval[aa], membrane_position[i])
    for i, aa in enumerate(seq_inverse):
        if i in membrane_span:
            sumer_inverse += np.polyval(hydrophobicity_polyval[aa], membrane_position[i])
        else:
            not_span_inverse += np.polyval(hydrophobicity_polyval[aa], membrane_position[i])
    if not_span >= LOOP_BYPASS:
        sumer += not_span
    if not_span_inverse >= LOOP_BYPASS:
        sumer_inverse += not_span_inverse
    if sumer_inverse > sumer:
        return sumer, 1
    else:
        return sumer_inverse, 0


def hydrophobicity_grade_increments(seq, mode):
    results = []
    dirs = []
    temp = grade_seq(seq[0:MIN_WIN], mode)
    results.append(temp[0])
    dirs.append(temp[1])
    for inc in range(1, 14):
        if inc+MIN_WIN <= len(seq):
            temp = grade_seq(seq[0:MIN_WIN+inc], mode)
            results.append(temp[0])
            dirs.append(temp[1])
    min_val = min(results)
    min_val_index = results.index(min_val)
    return min_val, min_val_index+MIN_WIN, dirs[min_val_index], results


def UniprotParser():
    refseq_re = re.compile('[A-Za-z]{2}_\d*\.\d*')
    trnasmem_re = re.compile('\d{1,}')
    seq_re = re.compile('[A-Z]*')

    # parse uniprot data into hash main_dict
    global main_dict, uniprot_entris
    uniprot_entris = 0
    uniprot_dat_file = open('/Users/jonathan/eden/membrane_prediciton/membrane_prediction/uniprot_sprot_10000.dat', 'r')
    for line0 in uniprot_dat_file:
        temp_dict = {}
        temp_dict['transmem'] = []
        seq = ''
        if line0[0:2] == '//':
           transmem_array = []
           for line1 in uniprot_dat_file:
                if line1[0:12] == 'DR   RefSeq;':
                    names = refseq_re.findall(line1)
                if line1[0:13] == 'FT   TRANSMEM':
                    transmem_array.append(tuple((trnasmem_re.findall(line1)[0:2])))
                if line1[0:2] == 'SQ':
                    for line2 in uniprot_dat_file:
                        if line2[0:2] == '//':
                            line1 = line2
                            break
                        seq += ''.join(seq_re.findall(line2))
                if line1[0:2] == '//':
                    if len(transmem_array) > 0 and len(names) > 0:
                        temp_dict['transmem'] = transmem_array
                        temp_dict['sequence'] = seq
                        main_dict[names[0]] = temp_dict
                        uniprot_entris += 1
                    break


def WriteDataSetToTSV(file):
    out = open(file, 'a')
    out.write('name\tsequence\ttransmem\tfull_window_grades\tthird_window_grades\ttwo_thirds_window_grades\tstarters\tfull_grades_inverse\tthird_grades_inverse\ttwo_thirds_grades_inverse\n')
    for protein in main_dict:
        message = str(protein + '\t' + main_dict[protein]['sequence'] + '\t' + str(main_dict[protein]['transmem']) + '\t')
        message += str(str(window_grades[protein]['full_grades']) + '\t' + str(window_grades[protein]['third_grades']) + '\t')
        message += str(str(window_grades[protein]['two_thirds_grades']) + '\t' + str(window_grades[protein]['starters']) + '\t')
        message += str(str(window_grades[protein]['full_grades_inverse']) + '\t' + str(window_grades[protein]['third_grades_inverse']) + '\t')
        message += str(str(window_grades[protein]['two_thirds_grades_inverse']) + '\n')
        out.write(message)


def WriteSingleSeqToTSV(grades, file):
    out = open(file, 'a')
    title = 'full_window_grades\tfull_window_win_len\tthird_window_grades\tthird_window_grades_win_len\ttwo_thirds_window_grades\ttwo_thirds_window_grades_win_len\tstarters\tfull_grades_inverse\tfull_grades_inverse_win_len\tthird_grades_inverse\tthird_grades_inverse_win_len\ttwo_thirds_grades_inverse\ttwo_thirds_grades_inverse_win_len\n'
    out.write(title)
    message = str(str(grades['full_grades']) + '\t' + str(grades['full_grades_win_len']) + '\t' + str(grades['third_grades']) + '\t' + str(grades['third_grades_win_len']) + '\t')
    message += str(str(grades['two_thirds_grades']) + str(grades['two_thirds_grades_win_len']) + '\t' + str(grades['starters']) + '\t')
    message += str(str(grades['full_grades_inverse']) + '\t' + str(grades['full_grades_inverse_win_len']) + '\t' + str(grades['third_grades_inverse']) + '\t')
    message += str(str(grades['third_grades_inverse_win_len']) + '\t' + str(grades['two_thirds_grades_inverse']) + '\t' + str(grades['two_thirds_grades_inverse_win_len'])+ '\n')
    out.write(message)
    print 'wrote to: ', file


def WriteSingleSeqToTSVusingCSV(grades, file):
    with open(file, 'wa+') as f:
        writer = csv.writer(f)
        writer.writerow(grades.keys())
        writer.writerow(grades.values())


def WindowGradesForSingleSequence(seq, name='default'):
    window_grades_single_chain = {}
    window_grades_single_chain['name'] = name
    full_grades = []
    full_grades_win_len = []
    full_dir = []
    third_grades = []
    third_grades_win_len = []
    third_dir = []
    two_thirds_grades = []
    two_thirds_grades_win_len = []
    two_thirds_dir = []
    starters = []
    for first in range(0, len(seq) - MIN_WIN):
        temp = hydrophobicity_grade_increments(seq[first:first+36], 'full')
        full_grades.append(temp[0])
        full_grades_win_len.append(temp[1])
        full_dir.append(temp[2])
        temp = hydrophobicity_grade_increments(seq[first:first+36], 'third')
        third_grades.append(temp[0])
        third_grades_win_len.append(temp[1])
        third_dir.append(temp[2])
        temp = hydrophobicity_grade_increments(seq[first:first+36], 'two_thirds')
        two_thirds_grades.append(temp[0])
        two_thirds_grades_win_len.append(temp[1])
        two_thirds_dir.append(temp[2])
        starters.append(first)
    window_grades_single_chain['full_grades'] = full_grades
    window_grades_single_chain['full_grades_win_len'] = full_grades_win_len
    window_grades_single_chain['full_dirs'] = full_dir
    window_grades_single_chain['third_grades'] = third_grades
    window_grades_single_chain['third_grades_win_len'] = third_grades_win_len
    window_grades_single_chain['third_dirs'] = third_dir
    window_grades_single_chain['two_thirds_grades'] = two_thirds_grades
    window_grades_single_chain['two_thirds_grades_win_len'] = two_thirds_grades_win_len
    window_grades_single_chain['two_thirds_dirs'] = two_thirds_dir
    window_grades_single_chain['starters'] = starters
    return window_grades_single_chain


def PymolMark(name, minima_tuples, sec_tuples=False):
    print '\n\nLocal (primary) minimas are ', minima_tuples
    print 'Local (secondary) minimas are ', sec_tuples, '\n\n'
    file = '/Users/jonathan/Desktop/test.pml'
    with open(file, 'wa+') as f:
        f.writelines('load ' + name + '.pdb,' + name + '\n')
        f.writelines('cmd.show_as("cartoon", "all")\n')
        # f.writelines('fetch ' + name + ',' + name + '\n')
        i = 1
        for TM in minima_tuples:
            f.writelines('select TM' + str(i) + ', ' + name + ' and resi ' + str(TM[0]) + '-' + str(TM[0]+TM[1]) + '\n')
            f.writelines('color red, TM' + str(i) + '\n')
            i += 1
        i = 1
        if sec_tuples:
            for TM in sec_tuples:
                f.writelines('select TM_S' + str(i) + ', ' + name + ' and resi ' + str(TM[0]) + '-' + str(TM[0]+TM[1]) + '\n')
                f.writelines('color blue, TM_S' + str(i) + '\n')
                i += 1
    subprocess.call(['/opt/local/bin/pymol', '-q', file, '&'])


def MinimaTuples(min_array, wins):
    result = []
    for i, val in enumerate(min_array.T):
        if not np.all(np.isnan(min_array[:, i])):
            result.append([i, int(wins[i])])
    return result


def SmoothenCurve(grd, win_len):
    grd_smooth = []
    win_len_smooth = []
    for tri in range(0, len(grd)):
        grd_smooth.append(np.mean(grd[tri:tri+SMOOTH_SIZE]))
        win_len_smooth.append(max(win_len[tri:tri+SMOOTH_SIZE]))
    return grd_smooth, win_len_smooth


def SmoothenCurveWithDirs(grd, win_len, dirs):
    grd_smooth_pos = []
    grd_smooth_inv = []
    win_len_smooth = []
    half = int(math.floor(SMOOTH_SIZE / 2))

    for i in range(0, half):
        if dirs[i] == 1:
            grd_smooth_pos.append(grd[i])
            grd_smooth_inv.append(1000)
        else:
            grd_smooth_inv.append(grd[i])
            grd_smooth_pos.append(1000)
        win_len_smooth.append(win_len[i])

    for start in range(half, len(grd)-half):
        if len(set(dirs[start:start+SMOOTH_SIZE])) == 1:
            if dirs[start] == 1:
                grd_smooth_pos.append(np.mean(grd[start:start+SMOOTH_SIZE]))
                grd_smooth_inv.append(1000)
            else:
                grd_smooth_inv.append(np.mean(grd[start:start+SMOOTH_SIZE]))
                grd_smooth_pos.append(1000)
            win_len_smooth.append(max(win_len[start:start+SMOOTH_SIZE]))
            # grd_smooth.append(np.mean(grd[start:start+SMOOTH_SIZE]))

        else:
            if dirs[start] == 1:
                grd_smooth_pos.append(grd[start])
                grd_smooth_inv.append(1000)
            else:
                grd_smooth_inv.append(grd[start])
                grd_smooth_pos.append(1000)
            win_len_smooth.append(win_len[start])

    for i in range(len(grd)-half, len(grd)):
        if dirs[i] == 1:
            grd_smooth_pos.append(grd[i])
            grd_smooth_inv.append(1000)
        else:
            grd_smooth_inv.append(grd[i])
            grd_smooth_pos.append(1000)
        win_len_smooth.append(win_len[i])
    # win_len_smooth.append(win_len[0:half])
    # grd_smooth.append(grd[-half:])
    return grd_smooth_pos, grd_smooth_inv, win_len_smooth


def MatrixByVecIndex(array, vec):
    result = []
    for i in range(0, len(array[0])):
        result.append(array[vec[i]][i])
    return result


def RecursiveMinima(grds, wins, point):
    if point+wins[point] >= len(grds[0]):
        min_end_vec = np.min(grds[:, point:len(grds[0])], axis=0)
        min_end_ind = np.argmin(min_end_vec)
        return min_end_ind + point
    if point >= len(wins):
        return point
    point_minima = np.min(grds[:, point:int(point+wins[point])])
    if point_minima == np.amin(grds[:, point]):
        return point
    elif point_minima < np.amin(grds[:, point]):
        next_point = int(np.argmin(grds[:, point:point+wins[point]]) % wins[point] + point)
        if next_point >= len(wins):
            return point
        point = RecursiveMinima(grds, wins, int(np.argmin(grds[:, point:point+wins[point]]) % wins[point] + point))
        return point


def LocalMinima(grds, wins):
    minimas = np.empty(shape=(6, len(grds[0])))
    minimas[:] = np.NAN
    start = 0
    while True:
        point = RecursiveMinima(grds, wins, start)
        if point > len(wins):
            break
        minimas[np.argmin(grds[:, point])][point] = grds[np.argmin(grds[:, point])][point]
        start = int(point+wins[point])
        if start >= len(grds[0]):
            break
    minimas[minimas > PRIMARY_MINIMA_THRESHOLD] = np.NAN
    return minimas


def SecondaryMinimas(grds, minimas, wins):
    # finds minimas that were not found by LocalMinima because they were 'blocked' by 'better' minimas
    # returns a 6 X len(wins) array of NANs everywhere but the secondary minimas, where the grade is found.
    sec_minimas = np.empty(shape=(6, len(grds[0])))
    sec_minimas[:] = np.NAN
    bound = []
    # collects all positions that are 'blocked' by primary minimas
    blocked = []
    for minima in minimas:
        [blocked.append(int(x)) for x in range(minima[0], minima[0]+minima[1]+1)]
    point = 0
    # iterates over all positions, collects non-blocked ones into bound, which is empty after the next blocked block
    # is reached.
    while point < len(grds[0]):
        if not any(x in blocked for x in range(point, point+int(wins[point]))):
            bound.append(point)
        else:
            if len(bound) != 0:
                if len(bound) == 1:
                    bound_min = np.min(grds[:, bound])
                else:
                    bound_min = np.min(grds[:, bound[0]:bound[-1]])
                if bound_min <= SECONDARY_MINIMA_THRESHOLD:
                    min_ind = np.where(grds == bound_min)
                    sec_minimas[min_ind[0], min_ind[1]] = bound_min
                bound = []
                for minima in minimas:
                    # jumps point to the end of the next blocked block
                    if minima[0] <= point+wins[point] < minima[0]+minima[1]:
                        point = minima[0]+minima[1]
        point += 1
    return sec_minimas


def PlotSinglePeptideWindows(grades):
    grd_smh = np.empty(shape=(6, len(grades['full_grades'])))
    win_smh = np.empty(shape=(6, len(grades['full_grades'])))
    (grd_smh[0], grd_smh[1], win_smh[0]) = SmoothenCurveWithDirs(grades['full_grades'], grades['full_grades_win_len'], grades['full_dirs'])
    win_smh[1] = win_smh[0]
    (grd_smh[2], grd_smh[3], win_smh[2]) = SmoothenCurveWithDirs(grades['third_grades'], grades['third_grades_win_len'], grades['third_dirs'])
    win_smh[3] = win_smh[2]
    (grd_smh[4], grd_smh[5], win_smh[4]) = SmoothenCurveWithDirs(grades['two_thirds_grades'], grades['two_thirds_grades_win_len'], grades['two_thirds_dirs'])
    win_smh[5] = win_smh[4]
    # grd_min = np.amin(grd_smh, axis=0)
    min_ind = np.argmin(grd_smh, axis=0)
    win_of_min = MatrixByVecIndex(win_smh, min_ind)

    minimas = LocalMinima(grd_smh, win_of_min)
    # np.set_printoptions(threshold=np.inf)

    minima_tuples = MinimaTuples(minimas, win_of_min)
    minimas_secondary = SecondaryMinimas(grd_smh, minima_tuples, win_of_min)
    sec_minima_tuples = MinimaTuples(minimas_secondary, win_of_min)
    PymolMark(grades['name'], minima_tuples, sec_minima_tuples)

    plot_array = np.empty(shape=[6, len(grades['starters'])])
    plot_array[:] = np.NAN
    for col, row in enumerate(min_ind):
        plot_array[row][col] = grd_smh[row][col]

    PlotParams()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    f = ax1.scatter(grades['starters'], plot_array[0], s=50, c='k')
    fi = ax1.scatter(grades['starters'], plot_array[1], s=50, facecolors='none', edgecolors='k')
    t = ax1.scatter(grades['starters'], plot_array[2], s=50, c='r')
    ti = ax1.scatter(grades['starters'], plot_array[3], s=50, facecolors='none', edgecolors='r')
    tw = ax1.scatter(grades['starters'], plot_array[4], s=50, c='b')
    twi = ax1.scatter(grades['starters'], plot_array[5], s=50, facecolors='none', edgecolors='b')

    fm = ax1.scatter(grades['starters'], minimas[0], s=120, c='k', marker='D')
    fmi = ax1.scatter(grades['starters'], minimas[1], s=120, c='k', marker='x')
    tm = ax1.scatter(grades['starters'], minimas[2], s=120, c='r', marker='D')
    tmi = ax1.scatter(grades['starters'], minimas[3], s=120, c='r', marker='x')
    twm = ax1.scatter(grades['starters'], minimas[4], s=120, c='b', marker='D')
    twmi = ax1.scatter(grades['starters'], minimas[5], s=120, c='b', marker='x')

    fsm = ax1.scatter(grades['starters'], minimas_secondary[0], s=200, c='k', marker='8')
    fsmi = ax1.scatter(grades['starters'], minimas_secondary[1], s=200, c='k', marker='s')
    tsm = ax1.scatter(grades['starters'], minimas_secondary[2], s=200, c='r', marker='8')
    tsmi = ax1.scatter(grades['starters'], minimas_secondary[3], s=200, c='r', marker='s')
    twsm = ax1.scatter(grades['starters'], minimas_secondary[4], s=200, c='b', marker='8')
    twsmi = ax1.scatter(grades['starters'], minimas_secondary[5], s=200, c='b', marker='s')

    ax1.plot(grades['starters'], [0]*len(grades['starters']), 'k--')
    for minima in minima_tuples:
        ax1.plot((minima[0], minima[0]), (-50, 50), 'k--')
        ax1.plot((minima[0]+minima[1], minima[0]+minima[1]), (-50, 50), 'r--')
    for minima in sec_minima_tuples:
        ax1.plot((minima[0], minima[0]), (-50, 50), 'b--')
        ax1.plot((minima[0]+minima[1], minima[0]+minima[1]), (-50, 50), 'g--')
    plt.ylim((-12, 20))
    # ax2 = ax1.twinx()
    # ax2.scatter(grades['starters'], win_of_min, 32, c='b', marker='.')
    plt.xlim((-5, len(grades['starters'])+5))
    # plt.ylim(np.amin(plot_array)-2, 15)
    plt.xticks(range(0, len(grades['starters']), 20), rotation=90)
    ax1.set_xlabel('Window Start Residue', fontsize=36)
    ax1.set_ylabel('$\Delta$G$_{transfer}$', fontsize=36)
    # ax2.set_ylabel('Window Length')
    ax1.legend((f, fi, t, ti, tw, twi, fm, fmi, tm, tmi, twm, twmi, fsm, fsmi, tsm, tsmi, twsm, twsmi),
               ('Full', 'Full -1', 'Third', 'Third -1', 'Two-Thirds', 'Two-Thirds -1', 'Full minima', 'Full -1 minima',
                'Third minima', 'Third -1 minima', 'Two Thirds minima', 'Two Thirds -1 minima', 'Full SM', 'Full -1 SM',
                'Third SM', 'Third -1 SM', 'Two Thirds SM', 'Two Thirds -1 SM'), scatterpoints=1, ncol=9,
               loc='lower center', prop={'size': 12})
    plt.suptitle(grades['name']+' plot with smooth '+str(SMOOTH_SIZE)+' and loop-bypass '+str(LOOP_BYPASS), fontsize=48)
    plt.show()


def PlotParams():
    # by: http://damon-is-a-geek.com/publication-ready-the-first-time-beautiful-reproducible-plots-with-matplotlib.html
    from matplotlib import rcParams
    rcParams['axes.labelsize'] = 9
    rcParams['xtick.labelsize'] = 9
    rcParams['ytick.labelsize'] = 9
    rcParams['legend.fontsize'] = 9

    # os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
    rcParams['font.family'] = 'serif'
    rcParams['font.serif'] = ['Computer Modern Roman']
    # rcParams['text.usetex'] = True


### MAIN ###
MakeHydrophobicityGrade()
# UniprotParser()
# CalculateWindowGrades()
# WriteDataSetToTSV('/Users/jonathan/Desktop/python_temp.txt')
# ss_grades = WindowGradesForSingleSequence('MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKN', '4j4q')
# ss_grades = WindowGradesForSingleSequence('TGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFSMRPEVATFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG', '1c3w')
# ss_grades = WindowGradesForSingleSequence('MELEEDLKGRADKNFSKMGKKSKKEKKEKKPAVSVLTMFRYAGWLDRLYMLVGTLAAIIHGVALPLMMLIFGDMTDSFASVGNVSKNSTNMSEADKRAMFAKLEEEMTTYAYYYTGIGAGVLIVAYIQVSFWCLAAGRQIHKIRQKFFHAIMNQEIGWFDVHDVGELNTRLTDDVSKINEGIGDKIGMFFQAMATFFGGFIIGFTRGWKLTLVILAISPVLGLSAGIWAKILSSFTDKELHAYAKAGAVAEEVLAAIRTVIAFGGQKKELERYNNNLEEAKRLGIKKAITANISMGAAFLLIYASYALAFWYGTSLVISKEYSIGQVLTVFFSVLIGAFSVGQASPNIEAFANARGAAYEVFKIIDNKPSIDSFSKSGHKPDNIQGNLEFKNIHFSYPSRKEVQILKGLNLKVKSGQTVALVGNSGCGKSTTVQLMQRLYDPLDGMVSIDGQDIRTINVRYLREIIGVVSQEPVLFATTIAENIRYGREDVTMDEIEKAVKEANAYDFIMKLPHQFDTLVGERGAQLSGGQKQRIAIARALVRNPKILLLDEATSALDTESEAVVQAALDKAREGRTTIVIAHRLSTVRNADVIAGFDGGVIVEQGNHDELMREKGIYFKLVMTQTAGNEIELGNEACKSKDEIDNLDMSSKDSGSSLIRRRSTRKSICGPHDQDRKLSTKEALDEDVPPASFWRILKLNSTEWPYFVVGIFCAIINGGLQPAFSVIFSKVVGVFTNGGPPETQRQNSNLFSLLFLILGIISFITFFLQGFTFGKAGEILTKRLRYMVFKSMLRQDVSWFDDPKNTTGALTTRLANDAAQVKGATGSRLAVIFQNIANLGTGIIISLIYGWQLTLLLLAIVPIIAIAGVVEMKMLSGQALKDKKELEGSGKIATEAIENFRTVVSLTREQKFETMYAQSLQIPYRNAMKKAHVFGITFSFTQAMMYFSYAACFRFGAYLVTQQLMTFENVLLVFSAIVFGAMAVGQVSSFAPDYAKATVSASHIIRIIEKTPEIDSYSTQGLKPNMLEGNVQFSGVVFNYPTRPSIPVLQGLSLEVKKGQTLALVGSSGCGKSTVVQLLERFYDPMAGSVFLDGKEIKQLNVQWLRAQLGIVSQEPILFDCSIAENIAYGDNSRVVSYEEIVRAAKEANIHQFIDSLPDKYNTRVGDKGTQLSGGQKQRIAIARALVRQPHILLLDEATSALDTESEKVVQEALDKAREGRTCIVIAHRLSTIQNADLIVVIQNGKVKEHGTHQQLLAQKGIYFSMVSVQAGAKRSYVHH', '3g61')
# ss_grades = WindowGradesForSingleSequence('GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRA', '1BRX')
# ss_grades = WindowGradesForSingleSequence('MDMLFAKTVVLAASAVGAGTAMIAGIGPGVGQGYAAGKAVESVARQPEAKGDIISTMVLGQAVAESTGIYSLVIALILLYANPFVGLLG', '1yce')
ss_grades = WindowGradesForSingleSequence('YAQALQSVPETQVSQLDNGVRVASEQSSQPTCTVGVWIDAGSRYESEKNNGAGYFLEHLAFKGTKNRPQNALEKEVESMGAHLNAYSSREHTAYYIKALSKDVPKAVELLADIVQNCSLEDSQIEKERDVIVRELQENDTSMREVVFNYLHATAFQGTGLAQSVEGPSENIRKLSRADLTEYLSTHYTAPRMVLAAAGGVEHQQLLELAQKHFGGVPFTYDDDAVPTLSKCRFTGSQIRHREDGLPLAHVAIAVEGPGWAHPDLVALQVANAIIGHYDRTYGGGLHSSSPLASIAVTNKLCQSFQTFSICYSETGLFGFYFVCDRMSIDDMMFVLQGQWMRLCTSISESEVLRGKNFLRNALVSHLDGTTPVCEDIGRELLTYGRRIPLEEWEERLAEVDARMVREVCSKYIYDQCPAVAGPGPIEQLPDYNRIRSGMFWLRPPHPQDLEITKLPNGLVIASLENYSPGSTIGVFIKAGSRYENSSNLGTSHLLRLASSLTTKGASSFKITRGIEAVGGKLSVESTRENMAYTVECLRDDVEILMEFLLNVTTAPEFRPWEVADLQPQLKIDKAVAFQNPQTHVIENLHAAAYRNALADSLYCPDYRIGKVTSVELHDFVQNHFTSARMALVGLGVSHPVLKNVAEQLLNIRGGLGLSGAKAKYRGGEIREQNGDSLVHAAIVAESAAIGGAEANAFSVLQHVLGANPHVKRGNPFDVSAFNASYSDSGLFGFYTISQAAYAGQVIKAAYNQVKTIAQGNVSNENVQAAKNKLKAKYLMSVESSEGFLEEVGSQALAAGSYNPPSTVLQQIDAVADADVIKAAKKFVSRQKSMAASGNLGHTPFVDELAPNIRKSHPLLKMINNSLIDLPAPSNISAWWNFGSLLAVCLMTQILTGLLLAMHYTADTSLAFSSVAHTCRNVQYGWLIRNLHANGASFFFICIFLHIGRGLYYGSYLYKETWNTGVILLLTLMATAFVGYVLPWGQMSFWGATVITNLFSAIPYIGHTLVEWAWGGFSVDNPTLTRFFALHFLLPFAIAGITIIHLTFLHESGSNNPLGISSDSDKIPFHPYYSFKDILGLTLMLTPFLTLALFSPNLLGDPENFTPANPLVTPPHIKPEWYFLFAYAILRSIPNKLGGVLALAASVLILFLIPFLHKSKQRTMTFRPLSQTLFWLLVANLLILTWIGSQPVEHPFIIIGQMASLSYFTILLILFPTIGTLENKMLNYSDLELHPPSYPWSHRGPLSSLDHTSIRRGFQVYKQVCSSCHSMDYVAYRHLVGVCYTEDEAKALAEEVEVQDGPNEDGEMFMRPGKLSDYFPKPYPNPEAARAANNGALPPDLSYIVRARHGGEDYVFSLLTGYCEPPTGVSVREGLYFNPYFPGQAIGMAPPIYNDVLEFDDGTPATMSQVAKDVCTFLRWAAEPEHDHRKRMGLKMLLMMGLLVPLVYYMKRHKWSVLKSRKLAYRPPKSHTDIKVPNFSDYRRPPDDYSTKSSRESDPSRKGFSYLVTAVTTLGVAYAAKNVVTQFVSSMSASADVLAMSKIEIKLSDIPEGKNMAFKWRGKPLFVRHRTKKEIDQEAAVEVSQLRDPQHDLERVKKPEWVILIGVCTHLGCVPIANAGDFGGYYCPCHGSHYDASGRIRKGPAPLNLEVPSYEFTSDDMVIVGSRWLEGIRKWYYNAAGFNKYGLMRDDTIYENDDVKEAIRRLPENLYDDRMFRIKRALDLNMRQQILPKEQWTKYEEDVPYLEPYLKEVIRERKEREEWDKRQFGHLTRVRHLITYSLSPFEQRPFPHYFSKGVPNVWRRLRACILRVAPPFLAFYLLYTWGTQEFEKSKRKNPAAYVNLVDPLTTVREQCEQLEKCVKARERLELCDERVSSRSQTEEDCTEELFDFLHARDHCVAHKLFNSLKTLTARLYSLLFRRTSTFALTIVVGALLFERAFDQGADAIYEHINEGKLWKHIKHKYENK', '1bcc')
# WriteSingleSeqToTSV(ss_grades, '/Users/jonathan/Desktop/bassaf_seq.txt')
# WriteSingleSeqToTSVusingCSV(ss_grades, '/Users/jonathan/Desktop/4j4q_win_gds.csv')

PlotSinglePeptideWindows(ss_grades)


# array = window_grades['YP_006514311.1']['full_grades']
# starters = window_grades['YP_006514311.1']['starters']
# plt.scatter(starters, array)
# print main_dict['YP_006514311.1']['transmem']
# plt.show()