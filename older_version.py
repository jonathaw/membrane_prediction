import re
import numpy as np
import matplotlib.pyplot as plt
import csv
import math


# MEMEBRANE_RANGE = range(-15, 15, 30 / 21)
MEMBRANE_FULL = range(0, 35)
MEMBRANE_THIRD = [0, 1, 4, 8, 11, 15, 18, 19, 22, 26, 29, 33, 36]
MEMBRANE_TWO_THIRDS = [0, 1, 2, 4, 5, 8, 9, 11, 12, 13, 15, 16, 18, 19, 20, 22, 23, 26, 27, 29, 30, 31, 33, 34, 36]
hydrophobicity_polyval = {}
window_grades = {}
uniprot_entris = 0
main_dict = {}
SMOOTH_SIZE = 1


def MakeHydrophobicityGrade():
    global hydrophobicity_polyval
    hydrophobicity_grade = open('/Users/jonathan/eden/membrane_prediciton/membrane_prediction/Poly_Values.txt', 'r')
    for line in hydrophobicity_grade:
        split = line.split(' ')
        hydrophobicity_polyval[split[0]] = [float(n) for n in split[1:6]]
    hydrophobicity_grade.close()


def grade_seq(seq, mode):
    sumer = 0
    not_span = 0
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
    if not_span >= 3:
        sumer += not_span
    return sumer


# def grade_seq(seq, mode):
#     sumer = 0
#     seq_inverse = seq[::-1]
#     sumer_inverse = 0
#     not_span = 0
#     not_span_inverse = 0
#     membrane_position = np.linspace(-15, 15, endpoint=True, num=len(seq))
#     if mode == 'full':
#         membrane_span = MEMBRANE_FULL
#     elif mode == 'third':
#         membrane_span = MEMBRANE_THIRD
#     elif mode == 'two_thirds':
#         membrane_span = MEMBRANE_TWO_THIRDS
#     for i, aa in enumerate(seq):
#         if i in membrane_span:
#             sumer += np.polyval(hydrophobicity_polyval[aa], membrane_position[i])
#         else:
#             not_span += np.polyval(hydrophobicity_polyval[aa], membrane_position[i])
#     for i, aa in enumerate(seq_inverse):
#         if i in membrane_span:
#             sumer_inverse += np.polyval(hydrophobicity_polyval[aa], membrane_position[i])
#         else:
#             not_span_inverse += np.polyval(hydrophobicity_polyval[aa], membrane_position[i])
#     if not_span >= 3:
#         sumer += not_span
#     if not_span_inverse >= 3:
#         sumer_inverse += not_span_inverse
#     if sumer_inverse > sumer:
#         return sumer, 1
#     else:
#         return sumer_inverse, 0


# def hydrophobicity_grade_increments(seq, mode):
#     results = []
#     dirs = []
#     temp = grade_seq(seq[0:22], mode)
#     results.append(temp[0])
#     dirs.append(temp[1])
#     for inc in range(1, 14):
#         if inc+22 <= len(seq):
#             temp = grade_seq(seq[0:22+inc], mode)
#             results.append(temp[0])
#             dirs.append(temp[1])
#     min_val = min(results)
#     min_val_index = results.index(min_val)
#     return min_val, min_val_index+22, dirs[min_val_index], results


def hydrophobicity_grade_increments(seq, mode):
    results = []
    results.append(grade_seq(seq[0:22], mode))
    for inc in range(1, 14):
        if inc+22 <= len(seq):
            results.append(grade_seq(seq[0:22+inc], mode))
    min_val = min(results)
    min_val_index = results.index(min_val)
    return min_val, min_val_index+22, results


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


# def CalculateWindowGrades():
#     ### go over all 22AA windows in all protein sequecnes in main_dict, and calculate relevant
#     ### hydrophobicity grades
#     global window_grades
#     for protein in main_dict:
#         full_grades = []
#         third_grades = []
#         two_thirds_grades = []
#         full_grades_inverse = []
#         third_grades_inverse = []
#         two_thirds_grades_inverse = []
#         starters = []
#         window_grades[protein] = {}
#         for first in range(0, len(main_dict[protein]['sequence']) - 22):
#             full_grades.append(hydrophobicity_grade_full(main_dict[protein]['sequence'][first:first+22]))
#             full_grades_inverse.append(hydrophobicity_grade_full(main_dict[protein]['sequence'][first:first+22][::-1]))
#             third_grades.append(hydrophobicity_grade_third(main_dict[protein]['sequence'][first:first+22]))
#             third_grades_inverse.append(hydrophobicity_grade_third(main_dict[protein]['sequence'][first:first+22][::-1]))
#             two_thirds_grades.append(hydrophobicity_grade_two_third((main_dict[protein]['sequence'][first:first+22])))
#             two_thirds_grades_inverse.append(hydrophobicity_grade_two_third((main_dict[protein]['sequence'][first:first+22][::-1])))
#             starters.append(first)
#         window_grades[protein]['full_grades'] = full_grades
#         window_grades[protein]['full_grades_inverse'] = full_grades_inverse
#         window_grades[protein]['third_grades'] = third_grades
#         window_grades[protein]['third_grades_inverse'] = third_grades_inverse
#         window_grades[protein]['two_thirds_grades'] = two_thirds_grades
#         window_grades[protein]['two_thirds_grades_inverse'] = two_thirds_grades_inverse
#         window_grades[protein]['starters'] = starters


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


# def WindowGradesForSingleSequence(seq):
#     window_grades_single_chain = {}
#     full_grades = []
#     full_grades_win_len = []
#     full_dir = []
#     third_grades = []
#     third_grades_win_len = []
#     third_dir = []
#     two_thirds_grades = []
#     two_thirds_grades_win_len = []
#     two_thirds_dir = []
#     starters = []
#     for first in range(0, len(seq) - 22):
#         temp = hydrophobicity_grade_increments(seq[first:first+36], 'full')
#         full_grades.append(temp[0])
#         full_grades_win_len.append(temp[1])
#         full_dir.append(temp[2])
#         temp = hydrophobicity_grade_increments(seq[first:first+36], 'third')
#         third_grades.append(temp[0])
#         third_grades_win_len.append(temp[1])
#         third_dir.append(temp[2])
#         temp = hydrophobicity_grade_increments(seq[first:first+36], 'two_thirds')
#         two_thirds_grades.append(temp[0])
#         two_thirds_grades_win_len.append(temp[1])
#         two_thirds_dir.append(temp[2])
#         starters.append(first)
#     window_grades_single_chain['full_grades'] = full_grades
#     window_grades_single_chain['full_grades_win_len'] = full_grades_win_len
#     window_grades_single_chain['full_dirs'] = full_dir
#     window_grades_single_chain['third_grades'] = third_grades
#     window_grades_single_chain['third_grades_win_len'] = third_grades_win_len
#     window_grades_single_chain['third_dirs'] = third_dir
#     window_grades_single_chain['two_thirds_grades'] = two_thirds_grades
#     window_grades_single_chain['two_thirds_grades_win_len'] = two_thirds_grades_win_len
#     window_grades_single_chain['two_thirds_dirs'] = two_thirds_dir
#     window_grades_single_chain['starters'] = starters
#     return window_grades_single_chain


def WindowGradesForSingleSequence(seq):
    window_grades_single_chain = {}
    full_grades = []
    full_grades_win_len = []
    third_grades = []
    third_grades_win_len = []
    two_thirds_grades = []
    two_thirds_grades_win_len = []
    full_grades_inverse = []
    full_grades_inverse_win_len = []
    third_grades_inverse = []
    third_grades_inverse_win_len = []
    two_thirds_grades_inverse = []
    two_thirds_grades_inverse_win_len = []
    starters = []
    seq_inverse = seq[::-1]
    for first in range(0, len(seq) - 22):
        temp = hydrophobicity_grade_increments(seq[first:first+36], 'full')
        full_grades.append(temp[0])
        full_grades_win_len.append(temp[1])
        temp = hydrophobicity_grade_increments(seq_inverse[first:first+36], 'full')
        full_grades_inverse.append(temp[0])
        full_grades_inverse_win_len.append(temp[1])
        temp = hydrophobicity_grade_increments(seq[first:first+36], 'third')
        third_grades.append(temp[0])
        third_grades_win_len.append(temp[1])
        temp = hydrophobicity_grade_increments(seq_inverse[first:first+36], 'third')
        third_grades_inverse.append(temp[0])
        third_grades_inverse_win_len.append(temp[1])
        temp = hydrophobicity_grade_increments(seq[first:first+36], 'two_thirds')
        two_thirds_grades.append(temp[0])
        two_thirds_grades_win_len.append(temp[1])
        temp = hydrophobicity_grade_increments(seq_inverse[first:first+36], 'two_thirds')
        two_thirds_grades_inverse.append(temp[0])
        two_thirds_grades_inverse_win_len.append(temp[1])
        starters.append(first)
    window_grades_single_chain['full_grades'] = full_grades
    window_grades_single_chain['full_grades_win_len'] = full_grades_win_len
    window_grades_single_chain['full_grades_inverse'] = full_grades_inverse[::-1]
    window_grades_single_chain['full_grades_inverse_win_len'] = full_grades_inverse_win_len[::-1]
    window_grades_single_chain['third_grades'] = third_grades
    window_grades_single_chain['third_grades_win_len'] = third_grades_win_len
    window_grades_single_chain['third_grades_inverse'] = third_grades_inverse[::-1]
    window_grades_single_chain['third_grades_inverse_win_len'] = third_grades_inverse_win_len
    window_grades_single_chain['two_thirds_grades'] = two_thirds_grades
    window_grades_single_chain['two_thirds_grades_win_len'] = two_thirds_grades_win_len
    window_grades_single_chain['two_thirds_grades_inverse'] = two_thirds_grades_inverse[::-1]
    window_grades_single_chain['two_thirds_grades_inverse_win_len'] = two_thirds_grades_inverse_win_len
    window_grades_single_chain['starters'] = starters
    return window_grades_single_chain


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

### make this function use the new smooth, and plot by the new system of grades
# def PlotSinglePeptideWindows(grades):
#     grd_smh = np.empty(shape=[6, len(grades['full_grades'])])
#     win_smh = np.empty(shape=[6, len(grades['full_grades'])])
#     (grd_smh[0], grd_smh[1], win_smh[0]) = SmoothenCurveWithDirs(grades['full_grades'], grades['full_grades_win_len'], grades['full_dirs'])
#     (grd_smh[2], grd_smh[3], win_smh[2]) = SmoothenCurveWithDirs(grades['third_grades'], grades['third_grades_win_len'], grades['third_dirs'])
#     win_smh[3] = win_smh[2]
#     (grd_smh[4], grd_smh[5], win_smh[3]) = SmoothenCurveWithDirs(grades['two_thirds_grades'], grades['two_thirds_grades_win_len'], grades['two_thirds_dirs'])
#     win_smh[5] = win_smh[4]
#     # grd_min = np.amin(grd_smh, axis=0)
#     min_ind = np.argmin(grd_smh, axis=0)
#     win_of_min = MatrixByVecIndex(win_smh, min_ind)
#
#     plot_array = np.empty(shape=[6, len(grades['starters'])])
#     plot_array[:] = np.NAN
#     for col, row in enumerate(min_ind):
#         plot_array[row][col] = grd_smh[row][col]
#
#     fig = plt.figure()
#     ax1 = fig.add_subplot(111)
#     f   = ax1.scatter(grades['starters'], plot_array[0], s=50, c='k')
#     fi  = ax1.scatter(grades['starters'], plot_array[1], s=50, facecolors='none', edgecolors='k')
#     t   = ax1.scatter(grades['starters'], plot_array[2], s=50, c='r')
#     ti  = ax1.scatter(grades['starters'], plot_array[3], s=50, facecolors='none', edgecolors='r')
#     tw  = ax1.scatter(grades['starters'], plot_array[4], s=50, c='b')
#     twi = ax1.scatter(grades['starters'], plot_array[5], s=50, facecolors='none', edgecolors='b')
#     ax2 = ax1.twinx()
#     ax2.scatter(grades['starters'], win_of_min, 32, c='b', marker='.')
#     plt.xlim((-5, len(grades['starters'])+5))
#     ax1.set_xlabel('Window')
#     ax1.set_ylabel('ddG')
#     ax2.set_ylabel('Window Length')
#     ax1.legend((f, fi, t, ti, tw, twi), ('Full', 'Full -1', 'Third', 'Third -1', 'Two-Thirds', 'Two-Thirds -1'),
#                scatterpoints=1, ncol=6, loc='lower left')
#     plt.show()


def PlotSinglePeptideWindows(grades):
    grd_smh = np.empty(shape=[6, len(grades['full_grades'])])
    win_smh = np.empty(shape=[6, len(grades['full_grades'])])
    (grd_smh[0], win_smh[0]) = SmoothenCurve(grades['full_grades'], grades['full_grades_win_len'])
    (grd_smh[1], win_smh[1]) = SmoothenCurve(grades['full_grades_inverse'], grades['full_grades_inverse_win_len'])
    (grd_smh[2], win_smh[2]) = SmoothenCurve(grades['third_grades'], grades['third_grades_win_len'])
    (grd_smh[3], win_smh[3]) = SmoothenCurve(grades['third_grades_inverse'], grades['third_grades_inverse_win_len'])
    (grd_smh[4], win_smh[4]) = SmoothenCurve(grades['two_thirds_grades'], grades['two_thirds_grades_win_len'])
    (grd_smh[5], win_smh[5]) = SmoothenCurve(grades['two_thirds_grades_inverse'], grades['two_thirds_grades_inverse_win_len'])
    grd_min = np.amin(grd_smh, axis=0)
    min_ind = np.argmin(grd_smh, axis=0)
    win_of_min = MatrixByVecIndex(win_smh, min_ind)

    plot_array = np.empty(shape=[6, len(grades['starters'])])
    plot_array[:] = np.NAN
    for col, row in enumerate(min_ind):
        plot_array[row][col] = grd_smh[row][col]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    f   = ax1.scatter(grades['starters'], plot_array[0], s=50, c='k')
    fi  = ax1.scatter(grades['starters'], plot_array[1], s=50, facecolors='none', edgecolors='k')
    t   = ax1.scatter(grades['starters'], plot_array[2], s=50, c='r')
    ti  = ax1.scatter(grades['starters'], plot_array[3], s=50, facecolors='none', edgecolors='r')
    tw  = ax1.scatter(grades['starters'], plot_array[4], s=50, c='b')
    twi = ax1.scatter(grades['starters'], plot_array[5], s=50, facecolors='none', edgecolors='b')
    ax2 = ax1.twinx()
    ax2.scatter(grades['starters'], win_of_min, 32, c='b', marker='.')
    plt.xlim((-5, len(grades['starters'])+5))
    ax1.set_xlabel('Window')
    ax1.set_ylabel('ddG')
    ax2.set_ylabel('Window Length')
    ax1.legend((f, fi, t, ti, tw, twi), ('Full', 'Full -1', 'Third', 'Third -1', 'Two-Thirds', 'Two-Thirds -1'),
               scatterpoints=1, ncol=6, loc='lower left')
    plt.show()


### MAIN ###
MakeHydrophobicityGrade()
# UniprotParser()
# CalculateWindowGrades()
# WriteDataSetToTSV('/Users/jonathan/Desktop/python_temp.txt')
ss_grades = WindowGradesForSingleSequence('MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA')
# ss_grades = WindowGradesForSingleSequence('TGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFSMRPEVATFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG')
# ss_grades = WindowGradesForSingleSequence('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
# ss_grades = WindowGradesForSingleSequence('MELEEDLKGRADKNFSKMGKKSKKEKKEKKPAVSVLTMFRYAGWLDRLYMLVGTLAAIIHGVALPLMMLIFGDMTDSFASVGNVSKNSTNMSEADKRAMFAKLEEEMTTYAYYYTGIGAGVLIVAYIQVSFWCLAAGRQIHKIRQKFFHAIMNQEIGWFDVHDVGELNTRLTDDVSKINEGIGDKIGMFFQAMATFFGGFIIGFTRGWKLTLVILAISPVLGLSAGIWAKILSSFTDKELHAYAKAGAVAEEVLAAIRTVIAFGGQKKELERYNNNLEEAKRLGIKKAITANISMGAAFLLIYASYALAFWYGTSLVISKEYSIGQVLTVFFSVLIGAFSVGQASPNIEAFANARGAAYEVFKIIDNKPSIDSFSKSGHKPDNIQGNLEFKNIHFSYPSRKEVQILKGLNLKVKSGQTVALVGNSGCGKSTTVQLMQRLYDPLDGMVSIDGQDIRTINVRYLREIIGVVSQEPVLFATTIAENIRYGREDVTMDEIEKAVKEANAYDFIMKLPHQFDTLVGERGAQLSGGQKQRIAIARALVRNPKILLLDEATSALDTESEAVVQAALDKAREGRTTIVIAHRLSTVRNADVIAGFDGGVIVEQGNHDELMREKGIYFKLVMTQTAGNEIELGNEACKSKDEIDNLDMSSKDSGSSLIRRRSTRKSICGPHDQDRKLSTKEALDEDVPPASFWRILKLNSTEWPYFVVGIFCAIINGGLQPAFSVIFSKVVGVFTNGGPPETQRQNSNLFSLLFLILGIISFITFFLQGFTFGKAGEILTKRLRYMVFKSMLRQDVSWFDDPKNTTGALTTRLANDAAQVKGATGSRLAVIFQNIANLGTGIIISLIYGWQLTLLLLAIVPIIAIAGVVEMKMLSGQALKDKKELEGSGKIATEAIENFRTVVSLTREQKFETMYAQSLQIPYRNAMKKAHVFGITFSFTQAMMYFSYAACFRFGAYLVTQQLMTFENVLLVFSAIVFGAMAVGQVSSFAPDYAKATVSASHIIRIIEKTPEIDSYSTQGLKPNMLEGNVQFSGVVFNYPTRPSIPVLQGLSLEVKKGQTLALVGSSGCGKSTVVQLLERFYDPMAGSVFLDGKEIKQLNVQWLRAQLGIVSQEPILFDCSIAENIAYGDNSRVVSYEEIVRAAKEANIHQFIDSLPDKYNTRVGDKGTQLSGGQKQRIAIARALVRQPHILLLDEATSALDTESEKVVQEALDKAREGRTCIVIAHRLSTIQNADLIVVIQNGKVKEHGTHQQLLAQKGIYFSMVSVQAGAKRSYVHH')
# WriteSingleSeqToTSV(ss_grades, '/Users/jonathan/Desktop/bassaf_seq.txt')
# WriteSingleSeqToTSVusingCSV(ss_grades, '/Users/jonathan/Desktop/4j4q_win_gds.csv')

PlotSinglePeptideWindows(ss_grades)


# array = window_grades['YP_006514311.1']['full_grades']
# starters = window_grades['YP_006514311.1']['starters']
# plt.scatter(starters, array)
# print main_dict['YP_006514311.1']['transmem']
# plt.show()