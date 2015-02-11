import re
import numpy as np
import matplotlib.pyplot as plt
import csv
import math
import subprocess
from collections import Counter
import database_parser
from other_functions import PsiReaderHelix

# import os

# MEMEBRANE_RANGE = range(-15, 15, 30 / 21)
MEMBRANE_FULL = range(0, 35)
# MEMBRANE_THIRD = [0, 1, 4, 8, 11, 15, 18, 19, 22, 26, 29, 33, 36]
# MEMBRANE_TWO_THIRDS = [0, 1, 2, 4, 5, 8, 9, 11, 12, 13, 15, 16, 18, 19, 20, 22, 23, 26, 27, 29, 30, 31, 33, 34, 36]
# MEMBRANE_THIRD = {'1': [0, 3, 4, 7, 11, 14, 15, 18, 22, 25, 26, 29, 33, 37, 38],
#                   '2': [0, 1, 4, 8, 11, 12, 15, 19, 22, 23, 26, 30, 33, 34, 37, 41],
#                   '3': [0, 4, 7, 8, 11, 15, 18, 19, 22, 26, 29, 30, 33, 37, 40, 41]}
# MEMBRANE_TWO_THIRDS = {'1': [0, 1, 2, 4, 5, 8, 9, 11, 12, 13, 15, 16, 19, 20, 22, 23, 24, 26, 27, 30, 31, 33, 34, 35,
#                              37, 38, 41, 42],
#                        '2': [0, 1, 4, 5, 7, 8, 9, 11, 12, 15, 16, 18, 19, 20, 22, 23, 26, 27, 29, 30, 31, 33, 34, 37,
#                              38, 40, 41, 42],
#                        '3': [0, 1, 3, 4, 5, 7, 8, 11, 12, 14, 15, 16, 18, 19, 22, 23, 25, 26, 27, 29, 30, 33, 34, 36,
#                              37, 38, 40, 41]}
MEMBRANE_THIRD = {'1': range(0, 35), '2': range(0, 35), '3': range(0, 35)}
MEMBRANE_TWO_THIRDS = {'1': range(0, 35), '2': range(0, 35), '3': range(0, 35)}
MEMBRANE_SPANS = {'full': {'1': MEMBRANE_FULL}, 'third': MEMBRANE_THIRD, 'two_thirds': MEMBRANE_TWO_THIRDS}
hydrophobicity_polyval = {}
window_grades = {}
uniprot_entris = 0
main_dict = {}
SMOOTH_SIZE = 1
LOOP_BYPASS = 300
MIN_WIN = 20
SECONDARY_MINIMA_THRESHOLD = 8
PRIMARY_MINIMA_THRESHOLD = 5
PSI_CUTOFF = 0.0
PSI_RES_PREC_CUTOFF = 0.1
HPHOBICITY_PSIPRED_BYPASS = -2 #SECONDARY_MINIMA_THRESHOLD - 2
PSI_EDGE = 4
PSI_CUTOFF_EDGE = 0.6
PSI_EDGE_PREC = 0.25


def MakeHydrophobicityGrade():
    global hydrophobicity_polyval
    # hydrophobicity_grade = open('Poly_Values.txt', 'r')
    hydrophobicity_grade = open('poly_value_11.2.txt', 'r')
    for line in hydrophobicity_grade:
        split = line.split()
        hydrophobicity_polyval[split[0]] = [float(n) for n in split[1:6]]
    hydrophobicity_grade.close()


def grade_seq_single_go(seq):
    seq_inverse = seq[::-1]
    membrane_position = np.linspace(-15, 15, endpoint=True, num=len(seq))
    best_score = 1000
    for mode, type in MEMBRANE_SPANS.iteritems():
        for typex, ranger in type.iteritems():
            sumer = 0
            not_span = 0
            sumer_inverse = 0
            not_span_inverse = 0
            for i, aa in enumerate(seq):
                if i in ranger:
                    sumer += np.polyval(hydrophobicity_polyval[aa], membrane_position[i])
                else:
                    not_span += np.polyval(hydrophobicity_polyval[aa], membrane_position[i])
            for i, aa in enumerate(seq_inverse):
                if i in ranger:
                    sumer_inverse += np.polyval(hydrophobicity_polyval[aa], membrane_position[i])
                else:
                    not_span_inverse += np.polyval(hydrophobicity_polyval[aa], membrane_position[i])
            if not_span >= LOOP_BYPASS:
                sumer += not_span
            if not_span_inverse >= LOOP_BYPASS:
                sumer_inverse += not_span_inverse
            if sumer < best_score:
                best_score = sumer
                dir = 1
                best_type = typex
                best_mode = mode
            if sumer_inverse < best_score:
                best_score = sumer_inverse
                dir = 0
                best_type = typex
                best_mode = mode
    return best_score, dir, best_mode, best_type


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


def hydrophobicity_grade_increments_ss_aware_with_type(seq, psi_pred):
    # a working version with psi-pred awareness. an increment which is lower than PSI_CUTOFF as a helix
    # will not be added
    results = []
    dirs = []
    modes = []
    types = []
    temp = grade_seq_single_go(seq[0:MIN_WIN])
    if IsHelical(psi_pred[0:MIN_WIN]) or temp[0] <= HPHOBICITY_PSIPRED_BYPASS:
        results.append(temp[0])
        dirs.append(temp[1])
        modes.append(temp[2])
        types.append(temp[3])
    else:
        results.append(1000)
        dirs.append(0)
        modes.append('full')
        types.append('1')
    for inc in range(1, 14):
        if inc+MIN_WIN <= len(seq):
            temp = grade_seq_single_go(seq[0:MIN_WIN+inc])
            if IsHelical(psi_pred[0:MIN_WIN+inc]) or temp[0] <= HPHOBICITY_PSIPRED_BYPASS:
                results.append(temp[0])
                dirs.append(temp[1])
                modes.append(temp[2])
                types.append(temp[3])
            else:
                results.append(1000)
                dirs.append(0)
                modes.append('full')
                types.append('1')
    min_val = min(results)
    min_val_index = results.index(min_val)
    return min_val, min_val_index+MIN_WIN, dirs[min_val_index], modes[min_val_index], types[min_val_index]


def hydrophobicity_grade_increments_ss_aware(seq, mode, psi_pred):
    # a working version with psi-pred awareness. an increment which is lower than PSI_CUTOFF as a helix
    # will not be added
    results = []
    dirs = []
    temp = grade_seq(seq[0:MIN_WIN], mode)
    if IsHelical(psi_pred[0:MIN_WIN]) or temp[0] <= HPHOBICITY_PSIPRED_BYPASS:
        results.append(temp[0])
        dirs.append(temp[1])
    else:
        results.append(1000)
        dirs.append(0)
    for inc in range(1, 14):
        if inc+MIN_WIN <= len(seq):
            temp = grade_seq(seq[0:MIN_WIN+inc], mode)
            if IsHelical(psi_pred[0:MIN_WIN+inc]) or temp[0] <= HPHOBICITY_PSIPRED_BYPASS:
                results.append(temp[0])
                dirs.append(temp[1])
            else:
                results.append(1000)
                dirs.append(0)
    min_val = min(results)
    min_val_index = results.index(min_val)
    return min_val, min_val_index+MIN_WIN, dirs[min_val_index], results


def IsHelical(psi_results):
    # tests whether a segment of psi-pred results has more non-helical residues
    # (measured by PSI_PRED helix score < PSI_PSICUTOFF) than PSR_RES_PREC_CUTOFF
    result = 0
    for aa in psi_results:
        if float(aa) < PSI_CUTOFF:
            result += 1
    precent = float(result) / len(psi_results)
    start = psi_results[0:PSI_EDGE]
    end = psi_results[-PSI_EDGE:]
    start_score = 0
    for aa in start:
        if float(aa) < PSI_CUTOFF_EDGE:
            start_score += 1
    start_score = float(start_score) / len(start)
    end_score = 0
    for aa in end:
        if float(aa) < PSI_CUTOFF_EDGE:
            end_score += 1
    end_score = float(end_score) / len(end)
    if precent >= PSI_RES_PREC_CUTOFF or start_score >= PSI_EDGE_PREC or end_score >= PSI_EDGE_PREC:
        return False
    else:
        return True


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


def WindowGradesForSingleSequence_no_ss(seq, name='default'):
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


def WindowGradesForSingleSequence(seq, name='default', uniprot='default'):
    # runs through an entire sequence, and sends segments of MIN_WIN+36 to hydrophobicity_grade_increments_ss_aware
    # doesn't do that when the WIN_MIN has a non helix residue, adds 1000 instead.
    window_grades_single_chain = {}
    window_grades_single_chain['name'] = name
    full_grades = [1000] * (len(seq) - MIN_WIN)
    full_grades_win_len = [100] * (len(seq) - MIN_WIN)
    full_dir = [100] * (len(seq) - MIN_WIN)
    third_grades = [1000] * (len(seq) - MIN_WIN)
    third_grades_win_len = [100] * (len(seq) - MIN_WIN)
    third_dir = [100] * (len(seq) - MIN_WIN)
    two_thirds_grades = [1000] * (len(seq) - MIN_WIN)
    two_thirds_grades_win_len = [100] * (len(seq) - MIN_WIN)
    two_thirds_dir = [100] * (len(seq) - MIN_WIN)
    types = [10] * (len(seq) - MIN_WIN)
    starters = []
    psi_pred = PsiReaderHelix(uniprot)
    for first in range(0, len(seq) - MIN_WIN):
        temp = hydrophobicity_grade_increments_ss_aware_with_type(seq[first:first+36], psi_pred[first:first+36])
        if temp[3] == 'full':
            full_grades[first] = temp[0]
            full_grades_win_len[first] = temp[1]
            full_dir[first] = temp[2]
        elif temp[3] == 'third':
            third_grades[first] = temp[0]
            third_grades_win_len[first] = temp[1]
            third_dir[first] = temp[2]
        elif temp[3] == 'two_thirds':
            two_thirds_grades[first] = temp[0]
            two_thirds_grades_win_len[first] = temp[1]
            two_thirds_dir[first] = temp[2]
        types[first] = temp[4]
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
    window_grades_single_chain['types'] = types
    window_grades_single_chain['starters'] = starters
    return window_grades_single_chain


def PymolMark(name, minima_tuples, sec_tuples=False, third_TM=False):
    # makes a pml file to describe the programs result, and initiates it.
    print '\n\nFirst minima list is (red): ', minima_tuples
    print 'Second minima list is (blue): ', sec_tuples
    if third_TM: 'Third minima List is (yellow) :', third_TM, '\n\n'
    file = 'test.pml'
    with open(file, 'wa+') as f:
        f.writelines('load ' + name.lower() + '.pdb,' + name + '\n')
        f.writelines('cmd.show_as("cartoon", "all")\n')
        # f.writelines('fetch ' + name + ',' + name + '\n')
        i = 1
        for TM in minima_tuples:
            f.writelines('select TM_1_' + str(i) + ', ' + name + ' and resi ' + str(TM[0]) + '-' + str(TM[0]+TM[1]) + '\n')
            f.writelines('color red, TM_1_' + str(i) + '\n')
            i += 1
        i = 1
        if sec_tuples:
            for TM in sec_tuples:
                f.writelines('select TM_2_' + str(i) + ', ' + name + ' and resi ' + str(TM[0]) + '-' + str(TM[0]+TM[1]) + '\n')
                f.writelines('color blue, TM_2_' + str(i) + '\n')
                i += 1
        if third_TM:
            for TM in third_TM:
                f.writelines('select TM_3_' + str(i) + ', ' + name + ' and resi ' + str(TM[0]) + '-' + str(TM[0]+TM[1]) + '\n')
                f.writelines('color yellow, TM_3_' + str(i) + '\n')
                i += 1
        f.writelines(['save ', name.lower()+'_TM_temp.pse\n'])
    subprocess.call(['/opt/local/bin/pymol', '-q', file, '&'])


def PymolMarkByRange(name, minima_tuples, sec_tuples=False, third_TM=False):
    # makes a pml file to describe the programs result, and initiates it.
    print '\n\nFirst minima list is (red): ', minima_tuples
    print 'Second minima list is (blue): ', sec_tuples
    if third_TM:
        print 'Third minima List is (purple) :', third_TM, '\n\n'
    file = 'test.pml'
    with open(file, 'wa+') as f:
        f.writelines('load ' + name.lower() + '.pdb,' + name + '\n')
        f.writelines('cmd.show_as("cartoon", "all")\n')
        # f.writelines('fetch ' + name + ',' + name + '\n')
        i = 1
        for TM in minima_tuples:
            f.writelines('select TM_1_' + str(i) + ', ' + name + ' and resi ' + str(TM[0]) + '-' + str(TM[1]) + '\n')
            f.writelines('color red, TM_1_' + str(i) + '\n')
            i += 1
        i = 1
        if sec_tuples:
            for TM in sec_tuples:
                f.writelines('select TM_2_' + str(i) + ', ' + name + ' and resi ' + str(TM[0]) + '-' + str(TM[1]) + '\n')
                f.writelines('color blue, TM_2_' + str(i) + '\n')
                i += 1
        if third_TM:
            for TM in third_TM:
                f.writelines('select TM_3_' + str(i) + ', ' + name + ' and resi ' + str(TM[0]) + '-' + str(TM[1]) + '\n')
                f.writelines('color purple, TM_3_' + str(i) + '\n')
                i += 1
        f.writelines(['save ', name.lower()+'_TM_temp.pse\n'])
    subprocess.call(['/opt/local/bin/pymol', '-q', file, '&'])


def PymolMarkExposed(name, grades, minima_tuples):
    # makes a pml file to describe the programs result, and initiates it.
    (grades_smh, win_of_min, min_ind) = Grades2Arry(grades)
    exposed = []
    for tup in minima_tuples:
        if min_ind[tup[0]] == 0 or min_ind[tup[0]] == 1:
            mode = 'full'
        elif min_ind[tup[0]] == 2 or min_ind[tup[0]] == 3:
            mode = 'third'
        elif min_ind[tup[0]] == 4 or min_ind[tup[0]] == 5:
            mode = 'two_thirds'
        if min_ind[tup[0]] == 0 or min_ind[tup[0]] == 2 or min_ind[tup[0]] == 4:
            dir = 0
        else:
            dir = 1
        typer = grades['types'][tup[0]]
        to_color = []
        if dir == 1:
            [to_color.append(x) for x in range(tup[0], tup[1]) if x-tup[0] in MEMBRANE_SPANS[mode][typer]]
        if dir == 0:
            membrane_span = MEMBRANE_SPANS[mode][typer]
            c = []
            [c.append(x in membrane_span) for x in range(0, 50)]
            c = c[::-1]
            [to_color.append(x) for x in range(tup[0], tup[1]) if c[x-tup[0]]]
        print 'found', to_color, typer, mode
        exposed.append(to_color)
    file = 'test.pml'
    with open(file, 'wa+') as f:
        f.writelines('load ' + name.lower() + '.pdb,' + name + '\n')
        f.writelines('cmd.show_as("cartoon", "all")\n')
        for i, colored in enumerate(exposed):
            f.writelines('select TM_1_' + str(i) + ', ' + name + ' and resi ' + '+'.join(str(x) for x in colored) + '\n')
            f.writelines('color red, TM_1_' + str(i) + '\n')
        f.writelines(['save ', name.lower()+'_TM_temp.pse\n'])
    subprocess.call(['/opt/local/bin/pymol', '-q', file])


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
    # smooths the data over SMOOTH_SIZE windows
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
    # a recursive function to find primary minimas. returns a point if it's equal to the minimum of the tested window.
    # if it's not it keeps lokking onward. misses minimas that in their window there's a lower point, even if that point
    # is not a local minima. that's what SecondaryMinimas is for
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
    # recursively calculates primary local minimas over the grades array. uses RecursiveMinima
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
    np.set_printoptions(threshold=np.inf)
    np.seterr(invalid='ignore')
    minimas[minimas > PRIMARY_MINIMA_THRESHOLD] = np.NAN
    np.seterr(invalid=None)
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
                    if point >= len(wins):
                        break
                    if minima[0] <= point+wins[point] < minima[0]+minima[1]:
                        point = minima[0]+minima[1]
        point += 1
    return sec_minimas


def Grades2Arry(grades):
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
    return grd_smh, win_of_min, min_ind


def MinimaOrganizer(grds, win_of_min):
    minimas = LocalMinima(grds, win_of_min)
    # np.set_printoptions(threshold=np.inf)

    minima_tuples = MinimaTuples(minimas, win_of_min)
    minimas_secondary = SecondaryMinimas(grds, minima_tuples, win_of_min)
    sec_minima_tuples = MinimaTuples(minimas_secondary, win_of_min)
    all_minimas_tuples = minima_tuples[:]
    [all_minimas_tuples.append(x) for x in sec_minima_tuples]
    ter_minimas = SecondaryMinimas(grds, all_minimas_tuples, win_of_min)
    for row in range(0, 6):
        for col, val in enumerate(ter_minimas[row]):
            if val != np.NAN:
                minimas_secondary[row][col] = val
    ter_minimas_tuples = MinimaTuples(ter_minimas, win_of_min)
    [sec_minima_tuples.append(x) for x in ter_minimas_tuples]
    return minimas, minima_tuples, minimas_secondary, sec_minima_tuples


def PlotSinglePeptide(grades):
    (grd_smh, win_of_min, min_ind) = Grades2Arry(grades)

    (minimas, minima_tuples, minimas_secondary, sec_minima_tuples) = MinimaOrganizer(grd_smh, win_of_min)
    # PymolMark(grades['name'], minima_tuples, sec_minima_tuples)
    plot_array = np.empty(shape=[6, len(grades['starters'])])
    plot_array[:] = np.NAN
    for col, row in enumerate(min_ind):
        plot_array[row][col] = grd_smh[row][col]

    PlotParams()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    f = ax1.scatter(grades['starters'], plot_array[0], s=100, c='k')
    fi = ax1.scatter(grades['starters'], plot_array[1], s=100, facecolors='none', edgecolors='k')
    t = ax1.scatter(grades['starters'], plot_array[2], s=100, c='r')
    ti = ax1.scatter(grades['starters'], plot_array[3], s=100, facecolors='none', edgecolors='r')
    tw = ax1.scatter(grades['starters'], plot_array[4], s=100, c='b')
    twi = ax1.scatter(grades['starters'], plot_array[5], s=100, facecolors='none', edgecolors='b')

    fm = ax1.scatter(grades['starters'], minimas[0], s=300, c='k', marker='D')
    fmi = ax1.scatter(grades['starters'], minimas[1], s=300, c='k', marker='*')
    tm = ax1.scatter(grades['starters'], minimas[2], s=300, c='r', marker='D')
    tmi = ax1.scatter(grades['starters'], minimas[3], s=300, c='r', marker='*')
    twm = ax1.scatter(grades['starters'], minimas[4], s=300, c='b', marker='D')
    twmi = ax1.scatter(grades['starters'], minimas[5], s=300, c='b', marker='*')

    fsm = ax1.scatter(grades['starters'], minimas_secondary[0], s=280, c='k', marker='8')
    fsmi = ax1.scatter(grades['starters'], minimas_secondary[1], s=280, c='k', marker='s')
    tsm = ax1.scatter(grades['starters'], minimas_secondary[2], s=280, c='r', marker='8')
    tsmi = ax1.scatter(grades['starters'], minimas_secondary[3], s=280, c='r', marker='s')
    twsm = ax1.scatter(grades['starters'], minimas_secondary[4], s=280, c='b', marker='8')
    twsmi = ax1.scatter(grades['starters'], minimas_secondary[5], s=280, c='b', marker='s')

    ax1.plot(grades['starters'], [0]*len(grades['starters']), 'k--', linewidth=3)
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
    ax1.set_xlabel('Window Start Residue', fontsize=40)
    ax1.set_ylabel('$\Delta$G$_{transfer}$', fontsize=40)
    # ax2.set_ylabel('Window Length')
    ax1.legend((f, fi, t, ti, tw, twi, fm, fmi, tm, tmi, twm, twmi, fsm, fsmi, tsm, tsmi, twsm, twsmi),
               ('Full', 'Full -1', 'Third', 'Third -1', 'Two-Thirds', 'Two-Thirds -1', 'Full minima', 'Full -1 minima',
                'Third minima', 'Third -1 minima', 'Two Thirds minima', 'Two Thirds -1 minima', 'Full SM', 'Full -1 SM',
                'Third SM', 'Third -1 SM', 'Two Thirds SM', 'Two Thirds -1 SM'), scatterpoints=1, ncol=9,
               loc='upper center', prop={'size': 12})
    ax1.tick_params(axis='both', which='major', labelsize=20)
    plt.suptitle(grades['name']+' plot with smooth '+str(SMOOTH_SIZE)+' and loop-bypass '+str(LOOP_BYPASS), fontsize=48)
    plt.show()


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


def Range2Tups(ranger):
    results = []
    temp = [ranger[0]]
    for i, val in enumerate(ranger):
        if val-1 > ranger[i-1]:
            temp.append(ranger[i-1])
            results.append(temp[:])
            temp = [val]
    temp.append(ranger[-1])
    results.append(temp)
    return results


def CompareOverlap(sw_tuples, our_tuples):
    SW_range = []
    [SW_range.append(range(int(tup[0]), int(tup[1]))) for tup in sw_tuples]
    SW_range = [y for x in SW_range for y in x]
    OUR_range = []
    [OUR_range.append(range(int(tup[0]), int(tup[1]))) for tup in our_tuples]
    OUR_range = [y for x in OUR_range for y in x]

    a_multiset = Counter(SW_range)
    b_multiset = Counter(OUR_range)
    overlap = list((a_multiset & b_multiset).elements())
    a_remainder = list((a_multiset - b_multiset).elements())
    b_remainder = list((b_multiset - a_multiset).elements())

    result = float((len(overlap) - len(a_remainder) - len(b_remainder))) / len(SW_range)
    return result, overlap, a_remainder, b_remainder


def AlignForPyMol(sw_seq, pdb_name, sw_tups, our_tups):
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    from pdb_sequences import pdb_seq_retrieve
    pdb_seq = pdb_seq_retrieve(pdb_name)
    matrix = matlist.blosum62
    gap_open = -10
    gap_extend = -0.5
    # top_aln = pairwise2.align.globalds(sw_seq, pdb_seq, matrix, gap_open, gap_extend)[0]
    top_aln = pairwise2.align.globalds(sw_seq, pdb_seq, matrix, gap_open, gap_extend)[0]
    (sw_aln, pdb_aln, score, begin, end) = top_aln
    print 'sw ', sw_aln
    print 'pdb ', pdb_aln
    sw_new_tups = []


def SW2PDB(sw_seq, pdb_name, sw_tups, our_tups):
    from pdb_sequences import pdb_seq_retrieve
    pdb_seq = pdb_seq_retrieve(pdb_name)
    sw_pdb_tups = []
    sw_color = []
    for tup in sw_tups:
        start = pdb_seq.find(sw_seq[int(tup[0]):int(tup[1])])
        i = 0
        while start == -1:
            i += 1
            start = pdb_seq.find(sw_seq[int(tup[0])+i:int(tup[1])-i])
        end = start + int(tup[1]) - int(tup[0]) - i
        sw_color.append([start, end-start])
        sw_pdb_tups.append([start, end])
    our_pdb_tups = []
    our_color = []
    for tup in our_tups:
        start = pdb_seq.find(sw_seq[int(tup[0]):int(tup[1])])
        i = 0
        while start == -1:
            i += 1
            start = pdb_seq.find(sw_seq[int(tup[0])+i:int(tup[1])-i])
        end = start + int(tup[1]) - int(tup[0]) - i
        our_pdb_tups.append([start, end])
        our_color.append([start, end-start])
    PymolMark(pdb_name, sw_color, our_color)
    # print sw_pdb_tups, sw_tups
    # print our_pdb_tups, our_tups


def TotalHydroPhobicity(seq):
    result = 0
    for aa in seq:
        result += np.polyval(hydrophobicity_polyval[aa], 0)
    return result


def WriteSW2CSV(num=False, name=False, uniprot=False):
    # reads num sequences out of a csv produced by database_parser.py (2nd funciton), and processes each sequences.
    if num:
        SW = database_parser.SWDB_parser_prediciton(num)
    elif name:
        SW = database_parser.SWDB_parser_prediciton_by_name(name)
    elif uniprot:
        SW = database_parser.SWDB_parser_prediciton_by_uniprot(uniprot)
    combined = open('/Users/jonathan/Documents/membrane_prediciton_data/combined_results.csv', 'wa+')
    csv_writer = csv.writer(combined)
    all_combined = {}

    results_dict = {'uniprot': '', 'seq_length': 0, 'total_hphobicity': 0, 'overlap': 0, 'sw_remainder': 0,
                    'our_remainder': 0, 'overlap_score': 0}
    csv_writer.writerow(results_dict.keys())

    for key, val in SW.iteritems():
        if key == '':
            continue
        ss_grades = WindowGradesForSingleSequence(val['seq'], key, val['uniprot'])
        PlotSinglePeptide(ss_grades)
        (grds_array, win_of_min, min_ind) = Grades2Arry(ss_grades)
        (pri_min, pri_min_tup, sec_min, sec_min_tup) = MinimaOrganizer(grds_array, win_of_min)

        all_combined[key] = {}
        all_combined[key]['seq'] = val['seq']
        all_combined[key]['pdb'] = val['pdb']
        all_combined[key]['terminal'] = val['term']
        SW_tuples = []
        for i, begin in enumerate(val['begin']):
            SW_tuples.append([begin, val['end'][i]])
        all_combined[key]['SW_tuples'] = SW_tuples
        all_combined[key]['OUR_primary_minima'] = []
        for tup in pri_min_tup:
            all_combined[key]['OUR_primary_minima'].append([tup[0], tup[0]+tup[1]])
        all_combined[key]['OUR_secondary_minima'] = []
        for tup in sec_min_tup:
            all_combined[key]['OUR_secondary_minima'].append([tup[0], tup[0]+tup[1]])
        all_our_tuples = [] # important re-initializisation
        all_our_tuples = all_combined[key]['OUR_primary_minima'][:]
        [all_our_tuples.append(x) for x in all_combined[key]['OUR_secondary_minima']]
        (overlap_score, overlap_list, sw_remainder, our_remainder) = CompareOverlap(SW_tuples, all_our_tuples)
        overlap_list.sort()
        sw_remainder.sort()
        our_remainder.sort()
        # AlignForPyMol(val['seq'], val['pdb'][0], SW_tuples, all_our_tuples)
        # PymolMarkByRange(val['pdb'][0], Range2Tups(overlap_list), Range2Tups(sw_remainder), Range2Tups(our_remainder))
        # for complexes with multiple chains separated into uniprot names pdbs
        PymolMarkByRange(uniprot, Range2Tups(overlap_list), Range2Tups(sw_remainder), Range2Tups(our_remainder))
        # SW2PDB(val['seq'], val['pdb'][0], SW_tuples, all_our_tuples)
        print overlap_score, ' for ', val['pdb']
        # print SW_tuples, all_our_tuples
        # csv_writer.writerow([val['uniprot'], overlap_score, len(overlap_list), len(sw_remainder), len(our_remainder)])
        # print 'aaa', overlap_score, overlap_list, sw_remainder, our_remainder
        results_dict = {'uniprot': val['uniprot'], 'seq_length': len(val['seq']), 'total_hphobicity':
            TotalHydroPhobicity(val['seq']), 'overlap': len(overlap_list), 'sw_remainder': len(sw_remainder),
                        'our_remainder': len(our_remainder), 'overlap_score': overlap_score}
        # csv_writer.writerow(results_dict.values())
        # PymolMarkExposed(uniprot, ss_grades, all_our_tuples)

    combined.close()


def PlotParams():
    # by: http://damon-is-a-geek.com/publication-ready-the-first-time-beautiful-reproducible-plots-with-matplotlib.html
    from matplotlib import rcParams
    rcParams['axes.labelsize'] = 9
    rcParams['xtick.labelsize'] = 9
    rcParams['ytick.labelsize'] = 9
    rcParams['legend.fontsize'] = 9

    # os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
    # rcParams['font.family'] = 'serif'
    # rcParams['font.serif'] = ['Computer Modern Roman']
    # rcParams['text.usetex'] = True


### MAIN ###
MakeHydrophobicityGrade()
# UniprotParser()
# CalculateWindowGrades()
# WriteDataSetToTSV('/Users/jonathan/Desktop/python_temp.txt')
# ss_grades = WindowGradesForSingleSequence('MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKN', '4j4q')
# ss_grades = WindowGradesForSingleSequence('TGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFSMRPEVATFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFG', '1c3w')
# ss_grades = WindowGradesForSingleSequence('MELEEDLKGRADKNFSKMGKKSKKEKKEKKPAVSVLTMFRYAGWLDRLYMLVGTLAAIIHGVALPLMMLIFGDMTDSFASVGNVSKNSTNMSEADKRAMFAKLEEEMTTYAYYYTGIGAGVLIVAYIQVSFWCLAAGRQIHKIRQKFFHAIMNQEIGWFDVHDVGELNTRLTDDVSKINEGIGDKIGMFFQAMATFFGGFIIGFTRGWKLTLVILAISPVLGLSAGIWAKILSSFTDKELHAYAKAGAVAEEVLAAIRTVIAFGGQKKELERYNNNLEEAKRLGIKKAITANISMGAAFLLIYASYALAFWYGTSLVISKEYSIGQVLTVFFSVLIGAFSVGQASPNIEAFANARGAAYEVFKIIDNKPSIDSFSKSGHKPDNIQGNLEFKNIHFSYPSRKEVQILKGLNLKVKSGQTVALVGNSGCGKSTTVQLMQRLYDPLDGMVSIDGQDIRTINVRYLREIIGVVSQEPVLFATTIAENIRYGREDVTMDEIEKAVKEANAYDFIMKLPHQFDTLVGERGAQLSGGQKQRIAIARALVRNPKILLLDEATSALDTESEAVVQAALDKAREGRTTIVIAHRLSTVRNADVIAGFDGGVIVEQGNHDELMREKGIYFKLVMTQTAGNEIELGNEACKSKDEIDNLDMSSKDSGSSLIRRRSTRKSICGPHDQDRKLSTKEALDEDVPPASFWRILKLNSTEWPYFVVGIFCAIINGGLQPAFSVIFSKVVGVFTNGGPPETQRQNSNLFSLLFLILGIISFITFFLQGFTFGKAGEILTKRLRYMVFKSMLRQDVSWFDDPKNTTGALTTRLANDAAQVKGATGSRLAVIFQNIANLGTGIIISLIYGWQLTLLLLAIVPIIAIAGVVEMKMLSGQALKDKKELEGSGKIATEAIENFRTVVSLTREQKFETMYAQSLQIPYRNAMKKAHVFGITFSFTQAMMYFSYAACFRFGAYLVTQQLMTFENVLLVFSAIVFGAMAVGQVSSFAPDYAKATVSASHIIRIIEKTPEIDSYSTQGLKPNMLEGNVQFSGVVFNYPTRPSIPVLQGLSLEVKKGQTLALVGSSGCGKSTVVQLLERFYDPMAGSVFLDGKEIKQLNVQWLRAQLGIVSQEPILFDCSIAENIAYGDNSRVVSYEEIVRAAKEANIHQFIDSLPDKYNTRVGDKGTQLSGGQKQRIAIARALVRQPHILLLDEATSALDTESEKVVQEALDKAREGRTCIVIAHRLSTIQNADLIVVIQNGKVKEHGTHQQLLAQKGIYFSMVSVQAGAKRSYVHH', '3g61')
# ss_grades = WindowGradesForSingleSequence('AGFDGGVIVEQGNHDELMREKGIYFKLVMTQTAGNEIELGNEACKSKDEIDNLDMSSKDSGSSLIRRRSTRKSICGPHDQDRKLSTKEALDEDVPPASFWRILKLNSTEWPYFVVGIFCAIINGGLQPAFSVIFSKVVGVFTNGGPPETQRQNSNLFSLLFLILGIISFITFFLQGFTFGKAGEILTKRLRYMVFKSMLRQDVSWFDDPKNTTGALTTRLANDAAQVKGATGSRLAVIFQNIANLGTGIIISLIYGWQLTLLLLAIVPIIAIAGVVEMKMLSGQALKDKKELEGSGKIATEAIENFRTVVSLTREQKFETMYAQSLQIPYRNAMKKAHVFGITFSFTQAMMYFSYAACFRFGAYLVTQQLMTFENVLLVFSAIVFGAMAVGQVSSFAPDYAKATVSASHIIRIIEKTPEIDSYSTQGLKPNMLEGNVQFSGVVFNYPTRPSIPVLQGLSLEVKKGQTLALVGSSGCGKSTVVQLLERFYDPMAGSVFLDGKEIKQLNVQWLRAQLGIVSQEPILFDCSIAENIAYGDNSRVVSYEEIVRAAKEANIHQFIDSLPDKYNTRVGDKGTQLSGGQKQRIAIARALVRQPHILLLDEATSALDTESEKVVQEALDKAREGRT', '3g61')
# ss_grades = WindowGradesForSingleSequence('VSVLTMFRYAGWLDRLYMLVGTLAAIIHGVALPLMMLIFGDMTDSFASVGNVSKNSTNMSEADKRAMFAKLEEEMTTYAYYYTGIGAGVLIVAYIQVSFWCLAAGRQIHKIRQKFFHAIMNQEIGWFDVHDVGELNTRLTDDVSKINEGIGDKIGMFFQAMATFFGGFIIGFTRGWKLTLVILAISPVLGLSAGIWAKILSSFTDKELHAYAKAGAVAEEVLAAIRTVIAFGGQKKELERYNNNLEEAKRLGIKKAITANISMGAAFLLIYASYALAFWYGTSLVISKEYSIGQVLTVFFSVLIGAFSVGQASPNIEAFANARGAAYEVFKIIDNKPSIDSFSKSGHKPDNIQGNLEFKNIHFSYPSRKEVQILKGLNLKVKSGQTVALVGNSGCGKSTTVQLMQRLYDPLDGMVSIDGQDIRTINVRYLREIIGVVSQEPVLFATTIAENIRYGREDVTMDEIEKAVKEANAYDFIMKLPHQFDTLVGERGAQLSGGQKQRIAIARALVRNPKILLLDEATSALDTESEAVVQAALDKAREGRTTIVIAHRLSTVRNADVIAGFDGGVIVEQGNHDELMREKGIYFKLVMTQTLDEDVPPASFWRILKLNSTEWPYFVVGIFCAIINGGLQPAFSVIFSKVVGVFTNGGPPETQRQNSNLFSLLFLILGIISFITFFLQGFTFGKAGEILTKRLRYMVFKSMLRQDVSWFDDPKNTTGALTTRLANDAAQVKGATGSRLAVIFQNIANLGTGIIISLIYGWQLTLLLLAIVPIIAIAGVVEMKMLSGQALKDKKELEGSGKIATEAIENFRTVVSLTREQKFETMYAQSLQIPYRNAMKKAHVFGITFSFTQAMMYFSYAACFRFGAYLVTQQLMTFENVLLVFSAIVFGAMAVGQVSSFAPDYAKATVSASHIIRIIEKTPEIDSYSTQGLKPNMLEGNVQFSGVVFNYPTRPSIPVLQGLSLEVKKGQTLALVGSSGCGKSTVVQLLERFYDPMAGSVFLDGKEIKQLNVQWLRAQLGIVSQEPILFDCSIAENIAYGDNSRVVSYEEIVRAAKEANIHQFIDSLPDKYNTRVGDKGTQLSGGQKQRIAIARALVRQPHILLLDEATSALDTESEKVVQEALDKAREGRTCIVIAHRLSTIQNADLIVVIQNGKVKEHGTHQQLLAQKGIYFSMVSVQA', '3g61')
# ss_grades = WindowGradesForSingleSequence('NGGLQPAFSVIFSKVVGVFTNGGPPETQRQNSNLFSLLFLILGIISFITFFLQGFTFGKAGEILTKRLR', 'part_of_3g61')
# ss_grades = WindowGradesForSingleSequence('GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRA', '1BRX')
# ss_grades = WindowGradesForSingleSequence('MDMLFAKTVVLAASAVGAGTAMIAGIGPGVGQGYAAGKAVESVARQPEAKGDIISTMVLGQAVAESTGIYSLVIALILLYANPFVGLLG', '1yce')
# ss_grades = WindowGradesForSingleSequence('YAQALQSVPETQVSQLDNGVRVASEQSSQPTCTVGVWIDAGSRYESEKNNGAGYFLEHLAFKGTKNRPQNALEKEVESMGAHLNAYSSREHTAYYIKALSKDVPKAVELLADIVQNCSLEDSQIEKERDVIVRELQENDTSMREVVFNYLHATAFQGTGLAQSVEGPSENIRKLSRADLTEYLSTHYTAPRMVLAAAGGVEHQQLLELAQKHFGGVPFTYDDDAVPTLSKCRFTGSQIRHREDGLPLAHVAIAVEGPGWAHPDLVALQVANAIIGHYDRTYGGGLHSSSPLASIAVTNKLCQSFQTFSICYSETGLFGFYFVCDRMSIDDMMFVLQGQWMRLCTSISESEVLRGKNFLRNALVSHLDGTTPVCEDIGRELLTYGRRIPLEEWEERLAEVDARMVREVCSKYIYDQCPAVAGPGPIEQLPDYNRIRSGMFWLRPPHPQDLEITKLPNGLVIASLENYSPGSTIGVFIKAGSRYENSSNLGTSHLLRLASSLTTKGASSFKITRGIEAVGGKLSVESTRENMAYTVECLRDDVEILMEFLLNVTTAPEFRPWEVADLQPQLKIDKAVAFQNPQTHVIENLHAAAYRNALADSLYCPDYRIGKVTSVELHDFVQNHFTSARMALVGLGVSHPVLKNVAEQLLNIRGGLGLSGAKAKYRGGEIREQNGDSLVHAAIVAESAAIGGAEANAFSVLQHVLGANPHVKRGNPFDVSAFNASYSDSGLFGFYTISQAAYAGQVIKAAYNQVKTIAQGNVSNENVQAAKNKLKAKYLMSVESSEGFLEEVGSQALAAGSYNPPSTVLQQIDAVADADVIKAAKKFVSRQKSMAASGNLGHTPFVDELAPNIRKSHPLLKMINNSLIDLPAPSNISAWWNFGSLLAVCLMTQILTGLLLAMHYTADTSLAFSSVAHTCRNVQYGWLIRNLHANGASFFFICIFLHIGRGLYYGSYLYKETWNTGVILLLTLMATAFVGYVLPWGQMSFWGATVITNLFSAIPYIGHTLVEWAWGGFSVDNPTLTRFFALHFLLPFAIAGITIIHLTFLHESGSNNPLGISSDSDKIPFHPYYSFKDILGLTLMLTPFLTLALFSPNLLGDPENFTPANPLVTPPHIKPEWYFLFAYAILRSIPNKLGGVLALAASVLILFLIPFLHKSKQRTMTFRPLSQTLFWLLVANLLILTWIGSQPVEHPFIIIGQMASLSYFTILLILFPTIGTLENKMLNYSDLELHPPSYPWSHRGPLSSLDHTSIRRGFQVYKQVCSSCHSMDYVAYRHLVGVCYTEDEAKALAEEVEVQDGPNEDGEMFMRPGKLSDYFPKPYPNPEAARAANNGALPPDLSYIVRARHGGEDYVFSLLTGYCEPPTGVSVREGLYFNPYFPGQAIGMAPPIYNDVLEFDDGTPATMSQVAKDVCTFLRWAAEPEHDHRKRMGLKMLLMMGLLVPLVYYMKRHKWSVLKSRKLAYRPPKSHTDIKVPNFSDYRRPPDDYSTKSSRESDPSRKGFSYLVTAVTTLGVAYAAKNVVTQFVSSMSASADVLAMSKIEIKLSDIPEGKNMAFKWRGKPLFVRHRTKKEIDQEAAVEVSQLRDPQHDLERVKKPEWVILIGVCTHLGCVPIANAGDFGGYYCPCHGSHYDASGRIRKGPAPLNLEVPSYEFTSDDMVIVGSRWLEGIRKWYYNAAGFNKYGLMRDDTIYENDDVKEAIRRLPENLYDDRMFRIKRALDLNMRQQILPKEQWTKYEEDVPYLEPYLKEVIRERKEREEWDKRQFGHLTRVRHLITYSLSPFEQRPFPHYFSKGVPNVWRRLRACILRVAPPFLAFYLLYTWGTQEFEKSKRKNPAAYVNLVDPLTTVREQCEQLEKCVKARERLELCDERVSSRSQTEEDCTEELFDFLHARDHCVAHKLFNSLKTLTARLYSLLFRRTSTFALTIVVGALLFERAFDQGADAIYEHINEGKLWKHIKHKYENK', '1bcc')
# ss_grades = WindowGradesForSingleSequence('DRPIFAWVIAIIIMLAGGLAILKLPVAQYPTIAPPAVTISASYPGADAKTVQDTVTQVIEQNMNGIDNLMYMSSNSDSTGTVQITLTFESGTDADIAQVQVQNKLQLAMPLLPQEVQQQGVSVEKSSSSFLMVVGVINTDGTMTQEDISDYVAANMKDAISRTSGVGDVQLFGSQYAMRIWMNPNELNKFQLTPVDVITAIKAQNAQVAAGQLGGTPPVKGQQLNASIIAQTRLTSTEEFGKILLKVNQDGSRVLLRDVAKIELGGENYDIIAEFNGQPASGLGIKLATGANALDTAAAIRAELAKMEPFFPSGLKIVYPYDTTPFVKISIHEVVKTLVEAIILVFLVMYLFLQNFRATLIPTIAVPVVLLGTFAVLAAFGFSINTLTMFGMVLAIGLLVDDAIVVVENVERVMAEEGLPPKEATRKSMGQIQGALVGIAMVLSAVFVPMAFFGGSTGAIYRQFSITIVSAMALSVLVALILTPALCATMLKFGWFNRMFEKSTHHYTDSVGGILRSTGRYLVLYLIIVVGMAYLFVRLPSSFLPDEDQGVFMTMVQLPAGATQERTQKVLNEVTHYYLTKEKNNVESVFAVNGFGFAGRGQNTGIAFVSLKDWADRPGEENKVEAITMRATRAFSQIKDAMVFAFNLPAIVELGTATGFDFELIDQAGLGHEKLTQARNQLLAEAAKHPMLTSVRPNGLEDTPQFKIDIDQEKAQALGVSINDINTTLGAAWGGSYVNDFIDRGRVKKVYVMSEAKYRMLPDDIGDWYVRAADGQMVPFSAFSSSRWEYGSPRLERYNGLPSMEILGQAAPGKSTGEAMELMEQLASKLPTGVGYDWSGNQAPSLYAISLIVVFLCLAALYESWSIPFSVMLVVPLGVIGALLAATFRGLTNDVYFQVGLLTTIGLSAKNAILIVEFAKDLMDKEGKGLIEATLDAVRMRLRPILMTSLAFILGVMPLVISTGAGSGAQNAVGTGVMGGMVTATVLAIFFVPVFFVVVRRRFSRK', '1iwg')
## 1brx with SW sequence
# ss_grades = WindowGradesForSingleSequence('GRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRA', '1BRX')
# new 3g61:
# ss_grades = WindowGradesForSingleSequence('VSVLTMFRYAGWLDRLYMLVGTLAAIIHGVALPLMMLIFGDMTDSFASVGNVSKNSTNMSEADKRAMFAKLEEEMTTYAYYYTGIGAGVLIVAYIQVSFWCLAAGRQIHKIRQKFFHAIMNQEIGWFDVHDVGELNTRLTDDVSKINEGIGDKIGMFFQAMATFFGGFIIGFTRGWKLTLVILAISPVLGLSAGIWAKILSSFTDKELHAYAKAGAVAEEVLAAIRTVIAFGGQKKELERYNNNLEEAKRLGIKKAITANISMGAAFLLIYASYALAFWYGTSLVISKEYSIGQVLTVFFSVLIGAFSVGQASPNIEAFANARGAAYEVFKIIDNKPSIDSFSKSGHKPDNIQGNLEFKNIHFSYPSRKEVQILKGLNLKVKSGQTVALVGNSGCGKSTTVQLMQRLYDPLDGMVSIDGQDIRTINVRYLREIIGVVSQEPVLFATTIAENIRYGREDVTMDEIEKAVKEANAYDFIMKLPHQFDTLVGERGAQLSGGQKQRIAIARALVRNPKILLLDEATSALDTESEAVVQAALDKAREGRTTIVIAHRLSTVRNADVIAGFDGGVIVEQGNHDELMREKGIYFKLVMTQTLDEDVPPASFWRILKLNSTEWPYFVVGIFCAIINGGLQPAFSVIFSKVVGVFTNGGPPETQRQNSNLFSLLFLILGIISFITFFLQGFTFGKAGEILTKRLRYMVFKSMLRQDVSWFDDPKNTTGALTTRLANDAAQVKGATGSRLAVIFQNIANLGTGIIISLIYGWQLTLLLLAIVPIIAIAGVVEMKMLSGQALKDKKELEGSGKIATEAIENFRTVVSLTREQKFETMYAQSLQIPYRNAMKKAHVFGITFSFTQAMMYFSYAACFRFGAYLVTQQLMTFENVLLVFSAIVFGAMAVGQVSSFAPDYAKATVSASHIIRIIEKTPEIDSYSTQGLKPNMLEGNVQFSGVVFNYPTRPSIPVLQGLSLEVKKGQTLALVGSSGCGKSTVVQLLERFYDPMAGSVFLDGKEIKQLNVQWLRAQLGIVSQEPILFDCSIAENIAYGDNSRVVSYEEIVRAAKEANIHQFIDSLPDKYNTRVGDKGTQLSGGQKQRIAIARALVRQPHILLLDEATSALDTESEKVVQEALDKAREGRTCIVIAHRLSTIQNADLIVVIQNGKVKEHGTHQQLLAQKGIYFSMVSVQA', '3g61')
WriteSW2CSV(uniprot='p02945')
# WriteSW2CSV(20)



# SW = SWDB_parser_prediciton(1)
# ss_grades = WindowGradesForSingleSequence(SW[SW.keys()[0]]['seq'], SW[SW.keys()[0]]['uniprot'])
# PlotSinglePeptide(ss_grades)

# WriteSingleSeqToTSV(ss_grades, '/Users/jonathan/Desktop/bassaf_seq.txt')
# WriteSingleSeqToTSVusingCSV(ss_grades, '/Users/jonathan/Desktop/4j4q_win_gds.csv')

# PlotSinglePeptideWindows(ss_grades)
# PlotSinglePeptide(ss_grades)