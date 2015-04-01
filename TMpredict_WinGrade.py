from WinGrade import *
from HphobicityScore import *


def main():
    import subprocess
    import re
    global hydrophobicity_polyval
    # import topdb_functions
    hydrophobicity_polyval = MakeHydrophobicityGrade()
    # temp = HphobicityScore('te
    # mp', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    # temp = HphobicityScore('temp', 'YSYRFVWWAISTAAMLYILY')
    # temp = HphobicityScore('1E12', 'MSITSVPGVVDAGVLGAQSAAAVRENALLSSSLWVNVALAGIAILVFVYMGRTIRPGRPRLIWGATLMIPLVSISSYLGLLSGLTVGMIEMPAGHALAGEMVRSQWGRYLTWALSTPMILLALGLLADVDLGSLFTVIAADIGMCVTGLAAAMTTSALLFRWAFYAISCAFFVVVLSALVTDWAASASSAGTAEIFDTLRVLTVVLWLGYPIVWAVGVEGLALVQSVGVTSWAYSVLDVFAKYVFAFILLRWVANNERTVAVAGQTLGTMSSDD', '../psipred/sw_fastas/P16102.ss2',hydrophobicity_polyval)
    # temp = HphobicityScore('1IWG', 'MPNFFIDRPIFAWVIAIIIMLAGGLAILKLPVAQYPTIAPPAVTISASYPGADAKTVQDTVTQVIEQNMNGIDNLMYMSSNSDSTGTVQITLTFESGTDADIAQVQVQNKLQLAMPLLPQEVQQQGVSVEKSSSSFLMVVGVINTDGTMTQEDISDYVAANMKDAISRTSGVGDVQLFGSQYAMRIWMNPNELNKFQLTPVDVITAIKAQNAQVAAGQLGGTPPVKGQQLNASIIAQTRLTSTEEFGKILLKVNQDGSRVLLRDVAKIELGGENYDIIAEFNGQPASGLGIKLATGANALDTAAAIRAELAKMEPFFPSGLKIVYPYDTTPFVKISIHEVVKTLVEAIILVFLVMYLFLQNFRATLIPTIAVPVVLLGTFAVLAAFGFSINTLTMFGMVLAIGLLVDDAIVVVENVERVMAEEGLPPKEATRKSMGQIQGALVGIAMVLSAVFVPMAFFGGSTGAIYRQFSITIVSAMALSVLVALILTPALCATMLKPIAKGDHGEGKKGFFGWFNRMFEKSTHHYTDSVGGILRSTGRYLVLYLIIVVGMAYLFVRLPSSFLPDEDQGVFMTMVQLPAGATQERTQKVLNEVTHYYLTKEKNNVESVFAVNGFGFAGRGQNTGIAFVSLKDWADRPGEENKVEAITMRATRAFSQIKDAMVFAFNLPAIVELGTATGFDFELIDQAGLGHEKLTQARNQLLAEAAKHPDMLTSVRPNGLEDTPQFKIDIDQEKAQALGVSINDINTTLGAAWGGSYVNDFIDRGRVKKVYVMSEAKYRMLPDDIGDWYVRAADGQMVPFSAFSSSRWEYGSPRLERYNGLPSMEILGQAAPGKSTGEAMELMEQLASKLPTGVGYDWTGMSYQERLSGNQAPSLYAISLIVVFLCLAALYESWSIPFSVMLVVPLGVIGALLAATFRGLTNDVYFQVGLLTTIGLSAKNAILIVEFAKDLMDKEGKGLIEATLDAVRMRLRPILMTSLAFILGVMPLVISTGAGSGAQNAVGTGVMGGMVTATVLAIFFVPVFFVVVRRRFSRKNEDIEHSHTVDHH', '../psipred/sw_fastas/P31224.ss2',hydrophobicity_polyval)
    # temp = HphobicityScore('1BRX', 'EAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSAGDGAAATS', '../psipred/sw_fastas/P02945.ss2',hydrophobicity_polyval)

    ### parameters and setup for running von heijne DB entries:
    # param_list = [0, 20, 0.2, 3]
    # vdb_dict = parse_v_db()
    # right = 0
    # wrong = 0
    # for v_entry in vdb_dict.values():
    #     if v_entry['name'].lower() == 'Atpe'.lower():
    #         print v_entry
    #         temp = HphobicityScore(v_entry['name'], v_entry['seq'], 'data_sets/VDB/'+v_entry['name']+'.ss2', hydrophobicity_polyval, param_list)
    #         temp.plot_win_grades()
            # if v_entry['cterm'] == temp.best_c_term:
            #     right += 1
            # else:
            #     wrong += 1
            # print temp.name, temp.best_c_term, temp.topo_best_val, temp.sec_best_c_term, temp.topo_sec_best_val, v_entry['cterm']


    # temp_topdb = topdb_functions.read_entries(False, 0, 1)
    # for temp_db in temp_topdb:
    #     temp = HphobicityScore(temp_db['name'], temp_db['seq'], temp_db['ss2'], hydrophobicity_polyval)
    #     topdb_functions.topo_compare(temp.topo_string, temp_db['topo'])
    #     print temp_db['name']
    #     print temp_db['seq']
        # temp.plot_win_grades()
        # print temp.topo
        # print temp.topo_string
        # print temp_db['topo']
        # pymol_mark_segments(temp.name, [[[i.begin, i.end] for i in temp.topo]])

    # db_entries = parsed_data_base_parser(1, 2)
    # print db_entries[0]
    #
    # temp = HphobicityScore(db_entries[0]['pdb'], db_entries[0]['seq'], '../psipred/sw_fastas/'+db_entries[0]['uniprot']+'.ss2', hydrophobicity_polyval)
    # print temp
    # temp.plot_win_grades()
    # temp.plot_energy_landscape()
    # print temp.make_topo_string()
    # pymol_mark_segments(temp.name, [[[i.begin, i.end] for i in temp.topo]])


    '''
    parses a range of SW entries, and prints the topology predcition reult
    '''
    # db_entries = parsed_data_base_parser(0,100)
    # topo_predict_score = {'good': 0, 'bad': 0}
    # for protein in db_entries:
    #     temp = HphobicityScore(protein['pdb'], protein['seq'], protein['uniprot'], hydrophobicity_polyval)
    #     if temp.n_term_orient == 'rev' and protein['orientation'] == 'out':
    #         topo_predict_score['good'] += 1
    #         print 'was correct', temp.n_term_orient, protein['orientation'], topo_predict_score
    #     else:
    #         topo_predict_score['bad'] += 1
    #         print 'was wrong', temp.n_term_orient, protein['orientation'], topo_predict_score
    # print 'prediction results:', topo_predict_score

    ### parsing and running rostlab_db entries:
    # param_list = [0, 20, 0.2, 3]
    param_list = [0, 20, 0.2, 5]
    rostlab_db_dict = parse_rostlab_db()
    # results = {'tm_num_correct': 0, 'tm_num_incorrect': 0, 'topo_correct': 0, 'topo_incorrect': 0}
    Q_ok_results = 0
    proteins = 0
    num = 0
    for name, entry in rostlab_db_dict.items():
        # if name != 'p11350':
        #     continue
        if float(entry['topo_string'].count('u') + entry['topo_string'].count('U'))/float(len(entry['topo_string'])) \
                > 0.2:
            continue
        if not -1 < num < 100:
            num += 1
            continue
        # print entry
        temp = HphobicityScore(name, entry['seq'], 'data_sets/rostlab_db/psipred/'+name+'.ss2', hydrophobicity_polyval, param_list)

        # print temp.topo_best
        topo_string = topo_string_rostlab_format(temp.topo_best, entry['seq'])

        ok_pred_tm = int(subprocess.Popen(['perl', './rostlab_evaluator.pl', entry['topo_string'], topo_string],
                                      stdout=subprocess.PIPE).stdout.read().split()[-1])
        pred_tm = len(temp.topo_best)
        obs_tm = len(re.findall('[1u20LU]h', entry['topo_string']))
        Q_ok_results += 1 if (ok_pred_tm/obs_tm == 1 and ok_pred_tm/pred_tm == 1) else 0
        proteins += 1
        print entry['topo_string']
        print topo_string
        dist = topo_string_distance(topo_string, entry['topo_string'])
        # print dist
        # print entry['topo_string']
        # print dist['aln']
        # print topo_string
        # results['tm_num_correct'] += 1 if dist['tm_agree'] else 0
        # results['tm_num_incorrect'] += 1 if not dist['tm_agree'] else 0
        # results['topo_correct'] += 1 if dist['topo_correct'] else 0
        # results['topo_incorrect'] += 1 if not dist['topo_correct'] else 0
        # print results
        num += 1
        # temp.plot_win_grades()
        results_writer(dist, entry, topo_string, temp, param_list, ok_pred_tm, pred_tm, obs_tm)
    # print results
    print Q_ok_results
    print 'Q_ok is', 100 * float(Q_ok_results)/float(proteins)


def results_writer(topo_score, entry_info, topo_string, hp_obj, param_list, ok_pred_tm, pred_tm, obs_tm):
    import matplotlib.pyplot as plt
    from time import strftime
    f = open('data_sets/rostlab_db/prediction/' + entry_info['name'] + '.prd', 'wr+')
    f.write('uniprot name\t%s\n' % entry_info['name'])
    f.write('PDB name    \t%s %s\n' % (entry_info['pdb'], entry_info['chain']))
    f.write('parameters: HP threshold %f MIN length %i psi helix %f psi res num %i\n' % (param_list[0], param_list[1],
                                                                                         param_list[2], param_list[3]))
    f.write('#TM observed %-3i\n#TM predicted %-3i\n' % (topo_score['tm_num_obs'], topo_score['tm_num_pre']))
    f.write('#TM agrees\n' if topo_score['tm_agree'] else '#TM DISAGREE\n')
    f.write('topo correct\n' if topo_score['topo_correct'] else 'topo INCORRECT\n')
    f.write('topo best energy        %f\n' % hp_obj.topo_best_val)
    f.write('topo second best energy %f\n' % hp_obj.topo_sec_best_val)
    f.write('topo confidence %f\n' % (hp_obj.topo_best_val / hp_obj.topo_sec_best_val))
    f.write('sequence %s\n' % entry_info['seq'])
    f.write('obs topo %s\n' % entry_info['topo_string'])
    f.write('aln topo %s\n' % topo_score['aln'])
    f.write('pre topo %s\n' % topo_string)
    f.write('\n\npredicted TM (Q) %i\n' % pred_tm)
    f.write('observed TM (Q) %i\n' % obs_tm)
    f.write('OK predicted TM %i\n' % ok_pred_tm)
    f.write('Qok or not ' + str(True) if (ok_pred_tm/obs_tm == 1 and ok_pred_tm/pred_tm == 1) else str(False))
    f.write('\n')
    f.write('prduced ' + strftime("%Y-%m-%d %H:%M:%S") + '\n')
    hp_obj.plot_win_grades()
    plt.savefig('data_sets/rostlab_db/prediction/' + entry_info['name'] + '.png')
    plt.close()
    f.close()


def topo_string_distance(tsp, tso):
    '''
    :param tsp: the predicted topo string
    :param tso: the observed topo string
    :return:number of disagreements (1/2/0 to h/H difference. disregarding U)
    '''
    score = {'tm_overlapp': 0, 'non_tm_overlapp': 0, 'pre_reminder': 0, 'obs_reminder': 0, 'tm_num_pre': 0,
             'tm_num_obs': 0, 'aln': ''}
    assert len(tsp) == len(tso), 'topo string length2 differ, ts1 %i, ts2 %i' % (len(tsp), len(tso))
    ### determine predicted N' and C' topology
    score['n_term_pre'] = 'in' if tsp[0] == '1' else 'out'
    score['c_term_pre'] = 'in' if tsp[-1] == '1' else 'out'
    ### determine observed N' and C' topology
    j = 0
    while tso[j].lower() == 'u': j += 1
    score['n_term_obs'] = 'in' if tso[j] == '1' else 'out' if tso[j] != '0' else 'unknown'
    j = len(tso)-1
    while tso[j].lower() == 'u': j -= 1
    score['c_term_obs'] = 'in' if tso[j] == '1' else 'out' if tso[j] != '0' else 'unknown'
    if tso[0] == 'H' or tso[0] == 'h':
        score['tm_num_obs'] += 1
    if tsp[0] == 'H' or tso[0] == 'h':
        score['tm_num_pre'] += 1
    ### compare topo-strings pos by pos
    for i in range(len(tsp)):
        if (tsp[i] == '1' or tsp[i] == '2' or tsp[i] == '0') and (tso[i] == 'H' or tso[i] == 'h'):
            score['obs_reminder'] += 1
            score['aln'] += '^'
        elif (tso[i] == '1' or tso[i] == '2' or tso[i] == '0') and (tsp[i] == 'H' or tsp[i] == 'h'):
            score['pre_reminder'] += 1
            score['aln'] += '#'
        elif (tsp[i] == 'h' or tsp[i] == 'H') and (tso[i] == 'h' or tso[i] == 'H'):
            score['tm_overlapp'] += 1
            score['aln'] += '|'
        elif (tsp[i] == '0' or tsp[i] == '1' or tsp[i] == '2') and \
                (tso[i] == '0' or tso[i] == '1' or tso[i] == '2' or tso[i] == 'u'):
            score['non_tm_overlapp'] += 1
            score['aln'] += '|' if tso[i] == tsp[i] else '\\'
        if i > 0:
            if (tsp[i] == 'H' or tsp[i] == 'h') and (tsp[i-1] == '0' or tsp[i-1] == '1' or tsp[i-1] == '2'):
                score['tm_num_pre'] += 1
            if (tso[i] == 'H' or tso[i] == 'h') and (tso[i-1] == '0' or tso[i-1] == '1' or tso[i-1] == '2'):
                score['tm_num_obs'] += 1
    score['tm_agree'] = True if score['tm_num_pre'] == score['tm_num_obs'] else False
    ### determine if topologies agree, or unknown
    score['topo_correct'] = True if score['n_term_obs'] == score['n_term_pre'] and \
        score['c_term_obs'] == score['c_term_pre'] else False \
        if (score['n_term_obs'] != 'unknown' and score['c_term_obs'] != 'unknown') else 'unknown'
    ### if either terminus of the observed topo is unknown, correctness is determined by the other side and TM num
    if score['topo_correct'] == 'unknown' and score['tm_num_obs'] == score['tm_num_pre']:
        if score['n_term_obs'] == score['n_term_pre'] or score['c_term_obs'] == score['c_term_pre']:
            score['topo_correct'] = True
        else:
            score['topo_correct'] = False
    return score


def topo_string_rostlab_format(topo, seq):
    '''
    :param topo:a topo (list of WinGrades describing a constructed topology)
    :param seq: the sequence
    :return:a string describing the topology in rostlab's format where 1:inside, 2: outdise H: TM helix
    '''
    global hydrophobicity_polyval
    topo_string = ''
    last_tm = WinGrade(0, 0, 'fwd', '', hydrophobicity_polyval)
    for tm in topo:
        topo_string += '1' * (tm.begin-last_tm.end) if tm.direction == 'fwd' else '2' * (tm.begin-last_tm.end)
        topo_string += 'H' * (tm.end - tm.begin)
        last_tm = tm
    topo_string += '2' * (len(seq)-last_tm.end) if last_tm.direction == 'fwd' else '1' * (len(seq)-last_tm.end)
    return topo_string


def MakeHydrophobicityGrade():
    '''
    :return: returns a dictionary of the polynom values for each residue
    '''
    global hydrophobicity_polyval
    # hydrophobicity_grade = open('Poly_Values.txt', 'r')
    # hydrophobicity_grade = open('poly_value_11.2.txt', 'r')
    # hydrophobicity_grade = open('poly_vals_23.2.txt', 'r')
    # hydrophobicity_grade = open('/home/labs/fleishman/jonathaw/membrane_prediciton/poly_vals_25.2.txt', 'r')
    hydrophobicity_grade = open('./poly_vals_25.2.txt', 'r')
    hydrophobicity_polyval = {}
    for line in hydrophobicity_grade:
        split = line.split()
        hydrophobicity_polyval[split[0]] = [float(n) for n in split[1:6]]
    hydrophobicity_grade.close()
    return hydrophobicity_polyval


def pymol_mark_segments(name, segments_set_set):
    """
    :param name: name of PDB file
    :param segments_set_set: a set of sets of two numbered lists, identifying the different types
    and ranges of segments to be colored
    :return: initiates a pymol session where the specified segments are colored
    """
    import subprocess
    from time import gmtime, strftime
    file_name = name + '_' + strftime("%H:%M", gmtime()) + '.pml'
    with open(file_name, 'wr+')as f:
        f.writelines('load ' + name.lower() + '.pdb,' + name + '\n')
        f.writelines('cmd.show("cartoon", "all")\n')
        seg_num = 1
        seg_set_num = 1
        colors = ['red', 'purple', 'blue', 'green', 'yellow', 'brown']
        for segment_set in segments_set_set:
            for segment in segment_set:
                f.writelines('select sele, %s and resi %i-%i' % (name, segment[0]+1, segment[1])+'\n')
                f.writelines('color ' + colors[seg_set_num] + ', sele\n')
                f.writelines('create seg%i, sele\n' % seg_num)
                seg_num += 1
            seg_set_num += 1
        f.writelines(['save ', name.lower()+'_TM_temp.pse\n'])
    subprocess.call(['/opt/local/bin/pymol', '-q', file_name])


def parsed_data_base_parser(num1=0, num2=1):
    """
    :param num1: start at entry #
    :param num2: finish at entry #
    :return: dictionary of entries information (PDB, uniprot, sequence, N' term orientation
    and starts and ends of TM helices
    """
    results = []
    i = 0
    with open('./database_new.txt', 'r') as f:
        for line in f.readlines()[0].split('\r'):
            if i < num1:
                i += 1
                continue
            line_split = line.split('\t')
            end_split = line_split[-1].split()
            starters = [int(a) for a in line_split[0].split()]
            enders = [int(a) for a in line_split[1].split()]
            assert type(end_split[3]) is str, "for %s PDB name is not string: %r" % \
                                              (results[-1]['uniprot'], results[-1]['pdb'])
            assert type(end_split[2]) is str, "for %s UNIPROT is not string: %r" % \
                                              (results[-1]['pdb'], results[-1]['uniprot'])
            assert type(line_split[2]) is str, "for %s seq is not string: %r" % \
                                               (results[-1]['uniprot'], results[-1]['seq'])
            assert type(end_split[0]) is str, "for %s orientation is not string: %r" % \
                                              (results[-1]['uniprot'], results[-1]['orientation'])
            assert type(starters) is list, "for %s begin list is not a list: %r" % \
                                           (results[-1]['uniprot'], results[-1]['begin'])
            assert type(enders) is list, "for %s begin list is not a list: %r" % \
                                         (results[-1]['uniprot'], results[-1]['end'])
            results.append({'pdb':end_split[3], 'seq': line_split[2], 'orientation': end_split[0], 'begin': starters,
                            'end': enders, 'uniprot': end_split[2]})
            if i == num2:
                break
            i += 1
    return results


def parse_v_db():
    # f = open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/V_Database.txt', 'r')
    f = open('./data_sets/V_Database.txt', 'r')
    resutls = {}
    i = 1
    for line in f:
        line_split = line.split()
        resutls[line_split[0]] = {'name': line_split[0], 'cterm': line_split[1], 'seq': line_split[4][:]}
        i += 1
    f.close()
    return resutls


def parse_rostlab_db():
    f = open('./data_sets/rostlab_db/rostlab_db.txt', 'r')
    cont_split = f.read().lower().split('>')
    results = {}
    for c in cont_split:
        if len(c.split()) > 3:
            continue
        split = c.split()
        name = split[0].split('|')[0]
        results[name] = {'name': name, 'seq': split[1].upper(), 'topo_string': split[2],
                         'pdb': split[0].split('|')[1].split(':')[0], 'chain': split[0].split('|')[1].split(':')[1][0]}
    f.close()
    return results


def ROC():
    import sys
    name = sys.argv[1]
    param_list = [float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])]
    vdb_dict = parse_v_db()
    global hydrophobicity_polyval
    hydrophobicity_polyval = MakeHydrophobicityGrade()
    for v_entry in vdb_dict.values():
        if v_entry['name'].lower() == name.lower():
            temp = HphobicityScore(v_entry['name'], v_entry['seq'], '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/VDB_noSRP/'+v_entry['name']+'.ss2', hydrophobicity_polyval, param_list)
            with open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/VDB_noSRP/ROC/'+name+'_'+
                    '_'.join(str(a) for a in param_list)+'.data', 'wr+') as f:
                f.writelines(v_entry['name']+'\n')
                f.writelines('database c_term\t'+v_entry['cterm']+'\n')
                f.writelines('predicted best c_term\t%s' % temp.best_c_term+'\n')
                f.writelines('predicted best grade\t%s' % temp.topo_best_val+'\n')
                f.writelines('predicted sec best c_term\t%s' % temp.sec_best_c_term+'\n')
                f.writelines('predicted sec best grade\t%s' % temp.topo_sec_best_val+'\n')


if __name__ == '__main__':
    main()
    # ROC()