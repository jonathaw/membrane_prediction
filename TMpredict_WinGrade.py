#!/usr/bin/env python2.7
# coding=utf-8
from WinGrade import *
from HphobicityScore import *


def main():
    # import subprocess
    # import re
    import os
    global hydrophobicity_polyval, args, param_list, tic
    import argparse
    import timeit

    parser = argparse.ArgumentParser()
    parser.add_argument('-hp_threshold', default=7.5, type=float, help='set the hp threshold for constructing the graph')#10
    parser.add_argument('-min_length', default=21, type=int, help='minimum window length') # 19
    # parser.add_argument('-psi_helix', default=0.2, type=float, help='no longer in use')#0.001
    # parser.add_argument('-psi_res_num', default=3, type=int, help='no longer in use')#4
    parser.add_argument('-mode', type=str, default='user', help='mode of run')
    parser.add_argument('-name', default=None, type=str, help='name of entry')
    # parser.add_argument('-known_tm_num', default=-100, type=int)
    # parser.add_argument('-c0', default=0.0, type=float)#0.5
    # parser.add_argument('-c1', default=0, type=float)#9.0 9.29
    # parser.add_argument('-c2', default=0, type=float)#-0.2 -0.645
    # parser.add_argument('-c3', default=0, type=float)#-0.006 0.00822
    parser.add_argument('-w', default=0, type=float, help='membrane deformation coeficent') ##0.011,  0.082 0.004
    parser.add_argument('-z_0', default=45, type=float, help='non-deformed membrane width') #  35 43
    parser.add_argument('-result_path', default=os.getcwd()+'/', help='path to write results to')
    parser.add_argument('-in_path', type=str, default=os.getcwd()+'/')
    parser.add_argument('-out_path', type=str, default=os.getcwd()+'/')
    parser.add_argument('-seq', default='', type=str, help='entry AA sequence')
    parser.add_argument('-with_msa', default=False, help='whether to use MSA or not')
    parser.add_argument('-msa_percentile', default=0, type=int, help='what MSA percentile to use')
    parser.add_argument('-with_cst', default=False, help='whether to use constraints')
    parser.add_argument('-cst_path', default=os.getcwd(), help='path to cst file')#+'/')
    parser.add_argument('-inc_max', default=10, type=int, help='maximal window increase')
    parser.add_argument('-fidelity', default=0, type=int, help='flanks on sides of tm_pos constraints')
    parser.add_argument('-msa_threshold', type=int, default=5)
    parser.add_argument('-db', default=None)
    parser.add_argument('-run_type', default='msa2plain')
    parser.add_argument('-ss2', default=None, type=str, help='name of ss2 file. if none is provided name.ss2 will be assumed')
    parser.add_argument('-create_html', default=True, help='whther to create an html')
    parser.add_argument('-with_sp', default=True, help='whether to check the TOPCONS output for spoctopus on signal peptide')
    args = vars(parser.parse_args())

    args['tic'] = timeit.default_timer()
    if args['ss2'] is None and args['mode'] != 'dG':
        args['ss2'] = args['name'].lower()+'.ss2'

    if args['with_msa'] == 'False':
        args['with_msa'] = False
    if args['create_html'] == 'True':
        args['create_html'] = True
    else:
        args['create_html'] = False
    # else:
    #     args['create_html'] = False
    if args['in_path'][-1] != '/':
        args['in_path'] += '/'

    if args['out_path'][-1] != '/':
        args['out_path'] += '/'

    if args['run_type'] in ['msa2plain', 'csts_msa2plain']:
        args['with_msa'] = True

    # import topdb_functions
    # hydrophobicity_polyval = MakeHydrophobicityGrade()
    if args['mode'] == 'ROC':
        rostlab_ROC(args)
    elif args['mode'] == 'single':
        args['name'] = args['name'].lower()
        process_single_protein(args['name'], args['result_path'])
    elif args['mode'] == 'dG':
        single_win_dG()
    elif args['mode'] == 'topo_VH':
        topo_VH()
    elif args['mode'] == 'ROC_by_single':
        ROC_rostlav_single_by_single()
    elif args['mode'] == 'new':
        if args['run_type'] == 'user_cst':
            args['with_cst'] = True
        args['with_msa'] = True
        args['original_name'] = args['name']
        args['name'] = args['name'].lower()
        process_single_new(args)
    elif args['mode'] == 'user':
        process_user(args)


def process_user(args):
    from ProcessEntry import create_topo_entry, process_entry
    import TMConstraint
    # args['db'] = None
    if args['with_cst'] or args['mode'] == 'csts_msa2plain':
        print args['with_cst']
        cst = TMConstraint.parse_cst(args['name'].lower(), args['in_path'])
    else:
        cst = TMConstraint.TMConstraint(args['name'])
    topo_entry = create_topo_entry(args['name'], args['seq'], args['in_path']+args['ss2'], args, cst, None,
                                   args['in_path'])
    process_entry(topo_entry, args['run_type'])


def process_single_new(args):
    from ProcessEntry import create_topo_entry, process_entry
    import TMConstraint
    import os

    if args['db'] == 'rost':
        rostlab_db_dict = parse_rostlab_db()
        entry = rostlab_db_dict[args['name'].lower()]
        ss2_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/psipred/'\
                   +args['name'].lower()+'.ss2'
        path_msa = '/home/labs/fleishman/jonathaw/membrane_prediction_DBs/BLASTs_9Aug/'
    elif args['db'] == 'vh':
        entry = topo_VH_parser(args['name'])
        ss2_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/VH_psipred/'\
                   +args['original_name']+'.ss2'
        if not os.path.isfile(ss2_path):
            for l in open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/VH_Data_Base_All_Sequences_No_SP_name_list.txt', 'r').read().split('\n'):
                if l.lower().rstrip() == args['name'].lower():
                    args['original_name'] = l.rstrip()
                    ss2_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/VH_psipred/'\
                               +args['original_name']+'.ss2'
                    break

        args['c_term_VH'] = entry['c_term_VH']
        # path_msa = '/home/labs/fleishman/jonathaw/membrane_prediction_DBs/BLAST_8Sep_VH/blast2fasta/'
        # path_msa = '/home/labs/fleishman/elazara/VH_MSA_60/blast2fasta/'
        path_msa = '/home/labs/fleishman/jonathaw/membrane_prediction_DBs/BLAST_8Sep_VH/blast2fasta/'
    else:
        assert "cant identify database (-db)"
    if args['with_cst']:
        entry_cst = TMConstraint.parse_cst(args['name'].lower(), args['cst_path'])
    else:
        entry_cst = TMConstraint.TMConstraint(args['name'])
    topo_entry = create_topo_entry(args['original_name'], entry['seq'], ss2_path, args, entry_cst, args['db'], path_msa)
    process_entry(topo_entry, args['run_type'])


def process_single_protein(name, path):
    import re
    import TMConstraint
    # path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/10overlap_uuu'
    # path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/production_28.4/'
    topc = spc_parser('/home/labs/fleishman/jonathaw/membrane_prediction_DBs/spoctopus_SPDB/'+name+'.spc')
    rostlab_db_dict = parse_rostlab_db()
    entry = rostlab_db_dict[name.lower()]
    if args['with_cst']:
        entry_cst = TMConstraint.parse_cst(name.lower(), args['cst_path'])
    else:
        entry_cst = TMConstraint.TMConstraint(args['name'])
    print entry
    if topc['topcons'].count('S') != 0:
        print 'topcons', topc['topcons']
        end_of_SP = [a for a in re.finditer('S*', topc['topcons']) if a != ''][0].end() - 1
        if end_of_SP == -1:
            end_of_SP = 0
        print 'end of SP', end_of_SP
        entry['seq_no_SP'] = 'u'*end_of_SP + entry['seq'][end_of_SP:]
        print 'no SP', entry['seq_no_SP']
        temp = HphobicityScore(name, entry['seq_no_SP'],
                        '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/psipred/'+name+'.ss2',
                        hydrophobicity_polyval, args, entry_cst)
        # topo_string = 'u'*end_of_SP + topo_string_rostlab_format(temp.topo_best, entry['seq_no_SP'])
        # sec_topo_string = 'u'*end_of_SP + topo_string_rostlab_format(temp.topo_sec_best, entry['seq_no_SP'])
        topo_string = topo_string_rostlab_format(temp.topo_best, entry['seq_no_SP'])
        sec_topo_string = topo_string_rostlab_format(temp.topo_sec_best, entry['seq_no_SP'])
    else:
        temp = HphobicityScore(name, entry['seq'],
                        '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/psipred/'+name+'.ss2',
                        hydrophobicity_polyval, args, entry_cst)
        topo_string = topo_string_rostlab_format(temp.topo_best, entry['seq'])
        sec_topo_string = topo_string_rostlab_format(temp.topo_sec_best, entry['seq'])
    print 'temp.topo_best', temp.topo_best
    print 'temp.topo_sec_best', temp.topo_sec_best
    print 'topo_string', topo_string
    # temp = HphobicityScore(name, entry['seq'],
    #                     '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/psipred/'+name+'.ss2',
    #                     hydrophobicity_polyval, args)
    # topo_string = topo_string_rostlab_format(temp.topo_best, entry['seq'])
    results_writer_skim(path, name, topo_string, sec_topo_string, temp.topo_best_val, temp.topo_sec_best_val, entry['seq'], entry_cst)
    # print 'Assaf, yoo my Boo!!! (ಠ‿ಠ)'
    # print entry['seq']
    # print entry['seq_no_SP']
    # print topc['seq']
    # print 'aaa', topo_string
    # print 'bbb', entry['pdbtm']
    # entry_results = {}
    # pred_tm = len(temp.topo_best)
    # entry_results['pdbtm'] = result_comparer(entry['pdbtm'], topo_string, pred_tm)
    # entry_results['opm'] = result_comparer(entry['opm'], topo_string, pred_tm)
    # overlap10 = {'pdbtm': result_comparer_10overlap(entry['pdbtm'], topo_string),
    #              'opm': result_comparer_10overlap(entry['pdbtm'], topo_string)}
    # print overlap10
    # results_writer(entry_results, entry, topo_string, temp, args, pred_tm, path, overlap10)
    # temp.plot_win_grades()


def ROC_rostlav_single_by_single():
    import os
    failed = []
    # path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/ROC_narrow/ROC_'\
    #        +str(args['c0'])+'_'+str(args['c1'])+'_'+str(args['c2'])\
    #        +'_'+str(args['c3'])+'/'
    path = os.getcwd()+'/ROC_'+str(args['c0'])+'_'+str(args['c1'])+'_'+str(args['c2'])+'_'+str(args['c3'])+'/'
    if not os.path.exists(path):
        os.makedirs(path)
    rostlab_db_dict = parse_rostlab_db()
    for name, entry in rostlab_db_dict.items():
        try:
            process_single_protein(name, path)
        except:
            failed.append(name)
            continue
    with open(path+'.roc') as o:
        o.writelines('num_failed %i\n' % len(failed))
        o.writelines('failed_list %r\n' % failed)


def results_writer_skim(path, name, pred_ts, sec_pred_ts, best_val, sec_best_val, seq, tmc):
    print 'writing.prd to', path+'/'+name+'.prd'
    with open(path+'/'+name+'.prd', 'wr+') as o:
        o.writelines('name %s\n' % name)
        try:
            o.writelines('pred_ts %s\n' % pred_ts)
            o.writelines('pred_sec_ts %s\n' % sec_pred_ts)
            o.writelines('seq %s\n' % seq)
            o.writelines('best_val %f\n' % best_val)
            o.writelines('sec_best_val %f\n' % sec_best_val)
        except:
            o.writelines('FAILED failed Failed !!!! :(\n')
            o.writelines('seq %s\n' % seq)
        for k, v in args.items():
            o.writelines('%s %r\n' % (k, v))
        o.writelines(str(tmc))


def archive_main():
    ### best results from 6.4 ROC: for pdbtm: hp_threshold=-3.0, min_length=19(use 18), psi_helix=0.2, psi_res_num=2 > ROC_-3.0_18_0.2_2
    ### best results from 6.4 ROC: for opm:   hp_threshold=-3.0, min_length=20(use 18), psi_helix=0.2, psi_res_num=1 > ROC_-3.0_18_0.2_1
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


def rostlab_ROC(param_list):
    import os
    # path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/ROC/ROC_'\
    #        +str(param_list['hp_threshold'])+'_'+str(param_list['min_length'])+'_'+str(param_list['psi_helix'])\
    #        +'_'+str(param_list['psi_res_num'])+'/'
    path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/ROC/ROC_'\
           +str(param_list['c0'])+'_'+str(param_list['c1'])+'_'+str(param_list['c2'])\
           +'_'+str(param_list['c3'])+'/'
    if not os.path.exists(path):
        os.makedirs(path)
    # path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/12.4Temp'
    rostlab_db_dict = parse_rostlab_db()
    Q_ok_results = {'opm': 0, 'pdbtm': 0}
    proteins, num, overlap10_ok = 0, 0, 0
    protein_names = []
    Q_ok, percent_topo_correct, overlap10 = {}, {}, {}
    num_topo_correct = {'opm': 0, 'pdbtm': 0}
    failed = []
    for name, entry in rostlab_db_dict.items():
        # if name != args['name'].lower():  # q99385
        #     continue
        # if float(entry['topo_string'].count('u') + entry['topo_string'].count('U'))/float(len(entry['topo_string'])) \
        #         > 0.2:
        #     print 'skipping entry due to UUUUUUU'
        #     continue
        # if not -1 < num < 150:
        #     num += 1
        #     continue
        print entry
        try:
            temp = HphobicityScore(name, entry['seq'], '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/psipred/'+name+'.ss2', hydrophobicity_polyval, param_list)
        except:
            failed.append(name)
            continue
        topo_string = topo_string_rostlab_format(temp.topo_best, entry['seq'])
        entry_results = {}
        pred_tm = len(temp.topo_best)
        entry_results['pdbtm'] = result_comparer(entry['pdbtm'], topo_string, pred_tm)
        entry_results['opm'] = result_comparer(entry['opm'], topo_string, pred_tm)
        num_topo_correct['pdbtm'] += 1 if entry_results['pdbtm']['topo_correct'] else 0
        num_topo_correct['opm'] += 1 if entry_results['opm']['topo_correct'] else 0
        Q_ok_results['pdbtm'] += 1 if entry_results['pdbtm']['Qok'] else 0
        Q_ok_results['opm'] += 1 if entry_results['opm']['Qok'] else 0
        overlap10['pdbtm'] = result_comparer_10overlap(entry['pdbtm'], topo_string)
        overlap10['opm'] = result_comparer_10overlap(entry['opm'], topo_string)
        overlap10['overlap10_both'] = any([overlap10['pdbtm']['10overlap'], overlap10['opm']['10overlap']])
        overlap10_ok += 1 if overlap10['overlap10_both'] else 0
        # temp.plot_win_grades()
        proteins += 1
        protein_names.append(entry['name'])
        results_writer(entry_results, entry, topo_string, temp, param_list, pred_tm, path, overlap10)
        num += 1
    Q_ok['pdbtm'] = 100.0 * float(Q_ok_results['pdbtm'])/float(proteins)
    Q_ok['opm'] = 100.0 * float(Q_ok_results['opm'])/float(proteins)
    percent_topo_correct['pdbtm'] = float(num_topo_correct['pdbtm']) * 100.0 / float(proteins)
    percent_topo_correct['opm'] = float(num_topo_correct['opm']) * 100.0 / float(proteins)
    percent_overlap10 = float(overlap10_ok)/float(proteins)
    # o = open('data_sets/rostlab_db/ROC/'+'_'.join(str(a) for a in param_list.values())+'.roc', 'wr+')
    ### for use when runing ROC
    # o = open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/ROC/roc'+
    #          str(param_list['hp_threshold'])+'_'+str(param_list['min_length'])+'_'+str(param_list['psi_helix'])
    #          +'_'+str(param_list['psi_res_num'])+'.roc', 'wr+')
    o = open(path+'ROC_results_'+str(param_list['c0'])+'_'+str(param_list['c1'])+'_'+str(param_list['c2'])+
             '_'+str(param_list['c3']), 'wr+')
    ### for use when NOT running ROC
    # o = open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/12.4Temp/temp.roc', 'wr+')
    o.write('hp_threshold %f\n' % param_list['hp_threshold'])
    o.write('min_length %i\n' % param_list['min_length'])
    o.write('psi_helix %f\n' % param_list['psi_helix'])
    o.write('psi_res_nume %i\n' % param_list['psi_res_num'])
    o.write('# proteins %i\n' % proteins)
    o.write('proteins: '+' '.join(a for a in protein_names)+'\n\n')
    for typer in ['pdbtm', 'opm']:
        o.write('Results for %s\n' % typer)
        o.write('Q_ok %f\n' % Q_ok[typer])
        o.write('# correct topo %i, precentage %f\n\n' % (num_topo_correct[typer], percent_topo_correct[typer]))
    o.write('failed %i proteins\n' % len(failed))
    o.write('failed: '+' '.join(failed)+'\n')
    o.write('Overlap10 results: %f\n' % percent_overlap10)
    o.close()


def topo_VH():
    from time import strftime
    import os
    vh_db = topo_VH_parser(args['name'])
    phobius = phobius_VH_parser(args['name'])
    ss2_path = '/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/VH_psipred/'+args['name']+'.ss2'
    hp_obj = HphobicityScore(vh_db['name'], vh_db['seq'], ss2_path, hydrophobicity_polyval, args)
    topo_string = topo_string_rostlab_format(hp_obj.topo_best, vh_db['seq'])
    pred_best_c_term = hp_obj.best_c_term
    pred_sec_best_c_term = hp_obj.sec_best_c_term
    # with open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/21July_VH_charges/'+args['name']+'.prd',
    #           'wr+') as o:
    print "I am printing here", os.getcwd()+'/'+args['name']+'.prd'
    with open(os.getcwd()+'/'+args['name']+'.prd',
              'wr+') as o:
        o.writelines('name %s\n' % args['name'])

        o.writelines('obs_c_term %s\n' % vh_db['c_term_VH'])
        o.writelines('phobius_c_term %s\n' % phobius['phobius_c_term'])
        o.writelines('best_c_term %s\n' % pred_best_c_term)
        o.writelines('best_val %f\n' % hp_obj.topo_best_val)
        o.writelines('best_tm_num %i\n' % len(hp_obj.topo_best))
        o.writelines('sec_best_c_term %s\n' % pred_sec_best_c_term)
        o.writelines('sec_best_val %f\n' % hp_obj.topo_sec_best_val)
        o.writelines('sec_best_tm_num %i\n' % len(hp_obj.topo_sec_best))
        o.writelines('best_sec_best_delta %f\n' % (hp_obj.topo_best_val-hp_obj.topo_sec_best_val))

        o.writelines('seq %s\n' % vh_db['seq'])
        o.writelines('pre %s\n' % topo_string)
        o.writelines('pred_correct %r\n' % (pred_best_c_term == vh_db['c_term_VH']))
        o.writelines('phobius_correct %r\n' % (phobius['phobius_c_term'] == vh_db['c_term_VH']))
        o.writelines('pred_phobius_agree %r\n' % (phobius['phobius_c_term'] == pred_best_c_term))
        o.write('produced ' + strftime("%Y-%m-%d %H:%M:%S") + '\n')


def phobius_VH_parser(name):
    """
    :param name: von-heijne database entry name
    :return: phobius results for the entry.
    """
    with open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/VH_Phobiuse_Topo.txt', 'r') as f:
        cont = f.read().split('\n')
    for line in cont:
        split = line.split()
        if len(split) < 4: continue
        if split[0].lower() == name.lower():
            return {'name': split[0].lower(), 'phobius_c_term': 'in' if split[3][-1] == 'i' else 'out'}


def topo_VH_parser(name):
    """
    :param name: von-heijne database entry name
    :return: dict of VH results for the name
    """
    with open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/VH_Data_Base_All_Sequences_No_SP.txt', 'r') as f:
        cont = f.read().split('\n')
    for line in cont:
        split = line.split()
        if len(split) < 5: continue
        if split[0].lower() == name.lower():
            return {'seq': split[4], 'c_term_VH': split[1], 'name': split[0].lower()}


def result_comparer_10overlap(obs_ts, pred_ts):
    import re
    result = {'10overlap': True}
    hhh = re.compile('[hH]*')
    obs_list = [(a.start(), a.end()) for a in hhh.finditer(obs_ts) if a.end()-a.start() > 1]
    pred_list = [(a.start(), a.end()) for a in hhh.finditer(pred_ts) if a.end()-a.start() > 1]
    result['obs_tm'] = len(obs_list)
    result['pred_tm'] = len(pred_list)

    for pred_seg in pred_list:
        if not any(segs_10overlap(pred_seg, a) for a in obs_list) and not seg_in_unknown(pred_seg, obs_ts):
            result['10overlap'] = False
            break
    return result


def seg_in_unknown(seg, obs_ts):
    return True if float(obs_ts[seg[0]:seg[1]+1].lower().count('u'))/float(seg[1]+1-seg[0]) > 0.98 else False


def segs_10overlap(seg1, seg2):
    res = 0
    seg2_list = [i for i in range(seg2[0], seg2[1]+1)]
    for i in range(seg1[0], seg1[1]+1):
        res += 1 if i in seg2_list else 0
    return True if res >= 10 else False


def result_comparer(obs_ts, pred_ts, pred_tm):
    import subprocess
    import re
    result = {}
    uuu_pred_ts = uuuu_pred_ts(obs_ts, pred_ts)  # wherever obs_ts has a u, uuu_pred_ts will too
    ok_pred_tm = int(subprocess.Popen(['perl', '/home/labs/fleishman/jonathaw/membrane_prediciton/rostlab_evaluator.pl',
                                       obs_ts, uuu_pred_ts], stdout=subprocess.PIPE).stdout.read().split()[-1])
    obs_tm = len(re.findall('[1u20LU]h', obs_ts))
    if pred_tm == ok_pred_tm and obs_tm == ok_pred_tm:
        result['Qok'] = True
    else:
        result['Qok'] = False
    result['tm_agree'] = True if obs_tm == pred_tm else False
    result['topo_comp'] = topo_string_distance(pred_ts, obs_ts)
    result['topo_correct'] = do_topos_agree_rostlab(pred_tm, obs_tm, obs_ts, pred_ts)
    result['obs_tm'] = obs_tm
    result['ok_pred_tm'] = ok_pred_tm
    return result


def uuuu_pred_ts(obs_ts, pred_ts):
    res = list(pred_ts)
    for i, aa in enumerate(obs_ts):
        if aa.lower() == 'u':
            res[i] = 'u'
    return ''.join(res)


def results_writer(entry_results, entry_info, topo_string, hp_obj, param_list, pred_tm, path, overlap10):
    import matplotlib.pyplot as plt
    from time import strftime
    # f = open('data_sets/rostlab_db/prediction/' + entry_info['name'] + '.prd', 'wr+')
    f = open(path+entry_info['name'] + '.prd', 'wr+')
    f.write('uniprot name\t%s\n' % entry_info['name'])
    f.write('PDB name    \t%s %s\n' % (entry_info['pdb'], entry_info['chain']))
    if param_list != None:
        f.write('hp_threshold %f\n' % param_list['hp_threshold'])
        f.write('min_length %i\n' % param_list['min_length'])
        f.write('psi_helix %f\n' % param_list['psi_helix'])
        f.write('psi_res_nume %i\n' % param_list['psi_res_num'])
    if hp_obj != None and hp_obj.topo_sec_best_val != 0:
        f.write('topo confidence %f\n' % (hp_obj.topo_best_val / hp_obj.topo_sec_best_val))
    for typer in ['pdbtm', 'opm']:
        f.write('Results for observed %s\n' % typer)
        f.write('#TM observed %-3i\n#TM predicted %-3i\n' % (entry_results[typer]['obs_tm'], pred_tm))
        f.write('#TM agrees\n' if entry_results[typer]['tm_agree'] else '#TM DISAGREE\n')
        f.write('topo correct\n' if entry_results[typer]['topo_correct'] else 'topo INCORRECT\n')
        if hp_obj != None:
            f.write('topo best energy        %f\n' % hp_obj.topo_best_val)
            f.write('topo second best energy %f\n' % hp_obj.topo_sec_best_val)
        else:
            f.write('topo best energy        %f\n' % 0.)
            f.write('topo second best energy %f\n' % 0.)
        f.write('sequence %s\n' % entry_info['seq'])
        f.write('obs topo %s\n' % entry_info[typer])
        f.write('aln topo %s\n' % entry_results[typer]['topo_comp']['aln'])
        f.write('pre topo %s\n' % topo_string)
        f.write('\n\npredicted TM (Q) %i\n' % pred_tm)
        f.write('observed TM (Q) %i\n' % entry_results[typer]['obs_tm'])
        f.write('OK predicted TM %i\n' % entry_results[typer]['ok_pred_tm'])
        if entry_results[typer]['obs_tm'] != 0 and pred_tm != 0:
            res = entry_results[typer]['ok_pred_tm']/entry_results[typer]['obs_tm'] == 1 and \
                  entry_results[typer]['ok_pred_tm']/pred_tm == 1
            f.write('Qok or not %r\n' % res)
        else:
            f.write('Qok cannot be calculated with pred TM %i and obs TM %i\n' % (pred_tm, entry_results[typer]['obs_tm']))
            if entry_results[typer]['obs_tm'] == pred_tm:
                f.write('Qok considered True, #TM equal\n')
        f.write('Topo is %r\n' % entry_results[typer]['topo_correct'])
        f.write('\n')
    for typer in ['pdbtm', 'opm']:
        f.write('by overlap 10 standard for %s:\n' % typer)
        f.write('observed tm %i\tpredicted tm %i\n' % (overlap10[typer]['obs_tm'], overlap10[typer]['pred_tm']))
        f.write('prediction OK by 10overlap: %r\n' % overlap10[typer]['10overlap'])
    over10_ok = any([overlap10['pdbtm']['10overlap'], overlap10['opm']['10overlap']])
    f.write('overlap10 is OK %r\n' % over10_ok)
    f.write('produced ' + strftime("%Y-%m-%d %H:%M:%S") + '\n')
    # hp_obj.plot_win_grades()
    # plt.savefig('data_sets/rostlab_db/prediction/' + entry_info['name'] + '.png')
    # plt.savefig('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/prediction/' + entry_info['name'] + '.png')
    plt.close()
    f.close()


def do_topos_agree_rostlab(pred_tm, obs_tm, obs_ts, pred_ts):
    j = 0
    while obs_ts[j].lower() == 'u': j += 1
    n_term_obs = 'in' if obs_ts[j] == '1' else 'out' if obs_ts[j] != '0' else 'unknown'
    j = len(obs_ts)-1
    while obs_ts[j].lower() == 'u': j -= 1
    c_term_obs = 'in' if obs_ts[j] == '1' else 'out' if obs_ts[j] != '0' else 'unknown'
    j = 0
    while pred_ts[j].lower() == 'u': j += 1
    n_term_pred = 'in' if pred_ts[j] == '1' else 'out' if (pred_ts[j] != '0' and pred_ts[j].lower() != 'u') else 'unknown'
    j = len(pred_ts)-1
    while pred_ts[j].lower() == 'u': j-=1
    c_term_pred = 'in' if pred_ts[j] == '1' else 'out' if (pred_ts[j] != '0' and pred_ts[j].lower() != 'u') else 'unknown'
    if pred_tm == obs_tm and n_term_obs == n_term_pred and c_term_obs == c_term_pred:
        return True
    if (n_term_pred == 'unknown' or c_term_pred == 'unknown' or n_term_obs == 'unknown' or c_term_obs == 'unknown'):
        if pred_tm == obs_tm and ((n_term_pred != 'unknown' and n_term_pred == c_term_pred) or
                                      (c_term_pred != 'unknown' and c_term_pred == c_term_obs)):
            return True
    return False


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
    last_tm = WinGrade(0, 0, 'fwd', '', hydrophobicity_polyval,
                       {k: v for k, v in args.items() if k in ['c0', 'c1', 'c2', 'c3', 'w', 'z_0']})
    for tm in topo:
        topo_string += '1' * (tm.begin-last_tm.end) if tm.direction == 'fwd' else '2' * (tm.begin-last_tm.end)
        topo_string += 'H' * (tm.end - tm.begin)
        last_tm = tm
    topo_string += '2' * (len(seq)-last_tm.end) if last_tm.direction == 'fwd' else '1' * (len(seq)-last_tm.end)
    return topo_string


# def MakeHydrophobicityGrade():
#     '''
#     :return: returns a dictionary of the polynom values for each residue
#     '''
#     global hydrophobicity_polyval
#     # hydrophobicity_grade = open('Poly_Values.txt', 'r')
#     # hydrophobicity_grade = open('poly_value_11.2.txt', 'r')
#     # hydrophobicity_grade = open('poly_vals_23.2.txt', 'r')
#     # hydrophobicity_grade = open('/home/labs/fleishman/jonathaw/membrane_prediciton/poly_vals_25.2.txt', 'r')
#     try:
#         # hydrophobicity_grade = open('/home/labs/fleishman/jonathaw/membrane_prediciton/polyval_21_5_15.txt', 'r')
#         hydrophobicity_grade = open('/home/labs/fleishman/jonathaw/membrane_prediciton/polyval_Kfix_14Jan2016.txt', 'r')
#     # else:
#     except:
#         # hydrophobicity_grade = open('/Volumes/labs/fleishman/jonathaw/membrane_prediciton/polyval_21_5_15.txt', 'r')
#         hydrophobicity_grade = open('/Volumes/labs/fleishman/jonathaw/membrane_prediciton/polyval_Kfix_14Jan2016.txt', 'r')
#     # hydrophobicity_grade = open('/Volumes/jonathaw-1/membrane_prediciton/poly_vals_25.2.txt', 'r')
#     # hydrophobicity_grade = open('./poly_vals_25.2.txt', 'r')
#     hydrophobicity_polyval = {}
#     for line in hydrophobicity_grade:
#         split = line.split()
#         hydrophobicity_polyval[split[0]] = [float(n) for n in split[1:6]]
#     hydrophobicity_grade.close()
#     return hydrophobicity_polyval


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
    # f = open('./data_sets/rostlab_db/rostlab_db.txt', 'r')
    f = open('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/rostlab_db.txt', 'r')
    cont_split = f.read().lower().split('>')
    results = {}
    for c in cont_split:
        if len(c.split()) > 3:
            continue
        split = c.split()
        name = split[0].split('|')[0]
        tech = split[0].split('|')[1].split('_')[1]
        if name not in results.keys():
            results[name] = {'name': name, 'seq': split[1].upper(), 'topo_string': split[2],
                             'pdb': split[0].split('|')[1].split(':')[0], 'chain': split[0].split('|')[1].split(':')[1][0],
                             tech: split[2]}
        else:
            results[name][tech] = split[2]
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


def single_win_dG():
    from ProcessEntry import MakeHydrophobicityGrade

    temp = WinGrade(0, len(args['seq']), 'fwd', args['seq'], MakeHydrophobicityGrade(),
                    {'w': args['w'], 'z_0': args['z_0']})
    print temp.grade


def topcons_parser(file_name):
    result = {}
    with open(file_name, 'r') as f:
        cont = f.read().split('\n')
    for i, line in enumerate(cont):
        split = line.split()
        if split == [] or split[0][0] == '#' or len(split) < 2: continue
        if split[1] == 'name:': result['name'] = split[2]
        elif split[0] == 'TOPCONS': result['topcons'] = cont[i+1]
        elif split[0] == 'OCTOPUS': result['octopus'] = cont[i+1]
        elif split[0] == 'Philius': result['philius'] = cont[i+1]
        elif split[0] == 'PolyPhobius': result['polyphobius'] = cont[i+1]
        elif split[0] == 'SCAMPI': result['scampi'] = cont[i+1]
        elif split[0] == 'SPOCTOPUS': result['spoctopus'] = cont[i+1]
    return result


def spc_parser(file_name):
    result = {}
    with open(file_name, 'r') as f:
        cont = f.read().split('\n')
    for line in cont:
        split = line.split()
        if split == []: continue
        result[split[0]] = split[1]
    return result


if __name__ == '__main__':
    main()
    # ROC()