from WinGrade import *
from HphobicityScore import *


def main():
    global hydrophobicity_polyval
    hydrophobicity_polyval = MakeHydrophobicityGrade()
    # temp = HphobicityScore('temp', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    # temp = HphobicityScore('temp', 'YSYRFVWWAISTAAMLYILY')
    # temp = HphobicityScore('1E12', 'MSITSVPGVVDAGVLGAQSAAAVRENALLSSSLWVNVALAGIAILVFVYMGRTIRPGRPRLIWGATLMIPLVSISSYLGLLSGLTVGMIEMPAGHALAGEMVRSQWGRYLTWALSTPMILLALGLLADVDLGSLFTVIAADIGMCVTGLAAAMTTSALLFRWAFYAISCAFFVVVLSALVTDWAASASSAGTAEIFDTLRVLTVVLWLGYPIVWAVGVEGLALVQSVGVTSWAYSVLDVFAKYVFAFILLRWVANNERTVAVAGQTLGTMSSDD', hydrophobicity_polyval)
    # temp = HphobicityScore('1IWG', 'MPNFFIDRPIFAWVIAIIIMLAGGLAILKLPVAQYPTIAPPAVTISASYPGADAKTVQDTVTQVIEQNMNGIDNLMYMSSNSDSTGTVQITLTFESGTDADIAQVQVQNKLQLAMPLLPQEVQQQGVSVEKSSSSFLMVVGVINTDGTMTQEDISDYVAANMKDAISRTSGVGDVQLFGSQYAMRIWMNPNELNKFQLTPVDVITAIKAQNAQVAAGQLGGTPPVKGQQLNASIIAQTRLTSTEEFGKILLKVNQDGSRVLLRDVAKIELGGENYDIIAEFNGQPASGLGIKLATGANALDTAAAIRAELAKMEPFFPSGLKIVYPYDTTPFVKISIHEVVKTLVEAIILVFLVMYLFLQNFRATLIPTIAVPVVLLGTFAVLAAFGFSINTLTMFGMVLAIGLLVDDAIVVVENVERVMAEEGLPPKEATRKSMGQIQGALVGIAMVLSAVFVPMAFFGGSTGAIYRQFSITIVSAMALSVLVALILTPALCATMLKPIAKGDHGEGKKGFFGWFNRMFEKSTHHYTDSVGGILRSTGRYLVLYLIIVVGMAYLFVRLPSSFLPDEDQGVFMTMVQLPAGATQERTQKVLNEVTHYYLTKEKNNVESVFAVNGFGFAGRGQNTGIAFVSLKDWADRPGEENKVEAITMRATRAFSQIKDAMVFAFNLPAIVELGTATGFDFELIDQAGLGHEKLTQARNQLLAEAAKHPDMLTSVRPNGLEDTPQFKIDIDQEKAQALGVSINDINTTLGAAWGGSYVNDFIDRGRVKKVYVMSEAKYRMLPDDIGDWYVRAADGQMVPFSAFSSSRWEYGSPRLERYNGLPSMEILGQAAPGKSTGEAMELMEQLASKLPTGVGYDWTGMSYQERLSGNQAPSLYAISLIVVFLCLAALYESWSIPFSVMLVVPLGVIGALLAATFRGLTNDVYFQVGLLTTIGLSAKNAILIVEFAKDLMDKEGKGLIEATLDAVRMRLRPILMTSLAFILGVMPLVISTGAGSGAQNAVGTGVMGGMVTATVLAIFFVPVFFVVVRRRFSRKNEDIEHSHTVDHH', hydrophobicity_polyval)
    # temp = HphobicityScore('1BRX', 'EAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSAGDGAAATS', hydrophobicity_polyval)

    # pymol_mark_segments(temp.name, [[[i.begin, i.end] for i in temp.topo_minimas]])
    # temp.plot_win_grades()

    # db_entries = parsed_data_base_parser(25,26)
    # temp = HphobicityScore(db_entries[1]['pdb'], db_entries[1]['seq'], db_entries[1]['uniprot'], hydrophobicity_polyval)
    # print temp
    # temp.plot_win_grades()
    # temp.plot_energy_landscape()

    '''
    parses a range of SW entries, and prints the topology predcition reult
    '''
    db_entries = parsed_data_base_parser(0,50)
    topo_predict_score = {'good': 0, 'bad': 0}
    for protein in db_entries:
        temp = HphobicityScore(protein['pdb'], protein['seq'], protein['uniprot'], hydrophobicity_polyval)
        if temp.n_term_orient == 'rev' and protein['orientation'] == 'out':
            topo_predict_score['good'] += 1
            print 'was correct', temp.n_term_orient, protein['orientation'], topo_predict_score
        else:
            topo_predict_score['bad'] += 1
            print 'was wrong', temp.n_term_orient, protein['orientation'], topo_predict_score
    print 'prediction results:', topo_predict_score


def MakeHydrophobicityGrade():
    '''
    :return: returns a dictionary of the polynom values for each residue
    '''
    global hydrophobicity_polyval
    # hydrophobicity_grade = open('Poly_Values.txt', 'r')
    hydrophobicity_polyval = {}
    hydrophobicity_grade = open('poly_value_11.2.txt', 'r')
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
            results.append({'pdb':end_split[3], 'seq': line_split[2], 'orientation': end_split[0], 'begin': starters,
                            'end': enders, 'uniprot': end_split[2]})
            assert type(results[-1]['pdb']) is str, "PDB name is not string: %r" % results[-1]['pdb']
            assert type(results[-1]['uniprot']) is str, "UNIPROT is not string: %r" % results[-1]['uniprot']
            assert type(results[-1]['seq']) is str, "seq is not string: %r" % results[-1]['seq']
            assert type(results[-1]['orientation']) is str, "orientation is not string: %r" % results[-1]['orientation']
            assert type(results[-1]['begin']) is list, "begin list is not a list: %r" % results[-1]['begin']
            assert type(results[-1]['end']) is list, "begin list is not a list: %r" % results[-1]['end']
            if i == num2:
                break
            i += 1
    return results


if __name__ == '__main__':
    main()