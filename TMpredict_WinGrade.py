from WinGrade import *
from HphobicityScore import *


def main():
    global hydrophobicity_polyval
    hydrophobicity_polyval = MakeHydrophobicityGrade()
    # temp = HphobicityScore('temp', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    # temp = HphobicityScore('temp', 'YSYRFVWWAISTAAMLYILY')
    # temp = HphobicityScore('1E12', 'MSITSVPGVVDAGVLGAQSAAAVRENALLSSSLWVNVALAGIAILVFVYMGRTIRPGRPRLIWGATLMIPLVSISSYLGLLSGLTVGMIEMPAGHALAGEMVRSQWGRYLTWALSTPMILLALGLLADVDLGSLFTVIAADIGMCVTGLAAAMTTSALLFRWAFYAISCAFFVVVLSALVTDWAASASSAGTAEIFDTLRVLTVVLWLGYPIVWAVGVEGLALVQSVGVTSWAYSVLDVFAKYVFAFILLRWVANNERTVAVAGQTLGTMSSDD', hydrophobicity_polyval)
    # temp = HphobicityScore('1IWG', 'MPNFFIDRPIFAWVIAIIIMLAGGLAILKLPVAQYPTIAPPAVTISASYPGADAKTVQDTVTQVIEQNMNGIDNLMYMSSNSDSTGTVQITLTFESGTDADIAQVQVQNKLQLAMPLLPQEVQQQGVSVEKSSSSFLMVVGVINTDGTMTQEDISDYVAANMKDAISRTSGVGDVQLFGSQYAMRIWMNPNELNKFQLTPVDVITAIKAQNAQVAAGQLGGTPPVKGQQLNASIIAQTRLTSTEEFGKILLKVNQDGSRVLLRDVAKIELGGENYDIIAEFNGQPASGLGIKLATGANALDTAAAIRAELAKMEPFFPSGLKIVYPYDTTPFVKISIHEVVKTLVEAIILVFLVMYLFLQNFRATLIPTIAVPVVLLGTFAVLAAFGFSINTLTMFGMVLAIGLLVDDAIVVVENVERVMAEEGLPPKEATRKSMGQIQGALVGIAMVLSAVFVPMAFFGGSTGAIYRQFSITIVSAMALSVLVALILTPALCATMLKPIAKGDHGEGKKGFFGWFNRMFEKSTHHYTDSVGGILRSTGRYLVLYLIIVVGMAYLFVRLPSSFLPDEDQGVFMTMVQLPAGATQERTQKVLNEVTHYYLTKEKNNVESVFAVNGFGFAGRGQNTGIAFVSLKDWADRPGEENKVEAITMRATRAFSQIKDAMVFAFNLPAIVELGTATGFDFELIDQAGLGHEKLTQARNQLLAEAAKHPDMLTSVRPNGLEDTPQFKIDIDQEKAQALGVSINDINTTLGAAWGGSYVNDFIDRGRVKKVYVMSEAKYRMLPDDIGDWYVRAADGQMVPFSAFSSSRWEYGSPRLERYNGLPSMEILGQAAPGKSTGEAMELMEQLASKLPTGVGYDWTGMSYQERLSGNQAPSLYAISLIVVFLCLAALYESWSIPFSVMLVVPLGVIGALLAATFRGLTNDVYFQVGLLTTIGLSAKNAILIVEFAKDLMDKEGKGLIEATLDAVRMRLRPILMTSLAFILGVMPLVISTGAGSGAQNAVGTGVMGGMVTATVLAIFFVPVFFVVVRRRFSRKNEDIEHSHTVDHH')
    temp = HphobicityScore('1BRX', 'EAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSAGDGAAATS', hydrophobicity_polyval)

    # pymol_mark_segments(temp.name, [[[i.begin, i.end] for i in temp.minimas]])
    temp.plot_win_grades()


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
    '''
    :param name: name of PDB file
    :param segments_set_set: a set of sets of two numbered lists, identifying the different types
    and ranges of segments to be colored
    :return: initiates a pymol session where the specified segments are colored
    '''
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
                f.writelines('create seg%i, %s and resi %i-%i' % (seg_num, name, segment[0]+1, segment[1])+'\n')
                f.writelines('color ' + colors[seg_set_num] + ', seg' + str(seg_num) + '\n')
                seg_num += 1
            seg_set_num += 1
        f.writelines(['save ', name.lower()+'_TM_temp.pse\n'])
    subprocess.call(['/opt/local/bin/pymol', '-q', file_name])


if __name__ == '__main__':
    main()