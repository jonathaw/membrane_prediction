#!/usr/bin/env python
def pymol_mark_segments(name, chain, pred_ts, obs_ts, seq, tech):
    """
    :param name: name of PDB file
    :param segments_set_set: a set of sets of two numbered lists, identifying the different types
    and ranges of segments to be colored
    :return: initiates a pymol session where the specified segments are colored
    """
    # import subprocess
    file_name = args['name'].lower()+'_'+tech+'.pml'
    obs_hp_seqs = ts2hp_seq(seq, obs_ts)
    pred_hp_seqs = ts2hp_seq(seq, pred_ts)
    overlap_hp_seqs = two_ts2hp_seq(seq, pred_ts, obs_ts)
    with open(file_name, 'wr+')as f:
        f.writelines('load ' + name.lower() + '.pdb,' + name + '\n')
        f.writelines('load ' + args['name'].lower() + '.pdb,' + name + '\n')
        f.writelines('import findseq\n')
        # f.writelines('remove %s and not chain %s' % (name, chain))
        f.writelines('cmd.show("cartoon", "all")\n')
        for obspredover, segs, col in zip(['obs', 'pred', 'overlap'], [obs_hp_seqs, pred_hp_seqs, overlap_hp_seqs],
                                          ['blue', 'red', 'purple']):
            for i, seg in enumerate(segs):
                f.writelines('findseq %s, %s, temp\n' % (seg, name))
                f.writelines('create %s_%i, temp\n' % (obspredover, i))
                f.write('delete temp\n')
                f.writelines('color %s, %s_%i\n' % (col, obspredover, i))
        # f.writelines(['save ', name.lower()+'_TM_temp.pse\n'])
    # subprocess.call(['/opt/local/bin/pymol', '-q', file_name])


def ts2hp_seq(seq, ts):
    '''
    :param seq: AA seq
    :param ts: topo string of prediction/observation
    :return: list of AA sequences of the TMs in ts
    '''
    import re
    hp_seqs = []
    hhh = re.compile('[hH]*')
    for it in hhh.finditer(ts):
        if it.end() - it.start() > 1:
            hp_seqs.append(seq[it.start():it.end()])
    return hp_seqs


def two_ts2hp_seq(seq, ts1, ts2):
    '''
    :param seq: AA seq
    :param ts1: topo string observation/prediction
    :param ts2: topo string observation/prediction
    :return: list of AA sequences of TMHs overlapping between ts1 and ts2
    '''
    import re
    hp_seqs = []
    segs = []
    hhh = re.compile('[hH]*')
    ts1_list = [(a.start(), a.end()) for a in hhh.finditer(ts1) if a.end()-a.start() > 1]
    ts2_list = [(a.start(), a.end()) for a in hhh.finditer(ts2) if a.end()-a.start() > 1]
    for s1 in ts1_list:
        for s2 in ts2_list:
            segs.append(tuple_tuple_overlap(s1, s2))
    segs = [a for a in segs if a is not None]
    for seg in segs:
        hp_seqs.append(''.join(seq[a] for a in seg))
    return hp_seqs


def tuple_tuple_overlap(tup1, tup2):
    '''
    :param tup1: tuple of start and end of segment
    :param tup2: tuple of start and end of segment
    :return: list of overlapping positions, or None
    '''
    result = []
    for i in range(tup1[0], tup1[1]):
        if tup2[0] <= i < tup2[1]:
          result.append(i)
    return result if result != [] else None


def TMpredict_reader(file_name):
    '''
    :param file_name: file name and path to a prediction result (.prd)
    :return: dict with files details
    '''
    result = {}
    print file_name
    with open(file_name, 'r') as f:
        cont = f.read().split('\n')
    for item in cont:
        split = item.split()
        if len(split) > 1:
            result[split[0]] = split[1]
    return result

if __name__ == '__main__':
    import argparse
    import os
    from TMpredict_WinGrade import parse_rostlab_db
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', type=str)
    parser.add_argument('-path', default=os.getcwd(), type=str)
    parser.add_argument('-tech', default='opm', type=str)
    args = vars(parser.parse_args())
    # pymol_mark_segments('4k1c', 'a', '222222222222222222222222222222222222222HHHHHHHHHHHHHHHHHHH1HHHHHHHHHHHHHHHHHHHHHH22222222222222HHHHHHHHHHHHHHHHHHHHH11HHHHHHHHHHHHHHHHHHHHHH22222222222222222HHHHHHHHHHHHHHHHHHHHHHHHH1111111111111111HHHHHHHHHHHHHHHHHHHHH2222222222222222222222222222HHHHHHHHHHHHHHHHHHHHHHHHH111111111111111111111111111111111111111111111111111111111111111111111111111111111111HHHHHHHHHHHHHHHHHH222222HHHHHHHHHHHHHHHHHH1111111111111',
    #                     'u111111111111111111111111111111111hhhhhhhhhhhhhhhhhh222222222222hhhhhhhhhhhhhhhhhhhhhh111111hhhhhhhhhhhhhhhhhhhh2222222222hhhhhhhhhhhhhhhhhhhhhh111111111111111hhhhhhhhhhhhhhhhh2222222222222222222222222hhhhhhhhhhhhhhhhhh1111111111111111111111111111111hhhhhhhhhhhhhhhhhhh222222222222222222222hhhhhhhhhhhhhhhhhh1111111111hhhhhhhhhhhhhhhhhhhhhh2222222222222222hhhhhhhhhhhhhhhhh111hhhhhhhhhhhhhhhhhh22222222222222222',
    #                     'MDATTPLLTVANSHPARNPKHTAWRAAVYDLQYILKASPLNFLLVFVPLGLIWGHFQLSHTLTFLFNFLAIIPLAAILANATEELADKAGNTIGGLLNATFGNAVELIVSIIALKKGQVRIVQASMLGSLLSNLLLVLGLCFIFGGYNRVQQTFNQTAAQTMSSLLAIACASLLIPAAFRATLPHGKEDHFIDGKILELSRGTSIVILIVYVLFLYFQLGSHHALFEQQEEETDEVMSTISRNPHHSLSVKSSLVILLGTTVIISFCADFLVGTIDNVVESTGLSKTFIGLIVIPIVGNAAEHVTSVLVAMKDKMDLALGVAIGSSLQVALFVTPFMVLVGWMIDVPMTLNFSTFETATLFIAVFLSNYLILDGESNWLEGVMSLAMYILIAMAFFYYPDEKTLDSIGNSL')
    # entry = TMpredict_reader('/home/labs/fleishman/jonathaw/membrane_prediction_DBs/ROC_6.4.2015/ROC_-3.0_18_0.2_2/p00423.prd')
    entry = TMpredict_reader(args['path']+'/'+args['name']+'.prd')
    # print entry
    rostlab_data = parse_rostlab_db()[args['name']]
    # print 'aaa', rostlab_data
    pymol_mark_segments(rostlab_data['pdb'], rostlab_data['chain'], entry['pred_ts'], rostlab_data[args['tech']],
                        rostlab_data['seq'], args['tech'])