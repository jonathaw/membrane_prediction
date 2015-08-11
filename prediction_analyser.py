'''
scripts to analyse and compare prediction results with the databases
'''


def pred_parser(name):
    with open(name, 'r') as f:
        cont = f.read().split('\n')
    result = {}
    for line in cont:
        spt = line.split()
        if len(spt) == 0: continue
        if spt[1][0] == "'":
            spt[1] = spt[1][1:-1]
        try:
            result[spt[0]] = float(spt[1])
        except:
            result[spt[0]] = spt[1]
    return result



if __name__ == '__main__':
    # compare_strings(pred[])
    pred_parser('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/23July_msa_memb_def/ROC_10/p40719.prd')
    pass