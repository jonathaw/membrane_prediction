def main():
    import os, re
    topcons_path = '/home/labs/fleishman/jonathaw/membrane_topcons/topo_VH_topcons/all_results/'
    # topcons_path = '/home/labs/fleishman/jonathaw/membrane_topcons/all_results/'
    file_list = [x for x in os.listdir(topcons_path)
                 if re.match('.*\.txt', x)]
    for file_i in file_list:
        topc = topcons_parser(topcons_path+file_i)
        with open('/home/labs/fleishman/jonathaw/membrane_prediction_DBs/spoctopus_SPDB/vh_db/'+topc['name']+'.spc', 'wr+') as o:
            for k, v in topc.items():
                o.writelines('%s %s\n' % (k, v))



def topcons_parser(file_name):
    result = {}
    with open(file_name, 'r') as f:
        cont = f.read().split('\n')
    for i, line in enumerate(cont):
        split = line.split()
        if split == []: continue
        if split[0] == 'Sequence:': result['seq'] = cont[i+1]
        if split[0][0] == '#' or len(split) < 2: continue
        if split[1] == 'name:': result['name'] = split[2]
        elif split[0] == 'TOPCONS': result['topcons'] = cont[i+1]
        elif split[0] == 'OCTOPUS': result['octopus'] = cont[i+1]
        elif split[0] == 'Philius': result['philius'] = cont[i+1]
        elif split[0] == 'PolyPhobius': result['polyphobius'] = cont[i+1]
        elif split[0] == 'SCAMPI': result['scampi'] = cont[i+1]
        elif split[0] == 'SPOCTOPUS': result['spoctopus'] = cont[i+1]
    return result



if __name__ == '__main__':
    main()