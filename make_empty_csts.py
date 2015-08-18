"""
a script to create empty (all None) .cst files for names with no .prd file
"""

def main():
    import os, re
    all_names = get_all_names()
    available_preds = [x[:-4] for x in os.listdir(os.getcwd()) if re.match('.*prd', x)]
    for name in all_names:
        if name not in available_preds:
            with open(name+'.cst', 'wr+') as fout:
                fout.write('name %s\n' % name)
                fout.write('tm_num None\n')
                fout.write('tm_pos None\n')
                fout.write('tm_pos_fidelity None\n')
                fout.write('c_term None\n')
                fout.write('n_term None\n')
                fout.write('non_tm_pos None\n')

def get_all_names():
    with open('/home/labs/fleishman/elazara/TM_benchmark/rostalb_uniprots_list.txt', 'r') as f:
        cont = f.read().split('\n')
    res = []
    for l in cont:
        if l != '':
            res.append(l.lower())
    return res

if __name__ == '__main__':
    main()