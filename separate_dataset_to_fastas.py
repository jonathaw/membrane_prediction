with open('data_sets/top_ap_TOPDB.txt', 'r') as f:
    cont = f.read().split('>')
db_dict = [{'name': i.split('\n')[0], 'seq': i.split('\n')[1], 'top': i.split('\n')[2]}
           for i in cont if len(i) > 1]

for entry in db_dict:
    f = open('data_sets/fastas/' + entry['name'] + '.fasta', 'wr+')
    f.writelines('>'+entry['name'] + '\n')
    f.writelines(entry['seq'] + '\n')
    f.close()