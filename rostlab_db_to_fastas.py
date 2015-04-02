f = open('./data_sets/rostlab_db/rostlab_db.txt', 'r')
cont_split = f.read().split('>')
for c in cont_split:
    if len(c.split()) > 3:
        continue
    name = c.split()[0].split('|')[0]
    with open('./data_sets/rostlab_db/fastas/' + name.lower() + '.fasta', 'wr+') as j:
        j.write('>' + name.lower() + '\n')
        j.write(c.split()[1] + '\n')