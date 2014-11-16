import re
import csv

uniprot_re = re.compile('^.*<uniprotNumber>(.*)</uniprotNumber>')
pdb_re = re.compile('^.*<pdbCode>(.*)</pdbCode>')
term_re = re.compile('^.*<nTerminal>(.*)</nTerminal>')
seq_re = re.compile('^.*<sequence>(.*)</sequence>')
begin_re = re.compile('^.*<beginIndex>(.*)</beginIndex>')
end_re = re.compile('^.*<endIndex>(.*)</endIndex>')
mptopp_re = re.compile('^.*</mptopoProtein>')

main_list = []

database = open('/Users/jonathan/Documents/membrane_prediciton/Topo_DATA_ALL_new.xml', 'r')
output = open('/Users/jonathan/Documents/membrane_prediciton/database.csv', 'wa+')
csv_writer = csv.writer(output)

temp_dict = {}
temp_dict['uniprot'] = ''
temp_dict['term'] = ''
temp_dict['seq'] = ''
temp_dict['pdb'] = []
temp_dict['begin'] = []
temp_dict['end'] = []

csv_writer.writerow(temp_dict.keys())

for line in database:
    if uniprot_re.search(line):
        temp_dict['uniprot'] = uniprot_re.search(line).group(1)
    if pdb_re.search(line):
        temp_dict['pdb'].append(pdb_re.search(line).group(1))
    if term_re.search(line):
        temp_dict['term'] = term_re.search(line).group(1)
    if seq_re.search(line):
        temp_dict['seq'] = seq_re.search(line).group(1)
    if begin_re.search(line):
        temp_dict['begin'].append(begin_re.search(line).group(1))
    if end_re.search(line):
        temp_dict['end'].append(end_re.search(line).group(1))
    if mptopp_re.search(line):
        print temp_dict, '\n\n'
        csv_writer.writerow(temp_dict.values())
        temp_dict = {}
        temp_dict['pdb'] = []
        temp_dict['begin'] = []
        temp_dict['end'] = []
