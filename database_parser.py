import re
import csv

pdb_re = re.compile('^.*<pdbCode>(.*)</pdbCode>')
uniprot_re = re.compile('^.*<uniprotNumber>(.*)</uniprotNumber>')
term_re = re.compile('^.*<nTerminal>(.*)</nTerminal>')
seq_re = re.compile('^.*<sequence>(.*)</sequence>')
begin_re = re.compile('^.*<beginIndex>(.*)</beginIndex>')
end_re = re.compile('^.*<endIndex>(.*)</endIndex>')
mptopp_re = re.compile('^.*</mptopoProtein>')

database = open('/Users/jonathan/Documents/membrane_prediciton_data/Topo_DATA_ALL_new.xml', 'r')
# database = open('/Users/jonathan/Documents/membrane_prediciton_data/temp_SW_single_seq.xml', 'r')

def SWDB_parser_prediciton(num):
    result_dict = {}
    temp_dict = {}
    temp_dict['uniprot'] = ''
    temp_dict['term'] = ''
    temp_dict['seq'] = ''
    temp_dict['pdb'] = []
    temp_dict['begin'] = []
    temp_dict['end'] = []
    temp_dict['seq_length'] = ''

    i = 0
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
            temp_dict['seq_length'] = len(temp_dict['seq'])
            result_dict[temp_dict['uniprot']] = temp_dict
            temp_dict = {}
            temp_dict['pdb'] = []
            temp_dict['begin'] = []
            temp_dict['end'] = []
            i += 1
        if i >= num:
            return result_dict


def SWDB_parser_prediciton_by_name(name):
    result_dict = {}
    temp_dict = {}
    temp_dict['uniprot'] = ''
    temp_dict['term'] = ''
    temp_dict['seq'] = ''
    temp_dict['pdb'] = []
    temp_dict['begin'] = []
    temp_dict['end'] = []
    temp_dict['seq_length'] = ''

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
            temp_dict['seq_length'] = len(temp_dict['seq'])
            result_dict[temp_dict['uniprot']] = temp_dict
            if [name.upper() is y.upper() for y in temp_dict['pdb']]:
                a = {}
                a[temp_dict['uniprot']] = temp_dict
                return a
            temp_dict = {}
            temp_dict['pdb'] = []
            temp_dict['begin'] = []
            temp_dict['end'] = []


def SW_parser_for_csv():
    output = open('/Users/jonathan/Documents/membrane_prediciton_data/database_new.csv', 'wa+')
    csv_writer = csv.writer(output)

    temp_dict = {}
    temp_dict['uniprot'] = ''
    temp_dict['term'] = ''
    temp_dict['seq'] = ''
    temp_dict['pdb'] = []
    temp_dict['begin'] = []
    temp_dict['end'] = []
    temp_dict['seq_length'] = ''

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
            temp_dict['seq_length'] = len(temp_dict['seq'])
            csv_writer.writerow(temp_dict.values())
            temp_dict = {}
            temp_dict['pdb'] = []
            temp_dict['begin'] = []
            temp_dict['end'] = []
