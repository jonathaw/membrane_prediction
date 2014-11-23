def PsiReader(uniprot):
    import re
    ss2_file = open('../psipred/sw_fastas/' + uniprot + '.ss2')
    line_re = re.compile('^\s*([0-9]*)\s*([A-Z]*)\s*([A-Z]*)\s*(0\.[0-9]*)\s*(0\.[0-9]*)\s*(0\.[0-9]*)')
    result = []
    for line in ss2_file:
        if line_re.search(line):
            if line_re.search(line).group(3) == 'C':
                temp = line_re.search(line).group(4)
            elif line_re.search(line).group(3) == 'H':
                temp = line_re.search(line).group(5)
            elif line_re.search(line).group(3) == 'E':
                temp = line_re.search(line).group(6)
            result.append((line_re.search(line).group(3), float(temp)))
    return result


