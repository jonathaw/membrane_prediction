from TMpredict_WinGrade import parse_rostlab_db
ttt = parse_rostlab_db()
# print ttt
print len(ttt)
for k, v in ttt.items():
    if 'pdbtm' not in v.keys() or 'opm' not in v.keys():
        print k