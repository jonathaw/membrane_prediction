def is_tmh(ind, ts):
    res = 0
    for i in range(ind, ind+20):
        if ts[i].lower() == 'h':
            res += 1
    return res / 20 >= 0.9

def is_not_helical(pos, psi):
    import numpy as np
    return False if (np.mean([psi[a]['e'] for a in range(pos[0], pos[1])]) <= 0.3 and
                    np.mean([psi[a]['c'] for a in range(pos[0], pos[1])]) <= 0.48 and
                    np.mean([psi[a]['h'] for a in range(pos[0], pos[1])]) >= 0.3) else True

from psipred_vs_mm_nomm import parse_psipred, psipred_avg
from TMpredict_WinGrade import parse_rostlab_db
import matplotlib.pyplot as plt
rost_db = parse_rostlab_db()
pdb_name = 'p02722'
psi = parse_psipred(pdb_name)
avg_b = []
avg_c = []
avg_h = []
indices = []
tmh = []
passed = []
for i in range(len(psi)-20):

    avg_b.append(psipred_avg(range(i+1, i+21), psi, 'e'))
    avg_c.append(psipred_avg(range(i+1, i+21), psi, 'c'))
    avg_h.append(psipred_avg(range(i+1, i+21), psi, 'h'))

    tmh.append(-0.1 if is_tmh(i, rost_db[pdb_name]['pdbtm']) else None)
    passed.append(-0.5 if not is_not_helical([i+1, i+21], psi) else None)

    indices.append(i)

plt.scatter(indices, avg_b, color='k')
plt.scatter(indices, avg_c, color='b')
plt.scatter(indices, avg_h, color='r')
plt.scatter(indices, tmh, color='m')
plt.scatter(indices, passed, color='g')
plt.show()

