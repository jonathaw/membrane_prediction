import numpy as np
import subprocess
from collections import Counter

def Range2Tups(ranger):
    results = []
    temp = [ranger[0]]
    for i, val in enumerate(ranger):
        if val-1 > ranger[i-1]:
            temp.append(ranger[i-1])
            results.append(temp[:])
            temp = [val]
    temp.append(ranger[-1])
    results.append(temp)
    return results

a = [1,2,3, 6,7,8,95, 22, 23, 24, 26]
print a.sort()
a.sort()
print a
# print Range2Tups([1,2,3, 6,7,8, 22, 23, 24, 26].sort())