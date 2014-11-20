import numpy as np
import subprocess
from collections import Counter

# list_a = [1, 2, 3, 4, 5]
# list_b = [3, 4, 6]
#
# def ListInList(list_a, list_b):
#     for a in list_a:
#         if a in list_b:
#             sumer += 1



a = [3,4,5,5,5,6]
b = [1,3,4,4,5,5,6,7]

a_multiset = Counter(a)
b_multiset = Counter(b)

overlap = list((a_multiset & b_multiset).elements())
a_remainder = list((a_multiset - b_multiset).elements())
b_remainder = list((b_multiset - a_multiset).elements())

print overlap, a_remainder, b_remainder