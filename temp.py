a = [0, 1, 4, 8, 11, 15, 18, 19, 22, 26, 29, 33, 36]

b = range(0, 36)

print a
print b
c = []
[c.append(x in a) for x in b]
print c