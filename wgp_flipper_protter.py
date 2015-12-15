#!/usr/bin/env python2.7
text = ['[3    to 24   in fwd =>  -3.405503 GNVAWILASTALVMLMVPGVG               0 0.000000', '35   to 58   in rev => -10.418108 GYFIWLLVTIILSIFSLAIMNVA             0 0.000000', '86   to 109  in fwd =>  -3.096447 LLFMMYQMMFAAVTIAILTSAIA             1 0.000000', '113  to 134  in rev =>  -6.461654 HAFPAYVFTLWLASLLIFSSV               0 0.000000', '152  to 173  in fwd =>   1.500724 GMVVHISSGFAALAVAMTIGK               1 0.000000', '188  to 210  in rev =>  -5.477971 LASGGNFGFWGFWLLAAGILTL              1 0.000000', '239  to 260  in fwd =>   1.195730 IKGKPGSLGIVSGAIAGLAAI               0 0.000000', '272  to 293  in rev =>  -4.704660 KKIRFDMALYCVIGAVLGIVI               1 0.000000', '300  to 321  in fwd =>   1.884355 AWAIHGIGGLWGSVAVGILAN               1 0.000000', '336  to 363  in rev =>  -6.846208 VAKALILTVLFAYATTSAVAILQSVLL         0 0.000000]']

for l in text:
    s = l.split()
    if s[4] == 'fwd':
        print s[7]
    else:
        print s[7][::-1]