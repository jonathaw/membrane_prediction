#!/usr/bin/env python2.7
text = ['[8    to 34   in rev =>  -4.755853 KRRVMGAYFFGVGPVMLMVLATSALI          0 0.000000', '39   to 61   in fwd =>  -9.471774 IALSFISLIITVLLWIFYGYSV              0 0.000000', '86   to 107  in rev =>  -7.407589 ASTLIAITVAAFMMQYMMFLL               1 0.000000', '110  to 131  in fwd =>  -7.600904 RAKVSSFILLSALWLTFVYAP               0 0.000000', '152  to 173  in rev =>  -1.749796 KGITMAVALAAFGSSIHVVMG               1 0.000000', '187  to 208  in fwd =>  -0.871715 PLTLIGAALLWFGWFGFNGGS               0 0.000000', '215  to 236  in rev =>   2.992121 VMWVFGAVAASTNTVVVANIA               1 0.000000', '239  to 260  in fwd =>   1.195730 IKGKPGSLGIVSGAIAGLAAI               0 0.000000', '272  to 293  in rev =>  -4.704660 KKIRFDMALYCVIGAVLGIVI               1 0.000000', '301  to 322  in fwd =>   3.096316 WAIHGIGGLWGSVAVGILANP               0 0.000000', '341  to 363  in rev =>  -5.644520 VAKALILTVLFAYATTSAVAIL              0 0.000000]']

for l in text:
    s = l.split()
    if s[4] == 'fwd':
        print s[7]
    else:
        print s[7][::-1]