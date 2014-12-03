load p02945.pdb,p02945
cmd.show_as("cartoon", "all")
select TM_1_0, p02945 and resi 20+21+23+24+25+27+28+31
color red, TM_1_0
select TM_1_1, p02945 and resi 53+54+58+62+65+66
color red, TM_1_1
select TM_1_2, p02945 and resi 89+93+96+97+100
color red, TM_1_2
select TM_1_3, p02945 and resi 143+144+147+148+150+151+152+154+155
color red, TM_1_3
select TM_1_4, p02945 and resi 180+181+184+185+187+188+189+191+192
color red, TM_1_4
select TM_1_5, p02945 and resi 201+202+204+205+206+208+209+212+213+215+216+217+219+220
color red, TM_1_5
select TM_1_6, p02945 and resi 112+113+116+120+123
color red, TM_1_6
save p02945_TM_temp.pse
