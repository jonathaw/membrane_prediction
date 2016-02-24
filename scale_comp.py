#!/usr/bin/env python2.7
from Equation import Expression
import matplotlib.pyplot as plt
import numpy as np
from math import e as e_
from scipy.constants import k as k_
from ProcessEntry import MakeHydrophobicityGrade

__author__ = 'jonathan'

T = 300.0
all_res = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
three_2_one = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
               'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
               'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}


def main_dinola():
    x = np.arange(-20., 20., 0.1)
    # degrado = get_degrade_scale()
    dinola = get_dinolad_scale()
    dinola_e_erf_zs = get_e_ref(x, dinola)
    ala_ref = {a: (-0.6*np.log(dinola['A'](a))-dinola_e_erf_zs[a])/0.6 for a in x}
    for i, res in enumerate(all_res):
        plt.subplot(4, 5, i+1)
        if res in ['C', 'S', 'T']:
            continue
        plt.plot(x, [(-0.6*np.log(dinola[res](a))-dinola_e_erf_zs[a])/0.6-ala_ref[a] for a in x])
        plt.title(res)
    plt.show()


def main_degrado():
    z = np.arange(-20., 20., 0.1)

    degrado = get_degrade_scale()
    degrado_ala = {z_: degrado['A'](z_) for z_ in z}

    elazar = get_elazar_scale()

    dinola = get_dinolad_scale()
    dinola_e_erf_zs = get_e_ref(z, dinola)
    dinola_ala = {a: (-0.6*np.log(dinola['A'](a))-dinola_e_erf_zs[a]) for a in z}

    senes = get_senes_scale()
    senes_ala = {z_: senes['A'](z_) for z_ in z}

    for i, res in enumerate(all_res):
        plt.subplot(4, 5, i+1)
        if res == 'P':
            plt.plot([z_ for z_ in z if z_ <= 0], [degrado[res][0](z_) - degrado_ala[z_] for z_ in z if z_ <= 0], c='k',
                     label='Degrado')
            plt.plot([z_ for z_ in z if z_ >= 0], [degrado[res][1](z_) - degrado_ala[z_] for z_ in z if z_ >= 0], c='k')
        else:
            plt.plot(z, [degrado[res](z_) - degrado_ala[z_] for z_ in z], c='k', label='Degrado')
        plt.plot(z, [elazar[res](z_) for z_ in z], c='b', label='Elazar')
        if res not in ['C', 'S', 'T']:
            plt.plot(z, [(-0.6*np.log(dinola[res](z_)))-dinola_e_erf_zs[z_] - dinola_ala[z_] for z_ in z], c='g',
                     label='DiNola')

        plt.plot(z, [senes[res](z_)-senes_ala[z_] for z_ in z], c='r', label='Senes')

        plt.ylim([-4., 5.])
        plt.title(res)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()


def prism_all():
    z = np.arange(-20., 20., 0.1)

    degrado = get_degrade_scale()
    degrado_ala = {z_: degrado['A'](z_) for z_ in z}

    elazar = get_elazar_scale()

    dinola = get_dinolad_scale()
    dinola_e_erf_zs = get_e_ref(z, dinola)
    dinola_ala = {a: (-0.6*np.log(dinola['A'](a))-dinola_e_erf_zs[a]) for a in z}

    senes = get_senes_scale()
    senes_ala = {z_: senes['A'](z_) for z_ in z}

    for res in all_res:
        with open('all_%s.txt' % res, 'w+') as fout:
            fout.write('Z\t%s\t%s\t%s\t%s\n' % ('Elazar', 'Degrado', 'Senes', 'DiNola'))
            for z_ in z:
                ela = elazar[res](z_)
                deg = degrado_score(res, z_, degrado, degrado_ala)
                sen = senes[res](z_) - senes_ala[z_]
                if res not in ['C', 'S', 'T']:
                    din = (-0.6*np.log(dinola[res](z_)))-dinola_e_erf_zs[z_] - dinola_ala[z_]
                else:
                    din = -100.
                fout.write('%f\t%f\t%f\t%f\t%f\n' % (z_, ela, deg, sen, din))


def get_degrade_scale():
    results = {}
    degrado = parse_degrado()
    for res, val in degrado.items():
        eqs = []
        # eq = Expression("sin(x+y^2)",["y","x"])
        if 'sig_right' in val['eq_type']:
            eqs.append(Expression("%f/(1+(z/%f)^%f)" % (val['e_0'], val['z_mid_left'], val['n_left']), argorder=['z']))
            eqs.append(Expression("%f/(1+(z/%f)^%f)" % (val['e_0'], val['z_mid_right'], val['n_right']),
                                  argorder=['z']))
            results[three_2_one[res]] = eqs
            continue

        for eq_type in val['eq_type']:
            if eq_type == 'gauss_left':
                eqs.append(Expression("%f*%f^((-((z-%f)^2))/(2*%f))" % (val['a_left'], e_, val['mu_left'],
                                                                     val['sigma_left']), argorder=['z']))
            if eq_type == 'gauss_right':
                eqs.append(Expression("%f*%f^((-((z-%f)^2))/(2*%f))" % (val['a_right'], e_, val['mu_right'],
                                                                     val['sigma_right']), argorder=['z']))
            if eq_type == 'sig_left' and 'sig_right' not in val['eq_type']:
                eqs.append(Expression("%f/(1+(abs(z/%f))^%f)" % (val['e_0'], val['z_mid_left'], val['n_left']),
                                      argorder=['z']))

        results[three_2_one[res]] = sum(eqs)
    return results


def parse_degrado():
    results = {}
    header = ['paq', 'a_left', 'mu_left', 'sigma_left', 'a_right', 'mu_right', 'sigma_right', 'e_0', 'z_mid_left',
              'n_left', 'z_mid_right', 'n_right', 'r_2']
    for l in open('/home/labs/fleishman/jonathaw/membrane_prediciton/degrade_data.tsv', 'r'):
    # for l in open('/Volumes/labs/fleishman/jonathaw/membrane_prediciton/degrade_data.tsv', 'r'):
        s = l.split()
        if s[0][:3] in three_2_one.keys():
            r = {}
            for i, t in enumerate(s[1:]):
                try:
                    r[header[i]] = float(t)
                except:
                    r[header[i]] = t

            results[s[0][:3]] = r.copy()
            results[s[0][:3]]['eq_type'] = []
            if r['a_left'] != '-':
                results[s[0][:3]]['eq_type'].append('gauss_left')
            if r['a_right'] != '-':
                results[s[0][:3]]['eq_type'].append('gauss_right')
            if r['n_left'] != '-':
                results[s[0][:3]]['eq_type'].append('sig_left')
            if r['n_right'] != '-':
                results[s[0][:3]]['eq_type'].append('sig_right')
    return results


def get_dinolad_scale():
    """
    :return: dict of all equations for all res
    """
    dinola = parse_dinola()
    results = {}
    for k, v in dinola.items():
        eq = Expression('%f+%f*%f^(-%f*(z-%f)^2)+%f*%f^(-%f*(z-%f)^2)' % (v[0], v[1], e_, v[2], v[3], v[4], e_, v[5], v[6]),
                        argorder=['z'])
        results[three_2_one[k]] = eq
    return results


def get_e_ref(zs, dinola):
    erefs = {}
    for z in zs:
        erefs[z] = -0.6*np.log(sum(dinola[res](z) for res in all_res if res not in ['C', 'S', 'T']))
    return erefs


def test():
    dinola = get_dinolad_scale()
    print dinola['A']
    print dinola['A'](-100)


def parse_dinola():
    results = {}
    for l in open('/home/labs/fleishman/jonathaw/membrane_prediciton/dinola_data.tsv', 'r'):
    # for l in open('/Volumes/labs/fleishman/jonathaw/membrane_prediciton/dinola_data.tsv', 'r'):
        s = l.split()
        if s[0] in three_2_one.keys():
            results[s[0]] = [float(a.replace('\xe2\x88\x92', '-')) for a in s[1:8]]
    return results


def get_elazar_scale():
    resutls = {}
    for res, val in MakeHydrophobicityGrade().items():
        resutls[res] = Expression("%f*z^4+%f*z^3+%f*z^2+%f*z+%f" % (val[0], val[1], val[2], val[3], val[4]),
                                  argorder=['z'])
    return resutls


def parse_senes():
    results = {}
    for l in open('/home/labs/fleishman/jonathaw/membrane_prediciton/senes_data.tsv', 'r'):
        s = l.split()
        if len(s) < 1:
            continue
        if s[0].upper() in three_2_one.keys():
            results[three_2_one[s[0].upper()]] = {'e0': float(s[1]), 'zmid': float(s[2]), 'n': float(s[3])}
    return results


def get_senes_scale():
    results = {}
    senes = parse_senes()
    for res, val in senes.items():
        if res in ['W', 'Y']:
            results[res] = Expression("%f*%f^(-((abs(z)-%f)^2)/(2*(%f^2)))" % (val['e0'], e_, val['zmid'], val['n']),
                                      argorder=['z'])
        else:
            results[res] = Expression("%f/(1+((abs(z)/%f)^%f))" % (val['e0'], val['zmid'], val['n']), argorder=['z'])
    results['C'] = Expression('0*z', argorder=['z'])
    return results


def prism_elazar():
    z = np.arange(-20., 20., 0.1)
    elazar = get_elazar_scale()
    print 'z\t' + '\t\t'.join([k for k in elazar.keys()])
    for z_ in z:
        print '%f\t' % z_ + '\t\t'.join(['%1.9f' % elazar[k](z_) for k in elazar.keys()])


def degrado_score(res, z_, degrado, degrado_ala):

    if res == 'P':
        if z_ >= 0:
            return degrado[res][1](z_) - degrado_ala[z_]
        else:
            return degrado[res][0](z_) - degrado_ala[z_]

    else:
        return degrado[res](z_) - degrado_ala[z_]


def prism_degrado():
    z = np.arange(-20., 20., 0.1)

    degrado = get_degrade_scale()
    degrado_ala = {z_: degrado['A'](z_) for z_ in z}

    print 'z\t' + '\t\t'.join([k for k in degrado.keys()])
    for z_ in z:
        print '%f\t' % z_ + '\t\t'.join(['%1.9f' % (degrado_score(k, z_, degrado, degrado_ala))
                                         for k in degrado.keys()])


if __name__ == '__main__':
    # main_degrado()
    # prism_elazar()
    # prism_degrado()
    prism_all()