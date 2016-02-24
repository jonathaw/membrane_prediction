#!/usr/bin/env python2.7
import matplotlib.pyplot as plt
import numpy as np

from ProcessEntry import MakeHydrophobicityGrade
from WinGrade import MEMBRANE_HALF_WIDTH


def main():
    polyval = MakeHydrophobicityGrade()
    z = np.linspace(-MEMBRANE_HALF_WIDTH, MEMBRANE_HALF_WIDTH, endpoint=True, num=1000)
    for i, aa in enumerate(list('ACDEFGHIKLMNPQRSTVWY')):
        plt.subplot(4, 5, i+1)
        plt.plot(z, [np.polyval(polyval[aa], z_) for z_ in z], c='g')
        plt.title(aa)
        plt.ylim([-3, 3])
    plt.show()


if __name__ == '__main__':
    main()