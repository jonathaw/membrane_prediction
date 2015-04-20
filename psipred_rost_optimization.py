"""
a script to find best psi_helix and psi_res_num values. examines ALL helices described in both pdbtm and opm
in the rostlab database. checks each observed helix if it would pass with certain psi_helix and psi_res_num
parameters. if so counts that for that parameter combination. print and shows a 3D scatter plot.
"""


def main():
    from TMpredict_WinGrade import parse_rostlab_db
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import operator
    psi_helix = [0.001, 0.005, 0.01, 0.1, 0.2, 0.3, 0.4]
    psi_res_num = [1, 2, 3, 4]
    total = 0
    results = {}
    for ph in psi_helix:
        for prn in psi_res_num:
            results[(ph, prn)] = 0

    rostlab_dict = parse_rostlab_db()
    for name, dict in rostlab_dict.items():
        psipred = PsiReaderHelix('/home/labs/fleishman/jonathaw/membrane_prediciton/data_sets/rostlab_db/psipred/'
                                 +name+'.ss2')
        assert len(psipred) == len(dict['seq']) == len(dict['pdbtm']) == len(dict['opm']), 'length unequal %s' % name
        for typ in ['pdbtm', 'opm']:
            helices = split_to_helices(dict['seq'], dict[typ], psipred)
            for h_seq, h_ss2 in helices:
                total += 1
                for ph in psi_helix:
                    for prn in psi_res_num:
                        if pass_helix(h_ss2, ph, prn):
                            results[(ph, prn)] += 1
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # for par, res in results.items():
    for par, res in sorted(results.items(), key=operator.itemgetter(1)):
        print 'psi_helix: %f psi_res_num %i result: %f' % (par[0], par[1], float(res)/float(total))
        ax.scatter(par[0], par[1], float(res)/float(total))
    ax.set_xlabel('psi_helix')
    ax.set_ylabel('psi_res_num')
    ax.set_zlabel('percent')
    plt.show()


def split_to_helices(seq, ts, psipred):
    """
    :param seq: AA sequence
    :param ts: topo-string (rostlab format) from either database
    :param psipred: psipred results as list of helix propensity
    :return: tuples of seq and psipred results list for that seq, corresponding to helices in the ts.
    """
    import re
    result = []
    hhh = re.compile('[hH]*')
    helices = [(a.start(), a.end()) for a in hhh.finditer(ts) if a.end()-a.start() > 1]
    for helix in helices:
        result.append((seq[helix[0]: helix[1]], psipred[helix[0]: helix[1]]))
    return result


def pass_helix(ss2, ph, prn):
    """
    :param ss2: list of psipred results
    :param ph: psi_helix parameter
    :param prn: psi_res_num parameter
    :return: whether this helix passes or not
    """
    didnt_pass = 0
    for grade in ss2:
        didnt_pass += 1 if grade < ph else 0
    return True if didnt_pass < prn else False


def PsiReaderHelix(ss2):
        """
        :return:reads the sequence's psipred results, as list of floats
        """
        ss2_file = open(ss2)
        result = []
        for line in ss2_file:
            split = line.split()
            if len(split) > 5:
                result.append(float(split[4]))
        ss2_file.close()
        return result


if __name__ == '__main__':
    main()