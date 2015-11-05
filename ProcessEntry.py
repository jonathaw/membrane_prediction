class TopoEntry():
    def __init__(self, name, seq, end_of_SP, psipred, hydro_polyval, param_list, csts, path_msa):
        self.name = name
        self.seq = seq
        self.seq_length = len(seq)
        self.psipred = psipred
        self.hydro_polyval = hydro_polyval
        self.param_list = param_list
        self.csts = csts
        self.path_msa = path_msa

    def __str__(self):
        msg = 'name %s\n' % self.name
        msg += 'seq %s\n' % self.seq
        msg += 'csts\n'
        msg += str(self.csts)
        msg += 'params\n'
        for k, v in self.param_list.items():
            msg += "%s %s\n" % (str(k), str(v))
        return msg


def create_topo_entry(name, seq, ss2, param_list, csts, db, msa_path):
    import re

    psipred = parse_psipred(ss2)
    h_polyval = MakeHydrophobicityGrade()
    if db == 'rost':
        topc = spc_parser(name)
        if topc['spoctopus'].count('S') != 0:
            end_of_SP = [a for a in re.finditer('S*', topc['spoctopus']) if a != ''][0].end() - 1
            if end_of_SP == -1:
                end_of_SP = 0
            seq_no_SP = 'u' * end_of_SP + seq[end_of_SP:]
        else:
            end_of_SP = 0
            seq_no_SP = seq
    else:
        end_of_SP = 0
        seq_no_SP = seq

    return TopoEntry(name=name, seq=seq_no_SP, end_of_SP=end_of_SP, psipred=psipred, hydro_polyval=h_polyval,
                     param_list=param_list, csts=csts, path_msa=msa_path)


def parse_psipred(ss2):
    results = {}
    with open(ss2, 'r') as f:
        cont = f.read().split('\n')
    for line in cont:
        split = line.split()
        if len(split) != 6: continue
        results[int(split[0]) - 1] = {'aa': split[1], 'c': float(split[3]), 'h': float(split[4]), 'e': float(split[5])}
    return results


def MakeHydrophobicityGrade():
    """
    :return: returns a dictionary of the polynom values for each residue
    """
    global hydrophobicity_polyval
    # hydrophobicity_grade = open('Poly_Values.txt', 'r')
    # hydrophobicity_grade = open('poly_value_11.2.txt', 'r')
    # hydrophobicity_grade = open('poly_vals_23.2.txt', 'r')
    # hydrophobicity_grade = open('/home/labs/fleishman/jonathaw/membrane_prediciton/poly_vals_25.2.txt', 'r')
    try:
        hydrophobicity_grade = open('/home/labs/fleishman/jonathaw/membrane_prediciton/polyval_21_5_15.txt', 'r')
    # else:
    except:
        hydrophobicity_grade = open('/Volumes/labs/fleishman/jonathaw/membrane_prediciton/polyval_21_5_15.txt', 'r')
    # hydrophobicity_grade = open('/Volumes/jonathaw-1/membrane_prediciton/poly_vals_25.2.txt', 'r')
    # hydrophobicity_grade = open('./poly_vals_25.2.txt', 'r')
    hydrophobicity_polyval = {}
    for line in hydrophobicity_grade:
        split = line.split()
        hydrophobicity_polyval[split[0]] = [float(n) for n in split[1:6]]
    hydrophobicity_grade.close()
    return hydrophobicity_polyval


def spc_parser(name):
    result = {}
    with open('/home/labs/fleishman/jonathaw/membrane_prediction_DBs/spoctopus_SPDB/' + name.lower() + '.spc',
              'r') as f:
        cont = f.read().split('\n')
    for line in cont:
        split = line.split()
        if split == []: continue
        result[split[0]] = split[1]
    return result


def write_results(topo_entry, run_type, best_path, sec_path, msa=False):
    from timeit import default_timer
    suffix = '.prd' if not msa else '_msa.prd'
    ts_best = topo_string_rostlab_format(topo_entry, best_path)
    ts_sec = topo_string_rostlab_format(topo_entry, sec_path)
    with open(topo_entry.param_list['out_path'] + topo_entry.name + suffix, 'wr+') as fout:
        fout.write(str(topo_entry))
        fout.write('run_type %s\n' % run_type)
        fout.write('best_path_ts %s\n' % ts_best)
        fout.write('sec_path_ts %s\n' % ts_sec)
        fout.write("best_path %s\n" % str(best_path))
        fout.write("sec_path %s\n" % str(sec_path))
        fout.write("ddG paths |%f|\n" % abs(best_path.total_grade - sec_path.total_grade))
        fout.write("time %f\n" % (default_timer()-topo_entry.param_list['tic']))


def process_entry(topo_entry, run_type, verbose=False):
    from output2html import create_html
    from WinGrade import wgp_msa_total_grade

    # this is a condition to prevent running as msa2plain when there's only one sequence in the MSA, which messes up
    if run_type == 'msa2plain':
        with open(topo_entry.path_msa+topo_entry.name+'.msa', 'r') as fin:
            cont = fin.read()
            if cont.count('>') == 1:
                run_type = 'plain'
                print 'found only 1 sequence in the MSA, so making this run as plain'

    if run_type == 'msa2plain':
        from TMConstraint import pred2cst

        topo_entry.param_list['with_msa'] = True
        wins = win_grade_generator(topo_entry)
        print 'found %i wins with msa' % len(wins)
        best_path, sec_path = topo_graph(topo_entry, wins)
        write_results(topo_entry, run_type, best_path, sec_path, msa=True)
        print best_path
        if best_path.total_grade >= 0. and wgp_msa_total_grade(best_path) >= 0.:
            print 'found no solution with MSA, will try as plain'
            run_type = 'plain'
        else:
            ts = topo_string_rostlab_format(topo_entry, best_path)
            topo_entry.param_list['with_msa'] = False
            topo_entry.csts = pred2cst(topo_entry.name, topo_entry.param_list['out_path'], ts, cst_mode=None,
                                       tm_pos_fidelity=topo_entry.param_list['fidelity'])
            wins = win_grade_generator(topo_entry)
            print 'found %i wins without msa' % len(wins)
            best_path, sec_path = topo_graph(topo_entry, wins)
            write_results(topo_entry, run_type, best_path, sec_path)
            create_html(topo_entry, best_path, sec_path, wins)

    if run_type == 'cst_only':
        wins = win_grade_generator(topo_entry, 'only')
        best_path, sec_path = topo_graph(topo_entry, wins)
        write_results(topo_entry, run_type, best_path, sec_path)
        create_html(topo_entry, best_path, sec_path, wins)

    if run_type == 'tm_num_cst':
        wins = win_grade_generator(topo_entry)
        topo_entry.csts.tm_pos = topo_graph_tm_num(topo_entry, wins)
        best_path, sec_path = topo_graph(topo_entry, wins)
        write_results(topo_entry, run_type, best_path, sec_path)
        create_html(topo_entry, best_path, sec_path, wins)

    if run_type == 'plain':
        topo_entry.param_list['with_msa'] = False
        topo_entry.param_list['with_cst'] = False
        wins = win_grade_generator(topo_entry)
        best_path, sec_path = topo_graph(topo_entry, wins)
        write_results(topo_entry, run_type, best_path, sec_path)
        create_html(topo_entry, best_path, sec_path, wins)

    if run_type == 'user_cst':
        topo_entry.param_list['with_msa'] = False
        topo_entry.param_list['with_cst'] = True
        # if topo_entry.csts.tm_pos
        wins = win_grade_generator(topo_entry, mode='only')
        try:
            best_path, sec_path = topo_graph(topo_entry, wins)
        except:
            print 'failed topo_graph, will try with all wins'
            wins = win_grade_generator(topo_entry, mode='only', tm_pos_mode='all')
            best_path, sec_path = topo_graph(topo_entry, wins)
        if topo_entry.csts.mode == 'only' and best_path.path == []:
            wins = win_grade_generator(topo_entry, mode='only', tm_pos_mode='all')
            for w in wins:
                w.set_grade(-100.0)
            best_path, sec_path = topo_graph(topo_entry, wins)
            best_path.add_to_all_wins(+100.0)
            sec_path.add_to_all_wins(+100.0)
        write_results(topo_entry, run_type, best_path, sec_path)
        create_html(topo_entry, best_path, sec_path, wins)

    if run_type not in ['msa2plain', 'cst_only', 'tm_num_cst', 'plain', 'user_cst']:
        print "unrecoginzed run-type"


def topo_string_rostlab_format(topo_entry, wgp):
    """
    :param topo:a topo (list of WinGrades describing a constructed topology)
    :param seq: the sequence
    :return:a string describing the topology in rostlab's format where 1:inside, 2: outdise H: TM helix
    """
    from WinGrade import WinGrade

    topo_string = ''
    last_tm = WinGrade(0, 0, 'fwd', '', topo_entry.hydro_polyval,
                       {k: v for k, v in topo_entry.param_list.items() if k in ['c0', 'c1', 'c2', 'c3', 'w', 'z_0']})
    for tm in wgp.path:
        topo_string += '1' * (tm.begin - last_tm.end) if tm.direction == 'fwd' else '2' * (tm.begin - last_tm.end)
        topo_string += 'H' * (tm.end - tm.begin)
        last_tm = tm
    topo_string += '2' * (len(topo_entry.seq) - last_tm.end) if \
        last_tm.direction == 'fwd' else '1' * (len(topo_entry.seq) - last_tm.end)
    return topo_string


def win_grade_generator(topo_entry, mode='all', tm_pos_mode='selective'):
    """
    :return:grades all segments of seq, and aggregates them as WinGrades
    """
    from WinGrade import WinGrade, count_charges
    from TMConstraint import seg_within_tm_pos, win_in_tm_pos
    import sys

    assert isinstance(topo_entry, TopoEntry)
    MIN_LENGTH = topo_entry.param_list['min_length']
    INC_MAX = topo_entry.param_list['inc_max']
    if topo_entry.param_list['with_msa']:
        from MSA_for_TMpredict import TMpredict_MSA, retrieve_seqs

        msa_obj = TMpredict_MSA(topo_entry.name, topo_entry.hydro_polyval, topo_entry.param_list,
                                topo_entry.param_list['msa_percentile'], path_msa=topo_entry.path_msa)
    psi = topo_entry.psipred
    grades = []
    pos1 = 0 if topo_entry.seq[0] != 'u' else max([i for i, aa in enumerate(topo_entry.seq) if aa == 'u']) + 1
    pos2 = topo_entry.seq_length + 1
    bar_width = len(range(pos1, pos2))
    sys.stdout.write("[%s]" % (" " * bar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (bar_width + 1))  # return to start of line, after '['
    if topo_entry.csts.tm_pos is not None:
        tm_pos_wins = {tmp: [] for tmp in topo_entry.csts.tm_pos}
    for i in range(pos1, pos2 + 1):
        for inc in range(min(INC_MAX, topo_entry.seq_length - MIN_LENGTH - i)):
            tm_pos, j = win_in_tm_pos((i, i + MIN_LENGTH + inc), topo_entry.csts.tm_pos,
                                      topo_entry.csts.tm_pos_fidelity)
            if mode == 'only':
                if tm_pos is False:
                    continue
            if tm_pos is not False:
                if topo_entry.param_list['with_msa']:
                    tm_pos_wins[tm_pos].append(retrieve_seqs(msa_obj, i, i + MIN_LENGTH + inc, 'fwd'))
                    tm_pos_wins[tm_pos].append(retrieve_seqs(msa_obj, i, i + MIN_LENGTH + inc, 'rev'))
                else:
                    tm_pos_wins[tm_pos].append(
                        WinGrade(i, i + MIN_LENGTH + inc, 'fwd', topo_entry.seq[i:i + MIN_LENGTH + inc],
                                 topo_entry.hydro_polyval, topo_entry.param_list))
                    tm_pos_wins[tm_pos].append(WinGrade(i, i + MIN_LENGTH + inc, 'rev',
                                                        topo_entry.seq[i:i + MIN_LENGTH + inc][::-1],
                                                        topo_entry.hydro_polyval, topo_entry.param_list))
            else:
                if not is_not_helical(topo_entry.seq, (i, i + MIN_LENGTH + inc), psi) and \
                                count_charges(topo_entry.seq[i:i + MIN_LENGTH + inc]) < 3 and \
                                len([a for a in topo_entry.seq[i:i + MIN_LENGTH + inc] if a in
                                        ['R', 'K', 'H', 'D', 'E', 'N', 'Q']]) < 5:
                    if topo_entry.param_list['with_msa']:
                        grades.append(retrieve_seqs(msa_obj, i, i + MIN_LENGTH + inc, 'fwd'))
                        grades.append(retrieve_seqs(msa_obj, i, i + MIN_LENGTH + inc, 'rev'))

                    else:
                        grades.append(WinGrade(i, i + MIN_LENGTH + inc, 'fwd', topo_entry.seq[i:i + MIN_LENGTH + inc],
                                               topo_entry.hydro_polyval, topo_entry.param_list))
                        grades.append(
                            WinGrade(i, i + MIN_LENGTH + inc, 'rev', topo_entry.seq[i:i + MIN_LENGTH + inc][::-1],
                                     topo_entry.hydro_polyval, topo_entry.param_list))

        sys.stdout.write("-")
        sys.stdout.flush()
    sys.stdout.write("]\n")
    if topo_entry.csts.tm_pos is not None:
        for k, cst_pos in tm_pos_wins.items():
            print 'at tm_pos %s found %i wins' % (k, len(cst_pos))
            if tm_pos_mode == 'selective':
                grades.extend(choose_wins_for_cst(topo_entry, cst_pos))
            elif tm_pos_mode == 'all':
                # print 'extending', '\n'.join([str(a.end) for a in cst_pos])
                grades.extend(cst_pos)
    return grades


def is_not_helical(seq, pos, psi, verbose=False):
    win_size = 6
    for i in range(pos[0], pos[1] - win_size + 2):
        if all(psi[j]['e'] >= 0.5 for j in range(i, i + win_size)) or \
                all(psi[j]['c'] >= 0.5 for j in range(i, i + win_size)) or \
                all(psi[j]['h'] <= 0.1 for j in range(i, i + win_size)):
            return True
    cs = []
    es = []
    hs = []
    for i in range(pos[0], pos[0] + 3):
        cs.append(psi[i]['c'] >= 0.5)
        es.append(psi[i]['e'] >= 0.5)
        hs.append(psi[i]['h'] <= 0.1)
    if all(cs) or all(es) or all(hs):
        return True
    cs = []
    es = []
    hs = []
    for i in range(pos[1] - 3, pos[1]):
        cs.append(psi[i]['c'] >= 0.5)
        es.append(psi[i]['e'] >= 0.5)
        hs.append(psi[i]['h'] <= 0.1)
    if all(cs) or all(es) or all(hs):
        return True


def choose_wins_for_cst(topo_entry, cst_pos, verbose=False):
    """
    :param topo_entry: a TopoEntry instance
    :param cst_pos: a bunch of wins from a bin of a certain tm_pos constraint
    :return: all wins that pass the thresholds, or the best one if none do
    """
    from WinGrade import count_charges, WinGrade

    passed = []
    best_fwd, best_rev = WinGrade(0, 0, 'fwd', '', grade=100.0), WinGrade(0, 0, 'fwd', '', grade=100.0)
    for win in cst_pos:
        if not is_not_helical(win.seq, [win.begin, win.end], topo_entry.psipred, verbose) and \
                        count_charges(win.seq) < 3 and \
                        len([a for a in win.seq if a in ['R', 'K', 'H', 'D', 'E', 'N', 'Q']]) < 5:
            passed.append(win)
        if win.direction == 'fwd':
            if win.grade <= best_fwd.grade:
                best_fwd = win
        if win.direction == 'rev':
            if win.grade <= best_rev:
                best_rev = win
    if passed is not []:
        return passed
    else:
        return [best_fwd, best_rev]


def find_cst_after(csts, win):
    """
    :param win: a WinGrade instance
    :return: a (beginning, end, direction) tuple from tm_pos that is after win
    """
    if csts.tm_pos is None:
        return None
    for cst in sorted(csts.tm_pos):
        if cst[0] - csts.tm_pos_fidelity > win.begin:
            return cst
    return None


def find_graph_path(pred, path, source_node):
    """
    :param pred:
    :param path:
    :param source_node:
    :return:
    """
    # find all wins in the minimal energy path from source to last win
    while path[-1].seq != source_node.seq:
        for k, v in pred.items():
            if v is None:
                continue
            if k.seq == path[-1].seq:
                path.append(v)
    path.pop(-1)  # get rid of source_node
    return path[::-1]  # revert the path to begin->end


def topo_graph(topo_entry, wins):
    """
    :return: a topology
    """
    import networkx as nx
    from WinGrade import WinGrade, WinGradePath, wgp_msa_total_grade
    import operator
    # print self.WinGrades
    G = nx.DiGraph()
    source_node = WinGrade(0, 0, 'fwd', '', topo_entry.hydro_polyval, topo_entry.param_list)  # define source win
    print "Constructing graph"
    if topo_entry.param_list['with_msa']:
        if topo_entry.csts.tm_pos is None:
            [G.add_edge(source_node, a, weight=a.msa_grade) for a in wins]  # add all win to source edges with msa
        else:
            [G.add_edge(source_node, a, weight=a.msa_grade) for a in wins if a.end <= topo_entry.csts.tm_pos[0][1] +
             topo_entry.csts.tm_pos_fidelity]  # add all win to source edges with msa
    else:
        if topo_entry.csts.tm_pos is None:
            [G.add_edge(source_node, a, weight=a.grade) for a in wins]  # add all win to source edges
        else:
            [G.add_edge(source_node, a, weight=a.grade) for a in wins if a.end <= topo_entry.csts.tm_pos[0][1] +
             topo_entry.csts.tm_pos_fidelity]  # add all win to source edges
    for win1 in wins:  # add all win to win edges where condition applies
        cst_after = find_cst_after(topo_entry.csts, win1)
        for win2 in wins:
            # make sure no edges skip a cst:
            if cst_after is None or win2.end <= cst_after[1] + topo_entry.csts.tm_pos_fidelity:
                # condition: non-overlapping, 2 is after 1, opposite directions
                if not win1.grade_grade_colliding(win2) and win2.begin > win1.end and win1.direction != win2.direction \
                        and win2.begin - win1.end >= 2:
                    if not topo_entry.param_list['with_msa']:
                        G.add_edge(win1, win2, weight=win2.grade)
                    elif topo_entry.param_list['with_msa']:
                        G.add_edge(win1, win2, weight=win2.msa_grade)
    print "About to Bellman-Ford"
    print 'found %i nodes' % G.number_of_nodes()
    if G.number_of_nodes() == 0:
        return WinGradePath([]), WinGradePath([])
    pred, dist = nx.bellman_ford(G, source_node)
    print "Finished Bellman-Fording"
    sorted_dist = sorted(dist.items(), key=operator.itemgetter(1))
    temp_direction = 'A'
    best_path, sec_best_path = [], []
    for last_path, total_grade in sorted_dist:
        temp_path = find_graph_path(pred, [last_path], source_node)
        if temp_path == [] or temp_path == [source_node]:
            continue
        if temp_direction != 'A' and topo_entry.csts.test_manager(temp_path) \
                and temp_path[-1].direction is not temp_direction:
            sec_best_path = temp_path[:]
            if sec_best_path[0].seq == 'SOURCE':
                sec_best_path.pop(0)
            break
        elif topo_entry.csts.test_manager(temp_path) and temp_path[-1].direction != temp_direction:
            best_path = temp_path[:]
            temp_direction = best_path[-1].direction
    best_wgp = WinGradePath(best_path)
    sec_wgp = WinGradePath(sec_best_path)
    return best_wgp, sec_wgp


def topo_graph_tm_num(topo_entry, wins):
    """
    a method to find minimal path with only tm_num nodes. a BFS algorithm that stops every path building once
    it gets to tm_num nodes. only saves the best solutions for c_term fwd & rev, and only if they pass all
    constraints
    """
    from collections import deque
    from copy import deepcopy
    import timeit
    import networkx as nx
    from WinGrade import WinGrade, WinGradePath

    G = nx.DiGraph()
    source_node = WinGrade(0, 0, 'fwd', '', topo_entry.hydro_polyval, topo_entry.param_list)  # define source win
    win_list_25 = [source_node] + [a for a in wins if a.length == 19]
    print "BFS building graph"
    if topo_entry.param_list['with_msa']:
        if topo_entry.csts.tm_pos is None:
            [G.add_edge(0, i, weight=a.msa_grade) for i, a in enumerate(win_list_25) if
             a != source_node]  # add all win to source edges with msa
        else:
            [G.add_edge(0, i, weight=a.msa_grade) for i, a in enumerate(win_list_25) if a.end <=
             topo_entry.csts.tm_pos[0][
                 1] + topo_entry.csts.tm_pos_fidelity and a != source_node]  # add all win to source edges with msa
    else:
        if topo_entry.csts.tm_pos is None:
            [G.add_edge(0, i, weight=a.grade) for i, a in enumerate(win_list_25) if
             a != source_node]  # add all win to source edges
        else:
            [G.add_edge(0, i, weight=a.grade) for i, a in enumerate(win_list_25) if a.end <=
             topo_entry.csts.tm_pos[0][
                 1] + topo_entry.csts.tm_pos_fidelity and a != source_node]  # add all win to source edges
    for win1 in win_list_25:  # add all win to win edges where condition applies
        cst_after = find_cst_after(topo_entry.csts, win1)
        for win2 in win_list_25:
            if cst_after is None or win2.begin < cst_after[1]:  # make sure no edges skip a cst
                # condition: non-overlapping, 2 is after 1, opposite directions
                if not win1.grade_grade_colliding(win2) and win2.begin > win1.end and win2.begin - win1.end >= 2 \
                        and win1.direction != win2.direction:
                    if not topo_entry.param_list['with_msa']:
                        G.add_edge(win_list_25.index(win1), win_list_25.index(win2), weight=win2.grade)
                    elif topo_entry.param_list['with_msa']:
                        G.add_edge(win_list_25.index(win1), win_list_25.index(win2), weight=win2.msa_grade)

    print "BFS finished building graph"
    tic = timeit.default_timer()
    print 'looking for %i wins' % topo_entry.csts.tm_num
    best_fwd = WinGradePath([])
    best_rev = WinGradePath([])
    Q = deque([[0]])
    num = 0
    mini_tic = timeit.default_timer()
    while Q:  # continue as long as there si something in Q
        path = Q.popleft()  # get the FIFO out
        if len(path) != num:
            print "reached %i in path, the Q is %12i, took %f seconds" % (len(path), len(Q),
                                                                          timeit.default_timer() - mini_tic)
            mini_tic = timeit.default_timer()
            num = len(path)
        if len(path) - 1 == topo_entry.csts.tm_num:  # the path has enough nodes
            path_wins = WinGradePath([win_list_25[i] for i in path])
            if topo_entry.csts.test_manager(path_wins.path):
                if path_wins.c_term == 'fwd' and path_wins.total_grade < best_fwd.total_grade:
                    best_fwd = deepcopy(path_wins)
                elif path_wins.c_term == 'rev' and path_wins.total_grade < best_rev.total_grade:
                    best_rev = deepcopy(path_wins)
            continue
        for neighbor in G.successors_iter(path[-1]):  # go over all successors of last win of path
            temp_path = deepcopy(path)
            temp_path.append(neighbor)
            Q.append(temp_path)
    print 'elapsed time for BFS', timeit.default_timer() - tic
    path_cst = best_fwd if best_fwd.total_grade < best_rev.total_grade else best_rev
    new_csts = []
    for win in path_cst.path:
        new_csts.append((max(0, win.begin), min(topo_entry.seq_length, win.end), None))
    return new_csts


def topo_graph_only(topo_entry, wins):
    from collections import deque
    from copy import deepcopy
    import timeit
    import networkx as nx
    from WinGrade import WinGrade, WinGradePath

    G = nx.DiGraph()
    source_node = WinGrade(0, 0, 'fwd', '', topo_entry.hydro_polyval, topo_entry.param_list)  # define source win
    print "BFS building graph"
    [G.add_edge(0, i, weight=a.grade) for i, a in enumerate(wins) if a.end <=
     topo_entry.csts.tm_pos[0][1] + topo_entry.csts.tm_pos_fidelity and a != source_node]  # add all win to source edges
    for win1 in wins:  # add all win to win edges where condition applies
        cst_after = find_cst_after(topo_entry.csts, win1)
        for win2 in wins:
            if cst_after is None or win2.begin < cst_after[1]:  # make sure no edges skip a cst
                # condition: non-overlapping, 2 is after 1, opposite directions
                if not win1.grade_grade_colliding(win2) and win2.begin > win1.end and win2.begin - win1.end >= 2 \
                        and win1.direction != win2.direction:
                    G.add_edge(wins.index(win1), wins.index(win2), weight=win2.grade)

    print "BFS finished building graph"
    tic = timeit.default_timer()
    print 'looking for %i wins' % len(topo_entry.csts.tm_pos)
    best_fwd, fwd_grade = WinGradePath([]), 1000
    best_rev, rev_grade = WinGradePath([]), 1000
    Q = deque([[0]])
    num = 0
    mini_tic = timeit.default_timer()
    while Q:  # continue as long as there si something in Q
        path = Q.popleft()  # get the FIFO out
        if len(path) != num:
            print "reached %i in path, the Q is %12i, took %f seconds" % (len(path), len(Q),
                                                                          timeit.default_timer() - mini_tic)
            mini_tic = timeit.default_timer()
            num = len(path)
        if len(path) - 1 == len(topo_entry.csts.tm_pos):
            # print 'found enpugh edges', path
            # the path has enough nodes
            path_wins = WinGradePath([wins[i] for i in path])
            if topo_entry.csts.test_manager(path_wins.path):
                if path_wins.c_term == 'fwd' and path_wins.total_grade < fwd_grade:
                    best_fwd = deepcopy(path_wins)
                    fwd_grade = best_fwd.total_grade
                elif path_wins.c_term == 'rev' and path_wins.total_grade < rev_grade:
                    best_rev = deepcopy(path_wins)
                    rev_grade = best_rev.total_grade
            continue
        for neighbor in G.successors_iter(path[-1]):  # go over all successors of last win of path
            temp_path = deepcopy(path)
            temp_path.append(neighbor)
            Q.append(temp_path)
    print 'elapsed time for BFS', timeit.default_timer() - tic
    print best_fwd
    print best_rev
    if best_fwd.total_grade < best_rev.total_grade:
        return best_fwd, best_rev
    else:
        return best_rev, best_fwd


if __name__ == '__main__':
    pass