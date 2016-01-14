#!/usr/bin/env python
# coding=utf-8
"""
a group of functions to create HTML output for a prediction
"""


def create_html(topo_entry, best_path, sec_path, wins):
    from ProcessEntry import topo_string_rostlab_format
    from timeit import default_timer
    import os
    # print 'making win plot'
    # make_win_plot(topo_entry, best_path, wins)
    best_ts = topo_string_rostlab_format(topo_entry, best_path)
    sec_ts = topo_string_rostlab_format(topo_entry, sec_path)
    print 'creating html at %s' % topo_entry.param_list['out_path']+topo_entry.name + '.html'
    with open(topo_entry.param_list['out_path']+topo_entry.name + '.html', 'wr+') as html:
        html.write('<html>\n')

        html.write('<style>\n')
        html.write('g { color: green }\n')
        html.write('r { color: red }\n')
        html.write('b { color: blue }\n')
        html.write('u { color: purple }\n')
        html.write('div { width: 80%; word-wrap: break-word; }\n')
        html.write('table, th, td { border: 1px solid black; border - collapse: collapse; }\n')
        html.write('th, td { padding: 5px; }\n')
        html.write('body{ font: 22px Arial, sans-serif; margin-left: 50px; }\n')

        html.write('</style>\n')
        html.write('<body>\n')
        html.write('<all>\n')

        html.write('<h1><center>TopoGraph results for job %s</center></h1>\n' % topo_entry.name)

        html.write('<h2><center>Best path results are:</center></h2>\n')
        html.write(path2html_table(best_path))
        html.write(topo_string2html(best_ts, topo_entry.seq) + '\n')
        # html.write('<img src = "%s" />\n' % (topo_entry.name+'.html'))

        html.write('<h2><center>Second best path results are:</center></h2>\n')
        html.write(path2html_table(sec_path))
        html.write(topo_string2html(sec_ts, topo_entry.seq) + '\n')

        html.write('<p>the difference in energy between the two best paths is %f</p>\n' %
                   abs(best_path.total_grade - sec_path.total_grade))

        html.write('<p>color index: <r> TMHs </r> <b> In </b> <g> Out </g> <u> Unknown </u></p>\n')

        html.write('<p>Protter view of the best path:</p>\n')
        protter_link, protter_down = protter_api(topo_entry, best_path)
        html.write('<img src="%s" alt="Protter view"">\n' % protter_down)# style="width:912px;height:684px;">\n' % protter_down)
        html.write('<p>purple is unknown, usually signal peptide</p>\n')
        html.write('<a href="%s">Protter link to design your figure</a>\n' % protter_link)
        # html.write('<p>Files to download</p>\n')
        html.write('<p>computation took %f seconds</p>\n' % (default_timer()-topo_entry.param_list['tic']))
        # what is the path to use for download from the browesr???
        # html.write('<p style="text-align:right">The Fleishman Lab &copy;</p>\n')

        html.write('</all>\n')
        html.write('</body>\n')
        html.write('</html>\n')
        html.write('protter text: %s\n' % protter_down)


def path2html_table(wgp):
    html = '<p>Found %i %s:</p>\n' % (len(wgp.path), 'TMHs' if len(wgp.path) > 1 else 'TMH')
    html += '<table style="width:80%">\n'
    # html += '<tr><th>#</th><th>Begin</th><th>End</th><th>N terminus orientation</th><th>&#916;G<sup>apparent</sup>' \
    #         '[kcal/mol]</th><th>Sequence</th><th>#charges</th><th>Membrane-deformation</th></tr>\n'
    html += '<tr><th>#</th><th>Begin</th><th>End</th><th>N terminus orientation</th><th>&#916;G<sup>app</sup>' \
            '[kcal/mol]</th><th>Sequence</th></tr>\n'
    for i, w in enumerate(wgp.path):
        html += '<p>%s</p>\n' % w.get_html(i+1)
    html += '</table>\n'
    html += 'Total grade for the path: %f\n' % wgp.total_grade
    return html


def make_win_plot(topo_entry, path, wins):
    import matplotlib.pyplot as plt

    plt.figure()
    for w in path.path:
        plt.plot((w.begin, w.end), (w.grade, w.grade), 'k--' if w.direction == 'fwd' else 'b', )

    for w in wins:
        plt.plot((w.begin, w.end), (w.grade, w.grade), 'k' if w.direction == 'fwd' else 'r', )

    plt.xlim([0, len(topo_entry.seq)])
    plt.ylim([min([a.grade for a in path.path]) - 2, 10])
    plt.xlabel('Sequence Position')
    plt.ylabel('$/Delta/$G')
    plt.title('Win Grades Energy Plot for %s' % topo_entry.name)
    print 'saving fig to', topo_entry.name + '.png'
    plt.savefig(topo_entry.name + '.png')


def topo_string2html(ts, seq):
    """
    :param ts: a topo string
    :return: an HTML representation of the topo string, where TMHs are Red (r), In is Blue b),
    Out is Green (g) and unknown is Purple (p)
    >>> ts = '1111HHHH2222HHHH1111UUUU'
    >>> se = 'AAAAHHHHBBBBHHHHBBBBUUUU'
    >>> topo_string2html(ts, se)
    '<pre><div><b>AAAA</b><r>HHHH</r><g>BBBB</g><r>HHHH</r><b>BBBB</b><u>UUUU</u></div></pre>'
    """
    html = '<pre><div>'
    if ts.count('H') == 0:
        return html
    ts2tag = {'1': ['<b>', '</b>'], '2': ['<g>', '</g>'], 'h': ['<r>', '</r>'], 'H': ['<r>', '</r>'],
              'U': ['<u>', '</u>']}
    split_points = [0] + [i + 1 for i in range(len(ts) - 1) if ts[i] != ts[i + 1]] + [len(seq)]
    split_segs = [[split_points[i], split_points[i + 1]] for i in range(len(split_points) - 1)]
    for seg in split_segs:
        html += ts2tag[ts[min(seg[0] + 1, len(ts)-1)]][0] + seq[seg[0]:seg[1]] + \
                ts2tag[ts[min(seg[0] + 1, len(ts)-1)]][1]
    html += '</div></pre>'
    return html


def protter_api(te, wgp):
    """
    :param te: topo_entry
    :param wgp: WinGradePath
    :return:
    """
    import urllib
    signal_peptide = [1, te.seq.count('u')] if 'u' in te.seq else [0, 0]
    wins = [w.seq if w.direction == 'fwd' else w.seq[::-1] for w in wgp.path]

    if wgp.path[0].direction == 'fwd' and signal_peptide == [0, 0]:
        nterm = 'intra'
    elif wgp.path[0].direction == 'rev' and signal_peptide == [0, 0]:
        nterm = 'extra'
    elif wgp.path[0].direction == 'fwd' and signal_peptide != [0, 0]:
        nterm = 'extra'
    elif wgp.path[0].direction == 'rev' and signal_peptide != [0, 0]:
        nterm = 'intra'

    if signal_peptide != [0, 0]:
        query = 'http://wlab.ethz.ch/protter/create?seq=%s&nterm=%s&tm=%s&mc=lightsalmon&lc=blue&tml=none&tex=;&' \
                'n:positives,s:circ,bc:cornflowerblue=R,K&n:signal peptide,cc:white,fc:mediumvioletred,bc:red=%i-%i&format=png' % \
                (te.original_seq, nterm, ','.join(wins), signal_peptide[0], signal_peptide[1])
    else:
        query = 'http://wlab.ethz.ch/protter/create?seq=%s&nterm=%s&tm=%s&mc=lightsalmon&lc=blue&tml=none&tex=;&' \
                'n:positives,s:circ,bc:cornflowerblue=R,K&format=png' % (te.original_seq, nterm, ','.join(wins))
    urllib.urlretrieve(query, "%s%s.png" % (te.param_list['out_path'], te.name))
    if signal_peptide != [0, 0]:
        link = 'http://wlab.ethz.ch/protter/#seq=%s&nterm=%s&tm=%s&mc=lightsalmon&lc=blue&tml=none&tex=;&' \
                'n:positives,s:circ,bc:cornflowerblue=R,K&n:signal peptide,cc:white,fc:mediumvioletred,bc:red=%i-%i&format=png' \
                % (te.original_seq, nterm, ','.join(wins), signal_peptide[0], signal_peptide[1])
    else:
        link = 'http://wlab.ethz.ch/protter/#seq=%s&nterm=%s&tm=%s&mc=lightsalmon&lc=blue&tml=none&tex=;&' \
                'n:positives,s:circ,bc:cornflowerblue=R,K&format=png' % (te.original_seq, nterm, ','.join(wins))
    return link, "%s.png" % te.name


if __name__ == '__main__':
    ts = '1111HHHH2222HHHH1111UUUU'
    se = 'AAAAHHHHBBBBHHHHBBBBHHHH'
    a = topo_string2html(ts, se)
    b = '<div><b>AAAA</b><r>HHHH</r><g>BBBB</g><r>HHHH</r><b>BBBB</b><p>UUUU</p></div>'
    print 'result', a
    print 'expect', b
    print a == b