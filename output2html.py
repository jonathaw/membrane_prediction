#!/usr/bin/env python
# coding=utf-8
"""
a group of functions to create HTML output for a prediction
"""


def create_html(topo_entry, best_path, sec_path, wins):
    from ProcessEntry import topo_string_rostlab_format
    import os
    # print 'making win plot'
    # make_win_plot(topo_entry, best_path, wins)
    best_ts = topo_string_rostlab_format(topo_entry, best_path)
    sec_ts = topo_string_rostlab_format(topo_entry, sec_path)
    print 'creating html at %s' % topo_entry.name + '.html'
    with open(topo_entry.name + '.html', 'wr+') as html:
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
        html.write('<p>Files to download</p>\n')
        # what is the path to use for download from the browesr???
        html.write('<a href="file:///%s">prediction file</a>\n' % (os.getcwd()+'/'+topo_entry.name+'.prd'))
        html.write('<p style="text-align:right">The Fleishman Lab &copy;</p>\n')

        html.write('</all>\n')
        html.write('</body>\n')
        html.write('</html>\n')


def path2html_table(wgp):
    html = '<p>Found %i %s:</p>\n' % (len(wgp.path), 'TMHs' if len(wgp.path) > 1 else 'TMH')
    html += '<table style="width:80%">\n'
    html += '<tr><th>#</th><th>Begin</th><th>End</th><th>N terminus orientation</th><th>&#916;G<sup>apparent</sup>' \
            '[kcal/mol]</th><th>Sequence</th><th>#charges</th><th>Membrane-deformation</th></tr>\n'
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
    ts2tag = {'1': ['<b>', '</b>'], '2': ['<g>', '</g>'], 'h': ['<r>', '</r>'], 'H': ['<r>', '</r>'],
              'U': ['<u>', '</u>']}
    split_points = [0] + [i + 1 for i in range(len(ts) - 1) if ts[i] != ts[i + 1]] + [len(seq)]
    split_segs = [[split_points[i], split_points[i + 1]] for i in range(len(split_points) - 1)]
    for seg in split_segs:
        html += ts2tag[ts[seg[0] + 1]][0] + seq[seg[0]:seg[1]] + ts2tag[ts[seg[0] + 1]][1]
    html += '</div></pre>'
    return html


if __name__ == '__main__':
    ts = '1111HHHH2222HHHH1111UUUU'
    se = 'AAAAHHHHBBBBHHHHBBBBHHHH'
    a = topo_string2html(ts, se)
    b = '<div><b>AAAA</b><r>HHHH</r><g>BBBB</g><r>HHHH</r><b>BBBB</b><p>UUUU</p></div>'
    print 'result', a
    print 'expect', b
    print a == b