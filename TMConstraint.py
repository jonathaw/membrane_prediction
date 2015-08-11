class TMConstraint():
    """
    a class to describe constraints on a protein's topology
    """
    def __init__(self, name, tm_num=None, tm_pos=None, c_term=None):
        self.name = name
        self.tm_num = tm_num
        self.tm_pos = tm_pos
        self.c_term = c_term


def parse_cst(file):
    with open(file, 'r') as f:
        cont = f.read().split('\n')
    for ln in cont:
        sp = ln.split()
        print sp