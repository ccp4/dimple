from dimple.contaminants.data import DATA
from dimple.cell import Cell, calculate_difference

class ContamCell(Cell):
    def __init__(self, entry):
        Cell.__init__(self, entry[2:8], symmetry=entry[1])
        self.pdb_id = entry[0]
        self.uniprot_src = entry[8]  # we may have a homolog of uniprot entry

def find_similar_rel(cell, rel_tol):
    rel_min = 1 - rel_tol
    rel_max = 1 + rel_tol
    similar_cells = []
    for entry in DATA:
        if entry[2] < rel_min * cell.a:
            continue
        if entry[2] > rel_max * cell.a:
            break  # because entries are sorted by a
        if not rel_min * cell.b < entry[3] < rel_max * cell.b:
            continue
        c = ContamCell(entry)
        if calculate_difference(cell, c) < rel_tol:
            similar_cells.append(c)
    similar_cells.sort(key=lambda x: calculate_difference(cell, x))
    return similar_cells

def check_and_print(cell, rel_tol=0.05):
    similar = find_similar_rel(cell, rel_tol)
    n = len(similar)
    if n == 0:
        return
    n_str = ('%d entries' % n if n > 1 else '1 entry')
    print(' Our little contaminant list has %s with similar unit cell:' % n_str)
    for c in similar:
        print('  %(pdb_id)s (%(uniprot_src)s): %(symmetry)s ' % c.__dict__ +
              '(%s)' % c.parameters_as_str())
