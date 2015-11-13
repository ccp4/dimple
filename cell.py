from math import radians, sin, cos, acos, sqrt, pi


class Cell(object):
    def __init__(self, parameters, symmetry):
        if parameters is None:
            self.cell = None
            return
        assert isinstance(parameters, tuple)
        assert len(parameters) == 6
        self.a, self.b, self.c = parameters[:3]
        self.alpha, self.beta, self.gamma = parameters[3:]
        self.cell = parameters
        self.symmetry = symmetry  # international SG symbol w/ spaces

    def get_volume(self):
        ca = cos(radians(self.alpha))
        cb = cos(radians(self.beta))
        cg = cos(radians(self.gamma))
        return self.a * self.b * self.c * sqrt((1 - ca*ca - cb*cb - cg*cg) +
                                               2 * ca*cb*cg)

    def asu_volume(self, z=None):
        if z is None:
            z = calculate_z_order(self.symmetry)
        return self.get_volume() / z

    def __str__(self):
        return str(self.cell)

    # The orthogonalization matrix we use is described in ITfC B p.262:
    # "An alternative mode of orthogonalization, used by the Protein
    # Data Bank and most programs, is to align the a1 axis of the unit
    # cell with the Cartesian X1 axis, and to align the a*3 aixs with the
    # Cartesian X3 axis."
    def get_orth_matrix(self):
        a, b, c = self.a, self.b, self.c
        alpha = radians(self.alpha)
        beta = radians(self.beta)
        gamma = radians(self.gamma)
        alpha_star = acos((cos(gamma)*cos(beta) - cos(alpha)) /
                          (sin(beta)*sin(gamma)))
        return Mat3(a, b*cos(gamma),  c*cos(beta),
                    0, b*sin(gamma), -c*sin(beta)*cos(alpha_star),
                    0, 0,             c*sin(beta)*sin(alpha_star))

    def get_frac_matrix(self):
        return self.get_orth_matrix().inverse()

    def max_shift_in_mapping(self, other):
        trans = self.get_frac_matrix().dot(other.get_orth_matrix())
        return (trans - Mat3.identity()).euclidean_norm()

    # This affects only primitive orthorhombic (P 2x 2x 2x).
    # Convert the "reference" (symmetry-based) settings to the "standard"
    # (cell-based) settings. See the SETTING keyword in POINTLESS.
    def to_standard(self):
        if self.alpha == self.beta == self.gamma == 90 and (
                self.a > self.b or self.b > self.c):
            sym_splitted = self.symmetry.split()
            if len(sym_splitted) == 4:  # is this condition redundant?
                reordered = sorted(zip(self.cell[:3], sym_splitted[1:]))
                new_cell, new_symm = zip(*reordered)
                return Cell(new_cell + (90., 90., 90.),
                            symmetry=' '.join((sym_splitted[0],) + new_symm))
        return self

    # returns symmetry with screw axes removed (changed to rotation axes)
    def unscrew_symmetry(self):
        return ' '.join(a[0] for a in self.symmetry.split())

class Mat3(object):
    "Matrix 3x3"
    def __init__(self, *args):
        if len(args) == 1:
            self.m = tuple(args[0])
        else:
            self.m = args
        assert len(self.m) == 9

    def __getitem__(self, index):
        return self.m[index]
    def __str__(self):
        return "[%g %g %g; %g %g %g; %g %g %g]" % self.m
    def __repr__(self):
        return "Mat3" + str(self.m)

    def __add__(self, other):
        assert isinstance(other, Mat3)
        return Mat3(a+b for a,b in zip(self.m, other.m))
    def __sub__(self, other):
        assert isinstance(other, Mat3)
        return Mat3(a-b for a,b in zip(self.m, other.m))
    # scalar must be float
    def __mul__(self, scalar):
        assert isinstance(scalar, float)
        return Mat3(a*scalar for a in self.m)

    @staticmethod
    def identity():
        return Mat3(1, 0, 0,
                    0, 1, 0,
                    0, 0, 1)


    def transpose(self):
        m = self.m
        return Mat3(m[0], m[3], m[6],
                    m[1], m[4], m[7],
                    m[2], m[5], m[8])

    def dot(self, other):
        a = self.m
        b = other.m
        return Mat3(sum(a[3*row+i] * b[3*i+col] for i in range(3))
                    for row in range(3)
                    for col in range(3))

    def det(self):
        m = self.m
        return (m[0] * (m[4] * m[8] - m[5] * m[7]) -
                m[1] * (m[3] * m[8] - m[5] * m[6]) +
                m[2] * (m[3] * m[7] - m[4] * m[6]))

    def trace(self):
        m = self.m
        return m[0] + m[4] + m[8]

    def inverse(self):
        d = self.det()
        if d == 0:
            raise ValueError("Matrix is not invertible")
        m = self.m
        return Mat3(( m[4] * m[8] - m[5] * m[7]) / d,
                    (-m[1] * m[8] + m[2] * m[7]) / d,
                    ( m[1] * m[5] - m[2] * m[4]) / d,
                    (-m[3] * m[8] + m[5] * m[6]) / d,
                    ( m[0] * m[8] - m[2] * m[6]) / d,
                    (-m[0] * m[5] + m[2] * m[3]) / d,
                    ( m[3] * m[7] - m[4] * m[6]) / d,
                    (-m[0] * m[7] + m[1] * m[6]) / d,
                    ( m[0] * m[4] - m[1] * m[3]) / d)

    def induced_1norm(self):  # aka 1-norm
        m = self.m
        return max(abs(m[0]) + abs(m[3]) + abs(m[6]),
                   abs(m[1]) + abs(m[4]) + abs(m[7]),
                   abs(m[2]) + abs(m[5]) + abs(m[8]))

    def euclidean_norm(self):  # aka induced 2-norm
        A = self.dot(self.transpose())
        # now get the largest eigenvalue of A (A is symmetric)
        # http://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
        p1 = A[1]**2 + A[2]**2 + A[5]**2
        if p1 == 0:
            return max(A[0], A[4], A[8])
        q = A.trace() / 3.
        p2 = (A[0] - q)**2 + (A[4] - q)**2 + (A[8] - q)**2 + 2 * p1
        p = sqrt(p2 / 6.)
        B = (A - Mat3.identity() * q) * (1. / p)
        r = B.det() / 2.
        if r <= -1:
            phi = pi / 3.
        elif r >= 1:
            phi = 0
        else:
            phi = acos(r) / 3.
        eig1 = q + 2 * p * cos(phi)
        return sqrt(eig1)

# space group
def match_symmetry(meta1, meta2):
    if not meta1 or not meta2:
        return None
    def sig(sym):
        first_chars = [a[0] for a in sym.split()]
        s = first_chars[0] + ''.join(sorted(first_chars[1:]))
        if s == 'I112': # I2 is equivalent to C2
            return 'C112'
        # note: I see names such as 'H 3' used in ccp4, never 'R 3'
        if s[0] == 'R':
            s = 'H' + s[1:]
        return s
    return sig(meta1.symmetry) == sig(meta2.symmetry)


_centering_n = {'P': 1, 'A': 2, 'B': 2, 'C': 2, 'I': 2,
                'R': 3, 'H': 3, 'S': 3, 'T': 3, 'F': 4}

_pg_symop = {'1': 1,
             '2': 2, '121': 2,
             '222': 4,
             '4': 4,
             '422': 8,
             '3': 3,
             '32': 6, '312': 6, '321': 6,
             '6': 6,
             '622': 12,
             '23': 12,
             '432': 24,
            }

def calculate_z_order(hm):
    pg = ''.join(a[0] for a in hm.split()[1:])
    return _centering_n[hm[0]] * _pg_symop[pg]


