from math import radians, sin, cos, acos, sqrt, pi


class Cell(object):
    def __init__(self, parameters):
        if parameters is None:
            self.cell = None
            return
        assert len(parameters) == 6
        self.a, self.b, self.c = parameters[:3]
        self.alpha, self.beta, self.gamma = parameters[3:]
        self.cell = parameters

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
        return Mat3( a, b*cos(gamma),  c*cos(beta),
                     0, b*sin(gamma), -c*sin(beta)*cos(alpha_star),
                     0, 0,             c*sin(beta)*sin(alpha_star))

    def get_frac_matrix(self):
        return self.get_orth_matrix().inverse()

    def max_shift_in_mapping(self, other):
        trans = self.get_frac_matrix().dot(other.get_orth_matrix())
        return (trans - Mat3.identity()).euclidean_norm()


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


