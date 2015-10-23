#!/usr/bin/env python2

import os, sys
import unittest
import numpy as np
from numpy.testing import assert_allclose
from numpy import linalg
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
                         os.path.abspath(__file__)))))
from dimple.cell import Cell, Mat3

def to_np(m):
    return np.array([[m[0], m[1], m[2]],
                     [m[3], m[4], m[5]],
                     [m[6], m[7], m[8]]])

class TestMat3(unittest.TestCase):
    def setUp(self):
        self.mat = Mat3(1.3, 2.3, 3.3, 4.3, 1.0, 6.3, 7.3, 8.3, 9.3)
        self.mat2 = Mat3(1.9, 2.9, 3.9, 4.9, 1.9, 6.9, 7.9, 8.9, 9.9)
        self.np_mat = to_np(self.mat)
        self.np_mat2 = to_np(self.mat2)
    def test_identity(self):
        assert_allclose(np.identity(3), to_np(Mat3.identity()))
    def test_det(self):
        #print self.mat.det()
        self.assertAlmostEqual(linalg.det(self.np_mat), self.mat.det())
    def test_transpose(self):
        assert_allclose(self.np_mat.transpose(), to_np(self.mat.transpose()))
    def test_dot(self):
        assert_allclose(np.dot(self.np_mat, self.np_mat2),
                        to_np(self.mat.dot(self.mat2)))
    def test_add(self):
        assert_allclose(self.np_mat+self.np_mat2, to_np(self.mat+self.mat2))
    def test_sub(self):
        assert_allclose(self.np_mat-self.np_mat2, to_np(self.mat-self.mat2))
    def test_norm1(self):
        self.assertAlmostEqual(linalg.norm(self.np_mat, 1),
                               self.mat.induced_1norm())
    def test_norm2(self):
        self.assertAlmostEqual(linalg.norm(self.np_mat, 2),
                               self.mat.euclidean_norm())

class TestCell(unittest.TestCase):
    def setUp(self):
        self.cell = Cell((22.84, 32.84, 42.84, 80.84, 90.84, 100.84))
    def test_orth(self):
        orth = self.cell.get_orth_matrix()
        expected = np.array([[22.84, -6.17612, -0.628045],
                             [0,     32.254,    6.82343],
                             [0,      0,       42.2884]])
        assert_allclose(to_np(orth), expected, rtol=1e-6)

    def test_frac(self):
        frac = self.cell.get_frac_matrix()
        orth = self.cell.get_orth_matrix()
        assert_allclose(to_np(frac.dot(orth)), np.identity(3), atol=1e-12)

if __name__ == '__main__':
    unittest.main()

