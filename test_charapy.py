import unittest
from test import support
import charapy as char
import numpy as np

class MyTestCaseFT(unittest.TestCase):

    def setUp(self):
        #borrar si no lo necesite
        ...

    def tearDown(self):
        #borrar si no lo necesite
        ... 

    def test_fitAB(self):
        foo_fit = char.Foo_fit()
        x = np.array([1,2])
        y = np.array([2,3])
        A = 2.4663
        B = -0.7095
        
        self.assertEqual(foo_fit.fit_AB(x,y)[0][0].round(4),A)
        self.assertEqual(foo_fit.fit_AB(x,y)[0][1].round(4),B)
        self.assertEqual(foo_fit.fit_LM(y,x)[0][0].round(4),A)
        self.assertEqual(foo_fit.fit_LM(y,x)[0][1].round(4),B)

class MyTestCaseDP(unittest.TestCase):

    def test_carbonnumber(self):
        dist_pedersen = char.Distribution_pedersen()
        x = 2
        c1 = 1
        c2 = 2
        cod = 2.6931
        self.assertEqual(dist_pedersen.carbon_number(x,c1,c2).round(4),cod)
        self.assertEqual(dist_pedersen.p_density(x,c1,c2).round(4),cod)

if __name__ == '__test__':
    unittest.main()