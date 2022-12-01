import unittest
from test import support
import charapy as char
import numpy as np

class MyTestCase1(unittest.TestCase):

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

if __name__ == '__test__':
    unittest.main()