"""
:module: UnitTest_overrp.py
:platform: Unix, Windows
:synopsis: Defines the unit test for the inverse potential classs

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, May2014
.. history:
..                Alex Upah <alexupah@iastate.edu> July 2022
..                  -rewrote test in proper unit test form
"""


import numpy as np
from hoodlt.Potentials.overrp_pot_cutoff import OverrpCut
import unittest


class Testoverrp(unittest.TestCase):
     
    def setUp(self):
        self.p_exp = [12, 10]
        self.sigma = [1.0, 3.3]
        self.eps = [1, 5.0]
        self.cut_off = [4.1, 5.1]
        self.eps_zero = 1.e-12
        
        # These values have been computed using maple
        self.rarg = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.5]
        self.rval = [[1.677721600e7, 4096.000000, 31.56929180, 1.000000000, 0.6871947673e-1, 0.7707346629e-2,
                      0.1212113856e-2, 0.2441406250e-3, 0.1677721600e-4, 0.1881676423e-5, 0.0],
                     [8.029884828e11, 7.841684413e8, 1.359868045e7, 7.657894934e5, 82226.02063, 13279.96141,
                      2842.686167, 747.8413009, 80.29884828, 12.96871230, .2248973118]]
            
        self.rg1 = [[-3.221225472e9, -1.966080000e5, -673.4782249, -12.00000000, -.5277655813, -0.4110584869e-1,
                     -0.4749507354e-2, -0.7324218750e-3, -0.3221225472e-4, -0.2508901897e-5, 0.0],
                    [-1.284781574e14, -3.136673765e10, -2.417543197e8, -7.657894934e6, -5.262465328e5,
                     -59022.05070, -9282.240557, -1869.603255, -128.4781574, -14.40968035, -.1110604013]]
           
        self.rg2 = [[4.50971566080e10, 2.75251200000e6, 9428.69514898, 168.000000000, 7.38871813865, .575481881651,
                    0.664931029530e-1, 0.102539062500e-1, 0.450971566080e-3, 0.351246265657e-4, 0.0],
                    [1.54173788882e15, 3.76400851763e11, 2.90105183608e9, 9.18947392000e7, 6.31495839262e6,
                     7.08264608420e5, 1.11386886670e5, 22435.2390625, 1541.73788882, 172.916164165, 1.33272481502]]
           
        self.rval_press = [[2.01326592000e8, 49152.0000000, 378.831501521, 12.0000000000, .824633720832,
                            0.924881595511e-1, 0.145453662709e-1, 0.292968750000e-2, 0.201326592000e-3,
                            0.225801170779e-4, 0.0],
                           [8.02988483831e12, 7.84168441242e9, 1.35986804828e8, 7.65789493400e6, 8.22260207443e5,
                            1.32799614090e5, 28426.8617046, 7478.41302148, 802.988483831, 129.687123135, 2.24897312554]]
                  
    # due to legacy issues these values have the opposite sign
        self.rg1_press = [[-4.50971566080e10, -2.75251200000e6, -9428.69514898, -168.000000000, -7.38871813865,
                           -.575481881651, -0.664931029530e-1, -0.102539062500e-1, -0.450971566080e-3,
                           -0.351246265657e-4, 0.0],
                          [-1.54173788882e15, -3.76400851763e11, -2.90105183608e9, -9.18947392000e7, -6.31495839262e6,
                           -7.08264608420e5, -1.11386886670e5, -22435.2390625, -1541.73788882,
                           -172.916164165, -1.33272481502]]
                 
        self.rg2_press = [[6.31360192512e11, 3.85351680000e7, 1.32001732086e5, 2352.00000000, 103.442053941,
                           8.05674634311, .930903441342, .143554687500, 0.631360192512e-2, 0.491744771919e-3, 0.0],
                          [1.85008546659e16, 4.51681022116e12, 3.48126220330e10, 1.10273687040e9, 7.57795007114e7,
                           8.49917530104e6, 1.33664264003e6, 2.69222868750e5, 18500.8546659, 2074.99396998,
                           15.9926977803]]
         
    def test_energy(self):
        for j in range(2):
            vp = OverrpCut(self.p_exp[j], self.sigma[j], self.eps[j], self.cut_off[j], self.eps_zero)
            for ind in range(len(self.rarg)):
                r_arg = np.array([self.rarg[ind], 0, 0])
                if np.abs(self.rval[j][ind]) > 1.0:
                    rat = np.abs(self.rval[j][ind])
                else:
                    rat = 1.0
                
                self.assertAlmostEqual((self.rval[j][ind] - vp.e1(r_arg))/rat, 0.0)
                
    def test_g_1(self):
        for j in range(2):
            vp = OverrpCut(self.p_exp[j], self.sigma[j], self.eps[j], self.cut_off[j], self.eps_zero)
            for ind in range(len(self.rarg)):
                r_arg = np.array([self.rarg[ind], 0, 0])
                if np.abs(self.rg1[j][ind]) > 1.0:
                    val = np.abs(self.rg1[j][ind] - vp.g1(r_arg))/np.abs(self.rg1[j][ind])
                else:
                    val = np.abs(self.rg1[j][ind] - vp.g1(r_arg))
                
                self.assertAlmostEqual(val, 0.0)
                
    def test_g_2(self):
        for j in range(2):
            vp = OverrpCut(self.p_exp[j], self.sigma[j], self.eps[j], self.cut_off[j], self.eps_zero)
            for ind in range(len(self.rarg)):
                if np.abs(self.rg2[j][ind]) > 1.0:
                    rat = np.abs(self.rg2[j][ind])
                else:
                    rat = 1.0
                
                val_mat = vp.g2(np.array([self.rarg[ind], 0, 0]))
                val = np.abs(self.rg2[j][ind] - val_mat) / rat
            
                self.assertAlmostEqual(val, 0.0)
                
    def test_pressure(self):
        for j in range(2):
            vp = OverrpCut(self.p_exp[j], self.sigma[j], self.eps[j], self.cut_off[j], self.eps_zero)
            for ind in range(len(self.rarg)):
                r_arg = np.array([self.rarg[ind], 0, 0])
                if np.abs(self.rval_press[j][ind]) > 1.0:
                    rat = np.abs(self.rval_press[j][ind])
                else:
                    rat = 1.0
                    
                val = np.abs(self.rval_press[j][ind] - vp.e2(r_arg)) / rat
                
                self.assertAlmostEqual(val, 0.0)
                
    def test_g1_press(self):
        for j in range(2):
            vp = OverrpCut(self.p_exp[j], self.sigma[j], self.eps[j], self.cut_off[j], self.eps_zero)
            for ind in range(len(self.rarg)):
                r_arg = np.array([self.rarg[ind], 0, 0])
                if np.abs(self.rg1_press[j][ind]) > 1.0:
                    rat = np.abs(self.rg1_press[j][ind])
                else:
                    rat = 1.0
                
                val = np.abs(-self.rg1_press[j][ind] - vp.g2(r_arg)) / rat
                
                self.assertAlmostEqual(val, 0.0)
            
    def test_g2_press(self):
        for j in range(2):
            vp = OverrpCut(self.p_exp[j], self.sigma[j], self.eps[j], self.cut_off[j], self.eps_zero)
            for ind in range(len(self.rarg)):
                val_mat = vp.g3(np.array([self.rarg[ind], 0, 0]))
                if np.abs(self.rg2_press[j][ind]) > 1.0:
                    rat = np.abs(self.rg2_press[j][ind])
                else:
                    rat = 1.0
                    
                val = np.abs(-self.rg2_press[j][ind] - val_mat) / rat
                
                self.assertAlmostEqual(val, 0.0)
                
                
if __name__ == '__main__':
    unittest.main()
