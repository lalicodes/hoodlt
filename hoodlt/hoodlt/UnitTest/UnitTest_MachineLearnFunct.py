"""
:module: UnitTest_MachineLearnFunct.py
:platform: Unix, Windows
:synopsis: Defines the unit test for force matching machine learning functions

.. moduleauthor:: Alex Travesset <trvsst@meslab.gov> , February 2024
"""

import unittest
import numpy as np
import numpy.linalg as la

import hoodlt.Utils.MachineLearnFunct as ML


class TestMachineLearnForceMatching(unittest.TestCase):

    def setUp(self):
        eta = 0.5
        r_c = 6.0
        self.param2 = np.array([eta, r_c])
        lamda = 1
        xi = 3
        mu = 0.12
        self.param3 = np.array([lamda, xi, mu])
        self.cut = 10.0

        self.ml_rot = ML.FRotG(self.param2, self.param3, self.cut)
        num = 5
        self.num = 5

        self.va = np.zeros([num, 3])
        self.vb = np.zeros_like(self.va)
        self.vc = np.zeros_like(self.va)

        self.va[:, 0] = np.arange(num) + np.sqrt(2)
        self.va[:, 1] = np.arange(num) - np.sqrt(3)
        self.va[:, 2] = np.arange(num) + np.exp(1)
        self.vb[:, 0] = np.arange(num) - np.pi
        self.vb[:, 1] = np.arange(num) - 1 / np.sqrt(7)
        self.vb[:, 2] = np.arange(num)
        self.vc[:, :] = self.va[:, :] - self.vb[:, :]

        self.ux = np.array([1.0, 0.0, 0.0])
        self.uy = np.array([0.0, 1.0, 0.0])
        self.uz = np.array([0.0, 0.0, 1.0])

    def test_fc(self):
        """
        test the cut-off function
        """

        r_vals = np.array([0.0, 1.0, 2.0, 4.0, 5.0, 9.999, 9.99999, 11.0, 12.00])
        f_vals = (np.tanh(1-r_vals/self.cut))**3
        f_vals[-1] = 0.0
        f_vals[-2] = 0.0

        cond = np.allclose(self.ml_rot.fc(r_vals), f_vals, atol=1e-14)

        self.assertTrue(cond)

    def test_hc(self):
        """
        test the cut-off function
        """

        r_vals = np.array([0.0, 1.0, 2.0, 4.0, 5.0, 9.999, 9.99999, 11.0, 12.00])
        f_vals = 3*(np.tanh(1 - r_vals / self.cut)/np.cosh(1-r_vals/self.cut)) ** 2 / self.cut
        f_vals[-1] = 0.0
        f_vals[-2] = 0.0

        cond = np.allclose(self.ml_rot.hc(r_vals), f_vals, atol=1e-14)

        self.assertTrue(cond)

        # compute the derivative approximately

        eps = 1e-6
        der_approx = (self.ml_rot.fc(r_vals+eps)-self.ml_rot.fc(r_vals))/eps

        cond = np.allclose(der_approx, -f_vals, atol=1e-8)
        self.assertTrue(cond)

    def test_cos(self):
        """
        test the cosinus of vectors
        """

        num = 6
        ang = np.array([0, np.pi/6, np.pi/3, np.pi/2, np.pi, np.pi+np.pi/6])
        arg_psi = 2*np.pi*np.arange(num)/num
        arg_theta = np.pi*np.arange(num)/num
        cos_comp = np.zeros([num, num])
        cos_result = np.zeros_like(cos_comp)
        vec_a = np.zeros(3)
        vec_b = np.zeros(3)
        ra = np.zeros([num, num, 3])
        rb = np.zeros_like(ra)

        for ind1 in range(num):
            for ind2 in range(num):
                vec_a[0] = np.sin(arg_theta[ind1])*np.cos(arg_psi[ind2])
                vec_a[1] = np.sin(arg_theta[ind1])*np.sin(arg_psi[ind2])
                vec_a[2] = np.cos(arg_theta[ind1])
                ra[ind1, ind2, :] = vec_a[:]
                vec_b[0] = np.sin(arg_theta[ind1]+ang[ind1])*np.cos(arg_psi[ind2])
                vec_b[1] = np.sin(arg_theta[ind1]+ang[ind1])*np.sin(arg_psi[ind2])
                vec_b[2] = np.cos(arg_theta[ind1]+ang[ind1])
                rb[ind1, ind2, :] = vec_b[:]

        cos_comp[:, :] = self.ml_rot.cos_theta(ra, rb)
        cos_result[:] = np.cos(ang[:num])

        cond = np.allclose(cos_comp, np.transpose(cos_result), atol=1e-9)
        self.assertTrue(cond)

        for ind1 in range(num):
            for ind2 in range(num):
                vec_a[0] = np.cos(arg_psi[ind2])
                vec_a[1] = np.sin(arg_psi[ind2])
                vec_a[2] = 0.0
                ra[ind1, ind2, :] = vec_a[:]
                vec_b[0] = np.cos(arg_psi[ind2]+ang[ind1])
                vec_b[1] = np.sin(arg_psi[ind2]+ang[ind1])
                vec_b[2] = 0.0
                rb[ind1, ind2, :] = vec_b[:]

        cos_comp[:, :] = self.ml_rot.cos_theta(ra, rb)
        cos_result[:] = np.cos(ang[:num])

        cond = np.allclose(cos_comp, np.transpose(cos_result), atol=1e-9)
        self.assertTrue(cond)

    def test_omega(self):
        """
        test omega function
        """

        vals_comp = np.zeros(self.num)
        vals_result = np.zeros_like(vals_comp)

        lamda, xi, mu = self.param3[:]

        term1 = (1+lamda*self.ml_rot.cos_theta(self.va, self.vb))**xi
        term2 = np.exp(-mu*(la.norm(self.va, axis=1)**2+la.norm(self.vb, axis=1)**2+la.norm(self.vc, axis=1)**2))
        term3_1 = self.ml_rot.fc(la.norm(self.va, axis=1))*self.ml_rot.fc(la.norm(self.vb, axis=1))
        term3 = term3_1*self.ml_rot.fc(la.norm(self.vc, axis=1))

        vals_result[:] = term1[:]*term2[:]*term3[:]
        vals_comp[:] = self.ml_rot.omega(self.va, self.vb)

        cond = np.allclose(vals_comp, vals_result, atol=1e-9)
        self.assertTrue(cond)

        # check other values of lambda, xi, mu
        lamda_1 = -1
        xi_1 = 4
        mu_1 = 0.01
        param3_1 = np.array([lamda_1, xi_1, mu_1])

        ter1 = 1.0*(1+lamda_1*self.ml_rot.cos_theta(self.va, self.vb))**xi_1
        ter2 = np.exp(-mu_1*(la.norm(self.va, axis=1)**2+la.norm(self.vb, axis=1)**2+la.norm(self.vc, axis=1)**2))
        ter3_1 = self.ml_rot.fc(la.norm(self.va, axis=1)) * self.ml_rot.fc(la.norm(self.vb, axis=1))
        ter3 = ter3_1 * self.ml_rot.fc(la.norm(self.vc, axis=1))

        vals_result[:] = ter1[:] * ter2[:] * ter3[:]
        frot = ML.FRotG(self.param2, param3_1, self.cut)
        vals_comp = frot.omega(self.va, self.vb)
        cond = np.allclose(vals_comp, vals_result, atol=1e-9)
        self.assertTrue(cond)

    def test_phi(self):
        """
        test phi function
        """

        vals_comp = np.zeros(self.num)
        vals_result = np.zeros_like(vals_comp)

        lamda, xi, mu = self.param3[:]

        term1 = 1/la.norm(self.vb, axis=1)-self.ml_rot.cos_theta(self.va, self.vb)/la.norm(self.va, axis=1)
        term2 = lamda*xi/(1+lamda*self.ml_rot.cos_theta(self.va, self.vb))

        vals_result[:] = term1[:]*term2[:]/la.norm(self.va, axis=1)
        vals_comp = self.ml_rot.phi(self.va, self.vb)

        cond = np.allclose(vals_comp, vals_result, atol=1e-9)
        self.assertTrue(cond)

        # check other values of lambda, xi, mu
        lamda_1 = -1
        xi_1 = 6
        mu_1 = 0.05
        param3_1 = np.array([lamda_1, xi_1, mu_1])

        term2 = lamda_1*xi_1/(1 + lamda_1*self.ml_rot.cos_theta(self.va, self.vb))

        vals_result[:] = term1[:]*term2[:]/la.norm(self.va, axis=1)

        frot = ML.FRotG(self.param2, param3_1, self.cut)
        vals_comp = frot.phi(self.va, self.vb)

        cond = np.allclose(vals_comp, vals_result, atol=1e-9)
        self.assertTrue(cond)

    def test_psi(self):
        """
        test psi function
        """

        vals_comp = np.zeros(self.num)
        vals_result = np.zeros_like(vals_comp)

        lamda, xi, mu = self.param3[:]

        term1 = 1/(la.norm(self.va, axis=1)*la.norm(self.vb, axis=1))
        term2 = lamda*xi/(1+lamda*self.ml_rot.cos_theta(self.va, self.vb))

        vals_result[:] = -term1[:]*term2[:]
        vals_comp = self.ml_rot.psi(self.va, self.vb)

        cond = np.allclose(vals_comp, vals_result, atol=1e-9)
        self.assertTrue(cond)

        # check other values of lambda, xi, mu
        lamda_1 = -1
        xi_1 = 1
        mu_1 = 0.005
        param3_1 = np.array([lamda_1, xi_1, mu_1])

        term2 = lamda_1*xi_1/(1+lamda_1*self.ml_rot.cos_theta(self.va, self.vb))

        vals_result[:] = -term1[:]*term2[:]

        frot = ML.FRotG(self.param2, param3_1, self.cut)
        vals_comp = frot.psi(self.va, self.vb)

        cond = np.allclose(vals_comp, vals_result, atol=1e-9)
        self.assertTrue(cond)

    def test_g2(self):
        """
        test g2 function
        """

        term1 = np.exp(-self.param2[0]*(la.norm(self.va, axis=1)-self.param2[1])**2)
        vals_comp = np.zeros(self.num)
        vals_result = np.zeros_like(vals_comp)

        vals_result[:] = term1[:]*self.ml_rot.fc(la.norm(self.va, axis=1))
        vals_comp = self.ml_rot.g2(self.va)

        cond = np.allclose(vals_comp, vals_result, atol=1e-9)
        self.assertTrue(cond)

    def test_der_g2(self):
        """
        test derivative of g2 function
        """

        vals_comp = self.ml_rot.der_g2_1(self.va)
        vals_result = np.zeros_like(vals_comp)

        eps = 1e-6

        vals_result[:, 0] = (self.ml_rot.g2(self.va + eps*self.ux) - self.ml_rot.g2(self.va))/eps
        vals_result[:, 1] = (self.ml_rot.g2(self.va + eps*self.uy) - self.ml_rot.g2(self.va))/eps
        vals_result[:, 2] = (self.ml_rot.g2(self.va + eps*self.uz) - self.ml_rot.g2(self.va))/eps

        cond = np.allclose(vals_comp, vals_result, atol=1e-9)
        self.assertTrue(cond)

    def test_g3(self):
        """
        test g3 function
        """

        vals_comp = np.zeros(self.num)
        vals_result = np.zeros_like(vals_comp)

        lamda, xi, mu = self.param3[:]

        term1 = 2**(1-xi)*(1+lamda*self.ml_rot.cos_theta(self.va, self.vb))**xi
        term2 = 1.0*np.exp(-mu*(la.norm(self.va, axis=1)**2+la.norm(self.vb, axis=1)**2+la.norm(self.vc, axis=1)**2))
        term3_1 = self.ml_rot.fc(la.norm(self.va, axis=1))*self.ml_rot.fc(la.norm(self.vb, axis=1))
        term3 = term3_1*self.ml_rot.fc(la.norm(self.vc, axis=1))

        vals_comp[:] = self.ml_rot.g3(self.va, self.vb)
        vals_result[:] = term1[:]*term2[:]*term3[:]

        cond = np.allclose(vals_comp, vals_result, atol=1e-9)
        self.assertTrue(cond)

    def test_der_g3_1(self):
        """
        derivative of g3
        """

        vals_comp = self.ml_rot.der_g3_1(self.va, self.vb)
        vals_result = np.zeros_like(vals_comp)

        eps = 1e-7

        d_va_x = self.va+eps*self.ux
        d_va_y = self.va+eps*self.uy
        d_va_z = self.va+eps*self.uz

        d_vb_x = self.vb+eps*self.ux
        d_vb_y = self.vb+eps*self.uy
        d_vb_z = self.vb+eps*self.uz

        vals_result[:, 0] = (self.ml_rot.g3(d_va_x, d_vb_x) - self.ml_rot.g3(self.va, self.vb))/eps
        vals_result[:, 1] = (self.ml_rot.g3(d_va_y, d_vb_y) - self.ml_rot.g3(self.va, self.vb))/eps
        vals_result[:, 2] = (self.ml_rot.g3(d_va_z, d_vb_z) - self.ml_rot.g3(self.va, self.vb))/eps

        cond = np.allclose(vals_comp, vals_result, atol=1e-9)
        self.assertTrue(cond)

    def test_der_g3_2(self):
        """
        derivative of g3
        """

        vals_comp = self.ml_rot.der_g3_2(self.va, self.vb)
        vals_result = np.zeros_like(vals_comp)

        eps = 1e-7

        vals_result[:, 0] = (self.ml_rot.g3(self.va-eps*self.ux, self.vb) - self.ml_rot.g3(self.va, self.vb)) / eps
        vals_result[:, 1] = (self.ml_rot.g3(self.va-eps*self.uy, self.vb) - self.ml_rot.g3(self.va, self.vb)) / eps
        vals_result[:, 2] = (self.ml_rot.g3(self.va-eps*self.uz, self.vb) - self.ml_rot.g3(self.va, self.vb)) / eps

        cond = np.allclose(vals_comp, vals_result, atol=1e-9)
        self.assertTrue(cond)

    def test_der_g3_3(self):
        """
        third derivative of g3
        """

        vals_comp = self.ml_rot.der_g3_3(self.va, self.vb)
        vals_result = np.zeros_like(vals_comp)

        eps = 1e-7

        vals_result[:, 0] = (self.ml_rot.g3(self.va, self.vb-eps*self.ux) - self.ml_rot.g3(self.va, self.vb)) / eps
        vals_result[:, 1] = (self.ml_rot.g3(self.va, self.vb-eps*self.uy) - self.ml_rot.g3(self.va, self.vb)) / eps
        vals_result[:, 2] = (self.ml_rot.g3(self.va, self.vb-eps*self.uz) - self.ml_rot.g3(self.va, self.vb)) / eps

        cond = np.allclose(vals_comp, vals_result, atol=1e-9)
        self.assertTrue(cond)

    def test_g3_all_derivatives(self):
        """
        Test that all derivatives add up to zero
        """

        vals_1 = self.ml_rot.der_g3_1(self.va, self.vb)
        vals_2 = self.ml_rot.der_g3_2(self.va, self.vb)
        vals_3 = self.ml_rot.der_g3_3(self.va, self.vb)

        vals_tot = vals_1+vals_2+vals_3
        vals_res = np.zeros_like(vals_tot)

        cond = np.allclose(vals_tot, vals_res, atol=1e-20)
        self.assertTrue(cond)

    def test_derivative_l2(self):
        """
        Test the derivative of l2
        """

        num_points = 3
        num_ncs = 3
        pos = np.zeros([num_points, num_ncs, 3])

        pos[0, 0, :] = np.array([0.0, 0.0, 2.0])
        pos[0, 1, :] = np.array([2.0, 0.0, 0.0])
        pos[0, 2, :] = np.array([0.0, -2.0, 0.0])

        pos[1, 0, :] = np.array([2.0, 2.0, 2.0])
        pos[1, 1, :] = np.array([2.0, -3.0, 3.0])
        pos[1, 2, :] = np.array([0.5, 0.0, -0.5])

        pos[2, 0, :] = np.array([0.0, -0.2, 2.0])
        pos[2, 1, :] = np.array([0.1, 0.1, 0.3])
        pos[2, 2, :] = np.array([0.2, 0.0, -0.15])

        eps = 1e-7
        der_comp = np.zeros([pos.shape[0], pos.shape[1], pos.shape[1], pos.shape[2]])
        id_p = np.identity(3)[np.newaxis, :, :]
        d_pos = np.zeros_like(pos)
        for ind1 in range(num_ncs):
            for ind2 in range(num_ncs):
                for ind3 in range(3):
                    d_pos[:, :, :] = pos[:, :, :]
                    d_pos[:, ind2, :] += eps*id_p[..., :, ind3]
                    der_comp[:, ind1, ind2, ind3] = (self.ml_rot.l2(ind1, d_pos)[:]-self.ml_rot.l2(ind1, pos)[:])/eps

        der_result = np.zeros_like(der_comp)
        for ind1 in range(num_ncs):
            for ind2 in range(num_ncs):
                der_result[:, ind1, ind2, :] = self.ml_rot.der_l2(ind2, ind1, pos)[:, :]

        cond = np.allclose(der_result, der_comp, atol=1e-8)
        self.assertTrue(cond)

    def test_derivative_psi_l2(self):
        """
        test the derivative of psi_l2
        """

        num_samples = 3
        num_ncs = 3
        pos = np.zeros([num_samples, num_ncs, 3])

        pos[0, 0, :] = np.array([0.0, 0.0, 2.0])
        pos[0, 1, :] = np.array([2.0, 0.0, 0.0])
        pos[0, 2, :] = np.array([0.0, -2.0, 0.0])

        pos[1, 0, :] = np.array([2.0, 2.0, 2.0])
        pos[1, 1, :] = np.array([2.0, -3.0, 3.0])
        pos[1, 2, :] = np.array([0.5, 0.0, -0.5])

        pos[2, 0, :] = np.array([0.0, -0.2, 2.0])
        pos[2, 1, :] = np.array([0.1, 0.1, 0.3])
        pos[2, 2, :] = np.array([0.2, 0.0, -0.15])

        eps = 1e-7
        der_comp = np.zeros([pos.shape[0], pos.shape[1], pos.shape[2]])
        id_p = np.identity(3)[np.newaxis, :, :]
        d_pos = np.zeros_like(pos)
        for ind1 in range(num_ncs):
            for ind2 in range(3):
                d_pos[:, :, :] = pos[:, :, :]
                d_pos[:, ind1, :] += eps*id_p[..., :, ind2]
                der_comp[:, ind1, ind2] = -(self.ml_rot.psi_l2(d_pos)[:]-self.ml_rot.psi_l2(pos)[:])/eps

        der_result = np.zeros_like(der_comp)
        der_result[:, :, :] = self.ml_rot.der_psi_l2(pos)[:, :, :]

        cond = np.allclose(der_result, der_comp, atol=1e-8)
        self.assertTrue(cond)

    def test_derivative_l3(self):
        """
        Test the derivative of l3
        """

        num_points = 4
        num_ncs = 4
        pos = np.zeros([num_points, num_ncs, 3])

        pos[0, 0, :] = np.array([0.0, 0.0, 2.0])
        pos[0, 1, :] = np.array([2.0, 0.0, 0.0])
        pos[0, 2, :] = np.array([0.0, -2.0, 0.0])
        pos[0, 3, :] = np.array([0.0, -2.0, 2.0])

        pos[1, 0, :] = np.array([2.0, 2.0, 2.0])
        pos[1, 1, :] = np.array([2.0, -3.0, 3.0])
        pos[1, 2, :] = np.array([0.5, 0.0, -0.5])
        pos[1, 3, :] = np.array([0.5, 0.5, -0.5])

        pos[2, 0, :] = np.array([0.0, -0.2, 2.0])
        pos[2, 1, :] = np.array([0.1, 0.1, 0.3])
        pos[2, 2, :] = np.array([0.2, 0.0, -0.15])
        pos[2, 3, :] = np.array([-0.2, 0.0, -0.15])

        pos[3, 0, :] = np.array([0.0, 0.0, 0.0])
        pos[3, 1, :] = np.array([0.0, 0.0, 4.0])
        pos[3, 2, :] = np.array([4.0, 0.0, 0.0])
        pos[3, 3, :] = np.array([0.0, 4.0, 0.0])

        eps = 1e-7
        der_comp = np.zeros([pos.shape[0], pos.shape[1], pos.shape[1], pos.shape[2]])
        id_p = np.identity(3)[np.newaxis, :, :]
        d_pos = np.zeros_like(pos)
        for ind1 in range(num_ncs):
            for ind2 in range(num_ncs):
                for ind3 in range(3):
                    d_pos[:, :, :] = pos[:, :, :]
                    d_pos[:, ind2, :] += eps * id_p[..., :, ind3]
                    der_comp[:, ind1, ind2, ind3] = (self.ml_rot.l3(ind1, d_pos)[:]-self.ml_rot.l3(ind1, pos)[:])/eps

        der_result = np.zeros_like(der_comp)
        for ind1 in range(num_ncs):
            for ind2 in range(num_ncs):
                der_result[:, ind1, ind2, :] = self.ml_rot.der_l3(ind2, ind1, pos)[:, :]

        cond = np.allclose(der_result, der_comp, atol=1e-8)
        self.assertTrue(cond)

    def test_derivative_psi_l3(self):
        """
        test the derivative of psi_l2
        """

        num_points = 4
        num_ncs = 4
        pos = np.zeros([num_points, num_ncs, 3])

        pos[0, 0, :] = np.array([0.0, 0.0, 2.0])
        pos[0, 1, :] = np.array([2.0, 0.0, 0.0])
        pos[0, 2, :] = np.array([0.0, -2.0, 0.0])
        pos[0, 3, :] = np.array([0.0, -2.0, 2.0])

        pos[1, 0, :] = np.array([2.0, 2.0, 2.0])
        pos[1, 1, :] = np.array([2.0, -3.0, 3.0])
        pos[1, 2, :] = np.array([0.5, 0.0, -0.5])
        pos[1, 3, :] = np.array([0.5, 0.5, -0.5])

        pos[2, 0, :] = np.array([0.0, -0.2, 2.0])
        pos[2, 1, :] = np.array([0.1, 0.1, 0.3])
        pos[2, 2, :] = np.array([0.2, 0.0, -0.15])
        pos[2, 3, :] = np.array([-0.2, 0.0, -0.15])

        pos[3, 0, :] = np.array([0.0, 0.0, 0.0])
        pos[3, 1, :] = np.array([0.0, 0.0, 4.0])
        pos[3, 2, :] = np.array([4.0, 0.0, 0.0])
        pos[3, 3, :] = np.array([0.0, 4.0, 0.0])

        eps = 1e-7
        der_comp = np.zeros([pos.shape[0], pos.shape[1], pos.shape[2]])
        id_p = np.identity(3)[np.newaxis, :, :]
        d_pos = np.zeros_like(pos)
        for ind1 in range(num_ncs):
            for ind2 in range(3):
                d_pos[:, :, :] = pos[:, :, :]
                d_pos[:, ind1, :] += eps*id_p[..., :, ind2]
                der_comp[:, ind1, ind2] = -(self.ml_rot.psi_l3(d_pos)[:]-self.ml_rot.psi_l3(pos)[:])/eps

        der_result = np.zeros_like(der_comp)
        der_result[:, :, :] = self.ml_rot.der_psi_l3(pos)[:, :, :]

        cond = np.allclose(der_result, der_comp, atol=1e-8)
        self.assertTrue(cond)


if __name__ == '__main__':
    unittest.main()
