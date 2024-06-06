"""
:module: MachineLearnFunct
:platform: Unix, Windows
:synopsis: Utility to calculate functions used for force matching in machine learning

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, February2024
"""

import numpy as np
import numpy.linalg as la


class FRotG(object):
    """
    Defines the functions necessary for machine learning calculations
    """

    def __init__(self, param2, param3, cut_off):
        """
        Constructor of the polyhedra

        :param param2: numpy array defining the parameters for `\\G^{(2)}` :math:`\\eta``, :math:`\\R_s``
        :param param3: numpy array defining the parameters for `\\G^{(3)}` :math:`\\lambda``, :math:`\\xi``, `\\mu`
        :param cut_off: value of the cut-off
        """

        # define parameters

        self.eta, self.r_s = param2
        self.lamda, self.xi, self.mu = param3
        self.r_cut = cut_off
        self.params = {'eta': self.eta, 'r_s': self.r_s, 'lambda': self.lamda,
                       'xi': self.xi, 'mu': self.mu, 'r_cut': self.r_cut}

    def fc(self, r):
        """
        Defines the cut-off function :math:`f_c(R)`

        :param r: distance

        return: value of the function
        """

        argval = 1-r/self.r_cut
        val = np.tanh(argval*np.heaviside(argval, 0))

        return val**3

    def hc(self, r):
        """
        Defines minus the derivative of the cut-off function :math:`h_c(R)=-\frac{d f_c}{dR}`

        :param r: distance

        :return: value of the function
        """

        argval = 1 - r / self.r_cut
        val = np.tanh(argval * np.heaviside(argval, 0))/np.cosh(argval)

        return 3*val**2/self.r_cut

    @staticmethod
    def cos_theta(r_a, r_b):
        """
        Defines the cosinus of the two vectors

        :param r_a: vector a
        :param r_b: vector b

        :return: value of the function
        """

        return np.sum(r_a*r_b, axis=-1)/(la.norm(r_a, axis=-1)*la.norm(r_b, axis=-1))

    def omega(self, r_a, r_b):
        """
        Defines the function :math:`\\Omega({\vec R}_a, {\vec R}_b)`

        :param r_a: vector a
        :param r_b: vector b

        :return: value of the function
        """

        l_a = la.norm(r_a, axis=-1)
        l_b = la.norm(r_b, axis=-1)
        l_c = la.norm(r_a-r_b, axis=-1)

        f_exp = np.exp(-self.mu*(l_a**2+l_b**2+l_c**2))
        f_cut = self.fc(l_a)*self.fc(l_b)*self.fc(l_c)

        return (1+self.lamda*self.cos_theta(r_a, r_b))**self.xi*f_exp*f_cut

    def phi(self, r_a, r_b):
        """
        Defines the function :math:`\\phi({\vec R}_a, {\vec R}_b)`

        :param r_a: vector a
        :param r_b: vector b

        :return: value of the function
        """

        l_a = la.norm(r_a, axis=-1)
        l_b = la.norm(r_b, axis=-1)
        cs = self.cos_theta(r_a, r_b)
        inv_a = 1/l_a

        return inv_a*(1/l_b-cs/l_a)*self.lamda*self.xi/(1+self.lamda*cs)

    def psi(self, r_a, r_b):
        """
        Defines the function :math:`\\psi({\\vec R}_a, {\\vec R}_b)`

        :param r_a: vector a
        :param r_b: vector b

        :return: value of the function
        """

        inv_a = 1/la.norm(r_a, axis=-1)
        inv_b = 1/la.norm(r_b, axis=-1)
        inv_cs = 1/(1+self.lamda*self.cos_theta(r_a, r_b))

        return -self.lamda*self.xi*inv_a*inv_b*inv_cs

    def g2(self, r_v):
        """
        Defines the function :math:`G^{(2)}({\vec R})`

        :param r_v: vector :math:`{\vec R}={\vec R}_{i}-{\vec R}_j`

        :return: value of the function

        """

        l_v = la.norm(r_v, axis=-1)

        return np.exp(-self.eta*(l_v-self.r_s)**2)*self.fc(l_v)

    def g3(self, r_v, r_w):
        """
        Defines the function :math:`G^{(2)}({\vec R}_1, {\vec R}_2)`

        Here,

        :param r_v: vector  :math:`{\vec R}_1={\vec R}_{i}-{\vec R}_j`
        :param r_w: vector :math:`{\vec R}_2={\vec R}_{i}-{\vec R}_k`

        :return: value of the function

        """

        l_v = la.norm(r_v, axis=-1)
        l_w = la.norm(r_w, axis=-1)
        l_u = la.norm(r_w-r_v, axis=-1)

        fac1 = 2**(1-self.xi)*(1+self.lamda*self.cos_theta(r_v, r_w))**self.xi
        fac2 = np.exp(-self.mu*(l_v**2+l_w**2+l_u**2))
        fac3 = self.fc(l_v)*self.fc(l_w)*self.fc(l_u)

        return fac1*fac2*fac3

    def der_g2_1(self, r_v):
        """
        Defines the function :math:`\\frac{\\partial}{\\partial {\vec R}_i}G^{(2)}({\\vec R})`

        :param r_v: vector :math:`{\\vec R}={\\vec R}_{i}-{\vec R}_j`

        :return: value of the function
        """

        l_v = la.norm(r_v, axis=-1)

        fac1 = np.exp(-self.eta*(self.r_s-l_v)**2)
        fac2 = self.hc(l_v)+2*self.eta*(l_v-self.r_s)*self.fc(l_v)

        fac = -fac1*fac2/l_v
        return r_v*fac[..., np.newaxis]

    def der_g2_2(self, r_v):
        """
        Defines the function :math:`\\frac{\\partial}{\\partial {\\vec R}_j}G^{(2)}({\\vec R})`

        :param r_v: vector :math:`{\\vec R}={\\vec R}_{i}-{\\vec R}_j`

        :return: value of the function
        """

        return -self.der_g2_1(r_v)

    def der_g3_1(self, r_v, r_w):
        """
        Defines the function :math:`\\frac{\\partial}{\\partial {\vec R}_i} G^{(3)}({\\vec R}_1, {\\vec R}_2)`

        Here,

        :param r_v: vector  :math:`{\\vec R}_1={\\vec R}_{i}-{\vec R}_j`
        :param r_w: vector  :math:`{\\vec R}_2={\\vec R}_{i}-{\vec R}_k`

        :return: value of the function

        """
        l_v = la.norm(r_v, axis=-1)
        l_w = la.norm(r_w, axis=-1)

        fac1 = 2**(1-self.xi)*self.omega(r_v, r_w)
        fac2_1 = fac1*(self.phi(r_v, r_w)-self.hc(l_v)/(l_v*self.fc(l_v))-2*self.mu)
        fac2_2 = fac1*(self.phi(r_w, r_v)-self.hc(l_w)/(l_w*self.fc(l_w))-2*self.mu)

        return fac2_1[..., np.newaxis]*r_v+fac2_2[..., np.newaxis]*r_w

    def der_g3_2(self, r_v, r_w):
        """
        Defines the function :math:`\\frac{\\partial}{\\partial {\\vec R}_{j}} G^{(3)}({\\vec R}_1, {\\vec R}_2)`

        :param r_v: vector  :math:`{\\vec R}_1={\\vec R}_{i}-{\\vec R}_j`
        :param r_w: vector  :math:`{\\vec R}_2={\\vec R}_{i}-{\\vec R}_k`

        :return: value of the function

        """

        l_v = la.norm(r_v, axis=-1)
        r_u = r_w-r_v
        l_u = la.norm(r_u, axis=-1)

        fac1 = 2**(1 - self.xi)*self.omega(r_v, r_w)
        fac2_1 = fac1*(self.phi(r_v, r_w)-self.hc(l_v)/(l_v * self.fc(l_v))-2*self.mu)
        fac2_2 = fac1*(self.psi(r_v, r_w)-self.hc(l_u)/(l_u * self.fc(l_u))-2*self.mu)

        return -fac2_1[..., np.newaxis]*r_v + fac2_2[..., np.newaxis]*r_u

    def der_g3_3(self, r_v, r_w):
        """
        Defines the function :math:`\\frac{\\partial}{\\partial {\\vec R}_{k}} G^{(2)}({\\vec R}_1, {\\vec R}_2)`

        :param r_v: vector  :math:`{\vec R}_1={\vec R}_{i}-{\vec R}_j`
        :param r_w: vector  :math:`{\vec R}_2={\vec R}_{i}-{\vec R}_k`

        :return: value of the function
        """

        l_w = la.norm(r_w, axis=-1)
        r_u = r_w-r_v
        l_u = la.norm(r_u, axis=-1)

        fac1 = 2**(1-self.xi)*self.omega(r_v, r_w)
        fac2_1 = fac1*(self.phi(r_w, r_v)-self.hc(l_w)/(l_w*self.fc(l_w))-2*self.mu)
        fac2_2 = fac1*(self.psi(r_v, r_w)-self.hc(l_u)/(l_u*self.fc(l_u))-2*self.mu)

        return -fac2_1[..., np.newaxis]*r_w - fac2_2[..., np.newaxis]*r_u

    def l2(self, index, pos_nc):
        """
        Defines the L2 function as defined in the notes

        index: particle index
        pos_nc: particle positions, numpy array with dimensions [number of configurations, number of particles, 3]
        """

        pos_rel = np.zeros([pos_nc.shape[0], pos_nc.shape[1]-1, pos_nc.shape[2]])
        pos_rel[:, :index, :] = pos_nc[:, index:(index+1), :] - pos_nc[:, :index, :]
        pos_rel[:, index:, :] = pos_nc[:, index:(index+1), :] - pos_nc[:, (index + 1):, :]

        return np.sum(self.g2(pos_rel), axis=1)

    def l3(self, index, pos_nc):
        """
        Defines the L3 function as defined in the notes

        index: particle index
        pos_nc: particle positions, numpy array with dimensions [number of configurations, number of particles, 3]
        """
        pos_rel = np.zeros([pos_nc.shape[0], pos_nc.shape[1]-1, pos_nc.shape[2]])
        pos_rel1 = np.zeros([pos_nc.shape[0], (pos_nc.shape[1]-1)*(pos_nc.shape[1]-2), pos_nc.shape[2]])
        pos_rel2 = np.zeros([pos_nc.shape[0], (pos_nc.shape[1]-1)*(pos_nc.shape[1]-2), pos_nc.shape[2]])
        pos_rel[:, :index, :] = pos_nc[:, index:(index+1), :] - pos_nc[:, :index, :]
        pos_rel[:, index:, :] = pos_nc[:, index:(index+1), :] - pos_nc[:, (index + 1):, :]

        ind_mat = pos_nc.shape[1]-1
        ind_mst = ind_mat-1

        for ind1 in range(ind_mat):
            pos_rel1[:, (ind1*ind_mst):(ind1+1)*ind_mst, :] = pos_rel[:, ind1:(ind1+1), :]
            pos_rel2[:, ind1*ind_mst:(ind1+ind1*ind_mst), :] = pos_rel[:, :ind1, :]
            pos_rel2[:, (ind1+ind1*ind_mst):(ind1+1)*ind_mst] = pos_rel[:, (ind1+1):, :]

        return np.sum(self.g3(pos_rel1, pos_rel2), axis=1)

    def der_l2(self, index_i, index_j, pos_nc):
        """
        Defines the derivative of the L2 function as defined in the notes

        index_i: derivative index
        index_j: particle index
        pos_nc: particle positions, numpy array with dimensions [number of configurations, number of particles, 3]
        """

        pos_rel = np.zeros([pos_nc.shape[0], pos_nc.shape[1]-1, pos_nc.shape[2]])
        pos_rel[:, :index_j, :] = pos_nc[:, index_j:(index_j+1), :] - pos_nc[:, :index_j, :]
        pos_rel[:, index_j:, :] = pos_nc[:, index_j:(index_j+1), :] - pos_nc[:, (index_j+1):, :]

        if index_i == index_j:
            return np.sum(self.der_g2_1(pos_rel), axis=1)

        return self.der_g2_2(pos_nc[:, index_j, :]-pos_nc[:, index_i, :])

    def der_l3(self, index_i, index_j, pos_nc):
        """
        Defines the derivative of the L3 function as defined in the notes

        index_i: derivative index
        index_j: particle index
        pos_nc: particle positions, numpy array with dimensions [number of configurations, number of particles, 3]
        """

        pos_rel = np.zeros([pos_nc.shape[0], pos_nc.shape[1] - 1, pos_nc.shape[2]])
        pos_rel1 = np.zeros([pos_nc.shape[0], (pos_nc.shape[1]-1)*(pos_nc.shape[1]-2), pos_nc.shape[2]])
        pos_rel2 = np.zeros([pos_nc.shape[0], (pos_nc.shape[1]-1)*(pos_nc.shape[1]-2), pos_nc.shape[2]])
        pos_rel[:, :index_j, :] = pos_nc[:, index_j:(index_j+1), :] - pos_nc[:, :index_j, :]
        pos_rel[:, index_j:, :] = pos_nc[:, index_j:(index_j+1), :] - pos_nc[:, (index_j + 1):, :]

        ind_mat = pos_nc.shape[1] - 1
        ind_mst = ind_mat-1

        for ind1 in range(ind_mat):
            pos_rel1[:, (ind1 * ind_mst):(ind1 + 1) * ind_mst, :] = pos_rel[:, ind1:(ind1 + 1), :]
            pos_rel2[:, ind1 * ind_mst:(ind1 + ind1 * ind_mst), :] = pos_rel[:, :ind1, :]
            pos_rel2[:, (ind1 + ind1 * ind_mst):(ind1 + 1) * ind_mst] = pos_rel[:, (ind1 + 1):, :]

        if index_i == index_j:
            return np.sum(self.der_g3_1(pos_rel1, pos_rel2), axis=1)

        ind_min = np.min([index_i, index_j])
        ind_max = np.max([index_i, index_j])
        indh = list(range(0, ind_min))+list(range(ind_min+1, ind_max))+list(range(ind_max+1, ind_mat+1))
        ind_a = np.array(indh, dtype=int)
        ind_b = index_i*np.ones_like(ind_a, dtype=int)

        s1 = np.sum(self.der_g3_2(pos_nc[:, index_j:(index_j+1), :]-pos_nc[:, ind_b, :],
                                  pos_nc[:, index_j:(index_j+1), :]-pos_nc[:, ind_a, :]), axis=1)
        s2 = np.sum(self.der_g3_3(pos_nc[:, index_j:(index_j+1), :]-pos_nc[:, ind_a, :],
                                  pos_nc[:, index_j:(index_j+1), :]-pos_nc[:, ind_b, :]), axis=1)

        return s1+s2

    def psi_l2(self, pos_nc):
        """
        Defines the entire potential as the sum of l2 functions

        pos_nc: particle positions, numpy array with dimensions [number of configurations, number of particles, 3]
        """

        res = np.zeros([pos_nc.shape[0]])
        for ind in range(pos_nc.shape[1]):
            res[:] += self.l2(ind, pos_nc)[:]

        return res

    def der_psi_l2(self, pos_nc):
        """
        Defines the derivative of the psi l2 function as defined in the notes

        index: particle index
        pos_nc: particle positions, numpy array with dimensions [number of configurations, number of particles, 3]
        """

        res = np.zeros_like(pos_nc)
        for ind in range(pos_nc.shape[1]):
            for index in range(pos_nc.shape[1]):
                res[:, ind, :] += self.der_l2(ind, index, pos_nc)

        return -res

    def psi_l3(self, pos_nc):
        """
        Defines the entire potential as the sum of l3 functions

        pos_nc: particle positions, numpy array with dimensions [number of configurations, number of particles, 3]
        """

        res = np.zeros([pos_nc.shape[0]])
        for ind in range(pos_nc.shape[1]):
            res[:] += self.l3(ind, pos_nc)[:]

        return res

    def der_psi_l3(self, pos_nc):
        """
        Defines the derivative of the psi l3 function as defined in the notes

        pos_nc: particle positions, numpy array with dimensions [number of configurations, number of particles, 3]
        """

        res = np.zeros_like(pos_nc)
        for ind in range(pos_nc.shape[1]):
            for index in range(pos_nc.shape[1]):
                res[:, ind, :] += self.der_l3(ind, index, pos_nc)

        return -res
