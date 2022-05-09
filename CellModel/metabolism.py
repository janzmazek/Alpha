import yaml
import pathlib
import numpy as np
from scipy.optimize import bisect
from scipy.interpolate import interp1d, interp2d, griddata, Rbf


class Cell(object):
    def __init__(self, type):
        data = pathlib.Path(__file__).parent / "parameters.yaml"
        with open(data, "r") as stream:
            try:
                par = yaml.safe_load(stream)[type]
            except yaml.YAMLError as exc:
                print(exc)
        self.par = par
        self.mitos = 1
        self.anaerobic = 1

    # GLYCOLYSIS
    def J_GK(self, G):
        J_max_GK, K_m_GK = self.par["J_max"], self.par["K_m"]
        return J_max_GK * G**2 / (K_m_GK**2 + G**2)

    def J_gly(self, G):
        return 2*self.J_GK(G)

    def J_lac(self, G):
        p_L = self.par["p_L"]
        return p_L*self.J_gly(G)

    def J_NADH_gly(self, G):
        return self.J_gly(G) - self.J_lac(G)

    def J_ATP_NADH_gly(self, G):
        return 1.5*self.J_NADH_gly(G)

    def J_ATP_gly(self, G):
        return self.anaerobic*self.J_gly(G)

    def J_pyr(self, G):
        return self.J_gly(G) - self.J_lac(G)

    # TCA CYCLE

    def J_anaplerosis(self, G):
        p_TCA = self.par["p_TCA"]
        return (1-p_TCA)*self.J_pyr(G)

    def J_NADH_pyr(self, G):
        return 5*(self.J_pyr(G)-self.J_anaplerosis(G))

    def J_ATP_NADH_pyr(self, G):
        return 2.5*self.J_NADH_pyr(G)

    # OXYGEN INPUT
    def J_O2_G(self, G):
        return 0.5*(self.J_NADH_pyr(G)+self.J_NADH_gly(G))

    def J_O2(self, G):
        J_O2_0, k_O2 = self.par["J_O2_0"], self.par["k_O2"]
        alpha = k_O2*self.J_GK(G) + J_O2_0

        J_O2_1, K_m_O2 = self.par["J_O2_1"], self.par["K_m_O2"]
        n_O2 = self.par["n_O2"]
        beta = J_O2_1*self.J_GK(G)**n_O2/(K_m_O2**n_O2+self.J_GK(G)**n_O2)
        return alpha + beta

    # BETA-OXIDATION
    def J_NADH_FFA(self, G):
        return 2*(self.J_O2(G)-self.J_O2_G(G))

    def J_ATP_NADH_FFA(self, G):
        return 2.5*self.J_NADH_FFA(G)  # !spremenjeno iz 2.5! 2.3

    # OXIDATIVE ATP PRODUCTION
    def J_ATP_ox(self, G):
        sum = self.J_ATP_NADH_gly(G)
        sum += self.J_ATP_NADH_pyr(G)
        sum += self.J_ATP_NADH_FFA(G)
        return self.mitos*sum

    # ATP PRODUCTION/HYDROLYSIS
    def J_ATP(self, G):
        return self.J_ATP_gly(G) + self.J_ATP_ox(G)

    def J_ATPase(self, ATP):
        k_ATPase, K_m_ATPase = self.par["k_ATPase"], self.par["K_m_ATPase"]
        return k_ATPase*ATP/(K_m_ATPase + ATP)

    # AXP CONCENTRATIONS
    def ATP(self, G):
        k_ATPase, K_m_ATPase = self.par["k_ATPase"], self.par["K_m_ATPase"]
        return K_m_ATPase*self.J_ATP(G)/(k_ATPase-self.J_ATP(G))

    def ADP(self, G):
        A_tot = self.par["A_tot"]
        return A_tot - self.ATP(G)

    def RAT(self, G):
        return self.ATP(G)/self.ADP(G)

    # KATP CHANNEL CONDUCTANCE
    def g_K_ATP(self, G):
        ag_K_ATP = self.par["ag_K_ATP"]
        MgADP = 0.165*self.ADP(G)
        ADP3 = 0.135*self.ADP(G)
        ATP4 = 0.05*self.ATP(G)

        up = 0.08*(1+2*MgADP/17)+0.89*(MgADP/17)**2
        down = (1+MgADP/17)**2*(1+ADP3/26+ATP4/1)
        return ag_K_ATP*up/down

    def g_K_ATP_vec(self):
        return np.vectorize(self.g_K_ATP)

    # Get glucose from gKATP
    def glucose(self, gkatp):
        f = lambda g: self.g_K_ATP(g)-gkatp
        try:
            result = bisect(f, 0, 100)
        except ValueError as e:
            if str(e) != "f(a) and f(b) must have different signs":
                raise ValueError(e)
            if gkatp > self.g_K_ATP(0):
                result = np.nan
            else:
                result = np.inf
        return result

    def glucose_vec(self, gKATP):
        return np.vectorize(self.glucose)(gKATP)

    # -------------------------- HORMONE SECRETION -------------------------- #
    def f_RS(self, gKATP):
        g_s_2, n_s = self.par["g_s_2"], self.par["n_s"]
        return (g_s_2**n_s)/(g_s_2**n_s + gKATP**n_s)

    def RS(self, gKATP):
        RS_0, g_K_ATP_ms = self.par["RS_0"], self.par["g_K_ATP_ms"]
        return (1 - RS_0)*self.f_RS(gKATP)/self.f_RS(g_K_ATP_ms) + RS_0


class Beta(Cell):
    def __init__(self):
        super(Beta, self).__init__("beta")


def norm_data(data):
    return (data-np.min(data))/(np.max(data)-np.min(data))


class Alpha(Cell):
    def __init__(self):
        super(Alpha, self).__init__("alpha")
        self.beta_mitos = 1

        cAMP_sAC_path = pathlib.Path(__file__).parent / "cAMP_sAC.txt"
        self.cAMP_sAC_data = np.loadtxt(cAMP_sAC_path).T
        self.cAMP_sAC_data[1] = norm_data(self.cAMP_sAC_data[1])

    # --------------------------------- cAMP -------------------------------- #
    def cAMP_sAC_interpolation(self, g):
        try:
            return interp1d(*self.cAMP_sAC_data)(g)
        except ValueError:
            return 0
