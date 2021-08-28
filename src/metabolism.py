import yaml
import pathlib
import numpy as np
from scipy.optimize import bisect
from scipy.interpolate import interp1d, interp2d, griddata, Rbf

PARACRINE = 0.23
INTRINSIC = 1-PARACRINE


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

    # GLYCOLYSIS
    def J_G6P(self, G):
        J_max, K_m = self.par["J_max"], self.par["K_m"]
        return J_max * G**2 / (K_m**2 + G**2)

    def J_ATP_gly(self, G):
        return 2*self.J_G6P(G)

    def J_NADH_gly(self, G):
        return 2*self.J_G6P(G)

    def J_pyr(self, G):
        return 2*self.J_G6P(G)

    def J_ATP_NADH_gly(self, G):
        p_L = self.par["p_L"]
        return 1.5*(self.J_NADH_gly(G)-p_L*self.J_pyr(G))

    # TCA CYCLE
    def J_NADH_pyr(self, G):
        p_L, p_TCA = self.par["p_L"], self.par["p_TCA"]
        return 5*p_TCA*(1-p_L)*self.J_pyr(G)

    def J_ATP_NADH_pyr(self, G):
        return 2.5*self.J_NADH_pyr(G)

    # OXYGEN INPUT
    def J_O2_G(self, G):
        p_L = self.par["p_L"]
        return 0.5*(self.J_NADH_pyr(G)+self.J_NADH_gly(G)-p_L*self.J_pyr(G))

    def J_O2(self, G):
        J_O2_0, k_O2 = self.par["J_O2_0"], self.par["k_O2"]
        alpha = k_O2*self.J_G6P(G) + J_O2_0

        J_O2_1, K_m_O2 = self.par["J_O2_1"], self.par["K_m_O2"]
        n_O2 = self.par["n_O2"]
        beta = J_O2_1*self.J_G6P(G)**n_O2/(K_m_O2**n_O2+self.J_G6P(G)**n_O2)
        return alpha + beta

    # BETA-OXIDATION
    def J_NADH_FFA(self, G):
        return 2*(self.J_O2(G)-self.J_O2_G(G))

    def J_ATP_NADH_FFA(self, G):
        return 2.3*self.J_NADH_FFA(G)  # !spremenjeno iz 2.5!

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


class Alpha(Cell):
    def __init__(self):
        super(Alpha, self).__init__("alpha")
        self.beta_mitos = 1
        camp_data = pathlib.Path(__file__).parent / "camp.txt"
        mesh_data = pathlib.Path(__file__).parent / "mesh.txt"
        self.cAMP_data = np.loadtxt(camp_data).T
        norm = self.cAMP_data[1]
        norm = (norm-np.min(norm))/(np.max(norm)-np.min(norm))
        self.cAMP_data[1] = norm
        self.mesh_data = np.loadtxt(mesh_data).T
        self.mesh_data[2] /= np.max(self.mesh_data[2])

    # --------------------------------- cAMP -------------------------------- #
    def cAMP_interpolation(self, gKATP):
        return interp1d(*self.cAMP_data)(gKATP)

    # ------------------------------- secretion ----------------------------- #
    def mesh_interpolation(self, gKATP, fcAMP):
        assert np.array(gKATP >= 0).all()  and np.array(gKATP <= 0.4).all()
        assert np.array(fcAMP >= 0).all() and np.array(fcAMP <= 1).all()
        x, y, z = self.mesh_data
        sparse_points = np.stack([x, y], -1)
        result = griddata(sparse_points, z, (gKATP, fcAMP))
        return result

    def RGS_Ca(self, gKATP):
        before = self.mitos
        self.mitos = 1
        G = self.glucose_vec(gKATP)
        self.mitos = before

        # if G == None:
        #     return 0 

        G = np.nan_to_num(G)

        Ca_max = 0.12
        G_05 = 17
        return Ca_max/(1 + np.exp(-0.3*(G-G_05)))

    def RGS_Ca_vec(self, gKATP):
        return np.vectorize(self.RGS_Ca)(gKATP)

    def RGS_intrinsic(self, gKATP, fcAMP):
        return self.mesh_interpolation(gKATP, fcAMP)+self.RGS_Ca_vec(gKATP)

    def RGS_paracrine(self, G):
        beta = Beta()
        beta.mitos = self.beta_mitos
        beta_gKATP = beta.g_K_ATP(G)
        RIS = beta.RS(beta_gKATP)
        return 1-RIS

    def RGS(self, G, gKATP, cAMP):
        return INTRINSIC*self.RGS_intrinsic(gKATP, cAMP)+PARACRINE*self.RGS_paracrine(G)
