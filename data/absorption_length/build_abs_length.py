#!/usr/bin/env python3
import sys
sys.dont_write_bytecode = True
import os
from scipy.interpolate import CubicSpline, PchipInterpolator, make_interp_spline
import matplotlib.pyplot as plt
import numpy as np
import math
from statistics import mean

# eV to nm
def ε_to_λ(ε):
    return 1239.841984 / ε

# nm to eV 
def λ_to_ε(λ):
    return 1239.841984 / λ

def Ti_from_Latt(L, L_att):
    return 1.0 - math.exp(- L_att / L)

def Latt_from_Ti(L, Ti):
    return -math.log(1.0-Ti) * L

""" Get approximate total decrease in transparency due to Fresnel reflections on two surfaces. """
def Tfactor_from_n(n):
    R_single_surface = ((n-1.0)/(1.0+n)) ** 2
    return (1 - R_single_surface) ** 2

""" Get approximate refraction index from transparency, neglecting absorption in material. """
def n_from_Tfactor(T):
    R_single_surface = (1 - math.sqrt(T))
    sqrtR = math.sqrt(R_single_surface)
    return (1 + sqrtR) / (1 - sqrtR)

class FSQ_RG715:
    """
    FSQ RG715 optical filter with 715 nm cutoff with exact parameters taken from
    https://www.us.schott.com/shop/advanced-optics/en/Matt-Filter-Plates/RG715/c/glass-RG715
    The parameters should be google-able even if link is dead.
    """

    def __init__(self, thickness_mm = 3, T_internal_fname = 'Filter_FSQ-RG715_T_internal.dat'):
        self.L = thickness_mm
        self.n_λs =   [486.1, 404.7, 435.8, 852.1, 546.1, 587.6]
        self.n_vals = [1.544, 1.554, 1.550, 1.528, 1.539, 1.537]
        self.n_λs, self.n_vals = zip(*sorted(zip(self.n_λs, self.n_vals)))
        self.n_vs_λ = make_interp_spline(self.n_λs, self.n_vals, k = 1)
        self.λs, self.Τis = np.loadtxt(T_internal_fname, comments="//", skiprows=1, unpack=True)
        self.Ti_vs_λ = make_interp_spline(self.λs, self.Τis, k = 1)

    def T_full(self, λ, thickness = None):
        n = self.n_vs_λ(λ)
        T_factor = Tfactor_from_n(n)
        return T_factor * self.T_in(λ, thickness=thickness)
    
    def T_in(self, λ, thickness = None):
        if not thickness:
            return self.Ti_vs_λ(λ)
        L_att = Latt_from_Ti(self.L, self.Ti_vs_λ(λ))
        return Ti_from_Latt(thickness, L_att)
    
    def get_n(self, λ):
        n_vs_λ = self.n_vs_λ
        if type(λ) in [list, np.ndarray, range]:
            return [n_vs_λ(λi) for λi in λ]
        return n_vs_λ(λ)
    
    def get_L_att(self, λ):
        L = self.L
        T_in = self.Ti_vs_λ
        if type(λ) in [list, np.ndarray, range]:
            return [Latt_from_Ti(L, T_in(λi)) for λi in λ]
        return Latt_from_Ti(L, T_in(λ))
    
    def get_T_full(self, λ, thickness = None):
        T_full = self.T_full
        if type(λ) in [list, np.ndarray, range]:
            return [T_full(λi, thickness=thickness) for λi in λ]
        return T_full(λ, thickness=thickness)

class FSQ_RG715_our:
    """
    FSQ RG715 optical filter with 715 nm cutoff with exact parameters calculated
    from tranparency measurements. The results of these measurements come on physical
    paper included with filters. The measurements were probably done at BINP.
    """

    def __init__(self, thickness_mm = 3, T_full_fname = 'Filter_FSQ-RG715_T_full.dat', λ_plateu_nm = 820.0):
        self.L = thickness_mm
        self.λs, self.Tfull = np.loadtxt(T_full_fname, comments="//", skiprows=1, unpack=True)
        T_plateu = mean([T for λ, T in zip(self.λs, self.Tfull) if λ > λ_plateu_nm])
        self.n = n_from_Tfactor(T_plateu)
        self.Tfull = [min(T, 0.9999 * T_plateu) for T in self.Tfull] # * 0.9999 is to avoid infinite attenuation length for λ = λ | T = T_plateu.
        self.Tfull_vs_λ = make_interp_spline(self.λs, self.Tfull, k = 1)

    def T_full(self, λ, thickness = None):
        if not thickness:
            return self.Tfull_vs_λ(λ)
        T_factor = Tfactor_from_n(self.n)
        Tinternal = self.Tfull_vs_λ(λ) / T_factor
        L_att = Latt_from_Ti(self.L, Tinternal)
        return T_factor * Ti_from_Latt(thickness, L_att)
    
    def T_in(self, λ, thickness = None):
        T_factor = Tfactor_from_n(self.n)
        Tinternal = self.Tfull_vs_λ(λ) / T_factor
        if not thickness:
            return Tinternal
        L_att = Latt_from_Ti(self.L, Tinternal)
        return Ti_from_Latt(thickness, L_att)
    
    def get_n(self, λ):
        n = self.n
        if type(λ) in [list, np.ndarray, range]:
            return [n for λi in λ]
        return n

    def get_L_att(self, λ):
        L = self.L
        T_in = self.T_in
        if type(λ) in [list, np.ndarray, range]:
            return [Latt_from_Ti(L, T_in(λi)) for λi in λ]
        return Latt_from_Ti(L, T_in(λ))
    
    def get_T_full(self, λ, thickness = None):
        T_full = self.T_full
        if type(λ) in [list, np.ndarray, range]:
            return [T_full(λi, thickness=thickness) for λi in λ]
        return T_full(λ, thickness=thickness)



"""
Export reflectivity index and attenuation length in Geant4 format.
If energies or λs for tabulation are not specified, input points are used.
"""
def save_to_files(filter, filename_n, filename_L_att, energies=None, λs=None):
    with open(filename_n, "w") as file:
        if hasattr(filter, 'n_λs'):
            out_λs = λs if λs else filter.n_λs
        else:
            out_λs = λs if λs else filter.λs
        out_εs = energies if energies else [λ_to_ε(λ) for λ in reversed(out_λs)]
        out_λs = [ε_to_λ(ε) for ε in out_εs]
        ns = filter.get_n(out_λs)
        file.write("//1\t2\t1\t1\t0\t0\t0\t0\n")
        file.write("//Energy [eV]\tRefraction index\n")
        for ε, n in zip(out_εs, ns):
            file.write(f"{ε}\t{n}\n")
    with open(filename_L_att, "w") as file:
        out_λs = λs if λs else filter.λs
        out_εs = energies if energies else [λ_to_ε(λ) for λ in reversed(out_λs)]
        out_λs = [ε_to_λ(ε) for ε in out_εs]
        L_atts = filter.get_L_att(out_λs)
        file.write("//1\t2\t1\t1\t0\t0\t0\t0\n")
        file.write("//Energy [eV]\tAttenuation length [mm]\n")
        for ε, L_att in zip(out_εs, L_atts):
            file.write(f"{ε}\t{L_att}\n")


if __name__ == "__main__":
    filter_715nm_cite = FSQ_RG715()
    filter_715nm_our = FSQ_RG715_our()

    fugure_sz = (14, 9)

    λs = range(124, 1000, 1)
    εs = [λ_to_ε(λ) for λ in λs]

    fig, L_att = plt.subplots(1, 1, squeeze=False, figsize=fugure_sz)
    L_att = L_att.flatten()[0]
    L_att.set_title(f"Attenuation length for optical filter FSQ RG715")
    L_att.set_xlabel(fr"Photon energy (eV)")
    L_att.set_ylabel(fr"Attenuation length (mm)")
    data_cite = filter_715nm_cite.get_L_att(λs)
    data_our = filter_715nm_our.get_L_att(λs)
    L_att.plot(λs, data_cite, '--', label=f"Data from www.us.schott.com manufacturer")
    L_att.plot(λs, data_our, '-', label=f"Data from our specifications")
    L_att.legend(loc='best')

    fig, L_att = plt.subplots(1, 1, squeeze=False, figsize=fugure_sz)
    L_att = L_att.flatten()[0]
    L_att.set_title(f"Total transparency for optical filter FSQ RG715")
    L_att.set_xlabel(fr"Photon energy (eV)")
    L_att.set_ylabel(fr"T")
    data_cite = filter_715nm_cite.get_T_full(λs)
    data_our = filter_715nm_our.get_T_full(λs)
    L_att.plot(λs, data_cite, '--', label=f"Data from www.us.schott.com manufacturer")
    L_att.plot(λs, data_our, '-', label=f"Data from our specifications")
    L_att.legend(loc='best')

    save_to_files(filter_715nm_cite, "../refractive_index/fsq_rg715_rindex_eV_manuf.dat", "fsq_rg715_absorption_length_eV_mm_manuf.dat", λs=λs)
    save_to_files(filter_715nm_our, "../refractive_index/fsq_rg715_rindex_eV_our.dat", "fsq_rg715_absorption_length_eV_mm_our.dat", λs=λs)
    
    plt.tight_layout()
    plt.show()
