#!/usr/bin/env python3
import sys
sys.dont_write_bytecode = True
import os
from scipy.interpolate import CubicSpline, PchipInterpolator, make_interp_spline
import matplotlib.pyplot as plt
import numpy as np
import math

# Convert energy in eVs to k in atomic units
def eV_to_k(e):
    return math.sqrt(2*e/27.211386) 

# Convert k in atomic units to energy in eVs
def k_to_eV(k):
    return k*k/2*27.211386

# Convert energy in eVs to k in atomic units
def eV_to_k_rel(e):
    # See Eq. (7) and below in McEachran17 (https://doi.org/10.1140/epjd/e2014-50166-7)
    e = e/27.211386
    # alpha ~ 1/137.03
    return math.sqrt(e*(e/pow(137.03, 2)+2)) 

# Convert k in atomic units to energy in eVs
def k_to_eV_rel(k):
    # See Eq. (7) and below in McEachran17 (https://doi.org/10.1140/epjd/e2014-50166-7)
    # alpha ~ 1/137.03
    return 27.211386 * (137.03*math.sqrt(k*k + pow(137.03, 2)) - pow(137.03, 2))

class CrossSectionBase:
    # Common interface for calculating cross sections for bith scalar and list values of energies
    def __init__(self, energy_shift=1.0, XS_shift=1.0):
        self.energy_shift = energy_shift
        self.XS_shift = XS_shift

    # Energy is in eV
    def get_elastic_XS(self, energy):
        energy_shift = self.energy_shift
        XS_shift = self.XS_shift
        if type(energy) in [list, np.ndarray]:
            return [XS_shift*self._get_elastic_XS(e/energy_shift) for e in energy]
        return XS_shift*self._get_elastic_XS(energy/energy_shift)

    # Energy is in eV
    def get_transport_XS(self, energy):
        energy_shift = self.energy_shift
        XS_shift = self.XS_shift
        if type(energy) in [list, np.ndarray]:
            return [XS_shift*self._get_transport_XS(e/energy_shift) for e in energy]
        return XS_shift*self._get_transport_XS(energy/energy_shift)

    def save_to_file(self, energies, filename_total, filename_transport):
        with open(filename_total, "w") as file:
            file.write("//1\t2\t1\t1\t0\t0\t0\t0\n")
            file.write("//Energy [eV]\tElastic XS[cm^2]\n")
            for x, y in zip(energies, self.get_elastic_XS(energies)):
                file.write(f"{x}\t{y}\n")
        with open(filename_transport, "w") as file:
            file.write("//1\t2\t1\t1\t0\t0\t0\t0\n")
            file.write("//Energy [eV]\tTrasnport XS[cm^2]\n")
            for x, y in zip(energies, self.get_transport_XS(energies)):
                file.write(f"{x}\t{y}\n")

class CrossSectionMcEachran14(CrossSectionBase):
    # Refer to https://doi.org/10.1140/epjd/e2014-50166-7
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.a_mt = [
            6.53,
            -146.99,
            198.67,
            -325.50,
            265.56,
            552.07,
            0.0,
        ]
        self.b_mt = [
            1.0,
            -10.68,
            -501.88,
            -187.37,
            76.11,
            -740.25,
            455.44,
        ]
        self.a_el = [
            6.53,
            -67.66,
            251.69,
            64.71,
            -64.89,
            -256.86,
            0.0,
        ]
        self.b_el = [
            1.0,
            -0.2964,
            4.424,
            -33.54,
            -57.92,
            -21.22,
            58.40,
        ]

    # Energy is in eV
    def _get_elastic_XS(self, energy):
        # See Eq. (7) and below in McEachran14 https://doi.org/10.1140/epjd/e2014-50166-7
        k = eV_to_k_rel(energy) # in atomic units
        k_terms = [1, k, pow(k, 2), pow(k, 2)*math.log(k), pow(k, 3), pow(k, 3)*math.log(k), pow(k, 4)]
        XS = sum([ai*ki for ai,ki in zip(self.a_el, k_terms)]) / sum([bi*ki for bi,ki in zip(self.b_el, k_terms)])
        XS *= 1e-20 * 1e4 # Eq. (7)'s LHS is in 1e-20 m^2. Converting it to cm^2
        return XS

    # Energy is in eV
    def _get_transport_XS(self, energy):
        # See Eq. (7) and below in McEachran14 https://doi.org/10.1140/epjd/e2014-50166-7
        k = eV_to_k_rel(energy) # in atomic units
        k_terms = [1, k, pow(k, 2), pow(k, 2)*math.log(k), pow(k, 3), pow(k, 3)*math.log(k), pow(k, 4)]
        XS = sum([ai*ki for ai,ki in zip(self.a_mt, k_terms)]) / sum([bi*ki for bi,ki in zip(self.b_mt, k_terms)])
        XS *= 1e-20 * 1e4 # Eq. (7)'s LHS is in 1e-20 m^2. Converting it to cm^2
        return XS


class CrossSectionMcEachran97(CrossSectionBase):
    def __init__(self, IntByEnergy=False, **kwargs):
        super().__init__(**kwargs)
        phase_shifts_data = "McEachran97_ArPhaseShifts.dat"
        self.l_max = 5
        self.int_by_energy = IntByEnergy
        self.phase_shifts_plus = {l: [] for l in range(0, self.l_max+1)}
        self.phase_shifts_minus = {l: [] for l in range(0, self.l_max+1)}
        line_n = 0
        reading_plus = True
        k = None
        ps_0 = None
        with open(phase_shifts_data) as file:
            while True:
                line_n += 1
                line = file.readline()
                if not line:
                    break
                line = line.lstrip()
                if line.startswith("//") or not line or line.isspace():
                    continue
                try:
                    words = line.split('\t')
                    if reading_plus:
                        k = float(words[0])
                        #eV = round(k*k/2*27.211386, 3)
                        ps_0 = float(words[1])
                        self.phase_shifts_plus[0].append((k_to_eV(k) if IntByEnergy else k, ps_0))
                        reading_plus = False
                        for l in range(1, self.l_max+1):
                            ps = None
                            try:
                                ps = float(words[l+1])
                            except Exception as e:
                                break
                            self.phase_shifts_plus[l].append((k_to_eV(k) if IntByEnergy else k, ps))
                    else:
                        self.phase_shifts_minus[0].append((k_to_eV(k) if IntByEnergy else k, ps_0))
                        reading_plus = True
                        for l in range(1, self.l_max+1):
                            ps = None
                            try:
                                ps = float(words[l-1])
                            except Exception as e:
                                break
                            self.phase_shifts_minus[l].append((k_to_eV(k) if IntByEnergy else k, ps))

                except Exception as e:
                    print(f"Error in file \"{phase_shifts_data}\" on line {line_n}:")
                    print(e)
                    continue
        #self.phase_shifts_plus_interp = {l: make_interp_spline(*zip(*k_ph_vals), k = 1) for l, k_ph_vals in self.phase_shifts_plus.items()}
        #self.phase_shifts_minus_interp = {l: make_interp_spline(*zip(*k_ph_vals), k = 1) for l, k_ph_vals in self.phase_shifts_minus.items()} 
        self.phase_shifts_plus_interp = {l: CubicSpline(*zip(*k_ph_vals)) for l, k_ph_vals in self.phase_shifts_plus.items()}
        self.phase_shifts_minus_interp = {l: CubicSpline(*zip(*k_ph_vals)) for l, k_ph_vals in self.phase_shifts_minus.items()}

    # Energy is in eV
    def _get_elastic_XS(self, energy):
        # See Eq. (11) in McEachran97 https://doi.org/10.1071/P96082
        k = eV_to_k(energy) # in atomic units
        phs_plus = {l: spline(energy if self.int_by_energy else k) for l, spline in self.phase_shifts_plus_interp.items()}
        phs_minus = {l: spline(energy if self.int_by_energy else k) for l, spline in self.phase_shifts_minus_interp.items()}
        XS = 0
        for l in range(0, self.l_max+1):
            XS += (l + 1) * pow(math.sin(phs_plus[l]), 2) + l * pow(math.sin(phs_minus[l]), 2)
        XS = XS * 4 * math.pi / k / k
        XS *= pow(5.291772e-11, 2) * 1e4 # a.u. to SI and then to cm^2
        return XS

    # Energy is in eV
    def _get_transport_XS(self, energy):
        # See Eq. (12) in McEachran97 https://doi.org/10.1071/P96082
        k = eV_to_k(energy) # in atomic units
        phs_plus = {l: spline(energy if self.int_by_energy else k) for l, spline in self.phase_shifts_plus_interp.items()}
        phs_minus = {l: spline(energy if self.int_by_energy else k) for l, spline in self.phase_shifts_minus_interp.items()}
        XS = 0
        for l in range(0, self.l_max):
            XS += (l + 1) * (l + 2) / (2*l + 3) * pow(math.sin(phs_plus[l] - phs_plus[l+1]), 2) \
                + l * (l + 1) / (2*l + 1) * pow(math.sin(phs_minus[l] - phs_minus[l+1]), 2) \
                + (l + 1) / (2*l + 1) / (2*l + 3) * pow(math.sin(phs_plus[l] - phs_minus[l+1]), 2)
        XS = XS * 4 * math.pi / k / k
        XS *= pow(5.291772e-11, 2) * 1e4 # a.u. to SI and then to cm^2
        return XS

class CrossSectionBiagi_v8_9(CrossSectionBase):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        data_fname = "e-Ar_Biagi-8.9_trasnport_XS.txt"
        self.XS_vs_E_data = []
        line_n = 0
        with open(data_fname) as file:
            while True:
                line_n += 1
                line = file.readline()
                if not line:
                    break
                if line_n < 76:
                    continue
                line = line.lstrip()
                if line.startswith("---") or line.startswith("xxx") or line.startswith("***") or not line or line.isspace():
                    continue
                try:
                    words = line.split('\t')
                    E = float(words[0])
                    XS = float(words[1]) * 1e4 # m^2 to cm^2
                    self.XS_vs_E_data.append((E, XS))
                except Exception as e:
                    print(f"Error in file \"{data_fname}\" on line {line_n}:")
                    print(e)
                    continue
        #self.XS_vs_E_interp = make_interp_spline(*zip(*self.XS_vs_E_data), k = 1)
        self.XS_vs_E_interp = CubicSpline(*zip(*self.XS_vs_E_data))

    # Energy is in eV
    def _get_elastic_XS(self, energy):
        return 0.0

    # Energy is in eV
    def _get_transport_XS(self, energy):
        return self.XS_vs_E_interp(energy)

class CrossSectionAvg(CrossSectionBase):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.XS_McEachran14 = CrossSectionMcEachran14(**kwargs)
        self.XS_Biagi89 = CrossSectionBiagi_v8_9(**kwargs)

    # Energy is in eV
    def _get_elastic_XS(self, energy):
        return 0.0

    # Energy is in eV
    def _get_transport_XS(self, energy):
        return 0.5 * (self.XS_McEachran14._get_transport_XS(energy) 
            + self.XS_Biagi89._get_transport_XS(energy))

if __name__ == "__main__":
    XS_McEachran97 = CrossSectionMcEachran97(IntByEnergy=True)
    l_max = XS_McEachran97.l_max

    # Plot McEachran97 phase shifts' interoplation
    fig, ax = plt.subplots(3, 2, squeeze=False, figsize=(10, 9))
    ax = ax.flatten()
    for l in range(0, l_max):
        ax[l].set_title(fr"$\delta_{{l={l}}}^+$")
        ax[l].set_xlabel("k (a.u.)")
        ax[l].plot(*zip(*XS_McEachran97.phase_shifts_plus[l]), 'o', label='McEachran97')
        ks = np.linspace(0, k_to_eV(1.0) if XS_McEachran97.int_by_energy else 1.0, num=301)
        ax[l].plot(ks, XS_McEachran97.phase_shifts_plus_interp[l](ks), '--', label='Cubic spline')
        ax[l].legend(loc='best')
    fig, ax = plt.subplots(3, 2, squeeze=False, figsize=(10, 9))
    ax = ax.flatten()
    for l in range(0, l_max):
        ax[l].set_title(fr"$\delta_{{l={l}}}^-$")
        ax[l].set_xlabel("k (a.u.)")
        ax[l].plot(*zip(*XS_McEachran97.phase_shifts_minus[l]), 'o', label='McEachran97')
        ks = np.linspace(0, k_to_eV(1.0) if XS_McEachran97.int_by_energy else 1.0, num=301)
        ax[l].plot(ks, XS_McEachran97.phase_shifts_minus_interp[l](ks), '--', label='Cubic spline')
        ax[l].legend(loc='best')

    # Plot McEachran97 elastic and momentum-trasnfer e-Ar cross section
    fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(6, 6))
    ax = ax.flatten()
    energies = np.logspace(math.log10(1e-2), math.log10(50.0), num=301, base=10)
    XS_total = XS_McEachran97.get_elastic_XS(energies)
    XS_transport = XS_McEachran97.get_transport_XS(energies)
    ax[0].set_title("McEachran97 e-Ar cross section in gas")
    ax[0].set_xlabel(fr"$\varepsilon$ (eV)")
    ax[0].set_ylabel(fr"$\sigma (\varepsilon)$ (cm$^2$)")
    ax[0].plot(energies, XS_total, '--', label="Elastic (total) cross section")
    ax[0].plot(energies, XS_transport, '--', label="Transport (momentum-transfer) cross section")
    ax[0].legend(loc='best')

    # Plot momentum-trasnfer e-Ar cross section for different sources
    fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(6, 6))
    ax = ax.flatten()
    energies = np.logspace(math.log10(1e-2), math.log10(50.0), num=321, base=10)
    XS_McEachran14 = CrossSectionMcEachran14()
    XS_Biagi89 = CrossSectionBiagi_v8_9()
    ax[0].set_title("e-Ar transport (momentum-transfer) cross section in gas")
    ax[0].set_xlabel(fr"$\varepsilon$ (eV)")
    ax[0].set_ylabel(fr"$\sigma_{{mt}} (\varepsilon)$ (cm$^2$)")
    ax[0].plot(energies, XS_McEachran97.get_transport_XS(energies), '--', label="McEachran97")
    ax[0].plot(energies, XS_McEachran14.get_transport_XS(energies), '--', label="McEachran14")
    ax[0].plot(energies, XS_Biagi89.get_transport_XS(energies), '--', label="BIAGI-v8.9")
    ax[0].legend(loc='best')

    XS_McEachran14.save_to_file(energies, "gas_elastic_XS_McEachran14.txt", "gas_transport_XS_McEachran14.txt")
    XS_Biagi89.save_to_file(energies, "gas_elastic_XS_BIAGI-v8.9.txt", "gas_transport_XS_BIAGI-v8.9.txt")

    XS_Biagi89_shifted = CrossSectionBiagi_v8_9(energy_shift = 0.73/11.6)
    XS_Biagi89_shifted.save_to_file(energies, "gas_elastic_XS_BIAGI-v8.9_shifted.txt", "gas_transport_XS_BIAGI-v8.9_shifted.txt")
    
    plt.tight_layout()
    plt.show()
