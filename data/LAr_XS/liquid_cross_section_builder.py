#!/usr/bin/env python3
import sys
sys.dont_write_bytecode = True
import gas_cross_section_builder as gas
from gas_cross_section_builder import eV_to_k, k_to_eV
from scipy.interpolate import CubicSpline, PchipInterpolator, make_interp_spline
import matplotlib.pyplot as plt
import numpy as np
import math

def read_data_vector_file(filename):
    line_n = 0
    xs = []
    ys = []
    with open(filename) as file:
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
                x = float(words[0])
                y = float(words[1])
                xs.append(x)
                ys.append(y)
            except Exception as e:
                print(f"Error in file \"{filename}\" on line {line_n}:")
                print(e)
                continue
    return np.array(xs), np.array(ys)
        

class CrossSectionBase:
    # Common interface for calculating cross sections for bith scalar and list values of energies
    def __init__(self, **kwargs):
        pass

    # Energy is in eV
    def get_energy_transfer_XS(self, energy):
        if type(energy) in [list, np.ndarray]:
            return [self._get_energy_transfer_XS(e) for e in energy]
        self._get_energy_transfer_XS(energy)

    # Energy is in eV
    def get_momentum_transfer_XS(self, energy):
        if type(energy) in [list, np.ndarray]:
            return [self._get_momentum_transfer_XS(e) for e in energy]
        self._get_momentum_transfer_XS(energy)

    def save_to_file(self, energies):
        filename_energy = self.default_ET_save_fname
        filename_momentum = self.default_MT_save_fname
        with open(filename_energy, "w") as file:
            file.write("//1\t2\t1\t1\t0\t0\t0\t0\n")
            file.write("//Energy [eV]\tET XS[cm^2]\n")
            for x, y in zip(energies, self.get_energy_transfer_XS(energies)):
                file.write(f"{x}\t{y}\n")
        with open(filename_momentum, "w") as file:
            file.write("//1\t2\t1\t1\t0\t0\t0\t0\n")
            file.write("//Energy [eV]\tMT XS[cm^2]\n")
            for x, y in zip(energies, self.get_momentum_transfer_XS(energies)):
                file.write(f"{x}\t{y}\n")

class CrossSectionFromFile(CrossSectionBase):
    def __init__(self, filename_energy, filename_momentum):
        super().__init__()
        self.XS_energy_trasnfer = read_data_vector_file(filename_energy)
        self.XS_momentum_trasnfer = read_data_vector_file(filename_momentum)

    # Energy is in eV
    def _get_energy_transfer_XS(self, energy):
        return np.interp(energy, *self.XS_energy_trasnfer)

    # Energy is in eV
    def _get_momentum_transfer_XS(self, energy):
        return np.interp(energy, *self.XS_momentum_trasnfer)

class CrossSectionAtrazhev(CrossSectionBase):
    # Building cross section in liquid Ar using Atrazhev
    # approach, but with modern data. See my notes.
    # S_at_0 is unitless, mobility_at_0 is in cm^2/V/s, XS_gas_momentum is in cm^2
    def __init__(self, mobility_at_0, S_at_0, XS_gas_momentum, **kwargs):
        super().__init__(**kwargs)
        self.En_min = 0.08
        self.XS_gas_momentum = XS_gas_momentum
        self.S_at_0 = S_at_0
        self.XS_at_0 = 1e8*1.227e-22/S_at_0/mobility_at_0
        structure_factor_fname = "S(eps)_Yarnell73.txt"
        S_vs_e = read_data_vector_file(structure_factor_fname)
        S_vs_e = filter(lambda x: x[1] > S_at_0, zip(*S_vs_e))
        es, vals = zip(*S_vs_e)
        self.S_energies = np.array(es)
        self.S_values = np.array(vals)

    # Energy is in eV
    def _get_energy_transfer_XS(self, energy):
        if (energy < self.En_min):
            return self.XS_at_0
        return max(self.XS_at_0, self.XS_gas_momentum.get_transport_XS(energy))

    # Energy is in eV
    def _get_momentum_transfer_XS(self, energy):
        return self._get_energy_transfer_XS(energy) * self._get_S(energy)

    def _get_S(self, energy):
        return max(self.S_at_0, np.interp(energy, self.S_energies, self.S_values, left=self.S_at_0, right=1.0))

class CrossSection_Boyle15(CrossSectionFromFile):
    def __init__(self, **kwargs):
        super().__init__("XS_energy_transfer_Boyle15.txt", "XS_momentum_transfer_Boyle15.txt", **kwargs)

class CrossSection_Atrazhev85(CrossSectionFromFile):
    def __init__(self, **kwargs):
        super().__init__("XS_energy_transfer_Atrazhev85.txt", "XS_momentum_transfer_Atrazhev85.txt", **kwargs)

class CrossSection_Atrazhev_AvgS_Avgmu_AvgXS(CrossSectionAtrazhev):
    # Refer to my notes.
    def __init__(self, **kwargs):
        super().__init__((518+567.8)/2, (0.0552+0.0537+0.0544)/3, gas.CrossSectionAvg(**kwargs), **kwargs)
        self.default_ET_save_fname = "XS_energy_transfer_Atrazhev_AvgS_Avgmu_AvgXS.txt"
        self.default_MT_save_fname = "XS_momentum_transfer_Atrazhev_AvgS_Avgmu_AvgXS.txt"

class CrossSection_Atrazhev_AvgS_Avgmu_B8_9(CrossSectionAtrazhev):
    # Refer to my notes.
    def __init__(self, **kwargs):
        super().__init__((518+567.8)/2, (0.0552+0.0537+0.0544)/3, gas.CrossSectionBiagi_v8_9(**kwargs), **kwargs)
        self.default_ET_save_fname = "XS_energy_transfer_Atrazhev_AvgS_Avgmu_B8.9.txt"
        self.default_MT_save_fname = "XS_momentum_transfer_Atrazhev_AvgS_Avgmu_B8.9.txt"

class CrossSection_Atrazhev_AvgS_Avgmu_B8_9_sh(CrossSectionAtrazhev):
    # Refer to my notes.
    def __init__(self, **kwargs):
        super().__init__((518+567.8)/2, (0.0552+0.0537+0.0544)/3, gas.CrossSectionBiagi_v8_9(energy_shift = 0.80/11.6, **kwargs), **kwargs)
        self.default_ET_save_fname = "XS_energy_transfer_Atrazhev_AvgS_Avgmu_B8.9_sh.txt"
        self.default_MT_save_fname = "XS_momentum_transfer_Atrazhev_AvgS_Avgmu_B8.9_sh.txt"

class CrossSection_Atrazhev_AvgS_Avgmu_McE14_sh(CrossSectionAtrazhev):
    # Refer to my notes.
    def __init__(self, **kwargs):
        super().__init__((518+567.8)/2, (0.0552+0.0537+0.0544)/3, gas.CrossSectionMcEachran14(energy_shift = 0.80/11.6, **kwargs), **kwargs)
        self.default_ET_save_fname = "XS_energy_transfer_Atrazhev_AvgS_Avgmu_McE14_sh.txt"
        self.default_MT_save_fname = "XS_momentum_transfer_Atrazhev_AvgS_Avgmu_McE14_sh.txt"

class CrossSection_Atrazhev_AvgS_Avgmu_McE14(CrossSectionAtrazhev):
    # Refer to my notes.
    def __init__(self, **kwargs):
        super().__init__((518+567.8)/2, (0.0552+0.0537+0.0544)/3, gas.CrossSectionMcEachran14(**kwargs), **kwargs)
        self.default_ET_save_fname = "XS_energy_transfer_Atrazhev_AvgS_Avgmu_McE14.txt"
        self.default_MT_save_fname = "XS_momentum_transfer_Atrazhev_AvgS_Avgmu_McE14.txt"

class CrossSection_Var10(CrossSectionAtrazhev):
    def __init__(self, **kwargs):
        super().__init__((518+567.8)/2, (0.0552+0.0537+0.0544)/3, gas.CrossSectionBiagi_v8_9(energy_shift = 0.64/11.6, XS_shift = 0.78, **kwargs), **kwargs)
        self.default_ET_save_fname = "XS_energy_transfer_Var10.txt"
        self.default_MT_save_fname = "XS_momentum_transfer_Var10.txt"


if __name__ == "__main__":
    XS_Atrazhev_AvgS_Avgmu_B8_9 = CrossSection_Atrazhev_AvgS_Avgmu_B8_9()
    XS_Atrazhev_AvgS_Avgmu_B8_9_sh = CrossSection_Atrazhev_AvgS_Avgmu_B8_9_sh()
    XS_Atrazhev_AvgS_Avgmu_McE14 = CrossSection_Atrazhev_AvgS_Avgmu_McE14()
    XS_Atrazhev_AvgS_Avgmu_McE14_sh = CrossSection_Atrazhev_AvgS_Avgmu_McE14_sh()
    XS_Atrazhev_AvgS_Avgmu_AvgXS = CrossSection_Atrazhev_AvgS_Avgmu_AvgXS()
    XS_Var10 = CrossSection_Var10()
    XS_Atrazhev85 = CrossSection_Atrazhev85()
    XS_Boyle15 = CrossSection_Boyle15()

    # Plot energy- and momentum-trasnfer e-LAr cross section for different sources
    fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(6, 6))
    ax = ax.flatten()
    energies = np.logspace(math.log10(1e-2), math.log10(30.0), num=351, base=10)
    ax[0].set_title("e-Ar energy-transfer effective cross sections in liquid")
    ax[0].set_xlabel(fr"$\varepsilon$ (eV)")
    ax[0].set_ylabel(fr"$q_{{et}} (\varepsilon)$ (cm$^2$)")
    ax[0].set_yscale('log')
    ax[0].set_xscale('log')
    ax[0].plot(energies, XS_Atrazhev_AvgS_Avgmu_AvgXS.get_energy_transfer_XS(energies), label='ET AvgS_Avgmu_AvgXS', linestyle='solid')
    ax[0].plot(energies, XS_Atrazhev_AvgS_Avgmu_B8_9.get_energy_transfer_XS(energies), label='ET AvgS_Avgmu_B8_9', linestyle='solid')
    ax[0].plot(energies, XS_Atrazhev_AvgS_Avgmu_B8_9_sh.get_energy_transfer_XS(energies), label='ET AvgS_Avgmu_B8_9_sh', linestyle='solid')
    ax[0].plot(energies, XS_Atrazhev_AvgS_Avgmu_McE14.get_energy_transfer_XS(energies), label='ET AvgS_Avgmu_McE14', linestyle='solid')
    ax[0].plot(energies, XS_Atrazhev_AvgS_Avgmu_McE14_sh.get_energy_transfer_XS(energies), label='ET AvgS_Avgmu_McE14_sh', linestyle='solid')
    ax[0].plot(energies, XS_Var10.get_energy_transfer_XS(energies), label='ET Var10', linestyle='solid')
    ax[0].plot(energies, XS_Atrazhev85.get_energy_transfer_XS(energies), label='ET Atrazhev85', linestyle='solid')
    ax[0].plot(energies, XS_Boyle15.get_energy_transfer_XS(energies), label='ET Boyle15', linestyle='solid')
    ax[0].legend(loc='best')

    fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(6, 6))
    ax = ax.flatten()
    ax[0].set_title("e-Ar momentum-transfer effective cross sections in liquid")
    ax[0].set_xlabel(fr"$\varepsilon$ (eV)")
    ax[0].set_ylabel(fr"$q_{{mt}} (\varepsilon)$ (cm$^2$)")
    ax[0].set_yscale('log')
    ax[0].set_xscale('log')
    ax[0].plot(energies, XS_Atrazhev_AvgS_Avgmu_AvgXS.get_momentum_transfer_XS(energies), '--', label='MT AvgS_Avgmu_AvgXS')
    ax[0].plot(energies, XS_Atrazhev_AvgS_Avgmu_B8_9.get_momentum_transfer_XS(energies), '--', label='MT AvgS_Avgmu_B8_9')
    ax[0].plot(energies, XS_Atrazhev_AvgS_Avgmu_B8_9_sh.get_momentum_transfer_XS(energies), '--', label='MT AvgS_Avgmu_B8_9_sh')
    ax[0].plot(energies, XS_Atrazhev_AvgS_Avgmu_McE14.get_momentum_transfer_XS(energies), '--', label='MT AvgS_Avgmu_McE14')
    ax[0].plot(energies, XS_Atrazhev_AvgS_Avgmu_McE14_sh.get_momentum_transfer_XS(energies), '--', label='MT AvgS_Avgmu_McE14_sh')
    ax[0].plot(energies, XS_Var10.get_momentum_transfer_XS(energies), '--', label='MT Var10')
    ax[0].plot(energies, XS_Atrazhev85.get_momentum_transfer_XS(energies), '--', label='MT Atrazhev85')
    ax[0].plot(energies, XS_Boyle15.get_momentum_transfer_XS(energies), '--', label='MT Boyle15')
    ax[0].legend(loc='best')
    
    XS_Atrazhev_AvgS_Avgmu_B8_9.save_to_file(energies)
    XS_Atrazhev_AvgS_Avgmu_McE14.save_to_file(energies)
    XS_Atrazhev_AvgS_Avgmu_AvgXS.save_to_file(energies)
    XS_Atrazhev_AvgS_Avgmu_B8_9_sh.save_to_file(energies)
    XS_Atrazhev_AvgS_Avgmu_McE14_sh.save_to_file(energies)
    XS_Var10.save_to_file(energies)
    plt.tight_layout()
    plt.show()
