import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

MECH = "chem.yaml"
SPECIES_LIST = ["AL2O3(L)", "AL2O3"]
Tmin, Tmax = 200, 4000
nT = 200

Tlist = np.linspace(Tmin, Tmax, nT)

gas = ct.Solution(MECH)

for sp in SPECIES_LIST:
    h = []
    s = []
    cp = []
    for T in Tlist:
        gas.TPX = T, ct.one_atm, f"{sp}:1.0"
        h.append(gas.enthalpy_mass)
        s.append(gas.entropy_mass)
        cp.append(gas.cp_mass)
    plt.figure(figsize=(12,4))
    plt.subplot(1,3,1)
    plt.plot(Tlist, h, label=sp)
    plt.xlabel("T [K]"); plt.ylabel("Enthalpy [J/kg]"); plt.title(f"{sp} Enthalpy")
    plt.grid(True)
    plt.subplot(1,3,2)
    plt.plot(Tlist, s, label=sp)
    plt.xlabel("T [K]"); plt.ylabel("Entropy [J/kg/K]"); plt.title(f"{sp} Entropy")
    plt.grid(True)
    plt.subplot(1,3,3)
    plt.plot(Tlist, cp, label=sp)
    plt.xlabel("T [K]"); plt.ylabel("Cp [J/kg/K]"); plt.title(f"{sp} Cp")
    plt.grid(True)
    plt.tight_layout()
    #plt.suptitle(sp)
    plt.show()