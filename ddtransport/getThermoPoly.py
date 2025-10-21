import cantera as ct 
import numpy as np



# ===== Setting ====== #
mech = "chem.yaml"
Tlow = 200
Thigh = 4000
grids = 200


#---------------------#

def transportTemplate(species_name,muPoly,kappaPoly,oneOverKappaPoly):
    muString    = "muCoeffs<8>       ("
    for i in range(len(muPoly) + 1):
        muString += f"{muPoly[i]:.5e}" + " "
    muString += ");\n"

    kappaString = "kappaCoeffs<8>    ("
    for i in range(len(kappaPoly) + 1)  :
        kappaString += f"{kappaPoly[i]:.5e}" + " "
    kappaString += ");\n"
    
    oneOverKappaString = "oneOverKappaCoeffs<8>    ("
    for i in range(len(kappaPoly) + 1)  :
        oneOverKappaString += f"{oneOverKappaPoly[i]:.5e}" + " "
    oneOverKappaString += ");\n"
    
    string = ""
    string += species_name + "\n"
    string += "{\n"
    string += "    transport\n"
    string += "    {\n"
    string += "        " + muString
    string += "        " + kappaString
    string += "        " + oneOverKappaString
    string += "    }\n"
    string += "}\n"
    return string

with open("thermoPolynomial","w+") as f:
    gas = ct.Solution(mech)
    for species in gas.species_names:
        reactants_species = "{}:1".format(species)
        TList = np.linspace(Tlow, Thigh, grids)
        muList = []
        kappaList = []
        oneOverKappaList = []
        for T in TList:
            gas.TPY = T, 101325,reactants_species
            gas.transport_model = "mixture-averaged"
            mu = gas.viscosity
            kappa = gas.thermal_conductivity
            muList.append(mu)
            kappaList.append(kappa)
            oneOverKappaList.append(1./kappa)
            
        pmu = np.polyfit(TList, muList, 7)
        pmu_poly = np.poly1d(pmu)
        pkappa = np.polyfit(TList, kappaList, 7)
        pkappa_poly = np.poly1d(pkappa)
        poneOverKappa = np.polyfit(TList, oneOverKappaList, 7)
        poneOverKappa_poly = np.poly1d(poneOverKappa)
        
        mu_rmse = np.sqrt(np.mean((pmu_poly(TList) - muList)**2))
        kappa_rmse = np.sqrt(np.mean((pkappa_poly(TList) - kappaList)**2))
        oneOverKappa_rmse = np.sqrt(np.mean((poneOverKappa_poly(TList) - oneOverKappaList)**2))
        
        mu_max_rel = np.max(np.abs((pmu_poly(TList) - muList)/muList))
        kappa_max_rel = np.max(np.abs((pkappa_poly(TList) - kappaList)/kappaList))
        oneOverKappa_max_rel = np.max(np.abs((poneOverKappa_poly(TList) - oneOverKappaList)/oneOverKappaList))
        
        print(f"{species}: mu_RMSE={mu_rmse:.3e}, mu_max_rel={mu_max_rel:.2%}; "
            f"kappa_RMSE={kappa_rmse:.3e}, kappa_max_rel={kappa_max_rel:.2%}; "
            f"1/kappa_RMSE={oneOverKappa_rmse:.3e}, 1/kappa_max_rel={oneOverKappa_max_rel:.2%}")

        f.write(transportTemplate(species, pmu_poly, pkappa_poly, poneOverKappa_poly))
    