import psi4
from psi4 import core

h2o = psi4.geometry("""
O
H 1 0.96
H 1 0.96 2 104.5
""")

psi4.set_options({"PCM":True})
psi4.set_options({ "pcm__input":  f"""
                                            Units = Angstrom
                                            Medium {{
                                                SolverType = 'CPCM'
                                                Solvent = 'water'
                                            }}

                                            Cavity {{
                                                RadiiSet = 'Bondi' # Bondi | UFF | Allinger
                                                Type = GePol
                                                Scaling = True # radii for spheres scaled by 1.2
                                                Area = 0.3
                                                Mode = "Implicit"
                                            }}
                                            """} )#without this calculation fails with psi4.driver.p4util.exceptions.ValidationError: Required option 'DDX_SOLVENT' is missing.

psi4.set_options({"PROPERTIES_ORIGIN":[0,0,0]})
core.print_global_options()
psi4.set_memory('500mb')
E, wfn = psi4.prop('scf/sto-3g', properties=["DIPOLE"], return_wfn = True)
# psi4.energy('scf/sto-3g')
print(wfn.variable('DIPOLE'))