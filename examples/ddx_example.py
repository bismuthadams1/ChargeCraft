import psi4
from psi4 import core

h2o = psi4.geometry("""
O
H 1 0.96
H 1 0.96 2 104.5
""")

# psi4.set_options({"ddx": "true", #supply a solvent here to see if epsilon then gets picked up
#                 "DDX_SOLVENT_EPSILON":'5',
#                 "DDX_RADII_SET": 'uff',
#                 "DDX_MODEL": 'PCM',
#                 "DDX_SOLVENT":""})  #without this calculation fails with psi4.driver.p4util.exceptions.ValidationError: Required option 'DDX_SOLVENT' is missing.

psi4.set_options({"PROPERTIES_ORIGIN":[0,0,0]})
core.print_global_options()
psi4.set_memory('500mb')
E, wfn = psi4.prop('scf/sto-3g', properties=["DIPOLE"], return_wfn = True)
# psi4.energy('scf/sto-3g')
print(wfn.variable('DIPOLE'))