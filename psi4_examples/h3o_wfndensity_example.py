import psi4
from psi4 import core
from psi4 import cubeprop
import os

h3o = psi4.geometry("""
nocom
noreorient
O       -2.5732456373      1.4108276299     -0.0272191596                 
H       -1.5824684287      1.4296837220     -0.0179470929                 
H       -2.9531674257      1.9953811385      0.6770802308                 
H       -2.9185085681      0.4847760989      0.0456030216                 
""")
h3o.set_molecular_charge(1)
h3o.set_multiplicity(1)


core.print_global_options()
core.print_options()
psi4.set_memory('500mb')
E, wfn = psi4.prop('HF/sto-3g', properties=["DIPOLE"], return_wfn = True)
print(wfn.Da().to_array())
print(wfn.Db().to_array())
