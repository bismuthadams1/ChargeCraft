from ase import Atoms
from ase.calculators.gaussian import Gaussian, GaussianOptimizer
import ase.io
import os

cwd = os.getcwd()

no_conformers = 10

converged_info = ''

for conf_no in range(no_conformers):
    conformer = ase.io.read(f'/mnt/storage/nobackup/nca121/test_jobs/QM_ESP_Psi4/QM_ESP_Psi4/conformer_{conf_no}.xyz')

    calc_opt = Gaussian(label = f'conformer_{conf_no}_opt',
                        method = 'hf',
                        basis = '6-31G*',
                        scf='qc'
                        )
    opt = GaussianOptimizer(conformer, calc_opt)
    converged = opt.run(fmax='tight', steps=100)

    converged_info += f'conformer: {conformer} converged? {converged}\n'

f = open(cwd+f"/gaussian_log", 'x')
f.write(converged_info)
f.close()



