import os
import json
import numpy as np
import pandas as pd
from pint import UnitRegistry
from pytunnel.pytunnel_script import run, particle_velocity

ureg = UnitRegistry()
Q = ureg.Quantity


# energies (relative to barrier height)
energies = [1.02, 1.05, 1.10, 1.15, 1.20, 1.40]

# desired tunneling time precision 
t_prec = 1E-15 #1fs


# load config
fil = open('config_atomic.json', 'r')
d = json.load(fil)
fil.close()

# propagation distance
dist_fm = Q(d['initial_sep']).m_as("fm") + Q(d['x_stop']).m_as("fm")

# barrier height and width
height_eV = Q(d['barrier_params']['height']).m_as("eV")
width_fm = Q(d['barrier_params']['width']).m_as("fm")

# output file
path_out = 'out_atomic.csv'
if os.path.isfile(path_out): os.remove(path_out)


# loop over energies
for i,e in enumerate(energies):

    print(f'{i+1}/{len(energies)}')

    # output directory
    d['output_dir'] = f'output_atomic_e{e}'

    # energy in eV
    e_eV = e * height_eV
    d['energy'] = f'{e_eV}eV'

    # velocity (fm/s)
    mass_MeV = Q(d['mass']).m_as("MeV")
    velo = particle_velocity(energy=e_eV*1e-6, mass=mass_MeV)

    # classical tunneling time
    t_clas = width_fm / velo #s

    # relative accuracy
    d['rel_accuracy'] = t_prec * velo / dist_fm
    
    run(**d) #run simulation

    # collect output
    df = pd.read_csv(os.path.join(d['output_dir'], 'log.csv'))
    delay = df['delay (s)'].values[-1]

    # quantum tunneling time
    t_quan = t_clas + delay

    fout = open(path_out, 'a')
    fout.write(f'{e}, {t_quan}\n')
    fout.close()

    
print(f'saved to {fout.name}')
fout.close()
