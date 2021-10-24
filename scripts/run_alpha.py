import os
import json
import numpy as np
import pandas as pd
from pint import UnitRegistry
from pytunnel.pytunnel_script import run, particle_velocity

ureg = UnitRegistry()
Q = ureg.Quantity


# energies (relative to barrier height)
energies = [0.1,0.5,0.95,1.0,1.02,1.05,1.1,1.15,1.2,1.25,1.30,1.35,1.4,1.5,1.8]

# desired tunneling time precision 
t_prec = 2E-22 #1fs


# load config
fil = open('config_alpha.json', 'r')
d = json.load(fil)
fil.close()

# propagation distance
dist_fm = 0.5 * Q(d['domain_size']).m_as("fm")

# barrier height and width
height_MeV = Q(d['barrier_params']['initial_height']).m_as("MeV")
width_fm = Q(d['barrier_params']['width']).m_as("fm")

# output file
path_out = 'out_alpha.csv'
if os.path.isfile(path_out): os.remove(path_out)


# loop over energies
for i,e in enumerate(energies):

    print(f'{i+1}/{len(energies)}')

    # output directory
    d['output_dir'] = f'output_alpha_e{e}'

    # energy in eV
    e_MeV = e * height_MeV
    d['energy'] = f'{e_MeV}MeV'

    # velocity (fm/s)
    mass_MeV = Q(d['mass']).m_as("MeV")
    velo = particle_velocity(energy=e_MeV, mass=mass_MeV)

    # classical tunneling time
    t_clas = width_fm / velo #s

    # relative accuracy
    d['rel_accuracy'] = t_prec * velo / dist_fm
    
    run(**d) #run simulation

    # collect output
    df = pd.read_csv(os.path.join(d['output_dir'], 'log.csv'))
    delay = df['delay_ev (s)'].values[-1]

    # quantum tunneling time
    t_quan = t_clas + delay

    fout = open(path_out, 'a')
    fout.write(f'{e}, {t_quan}\n')
    fout.close()

    
print(f'saved to {fout.name}')
fout.close()
