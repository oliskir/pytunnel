import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
from pint import UnitRegistry
from pytunnel.pytunnel_script import particle_velocity

ureg = UnitRegistry()
Q = ureg.Quantity


# quantum
df = pd.read_csv('out_atomic.csv', names=['E/V','time'])
plt.plot(df['E/V'], df['time']/1e-15, '-o', label='Quantum')


fil = open('config_atomic.json', 'r')
d = json.load(fil)
fil.close()

# classical
mass_MeV = Q(d['mass']).m_as("MeV")
energies_eV = np.linspace(0.01,0.6,100) 
velo = particle_velocity(energy=energies_eV*1e-6, mass=mass_MeV)
height_eV = Q(d['barrier_params']['height']).m_as("eV")
width_fm = Q(d['barrier_params']['width']).m_as("fm")
time_class = width_fm / velo

plt.plot(energies_eV/height_eV, time_class/1e-15, '-', label='Classical')


plt.xlim(0,2.0)
plt.ylim(0, 80)

plt.xlabel(r'$E/V_{max}$')
plt.ylabel(r'Time (fs)')

plt.grid()
plt.legend(prop={'size': 11})

plt.title(r'Tunneling time (atomic)')

plt.show()
