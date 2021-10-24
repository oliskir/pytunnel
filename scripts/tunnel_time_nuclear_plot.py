import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
from pint import UnitRegistry
from pytunnel.pytunnel_script import particle_velocity

ureg = UnitRegistry()
Q = ureg.Quantity


# quantum
#df = pd.read_csv('out_nuclear_00.csv', names=['E/V','time'])
#plt.plot(df['E/V'], df['time']/1e-21, '-o', label='Quantum')


fil = open('config_nuclear.json', 'r')
d = json.load(fil)
fil.close()

# classical
mass_MeV = Q(d['mass']).m_as("MeV")
energies_MeV = np.linspace(0.02,2.0,100) 
velo = particle_velocity(energy=energies_MeV, mass=mass_MeV)
height_MeV = Q(d['barrier_params']['height']).m_as("MeV")
width_fm = Q(d['barrier_params']['width']).m_as("fm")
time_class = width_fm / velo

plt.plot(energies_MeV/height_MeV, time_class/1e-21, '-', label='Classical')


plt.xlim(0,2.0)
plt.ylim(0, 20)

plt.xlabel(r'$E/V_{max}$')
plt.ylabel(r'Time ($10^{-21}$ s)')

plt.grid()
plt.legend(prop={'size': 11})

plt.title(r'Tunneling time')

plt.show()
