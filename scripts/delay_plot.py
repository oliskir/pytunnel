import matplotlib.pyplot as plt
import pandas as pd


dfs, labels = [], []

#dfs.append(pd.read_csv('output_atomic_e0.1/log.csv'))
#labels.append(r'$E/V=0.1$')
#dfs.append(pd.read_csv('output_atomic_e0.5/log.csv'))
#labels.append(r'$E/V=0.5$')
#dfs.append(pd.read_csv('output_atomic_e0.95/log.csv'))
#labels.append(r'$E/V=0.95$')
#dfs.append(pd.read_csv('output_atomic_e1.0/log.csv'))
#labels.append(r'$E/V=1.0$')

dfs.append(pd.read_csv('output_nuclear_s600fm_e0.1/log.csv'))
labels.append(r'$\sigma=600$ fm')
dfs.append(pd.read_csv('output_nuclear_s300fm_e0.1/log.csv'))
labels.append(r'$\sigma=300$ fm')


for label,df in zip(labels,dfs):
#    plt.plot(df['pos_T (fm)']*1e-9, df['delay (s)']/1e-15, '-', label=label)
    plt.plot(df['pos_T (fm)'], df['delay (s)']/1e-21, '-', label=label)

#plt.xlim(0, 1.0)
plt.xlim(0, 1000)
plt.ylim(-20, 5)

#plt.xlabel(r'$x$ ($\mu$m)')
plt.xlabel(r'$x$ (fm)')
#plt.ylabel(r'$\delta t$ $(10^{-15}$ s)')
plt.ylabel(r'$\delta t$ $(10^{-21}$ s)')

plt.grid()
plt.legend(prop={'size': 6})

plt.title(r'Delay wrt free prop.') # ($w=20$ nm, $V=0.3$ eV, $E=0.15$ eV, $m=32.2$ keV)')

plt.show()
