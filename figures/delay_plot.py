import matplotlib.pyplot as plt
import pandas as pd


dfs, labels = [], []

dfs.append(pd.read_csv('log_nuclear.csv'))
labels.append(r'$w=20$ fm, $V=1$ MeV, $E=0.5$ MeV')

#dfs.append(pd.read_csv('log_D3um_s0.2um_FINE.csv'))
#labels.append(r'd=3$\mu$m, $\sigma=0.2 \mu$m FINE')

#dfs.append(pd.read_csv('log_D3um_s0.3um_FINE.csv'))
#labels.append(r'd=3$\mu$m, $\sigma=0.3 \mu$m FINE')

#dfs.append(pd.read_csv('log_D3um_s0.4um_FINE.csv'))
#labels.append(r'd=3$\mu$m, $\sigma=0.4 \mu$m FINE')

#dfs.append(pd.read_csv('log_D3um_s0.5um_FINE.csv'))
#labels.append(r'd=3$\mu$m, $\sigma=0.5 \mu$m FINE')

#dfs.append(pd.read_csv('log_D2um_s0.3um.csv'))
#labels.append(r'd=2$\mu$m, $\sigma=0.3 \mu$m')

#dfs.append(pd.read_csv('log_D2um_s0.2um.csv'))
#labels.append(r'd=2$\mu$m, $\sigma=0.2 \mu$m')

for label,df in zip(labels,dfs):
#    plt.plot(df['pos_T (fm)']*1e-9, df['delay (s)']/1e-15, '-', label=label)
    plt.plot(df['pos_T (fm)'], df['delay (s)']/1e-22, '-o', label=label)

plt.xlim(0,100)
#plt.xlim(0,100)
plt.ylim(-50, 10)

plt.xlabel(r'$x$ (fm)')
plt.ylabel(r'$\delta t$ $(10^{-22}$ s)')

plt.grid()
plt.legend(prop={'size': 6})

plt.title(r'Delay wrt free prop.') # ($w=20$ nm, $V=0.3$ eV, $E=0.15$ eV, $m=32.2$ keV)')

plt.show()
