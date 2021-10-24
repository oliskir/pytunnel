import os
import json
import numpy as np
from pytunnel.pytunnel_script import run


pos = np.linspace(-20, 40, num=30) 

fil = open('config_alpha.json', 'r')
d = json.load(fil)
fil.close()

os.remove('out.csv')

for i,p in enumerate(pos):

    print(f'{i+1}/{len(pos)}')

    d['barrier_params']['transition_pos'] = f'{p}fm' #replace barrier parameter value

    run(**d) #run simulation

    f = open('output/output.txt', 'r')
    for i in range(5): l = f.readline()
    T = float(l[5:])
    f.close()

    fout = open('out.csv', 'a')
    fout.write(f'{p}, {T}\n')
    fout.close()
    
print(f'saved to {fout.name}')
fout.close()
