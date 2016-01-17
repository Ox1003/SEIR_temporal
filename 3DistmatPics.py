# -*- coding: utf-8 -*-
"""
heatmaps.py

creates heatmaps based on input distance matrices
@author: 
EnCt28648cf14f7c30972de440f01c53b7544e83c9ee78648cf14f7c30972de440f01puAufN5B6QB
rpWuQm1ZPfSb3bcvooZYGIwEmS

Decrypt it at https://encipher.it
"""

import matplotlib.pyplot as plt
import numpy as np

Deff1=np.load('D.npy')
Deff2=np.load('Deff.npy')
spars=1

if spars==1:    
    
    x    = range(242,-1,-1)
    y    = range(243)
    x, y = np.meshgrid(x, y)
    
    fig=plt.figure(figsize=(9, 4.5))
    plt.rcParams.update({'font.size': 9})
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.subplot(1,2,1)
    z=Deff1
    
    z[z==0]=1.1*np.max(z)
    z[np.diag_indices(242)]=0.0
    plt.pcolormesh(x, y, np.rot90(z),clim=(243,243))
    plt.colorbar() #need a colorbar to show the intensity scale
    plt.title('First order contacts (inf distance scaled to value %s)'%np.round(np.max(z),1))
    plt.xlim(xmax=243) #or xl
    plt.ylim(ymax=243) #or yl
    
    plt.subplot(1,2,2)
    z=Deff2
    z[np.diag_indices(242)]=0.0
    plt.pcolormesh(x,y,np.rot90(z), clim=(243,243))
    plt.colorbar() #need a colorbar to show the intensity scale
    plt.title('Up to third order contacts')
    plt.xlim(xmax=243) #or xl
    plt.ylim(ymax=243) #or yl
    plt.show()  
    plt.subplots_adjust(hspace=0.6)
    plt.suptitle('Effective distances (1-log(P_{mn}))',fontsize=12)     
    fig.savefig('effDistMatricesDeff.png', dpi=500) 

