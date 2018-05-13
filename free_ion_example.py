import numpy as np
import matplotlib.pyplot as plt
import spectrapy

## Example of free ion calculation using spectrapy
##
## This makes Figure 3.1 in Mike Reids notes:
## which shows the energy levesl of Pr3+ as the strength
## of the spin orbit interaction is increased from zero to
## about 150% of it's actual value.
##
## http://www2.phys.canterbury.ac.nz/~mfr24/electronicstructure/00electronic.pdf


nf=2 # the number of f electrions

#this reads in the crystal field parameters for 
#cfparams = spectrapy.readLaF3params(nf)

(LSJlevels,fi_mat,LSterms,Uk,V)= spectrapy.read_crosswhite(nf)


cfparams = spectrapy.readLaF3params(nf)
zeta0 = cfparams['ZETA']
zetavals = np.linspace(0,1.5*zeta0,100)

numLSJ=len(LSJlevels)
nrglevels = np.zeros([len(zetavals),numLSJ])

H0 = np.zeros([numLSJ,numLSJ])
for k in cfparams.keys():
    if k in fi_mat:
        print("using parameter ",k)
        H0 = H0+cfparams[k]*fi_mat[k]
(evals,evects) = np.linalg.eig(H0)
E0 = np.min(evals)

l=0
cfparams2 = cfparams
for zeta in zetavals:
    cfparams2['ZETA']=zeta
    H = np.zeros([numLSJ,numLSJ])
    for k in cfparams2.keys():
        if k in fi_mat:
            H = H+cfparams2[k]*fi_mat[k]
    (evals,evects) = np.linalg.eig(H)
    nrglevels[l,:]=np.sort(evals)-E0
    l=l+1
plt.axis([0,max(zetavals),-2000,25000])
plt.plot(zetavals,nrglevels,[zeta0,zeta0],[-10000,+50000])
plt.show()
    
            
        
        
    
    
