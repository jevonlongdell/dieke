import numpy as np
import matplotlib.pyplot as plt
import spectrapy


#nf = input('Enter number of f electrons: ')
nf=2
cfparams = spectrapy.readLaF3params(nf)
(LSJlevels,fi_mat,LSterms,Uk,V)= spectrapy.read_crosswhite(nf)

#make figure 3.1 in mikes notes

cfparams = spectrapy.readLaF3params(nf)
zeta0 = cfparams['ZETA']
zetavals = np.linspace(0,1.5*zeta0,100)

numLSJ=len(LSJlevels)
nrglevels = np.zeros([len(zetavals),numLSJ])

H0 = np.zeros([numLSJ,numLSJ])
for k in cfparams.keys():
    if fi_mat.has_key(k):
        print "using parameter ",k
        H0 = H0+cfparams[k]*fi_mat[k]
(evals,evects) = np.linalg.eig(H0)
E0 = np.min(evals)

l=0
cfparams2 = cfparams
for zeta in zetavals:
    cfparams2['ZETA']=zeta
    H = np.zeros([numLSJ,numLSJ])
    for k in cfparams2.keys():
        if fi_mat.has_key(k):
            H = H+cfparams2[k]*fi_mat[k]
    (evals,evects) = np.linalg.eig(H)
    nrglevels[l,:]=np.sort(evals)-E0
    l=l+1
plt.axis([0,max(zetavals),-2000,25000])
plt.plot(zetavals,nrglevels,[zeta0,zeta0],[-10000,+50000])
plt.show()
    
            
        
        
    
    
