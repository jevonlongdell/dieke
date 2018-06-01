import numpy as np
import matplotlib.pyplot as plt
import spectrapy


#nf = input('Enter number of f electrons: ')
nf=12


(LSterms,Uk,LSJlevels,freeion_mat,LSJmJstates,full_freeion_mat,Ckq) = spectrapy.makeMatricies(nf)
#LSterms     - list of LSterm labels
#Uk          - Uk in terms of these terms
#LSJlevels    - list of LSJ levels
#freeion_mat - dictionary of free ion matricies in terms of those levels
#LSJmJstates - list of LSJmJ states
#full_freeion_mat - the free ion matricies again but now in terms of LSJmJ
#Ckq         - Ckq matricies 


#cfparams = spectrapy.readLaF3params(nf)
#reads Carnall 89 crystal field parameters from spreadsheet

cfparams = {}
cfparams['F2'] = 101228
cfparams['F4'] = 71386
cfparams['F6'] = 51369
cfparams['ZETA'] = 2630
#cfparams['B20'] = -194
#cfparams['B40'] = 360
#cfparams['B44'] = 850
#cfparams['B60'] = -622
#cfparams['B64'] = -41



numLSJ=len(LSJlevels)   
numLSJmJ=len(LSJmJstates)
cfstrvals =  np.linspace(0,1.5,15)


H0 = np.zeros([numLSJmJ,numLSJmJ])
for k in cfparams.keys():
    if full_freeion_mat.has_key(k):
        print("Adding free ion parameter %s\n"%(k))
        H0 = H0+cfparams[k]*full_freeion_mat[k]

H=H0
for k in [2,4,6]:
    for q in range(0,k+1):
        if cfparams.has_key('B%d%d'%(k,q)):
            if q==0:
                H = H+cfparams['B%d%d'%(k,q)]*Ckq[(k,q)]
            else:
                H = H+cfparams['B%d%d'%(k,q)]*Ckq[(k,q)]
                H = H+((-1)**q)*np.conj(cfparams['B%d%d'%(k,q)])*Ckq[(k,-q)]
(evals,evects) = np.linalg.eig(H)  
E0 = np.min(evals)
calc_nrg_levels = np.sort(evals-E0)

numLSJmJ = len(LSJmJstates)
nrglevels = np.zeros([len(cfstrvals),numLSJmJ])





l=0
cfparams2 = cfparams
for cfstr in cfstrvals:
    H = H0
    for k in [2,4,6]:
        for q in range(0,k+1):
            if cfparams.has_key('B%d%d'%(k,q)):
                #print("adding in B %d %d"%(k,q))
                if q==0:
                    H = H+cfstr*cfparams['B%d%d'%(k,q)]*Ckq[(k,q)]
                else:
                    H = H+cfstr*cfparams['B%d%d'%(k,q)]*Ckq[(k,q)]
                    H = H+cfstr*((-1)**q)*np.conj(cfparams['B%d%d'%(k,q)])*Ckq[(k,-q)]
    (evals,evects) = np.linalg.eig(H)
    nrglevels[l,:]=np.sort(np.real(evals))-E0
    l=l+1

plt.ion()
plt.axis([0,max(cfstrvals),-2000,50000])
plt.plot(cfstrvals,nrglevels,'x-',[1,1],[-10000,+50000])
plt.draw()

for k in range(91):
    print calc_nrg_levels[k]


    
            
        
        
    
    

        
    
    
