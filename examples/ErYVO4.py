import numpy as np
import dieke
import pickle
import gzip
import matplotlib.pyplot as plt

nf = 11  # 11 f-electrons means we're dealing with Er

# This object contains the matricies we need for the
# calculations all in the dictionary "FreeIonMatrix"

try:
    f = gzip.open('Ermatricies.dat.gz', 'rb')
    print("Loading matricies")
    Er = pickle.load(f)
except FileNotFoundError:
    print("Making matricies")
    Er = dieke.RareEarthIon(nf)
    pickle.dump(Er, gzip.open('Ermatricies.dat.gz', 'wb'))


# Read in a a set of crystal field parameters from Er:LaF3
# dieke reads these from (incomplete) carnall89params.xls
cfparams = dieke.readLaF3params(nf)


# Number of levels and states
numLSJ = Er.numlevels()
numLSJmJ = Er.numstates()


# Make a free-ion Hamiltonian use LaF3 parameters
H0 = np.zeros([numLSJmJ, numLSJmJ])
for k in cfparams.keys():
    if k in Er.FreeIonMatrix:
        print("Adding free ion parameter %s\n" % (k))
        H0 = H0+cfparams[k]*Er.FreeIonMatrix[k]



# cfparams['B20'] = 342
# cfparams['B40'] = 109
# cfparams['B60'] = -733
# cfparams['B44'] = -695
# cfparams['B64'] = -105

# cfparams['B20'] = 308.2
# cfparams['B40'] = 102.2
# cfparams['B60'] = -735.7
# cfparams['B44'] = -716.3
# cfparams['B64'] = -25.5

# cfparams['B20'] = 151.7
# cfparams['B40'] = 19.9
# cfparams['B60'] = -44.2
# cfparams['B44'] = -921
# cfparams['B64'] = -97


cfparams['B20'] = -206
cfparams['B40'] = 362
cfparams['B44'] = 926
cfparams['B60'] = -688
cfparams['B64'] = 31.5



# zeroify all cf params that are incompatible with YVOs symmetry
H = H0
for k, q in [(2, 0), (4, 0), (6, 0), (4, 4), (6, 4)]:
        if 'B%d%d' % (k, q) in cfparams:
            
            if q == 0:
                H = H+cfparams['B%d%d' % (k, q)]*Er.Cmatrix(k, q)
            else:
                Bkq = cfparams['B%d%d' % (k, q)]
                Bkmq = (-1)**q*np.conj(Bkq)
                Ckq = Er.Cmatrix(k, q)
                Ckmq = Er.Cmatrix(k, -q)
                H = H + Bkq*Ckq + Bkmq*Ckmq


# we now have the hamiltonian at zero magnetic field
H0 = H               
(evals, vec) = np.linalg.eig(H0)

idx = np.argsort(evals)
evals = evals[idx]
E0 = evals[0] #ground state energy
evals = evals - E0 #make everything relative to GS
vec = vec[:, idx]


Bvals = np.linspace(0,0.1,10)
nrglevels = np.zeros([len(Bvals), 4+0*len(evals)])

muB = 0.46686 # Bohr magneton in wavenumbers per testla

for i, B in enumerate(Bvals):
#    cfparams2['ZETA'] = zeta
    H = H0 +B*muB*(Er.FreeIonMatrix['Lz']+ 2*Er.FreeIonMatrix['Sz'])
    (evals, evects) = np.linalg.eig(H)
    evals = np.sort(np.real(evals-E0))
    evals = np.take(evals,[0,1,18,19])
    nrglevels[i, :] = np.sort(np.real(evals-E0))


nrglevels = nrglevels*29.98 #  to GHz
F0 = 6589.0*29.98



    
plt.ion()  # interactive plotting, so plt.show() doesn't block
#plt.axis([0, max(Bvals), -10, 7000])
plt.plot(Bvals, nrglevels[:,3]-nrglevels[:,0] - F0,
         Bvals, nrglevels[:,3]-nrglevels[:,1] - F0,
         Bvals, nrglevels[:,2]-nrglevels[:,0] - F0,
         Bvals, nrglevels[:,2]-nrglevels[:,1] - F0)  #  , [zeta0, zeta0], [-10000, +50000])
plt.show()






