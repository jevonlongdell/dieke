import numpy as np
import dieke
import pickle
import gzip
import matplotlib.pyplot as plt
import sys

nf = 11  # 11 f-electrons means we're dealing with Er


#only compute the matrix elements the first run through
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



cfparams['B22'] = 0
cfparams['B42'] = 0
cfparams['B62'] = 0
cfparams['B66'] = 0

     

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
vec = vec[:, idx] #sort eigenvalues


# Add a fixed energy to the hamiltonian to make the ground state
# zero energy
H = H - np.identity(len(evals))*E0





Bvals = np.linspace(0,2.1,20)
nrglevels = np.zeros([len(Bvals), len(evals)])

muB = 0.46686 # Bohr magneton in wavenumbers per testla

# angle magnetic field makes with Z access
th = 0



print('Calculating levels')
for i, B in enumerate(Bvals):
    print(i,end=" ")
    sys.stdout.flush()
#    cfparams2['ZETA'] = zeta
    H = H0 +B*muB*np.cos(th)*(Er.FreeIonMatrix['Lz']+ 2*Er.FreeIonMatrix['Sz'])+\
        B*muB*np.sin(th)*(Er.FreeIonMatrix['Lx']+ 2*Er.FreeIonMatrix['Sx'])
    (evals, evects) = np.linalg.eig(H)
    evals = np.sort(np.real(evals))
#    evals = np.take(evals,[0,1,18,19])
    nrglevels[i, :] = np.sort(np.real(evals-E0))
print()


# nrglevels = nrglevels*29.98 #  to GHz
# F0 = 6589.0*29.98


nlevels = 30

lines = np.zeros([len(Bvals),nlevels*(nlevels-1)])


print('Calculating lines')
count = 0
for i in [0,1]: #range(nlevels-1): #lines out of the ground two states
	for j in range(i+1,nlevels):
		lines[:,count] = np.abs(nrglevels[:,i]-nrglevels[:,j])
		count=count+1


   
#plt.ion()  # interactive plotting, so plt.show() doesn't block

f1 = plt.figure(1)
plt.title('Spectral lines around 1520nm')
plt.axis([0, max(Bvals), 6584, 6594])
plt.plot(Bvals, (lines))  #  , [zeta0, zeta0], [-10000, +50000])
plt.xlabel('B field (T)')
plt.ylabel('Energy (cm^-1)')
plt.show()



f1 = plt.figure(2)
plt.title('Energy Levels for Er:YVO4')
plt.subplot(1,2,1)
plt.axis([0, max(Bvals), -3,3])
plt.plot(Bvals, nrglevels)  #  , [zeta0, zeta0], [-10000, +50000])
plt.xlabel('B field (T)')
plt.ylabel('Energy (cm^-1)')

plt.subplot(1,2,2)
plt.axis([0, max(Bvals), 6584,6594])
plt.plot(Bvals, nrglevels)  #  , [zeta0, zeta0], [-10000, +50000])
plt.xlabel('B field (T)')
plt.ylabel('Energy (cm^-1)')

plt.show()






