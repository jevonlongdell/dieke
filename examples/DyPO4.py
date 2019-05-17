import numpy as np
import dieke
import pickle
import gzip

# Calculates some approximation to the energy levels in DyPO4
# uses carnal free ion parameters
# uses crystal field parameters from Loong et al Jounal of Alloys and Compuonds
# 207/208 1994 165-169
# https://doi.org/10.1016/0925-8388(94)90198-8

nf = 9  # 9 f-electrons means we're dealing with Dy

# This object contains the matricies we need for the
# calculations all in the dictionary "FreeIonMatrix"

try:
    f = gzip.open('Dymatricies.dat.gz', 'rb')
    print("Loading matricies")
    Dy = pickle.load(f)
except FileNotFoundError:
    print("Making matricies")
    Dy = dieke.RareEarthIon(nf)
    pickle.dump(Dy, gzip.open('Dymatricies.dat.gz', 'wb'))


# Read in a a set of crystal field parameters from Dy:LaF3
# dieke reads these from (incomplete) carnall89params.xls
cfparams = dieke.readLaF3params(nf)


# Number of levels and states
numLSJ = Dy.numlevels()
numLSJmJ = Dy.numstates()

# Variable multiplier for the crystal field strength
# this will be the x-axis of the graph we will draw
cfstrvals = np.linspace(0, 1.5, 15)


# Make a free-ion Hamiltonian use LaF3 parameters
H0 = np.zeros([numLSJmJ, numLSJmJ])
for k in cfparams.keys():
    if k in Dy.FreeIonMatrix:
        print("Adding free ion parameter %s\n" % (k))
        H0 = H0+cfparams[k]*Dy.FreeIonMatrix[k]



# cfparams['B20'] = 342
# cfparams['B40'] = 109
# cfparams['B60'] = -733
# cfparams['B44'] = -695
# cfparams['B64'] = -105

cfparams['B20'] = 308.2
cfparams['B40'] = 102.2
cfparams['B60'] = -735.7
cfparams['B44'] = -716.3
cfparams['B64'] = -25.5

# cfparams['B20'] = 151.7
# cfparams['B40'] = 19.9
# cfparams['B60'] = -44.2
# cfparams['B44'] = -921
# cfparams['B64'] = -97


# zeroify all cf params that are incompatible with YPO4s symmetry
H = H0
for k, q in [(2, 0), (4, 0), (6, 0), (4, 4), (6, 4)]:
        if 'B%d%d' % (k, q) in cfparams:
            
            if q == 0:
                H = H+cfparams['B%d%d' % (k, q)]*Dy.Cmatrix(k, q)
            else:
                Bkq = cfparams['B%d%d' % (k, q)]
                Bkmq = (-1)**q*np.conj(Bkq)
                Ckq = Dy.Cmatrix(k, q)
                Ckmq = Dy.Cmatrix(k, -q)
                H = H + Bkq*Ckq + Bkmq*Ckmq

(evals, vec) = np.linalg.eig(H)

idx = np.argsort(evals)

evals = evals[idx]
evals = evals - evals[0]
vec = vec[:, idx]


print("Makeup of first 8 states")
for st in range(8):
    for i in range(numLSJmJ):
        if abs(vec[i, st]) > 1e-2:
            print("%5.5f %+5.5f  %s" % (vec[i, st].real,
                                        vec[i, st].imag,
                                        Dy.LSJmJstateLabels[i]))
    print('\n\n')


for i in range(98):
    print(evals[i])
