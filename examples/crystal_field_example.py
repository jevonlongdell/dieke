import numpy as np
import matplotlib.pyplot as plt
import dieke


# Calculates the energy levels of Pr:LaF3 as the strength of the
# crystal field terms are varied from zero to 150% or their actual
# values.


nf = 2  # 2 f-electrons means we're dealing with Pr

# This object contains the matricies we need for the
# calculations all in the dictionary "FreeIonMatrix"
Pr = dieke.RareEarthIon(nf)


# Read in a a set of crystal field parameters from Pr:LaF3
# dieke reads these from (incomplete) carnall89params.xls
cfparams = dieke.readLaF3params(nf)


# Number of levels and states
numLSJ = Pr.numlevels()
numLSJmJ = Pr.numstates()

# Variable multiplier for the crystal field strength
# this will be the x-axis of the graph we will draw
cfstrvals = np.linspace(0, 1.5, 15)


# Make a free-ion Hamiltonian
H0 = np.zeros([numLSJmJ, numLSJmJ])
for k in cfparams.keys():
    if k in Pr.FreeIonMatrix:
        print("Adding free ion parameter %s\n" % (k))
        H0 = H0+cfparams[k]*Pr.FreeIonMatrix[k]

## Add in the crystal field terms and diagonalise the result
#H = H0
#for k in [2, 4, 6]:
#    for q in range(0, k+1):
#        if 'B%d%d' % (k, q) in cfparams:
#            if q == 0:
#                H = H+cfparams['B%d%d' % (k, q)]*Pr.Cmatrix(k, q)
#                print("adding in B %d %d \t = %d"%(k, q, cfparams['B%d%d' % (k, q)]))
#            else:
#                H = H+cfparams['B%d%d' % (k, q)]*Pr.Cmatrix(k, q)
#                print("adding in B %d %d \t = %d"%(k, q, cfparams['B%d%d' % (k, q)]))
#                H = H+((-1)**q)*np.conj(cfparams['B%d%d' % (k, q)]) \
#                    * Pr.Cmatrix(k, -q)
#                print("adding in B %d %d \t = %d"%(k, -q, ((-1)**q)*np.conj(cfparams['B%d%d' % (k, q)])))
#(evals, evects) = np.linalg.eig(H)
#E0 = np.min(evals)
#calc_nrg_levels = np.sort(evals-E0)

# TESTING Using Hermitian matrices
Omega = {}
for k in [2, 4, 6]:
    for q in range(0, k+1):
        if 'B%d%d' % (k, q) in cfparams:
            if q == 0:
                label = "O%d%d" % (k,q)
                Omega[label] = Pr.Cmatrix(k, q) 
            else:
                label_pos = "O%d%d" % (k,q)
                label_neg = "O%d%d" % (k,-q)
                Omega[label_pos] = Pr.Cmatrix(k, q) + ((-1)**q)*Pr.Cmatrix(k,-q)
                Omega[label_neg] = complex(0,1)*(Pr.Cmatrix(k, -q) - ((-1)**q)*Pr.Cmatrix(k,q))

H = H0
for k in [2, 4, 6]:
    for q in range(0, k+1):
        if 'B%d%d' % (k, q) in cfparams:
            if q == 0:
                H = H+cfparams['B%d%d' % (k, q)]*Omega['O%d%d' % (k, q)]
                print("adding in B %d %d \t = %d"%(k, q, cfparams['B%d%d' % (k, q)]))
            else:
                H = H+cfparams['B%d%d' % (k, q)]*Omega['O%d%d' % (k, q)]
                print("adding in B %d %d \t = %d"%(k, q, cfparams['B%d%d' % (k, q)]))
                
                H = H+((-1)**q)*np.conj(cfparams['B%d%d' % (k, q)]) \
                    * Omega['O%d%d' % (k, -q)]
                print("adding in B %d %d \t = %d"%(k, -q, ((-1)**q)*np.conj(cfparams['B%d%d' % (k, q)])))
(evals, evects) = np.linalg.eig(H)
E0 = np.min(evals)
calc_nrg_levels = np.sort(evals-E0)


## TESTING change in cf params
#
#H = H0
#for k in [2, 4, 6]:
#    for q in range(0, k+1):
#        if 'B%d%d' % (k, q) in cfparams:
#            if q == 0:
#                H = H+cfparams['B%d%d' % (k, q)]*Pr.Cmatrix(k, q)
#            else:
#                Bkq = cfparams['B%d%d' % (k, q)]
#                Bknegq = 1/((-1)**q)*Bkq
#                Bkq_hat = complex(Bkq,Bknegq)
#                Bknegq_hat = 1/((-1)**q)*np.conj(Bkq_hat)
#                print('B%d%d hat = %d + j%d' % (k, q, np.real(Bkq_hat),np.imag(Bkq_hat)))
#                H = H+ Bkq_hat*Pr.Cmatrix(k, q)
#                H = H+ Bknegq_hat*Pr.Cmatrix(k,-q)
#(evals, evects) = np.linalg.eig(H)
#E0 = np.min(evals)
#calc_nrg_levels = np.sort(evals-E0)
#
#
# Print out the 20 lowest enery levels
for k in range(20):
    print(np.real(calc_nrg_levels[k]))


## Hermition tests
#    
#EPS=3e-7
#Hcf = np.zeros([numLSJmJ, numLSJmJ])
#for k in [2, 4, 6]:
#    for q in range(0, k+1):
#        if 'B%d%d' % (k, q) in cfparams:
#            if q == 0:
#                matrix_add = cfparams['B%d%d' % (k, q)]*Pr.Cmatrix(k, q)
#                
#                # Checking if matrix is hermitian
#                if np.linalg.norm(matrix_add-matrix_add.H) < EPS:
#                    print("B %d %d \t is Hermitian"%(k, q))
#                else:
#                    print("B %d %d \t not Hermitian"%(k, q))
#                
#                Hcf = Hcf + matrix_add
#            else:
#                matrix_add_1 = cfparams['B%d%d' % (k, q)]*Pr.Cmatrix(k, q)
#                matrix_add_2 = ((-1)**q)*np.conj(cfparams['B%d%d' % (k, q)]) \
#                    * Pr.Cmatrix(k, -q)
#                # Checking if Hermitian
#                if np.linalg.norm(matrix_add_1 - matrix_add_1.H) < EPS:
#                    print("B %d %d \t is Hermitian"%(k, q))
#                else:
#                    print("B %d %d \t not Hermitian"%(k, q))
#                
#                if np.linalg.norm(matrix_add_2 - matrix_add_2.H) < EPS:
#                    print("B %d %d \t is Hermitian"%(k, -q))
#                else:
#                    print("B %d %d \t not Hermitian"%(k, -q))
#                
#                Hcf = Hcf+matrix_add_1 + matrix_add_2
#                
#                
## Checking Hcf hermitian
#                
#if np.linalg.norm(Hcf - Hcf.H) < EPS:
#    print("Hcf \t is Hermitian")
#else:
#    print("Hcf \t not Hermitian")
    