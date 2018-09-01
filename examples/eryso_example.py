import numpy as np
#mport matplotlib.pyplot as plt
import dieke



nf = 11  # 2 f-electrons means we're dealing with Er

# This object contains the matricies we need for the
# calculations all in the dictionary "FreeIonMatrix"
Er = dieke.RareEarthIon(nf)

# Read in a a set of crystal field parameters from Er:LaF3
# dieke reads these from (incomplete) carnall89params.xls
cfparams = dieke.readLaF3params(nf)


#change most of the parameters to those in eryso_cf paper
#leave some from carnall89

cfparams['E0'] = 35503.5 #this is ignored
cfparams['ZETA'] = 2362.9
cfparams['F2'] = 96029.6
cfparams['F4'] = 67670.6
cfparams['F6'] = 53167.1

cfparams['B20'] = -149.8
cfparams['B21'] = 420.6+396.0j
cfparams['B22'] = -228.5+27.6j
cfparams['B40'] = 1131.2
cfparams['B41'] = 985.7+34.2j
cfparams['B42'] = 296.8+145.0j
cfparams['B43'] = -402.3-381.7j
cfparams['B44'] = -282.3+1114.3j
cfparams['B60'] = -263.2
cfparams['B61'] = 111.9+222.9j
cfparams['B62'] = 124.7+195.9j
cfparams['B63'] = -97.9+139.7j
cfparams['B64'] = -93.7-145.0j
cfparams['B65'] = 13.9+109.5j
cfparams['B66'] = 3.0-108.6j
cfparams['A'] = 0.005466 #this is ignored
cfparams['Q'] = 0.0716 #this is ignored


# Number of levels and states
numLSJ = Er.numlevels()
numLSJmJ = Er.numstates()


# Make a free-ion Hamiltonian
H0 = np.zeros([numLSJmJ, numLSJmJ])
for k in cfparams.keys():
    if k in Er.FreeIonMatrix:
        print("Adding free ion parameter %s" % (k))
        H0 = H0+cfparams[k]*Er.FreeIonMatrix[k]

# Add in the crystal field terms and diagonalise the result
H = H0
for k in [2, 4, 6]:
    for q in range(0, k+1):
        if 'B%d%d' % (k, q) in cfparams:
            if q == 0:
                H = H+cfparams['B%d%d' % (k, q)]*Er.Cmatrix(k, q)
            else:
                H = H+cfparams['B%d%d' % (k, q)]*Er.Cmatrix(k, q)
                H = H+((-1)**q)*np.conj(cfparams['B%d%d' % (k, q)]) \
                    * Er.Cmatrix(k, -q)
(evals, evects) = np.linalg.eig(H)



E0 = np.min(evals)
nrg_levels = np.sort(np.real(evals-E0))#cfparams['E0'])
nrg_levels = nrg_levels[0::2] #take every second energy level

# # Plot the results
# plt.ion()
# plt.axis([0, max(cfstrvals), -2000, 25000])
# plt.plot(cfstrvals, nrglevels, 'x-', [1, 1], [-10000, +50000])
# plt.draw()

nrg_levels = np.array(list(map(int,nrg_levels)))
levels_seb = np.array([15, 47, 75, 130,  199, 388, 462, 508,
                    6522, 6560, 6583, 6640, 6777, 6833, 6867])

# Print out the 15 lowest enery level doublets
print('-------\n Calculated energy levels\n comparison with Sebastian \n ---------')
print ("jevon, sebastian, diff")
for k in range(15):
    print(nrg_levels[k],levels_seb[k],nrg_levels[k]-levels_seb[k])

print()
print('---------------------')
print('Parameters I\'m using')
print('---------------------')
print(cfparams)


    
# comdiff_seb = np.mean(levels_seb[8:14]) - np.mean(levels_seb[0:7])

# comdiff = np.mean(nrg_levels[8:14]) - np.mean(nrg_levels[0:7])


# print(comdiff, comdiff_seb, comdiff-comdiff_seb)
