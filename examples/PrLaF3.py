import numpy as np
import dieke

###########################################
# Calculates the energy levels of Pr:LaF3 
###########################################


nf = 2  # 2 f-electrons means we're dealing with Pr

# This object contains the matricies we need for the
# calculations all in the dictionaries "FreeIonMatrix"
# and "Cmatrix"
print('Populating matricies')
Pr = dieke.RareEarthIon(nf)


# Creating the rare earth ion object can be timeconsuming,
# especially as the number of f-electrons gets closer to 14/2=7.
# One practical way to get around this is using code like the following
# which caches the matricies.
#
# import gzip
# import pickle
# try:
#     f = gzip.open('pr_matricies.dat.gz', 'rb')
#     print("Loading matricies")
#     Pr = pickle.load(f)
# except FileNotFoundError:
#     print("Making matricies")
#     Pr = dieke.RareEarthIon(nf)
#     pickle.dump(Er, gzip.open('pr_matricies.dat.gz', 'wb'))




# Read in a a set of crystal field parameters for Pr:LaF3
# dieke reads these from carnall89params.xls included in 
# the dieke package. (This spreadsheet isn't fully populated
# and this is an option if you  want to help the dieke project.)
#
# Data from:
# A systematic analysis of the spectra of the lanthanides doped into single crystal LaF3
# W. T. Carnall, J. Chem. Phys. 90, 3443-3457(1989)
# https://doi.org/10.1063/1.455853
cfparams = dieke.readLaF3params(nf)


# Number of levels and states
numLSJ = Pr.numlevels()
numLSJmJ = Pr.numstates()

# Make a free-ion Hamiltonian
H0 = np.zeros([numLSJmJ, numLSJmJ])
for k in cfparams.keys():
    if k in Pr.FreeIonMatrix:
        print("Adding free ion parameter %s\n" % (k))
        H0 = H0 + cfparams[k]*Pr.FreeIonMatrix[k]

# Add in the crystal field terms
H = H0
for k in [2, 4, 6]:
    for q in range(0, k+1):
        if 'B%d%d' % (k, q) in cfparams:
            if q == 0:
                H = H+cfparams['B%d%d' % (k, q)]*Pr.Cmatrix(k, q)
            else:
                Bkq = cfparams['B%d%d' % (k, q)]
                Bkmq = (-1)**q*np.conj(Bkq)  # B_{k,-q}
                Ckq = Pr.Cmatrix(k, q)
                Ckmq = Pr.Cmatrix(k, -q)
                # See Page 44, Eq 3.1 of the "Crystal Field Handbook"
                H = H + Bkq*Ckq + Bkmq*Ckmq



# Diagonalise the result and print out the energy 
# of the lowest 20 levels relative to the ground state.
(evals, evects) = np.linalg.eig(H)
E0 = np.min(evals)
calc_nrg_levels = np.sort(evals-E0)

for k in range(20):
    print(np.real(calc_nrg_levels[k]))
