import numpy as np
import matplotlib.pyplot as plt
import dieke

# Example of free ion calculation using dieke
#
# This makes Figure 3.1 in Mike Reids notes which shows
# the energy levesl of Pr3+ as the strength of the spin
# orbit interaction is increased from zero to about 150%
# of it's actual value.
#
# www2.phys.canterbury.ac.nz/~mfr24/electronicstructure/00electronic.pdf


nf = 2  # the number of f electrions

# This object contains the matricies we need for the
# calculations all in the dictionary "FreeIonMatrix"
Pr = dieke.IsotropicRareEarthIon(nf)


# Read the crystal field parameters
cfparams = dieke.readLaF3params(nf)

# Get the spin orbit coupling parameter
zeta0 = cfparams['ZETA']

# Make a list of spin orbit parameters for the plot
zetavals = np.linspace(0, 1.5*zeta0, 100)

# Make an empy matrix to put the results in as well as an
# empty matrix for the Hamiltonian
numLSJ = Pr.numlevels()

nrglevels = np.zeros([len(zetavals), numLSJ])
H0 = np.zeros([numLSJ, numLSJ])

# Make the Hamiltonian, diagonalise it and work out the
# energy of the ground state
for k in cfparams.keys():
    print("using parameter ", k)
    if k in Pr.FreeIonMatrix:
        H0 = H0+cfparams[k]*Pr.FreeIonMatrix[k]

(evals, evects) = np.linalg.eig(H0)
E0 = np.min(evals)


# Loop over each of our zetavals, and calulate the Hamiltonian
# for that value of zeta, diagonalise and store energy levels.
cfparams2 = cfparams
for i, zeta in enumerate(zetavals):
    cfparams2['ZETA'] = zeta
    H = np.zeros([numLSJ, numLSJ])
    for k in cfparams2.keys():
        if k in Pr.FreeIonMatrix:
            H = H + cfparams2[k]*Pr.FreeIonMatrix[k]
    (evals, evects) = np.linalg.eig(H)
    nrglevels[i, :] = np.sort(evals)-E0

plt.ion()  # interactive plotting, so plt.show() doesn't block
plt.axis([0, max(zetavals), -2000, 25000])
plt.plot(zetavals, nrglevels, [zeta0, zeta0], [-10000, +50000])
plt.show()
