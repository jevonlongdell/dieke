import numpy as np
import matplotlib.pyplot as plt
import dieke

# does example 3.5.2 of mikes notes

nf = 2  # the number of f electrions

# This object contains the matricies we need for the
# calculations all in the dictionary "FreeIonMatrix"
Pr = dieke.IsotropicRareEarthIon(nf)


# Read the crystal field parameters
#cfparams = dieke.readLaF3params(nf)

cfparams ={'nf': 2.0,
 'F2': 68878.0,
 'F4': 50347.0,
 'F6': 32901.0}


numLSJ = Pr.numlevels()

for operator in ['F2','F4','F6','ZETA']:
    mat = Pr.FreeIonMatrix[operator]
    print(operator)
    for i in range(numLSJ):
        for j in range(i, numLSJ):
            if np.abs(mat[i, j])>1e-6:
                print("<%s | %s | %s > = %f" % (Pr.LSJlevelLabels[i],operator,
                                           Pr.LSJlevelLabels[i],
                                           mat[i, j]))


Ce = dieke.RareEarthIon(1)
mat = Ce.Cmatrix(2,0)
for i in range(Ce.numstates()):
    for j in range(Ce.numstates()):
        if np.abs(mat[i,j])>1e-6:
            print("<%s | C20 | %s > = %f" % (Ce.LSJmJstateLabels[i],
                                           Ce.LSJmJstateLabels[j],
                                           mat[i, j]))


                
# # Make an empy matrix to put the results in as well as an
# # empty matrix for the Hamiltonian
# numLSJ = Pr.numlevels()

H0 = np.zeros([numLSJ, numLSJ])

# Make the Hamiltonian, diagonalise it and work out the
# energy of the ground state
for k in cfparams.keys():
    if k in Pr.FreeIonMatrix:
        print("using parameter ", k)
        H0 = H0+cfparams[k]*Pr.FreeIonMatrix[k]
(evals, evects) = np.linalg.eig(H0)
E0 = np.min(evals)

print(np.sort(evals))



# # Loop over each of our zetavals, and calulate the Hamiltonian
# # for that value of zeta, diagonalise and store energy levels.
# cfparams2 = cfparams
# for i, zeta in enumerate(zetavals):
#     cfparams2['ZETA'] = zeta
#     H = np.zeros([numLSJ, numLSJ])
#     for k in cfparams2.keys():
#         if k in Pr.FreeIonMatrix:
#             H = H + cfparams2[k]*Pr.FreeIonMatrix[k]
#     (evals, evects) = np.linalg.eig(H)
#     nrglevels[i, :] = np.sort(evals)-E0

# plt.ion()  # interactive plotting, so plt.show() doesn't block
# plt.axis([0, max(zetavals), -2000, 25000])
# plt.plot(zetavals, nrglevels, [zeta0, zeta0], [-10000, +50000])
# plt.show()

