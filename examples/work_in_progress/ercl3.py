import numpy as np
import matplotlib.pyplot as plt
import dieke


# Calculates the energy levels of Pr:LaF3 as the strength of the
# crystal field terms are varied from zero to 150% or their actual
# values.


nf = 11  # 11 f-electrons means we're dealing with Er


# all the matricies we need
(LSterms, Uk, LSJlevels, freeion_mat,
 LSJmJstates, full_freeion_mat, Ckq) = dieke.makeMatricies(nf)
# # LSterms          - list of LSterm labels
# # Uk               - Uk in terms of these terms
# # LSJlevels        - list of LSJ levels
# # freeion_mat      - dictionary of free ion matricies in terms of those levels
# # LSJmJstates      - list of LSJmJ states
# # full_freeion_mat - the free ion matricies again but now in terms of LSJmJ
# # Ckq              - Ckq matricies



# # Read in a a set of crystal field parameters from Pr:LaF3
# # dieke reads these from (incomplete) carnall89params.xls
# cfparams = dieke.readLaF3params(nf)


# # Number of levels and states
# numLSJ = len(LSJlevels)
# numLSJmJ = len(LSJmJstates)
# # Variable multiplier for the crystal field strength
# # this will be the x-axis of the graph we will draw
# cfstrvals = np.linspace(0, 1.5, 15)


# # Make a free-ion Hamiltonian
# H0 = np.zeros([numLSJmJ, numLSJmJ])
# for k in cfparams.keys():
#     if k in full_freeion_mat:
#         print("Adding free ion parameter %s\n" % (k))
#         H0 = H0+cfparams[k]*full_freeion_mat[k]

# # Add in the crystal field terms and diagonalise the result
# H = H0
# for k in [2, 4, 6]:
#     for q in range(0, k+1):
#         if 'B%d%d' % (k, q) in cfparams:
#             if q == 0:
#                 H = H+cfparams['B%d%d' % (k, q)]*Ckq[(k, q)]
#             else:
#                 H = H+cfparams['B%d%d' % (k, q)]*Ckq[(k, q)]
#                 H = H+((-1)**q)*np.conj(cfparams['B%d%d' % (k, q)]) \
#                     * Ckq[(k, -q)]
# (evals, evects) = np.linalg.eig(H)
# E0 = np.min(evals)
# calc_nrg_levels = np.sort(evals-E0)


# # Make a matrix for the energy levels vs crystal field strength
# # (this is what we will plot)
# numLSJmJ = len(LSJmJstates)
# nrglevels = np.zeros([len(cfstrvals), numLSJmJ])


# # Work out the energy levels as the crystal field strengths are varied.
# level_index = 0
# cfparams2 = cfparams
# for cfstr in cfstrvals:
#     H = H0
#     for k in [2, 4, 6]:
#         for q in range(0, k+1):
#             if 'B%d%d' % (k, q) in cfparams:
#                 # print("adding in B %d %d"%(k, q))
#                 if q == 0:
#                     H = H+cfstr*cfparams['B%d%d' % (k, q)] * Ckq[(k, q)]
#                 else:
#                     H = H+cfstr*cfparams['B%d%d' % (k, q)] * Ckq[(k, q)]
#                     H = H+cfstr*((-1)**q) \
#                         * np.conj(cfparams['B%d%d' % (k, q)]) * Ckq[(k, -q)]
#     (evals, evects) = np.linalg.eig(H)
#     nrglevels[level_index, :] = np.sort(np.real(evals))-E0
#     level_index = level_index+1


# # Plot the results
# plt.ion()
# plt.axis([0, max(cfstrvals), -2000, 25000])
# plt.plot(cfstrvals, nrglevels, 'x-', [1, 1], [-10000, +50000])
# plt.draw()


# # Print out the 20 lowest enery levels
# for k in range(20):
#     print(calc_nrg_levels[k])
