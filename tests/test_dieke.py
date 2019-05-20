#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `dieke` package."""

import pytest
# from multiprocessing import Pool
import numpy as np
import dieke
from numpy.linalg import norm


def test_LS_commutation_relations():
  
    eps = 1e-7
    
    Re = dieke.RareEarthIon(2)  # Pr

    L = Re.FreeIonMatrix['L']
    Lx = Re.FreeIonMatrix['Lx']
    Ly = Re.FreeIonMatrix['Ly']
    Lz = Re.FreeIonMatrix['Lz']
    S = Re.FreeIonMatrix['S']
    Sx = Re.FreeIonMatrix['Sx']
    Sy = Re.FreeIonMatrix['Sy']
    Sz = Re.FreeIonMatrix['Sz']
    
    J = Re.FreeIonMatrix['J']
    Jx = Lx+Sx
    Jy = Ly+Sy
    Jz = Lz+Sz
    
    Lvec = [Lx,Ly,Lz]
    Svec = [Sx,Sy,Sz]

    for k in range(3):
        for j in range(3):
            assert(norm(Lvec[k]*Svec[j]-Svec[j]*Lvec[k])<eps)
        
    assert(norm(Lx*Ly-Ly*Lx-1j*Lz)<eps)
    assert(norm(Sx*Sy-Sy*Sx-1j*Sz)<eps)
    assert(norm(Jx*Jy-Jy*Jx-1j*Jz)<eps)

    assert(norm(Lx*Lx + Ly*Ly + Lz* Lz - L*(L+np.eye(len(L)))) < eps)
    assert(norm(Sx*Sx + Sy*Sy + Sz* Sz - S*(S+np.eye(len(S)))) < eps)
    assert(norm(Jx*Jx + Jy*Jy + Jz* Jz - J*(J+np.eye(len(J)))) < eps)



def test_Pr_LaF3():
    """
    Tests that we give the same results as Carnall 1989 for Pr:LaF3
    """

    #energy levels from Carnall's paper
    carnall_levels = [0.2, 71, 95, 138, 183, 221, 333, 444, 463,
                      2126, 2158, 2191, 2284, 2290, 2295, 2318,
                      2399, 2412,2438, 2540, 4179,4200, 4283,
                      4321, 4384, 4467, 4478, 4496]
    carnall_levels = np.array(carnall_levels)
    
    nf = 2  # 2 f-electrons means we're dealing with Pr
    Pr = dieke.RareEarthIon(nf)
    cfparams = dieke.readLaF3params(nf)
    
    # Number of levels and states
    numLSJ = Pr.numlevels()
    numLSJmJ = Pr.numstates()

    

    # Make a free-ion Hamiltonian
    H0 = np.zeros([numLSJmJ, numLSJmJ])
    for k in cfparams.keys():
        if k in Pr.FreeIonMatrix:
            H0 = H0+cfparams[k]*Pr.FreeIonMatrix[k]

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
                    # See page 44, eq 3.1 of the crystal field handbook
                    H = H + Bkq*Ckq + Bkmq*Ckmq



    # Diagonalise the result
    (evals, evects) = np.linalg.eig(H)
    E0 = np.min(evals)
    calc_nrg_levels = np.sort(evals-E0)

    # Assert that they all agree within one wave number
    assert(np.max(np.abs(
        calc_nrg_levels[0:len(carnall_levels)]-carnall_levels))
           < 1.0 )
           
    


def test_ErYSO():
    """
     Tests that we give the same as Sebastians calcs for ErYSO
     (arXiv:1809.01058)
    """
    
    #energy levels from Paper 
    seb_levels = [15, 47, 75, 130, 199, 388, 462, 508,
                  6522, 6560, 6583, 6640, 6777, 6833, 6867,
                  10206, 10236, 10267, 10339, 10381, 10398]
    seb_levels = np.array(seb_levels)
    

    nf = 11  # 11 f-electrons means we're dealing with Er

    Er = dieke.RareEarthIon(nf)


    
    # Number of levels and states
    numLSJ = Er.numlevels()
    numLSJmJ = Er.numstates()

    cfparams = dieke.readLaF3params(nf)
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



    # Make a free-ion Hamiltonian
    H0 = np.zeros([numLSJmJ, numLSJmJ])
    for k in sorted(cfparams.keys()):
        if k in Er.FreeIonMatrix:
            print("Adding free ion parameter \'%s\' = %g" % (k, cfparams[k]))
            H0 = H0+cfparams[k]*Er.FreeIonMatrix[k]

    # Add in the crystal field terms and diagonalise the result

    H = H0
    for k in [2, 4, 6]:
        for q in range(0, k+1):
            if 'B%d%d' % (k, q) in cfparams:
                if q == 0:
                    H = H+cfparams['B%d%d' % (k, q)]*Er.Cmatrix(k, q)
                else:
                    Bkq = cfparams['B%d%d' % (k, q)]
                    Bkmq = (-1)**q*np.conj(Bkq)
                    Ckq = Er.Cmatrix(k, q)
                    Ckmq = Er.Cmatrix(k, -q)
                    #See page 44, eq 3.1 of the crystal field handbook
                    H = H + Bkq*Ckq + Bkmq*Ckmq


    (evals, evects) = np.linalg.eig(H)
    E0 = np.min(evals)
    calc_nrg_levels = np.sort(evals-E0)
    calc_nrg_levels = calc_nrg_levels[::2]  # ignore every second element
    E0seb = np.min(seb_levels)
    seb_levels = seb_levels-E0seb

    print('    Jevon  Sebastian Difference')
    for k in range(len(seb_levels)):
        print("%9.1f %9.1f %6.1f"%(np.real(calc_nrg_levels[k]),seb_levels[k],np.real(calc_nrg_levels[k]) - seb_levels[k]))
    

    
    # This isn't much of a test at the moment because there seems to be
    # an issue with the free ion parameters
    
    assert(np.max(np.abs(
        calc_nrg_levels[0:len(seb_levels)]-seb_levels))
           < 26.0 )

    # This is more of a test because it it only looks at the ground multiplet
    assert(np.max(np.abs(
        calc_nrg_levels[0:8]-seb_levels[0:8]))
           < 0.6 )


    
