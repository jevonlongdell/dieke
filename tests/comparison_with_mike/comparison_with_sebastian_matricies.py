#!/usr/bin/python2

# note: python3 pickle doesn't seem to like python2 pickles

import pickle
import gzip
import dieke
import numpy as np




try:
    FileNotFoundError
except NameError:
    #py2
    FileNotFoundError = IOError


EPS=1e-5


#get sebastians matricies
sebmats = pickle.load(gzip.open('f11cf_matel.p.gz', 'rb'))

nf = 11  # 2 f-electrons means we're dealing with Er

# This object contains the matricies we need for the
# calculations all in the dictionary "FreeIonMatrix"
# or via the function "Cmatrix"



#cache my version of the matrix elements to make stuff faster
try:
    f = gzip.open('jevmats.dat.gz', 'rb')
    print("Loading Jevon's matricies")
    jevmats = pickle.load(f)
except FileNotFoundError:
    print("Making Jevon's matricies")
    Er = dieke.RareEarthIon(nf)
    jevmats = Er.FreeIonMatrix
    # Add the crystal field matricies to my dict
    for k in [2, 4, 6]:
        for q in range(k+1):
            jevmats['C%d%d' % (k, q)] = Er.Cmatrix(k, q)   
    pickle.dump(jevmats, gzip.open('jevmats.dat.gz', 'wb'))



#scale cf matricies

# jevmats['C20'] = -jevmats['C20']
# jevmats['C40'] = -jevmats['C40']
# jevmats['C60'] = -jevmats['C60']
# jevmats['C21'] = -np.sqrt(2)*jevmats['C21']


matnames = set(jevmats.keys())
matnames.update(sebmats.keys())
matnames = list(matnames)
matnames.sort()

for m in matnames:
    if not(m in jevmats):
        print("Jev missing %s"%(m))
    if not(m in sebmats):
        print("Seb missing %s"%(m))
    if (m in jevmats) and (m in sebmats):
        ms = np.matrix(sebmats[m])
        mj = np.matrix(jevmats[m])
        assert(np.linalg.norm(ms-ms.H) < EPS)
        assert(np.linalg.norm(mj-mj.H) < EPS)
        normofdiff = np.linalg.norm(ms-mj)
        if (normofdiff > EPS):
            # print("Matricies differ for %s, norm of diff = %g" % (m,
            #                                                      normofdiff))
            (evalss, _) = np.linalg.eig(ms)
            (evalsj, _) = np.linalg.eig(mj)
            #make sure evals are real
            assert(np.linalg.norm(np.imag(evalss)) < EPS)
            assert(np.linalg.norm(np.imag(evalsj)) < EPS)
            evalss = np.real(evalss)
            evalsj = np.real(evalsj)
            evalss.sort()
            evalsj.sort()
            normofdiff = np.linalg.norm(evalss-evalsj)
            if (normofdiff > EPS):
                print("Eigenvalues differ for %s norm of diff = %g" % (m,
                                                                       normofdiff))
            else:
                print("Eigenvalues agree for %s, norm of diff  < %g" %(m,
                                                                      EPS))
        else:
            print("!!!! Matricies agree for %s, norm of diff = %g !!!!" % (m,
                                                                 normofdiff))


            
# Read in a a set of crystal field parameters from Er:LaF3
# dieke reads these from (incomplete) carnall89params.xls
cfparams = dieke.readLaF3params(nf)



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
