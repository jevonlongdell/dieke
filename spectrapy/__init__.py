import numpy as np
from .wigner import Wigner6j, Wigner3j
from scipy.misc import factorial
from fractions import Fraction
import pandas
import os


# LSterms     - list of LSterm labels
# Uk          - Uk in terms of these terms
# LSJleves    - list of LSJ levels
# freeion_mat - dictionary of free ion matricies in terms of those levels
# LSJmJstates - list of LSJmJ states
# full_freeion_mat - the free ion matricies again but now in terms of LSJmJ
# Ckq         - Ckq matricies

def makeMatricies(nf):
    (LSJlevels, freeion_mat, LSterms, Uk, V) = read_crosswhite(nf)
    (LSJmJstates, full_freeion_mat) = makeFullFreeIonOperators(
                                              nf, LSJlevels, freeion_mat)
    Ckq = makeCkq(LSJmJstates, LSJlevels, LSterms, Uk)
    return (LSterms, Uk, LSJlevels, freeion_mat, LSJmJstates,
            full_freeion_mat, Ckq)


def readLaF3params(nf):
    print(__file__)
    pd = pandas.read_excel(os.path.join(__path__[0],'carnall89params.xls'),
                           skiprows=2).set_index('param')
    RareEarths = ['La', 'Ce', 'Pr', 'Nd', 'Pm',
                  'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho',
                  'Er', 'Tm', 'Yb']
    p = {}
    re = RareEarths[nf]
    for k in pd[re].keys():
        if not np.isnan(pd[re][k]):
            p[k] = pd[re][k]
    if 'M0' in p:
        p['M2'] = 0.56*p['M0']
        p['M4'] = 0.31*p['M0']
    if 'P2' in p:
        p['P4'] = 0.5*p['P2']
        p['P6'] = 0.1*p['P2']
    return p

#######################################################
# Functions to read state labels and return quantum numbers
#######################################################


def LfromTermLabel(termlabel):
    letters = 'SPDFGHIKLMNOQRT'
    L = letters.find(termlabel[3])
    return L


def SfromTermLabel(termlabel):
    mult = int(termlabel[2])
    return (mult-1)/2.0


def LfromLevelLabel(levellabel):
    letters = 'SPDFGHIKLMNOQRT'
    L = letters.find(levellabel[3])
    return L


def SfromLevelLabel(levellabel):
    mult = int(levellabel[2])
    return (mult-1)/2.0


def termFromLevelLabel(level):
    return level[0:4]


def levelFromStateLabel(state):
    return state[0:9]


def JfromLevelLabel(levellabel):
    return float(Fraction(levellabel[5:9]))


def LfromStateLabel(levellabel):
    letters = 'SPDFGHIKLMNOQRT'
    L = letters.find(levellabel[3])
    return L


def SfromStateLabel(levellabel):
    mult = int(levellabel[2])
    return (mult-1)/2.0


def JfromStateLabel(levellabel):
    if levellabel[7] == '/':
        return int(levellabel[5:7])/2.0
    else:
        return int(levellabel[5:7])


def mJfromStateLabel(levellabel):
    if levellabel[13] == '/':
        return int(levellabel[9:13])/2.0
    else:
        return int(levellabel[9:13])


#####################################
# Convert Ek parametres to F^k params
######################################

# <l||C^k||l'> eq 1.20 fro Guokui and Liu
def reducedCk(l, k, lprime):
    return (-1)**l*np.sqrt((2*l+1)*(2*lprime+1))*Wigner3j(
        l, k, lprime, 0, 0, 0)


# <TLSJ||Uk||t'L'S'J'>" eq 1.38 from Guokui and Liu
def makesinglyreducedUk(doublyReducedUk, LSterms, LSJlevels):
    LStermdict = {}
    for k in range(len(LSterms)):
        LStermdict[LSterms[k]] = k
    kvals = [2, 4, 6]  # values of k for which we need to worry about
    singlyreducedUk = np.zeros([3, len(LSJlevels), len(LSJlevels)])
    for i in range(len(LSJlevels)):
        L = LfromLevelLabel(LSJlevels[i])
        S = SfromLevelLabel(LSJlevels[i])
        J = JfromLevelLabel(LSJlevels[i])
        iterm = termFromLevelLabel(LSJlevels[i])
        for j in range(len(LSJlevels)):
            Lprime = LfromLevelLabel(LSJlevels[j])
            # Sprime = SfromLevelLabel(LSJlevels[j])
            Jprime = JfromLevelLabel(LSJlevels[j])
            jterm = termFromLevelLabel(LSJlevels[j])
            for k_idx in range(len(kvals)):
                k = kvals[k_idx]
                # Equation1.37 from Guokui and Liu
                singlyreducedUk[k_idx, i, j] = \
                    (-1)**(S+Lprime+J+k)*np.sqrt((2*J+1)*(2*Jprime+1)) * \
                    Wigner6j(J, Jprime, k, Lprime, L, S) * \
                    doublyReducedUk[k_idx, LStermdict[iterm],
                                    LStermdict[jterm]]
    return singlyreducedUk


class WignerDict:
    def __init__(self):
        self.wdict = {}

    def w3j(self, twicej1, twicej2, twicej3, twicem1, twicem2, twicem3):
        wargs = (twicej1, twicej2, twicej3, twicem1, twicem2, twicem3)
        if wargs in self.wdict:
                return self.wdict[wargs]
        else:
            w3jtemp = Wigner3j(twicej1/2.0, twicej2/2.0, twicej3/2.0,
                               twicem1/2.0, twicem2/2.0, twicem3/2.0)
            self.wdict[(wargs)] = w3jtemp
            return w3jtemp


# Equation 1.37 from Guokui and Liu
def makeCkq(LSJmJstates, LSJlevels, LSterms, doublyReducedUk):
    wignerlookup = WignerDict()
    numstates = len(LSJmJstates)
    leveldict = {}
    for k in range(len(LSJlevels)):
        leveldict[LSJlevels[k]] = k
    print("Making singly reduced Uk matricies")
    singlyreducedUk = makesinglyreducedUk(doublyReducedUk, LSterms, LSJlevels)
    multiplet_size = []
    multiplet_start = []

    count = 0
    for lvl in LSJlevels:
        twiceJ = int(JfromLevelLabel(lvl)*2+0.5)
        twicemJvals = range(-twiceJ, twiceJ+1, 2)
        # the +1 in line above is only to make sure mJ
        # goes between -J and J inclusive
        assert(len(twicemJvals) == twiceJ+1)
        multiplet_start.append(count)
        multiplet_size.append(twiceJ+1)
        for twicemJ in twicemJvals:
            if (twiceJ % 2) == 0:
                assert(LSJmJstates[count] == '%s %3d  ' % (lvl, twicemJ//2))
            else:
                assert(LSJmJstates[count] == '%s %3d/2' % (lvl, twicemJ))
            count = count+1

    Ckq = {}
    for k in [2, 4, 6]:
        lCkl = reducedCk(3, k, 3)
        for q in range(-k, k+1):
            print("Making C%d%d matrix." % (k, q))
            Ckq[(k, q)] = np.zeros([numstates, numstates])

            for i in range(len(LSJlevels)):
                istart = multiplet_start[i]
                isize = multiplet_size[i]
                # istop = istart+isize # commented this out because never used?
                for j in range(len(LSJlevels)):
                    if abs(singlyreducedUk[k//2-1, i, j]) < 1e-10:
                        continue
                    jstart = multiplet_start[j]
                    jsize = multiplet_size[j]
                    # jstop = jstart+jsize #commented out because never used?
                    twiceJ = isize-1
                    J = twiceJ/2.0
                    twiceJprime = jsize-1
                    # Jprime = twiceJprime/2.0 #never used?
                    for ii in range(isize):  # ii = inner i
                        twicemJ = -twiceJ+2*ii
                        mJ = -J + ii
                        for ij in range(jsize):
                            twicemJprime = -twiceJprime + 2*ij
                            # mJprime=-Jprime+ij#commented because never used?
                            threejtemp = wignerlookup.w3j(twiceJ, 2*k,
                                                          twiceJprime,
                                                          -twicemJ, 2*q,
                                                          twicemJprime)
                            if(threejtemp != 0):
                                Ckq[(k, q)][istart+ii, jstart+ij] = \
                                    (-1)**(J-mJ)*threejtemp * \
                                    singlyreducedUk[k//2-1, i, j]*lCkl
    return Ckq


# takes free ion operators defined over LSJ levels
# and epands them to those defined in terms of LSJmJ states
def makeFullFreeIonOperators(nf, LSJlevels, fi_mat):
    numstates = factorial(14)//(factorial(nf)*factorial(14-nf))
    numstates = int(numstates+0.5)
    full_fi_mat = {}

    for key in fi_mat.keys():
        full_fi_mat[key] = np.zeros([numstates, numstates])
    multiplet_size = []
    multiplet_start = []

    LSJmJstates = []
    count = 0
    for lvl in LSJlevels:
        twiceJ = int(JfromLevelLabel(lvl)*2+0.5)
        twicemJvals = range(-twiceJ, twiceJ+1, 2)
        # the in the line above is +1 is only to make
        # sure mJ goes between -J and J inclusive
        assert(len(twicemJvals) == twiceJ+1)
        multiplet_start.append(count)
        count = count+twiceJ+1
        multiplet_size.append(twiceJ+1)

        for twicemJ in twicemJvals:
            if (twiceJ % 2) == 0:
                LSJmJstates.append('%s %3d  ' % (lvl, twicemJ/2))
            else:
                LSJmJstates.append('%s %3d/2' % (lvl, twicemJ))

    for i in range(len(LSJlevels)):
        istart = multiplet_start[i]
        isize = multiplet_size[i]
        istop = istart+isize
        for j in range(len(LSJlevels)):
            if isize != multiplet_size[j]:
                continue
            jstart = multiplet_start[j]
            jstop = jstart+isize
            for key in full_fi_mat.keys():
                full_fi_mat[key][istart:istop, jstart:jstop] = \
                                np.eye(isize)*fi_mat[key][i, j]
    return (LSJmJstates, full_fi_mat)


def read_crosswhite(nf):
    # Use 14-nf for nf>7
    reduced_tensor_file = 'data/f%dnm.dat' % (7-abs(7-nf))
    f = open(reduced_tensor_file, 'r')

    # Read first line
    line = f.readline().split()
    line = list(map(int, line))
    assert(7-abs(7-nf) == line[0])
    numLS = line[1]  # number of LS states
    nJsub = line[2]  # number of different J subspaces
    ndim = line[3:]  # dimensions for each of the J subspaces
    assert(len(ndim) == nJsub)

    # Read second line
    LSterms = f.readline().split()
    while(len(LSterms) != numLS):
        line = f.readline()
        assert(line != '')
        line = line.split()
        map(LSterms.append, line)

    # Read third line
    x = list(map(int, f.readline().split()))
    while(len(x) < 2*numLS):
        line = f.readline()
        assert(line != '')
        for k in map(int, line.split()):
            if k > 10:  # the numbers have run together
                x.append(k//100)
                x.append(k % 100)
            else:
                x.append(k)
    mult = x[0::2]  # 2S+1
    Lvalue = x[1::2]  # L

    # check state labels
    assert(len(mult) == len(LSterms))
    assert(len(Lvalue) == len(LSterms))
    for k in range(len(mult)):
        letters = 'SPDFGHIKLMNOQRT'
        assert(mult[k] == int(LSterms[k][0]))
        assert(letters[Lvalue[k]] == LSterms[k][1])

    # Rewrite state labels so that they are of the form
    # "2 1D" rather than "1D2"
    for k in range(len(LSterms)):
        if len(LSterms[k]) == 2:
            LSterms[k] = "1 %s" % (LSterms[k])
        else:
            LSterms[k] = "%s %s" % (LSterms[k][2], LSterms[k][0:2])

    # read next block
    # a dictionary containing all the LS states that have
    # a particular value for  2J
    statesby2J = {}

    Jmin = (nf % 2)
    for k in range(nJsub):
        statesby2J[2*k+Jmin] = f.readline().split()
        # this bit is beause some of the lines wrap
        while (len(statesby2J[2*k+Jmin]) != ndim[k]):
            line = f.readline()
            assert(line != '')
            list(map(statesby2J[2*k+Jmin].append, line.split()))

    Uk = np.zeros([3, numLS, numLS])
    V = np.zeros([3, numLS, numLS])

    lines = f.readlines()
    for line in lines:
        entries = line.split()
        (i, j) = list(map(int, entries[0:2]))
        Uk[:, i-1, j-1] = list(map(float, entries[2:5]))
        V[:, i-1, j-1] = list(map(float, entries[5:8]))
    f.close()

    ##############################
    # read in free ion matricies
    #############################
    freeionfilename = 'data/f%dmp.dat' % (nf,)
    f = open(freeionfilename, 'r')

    fi_mat = {}  # a dictionary to hold our free ion matricies
    # the key will be the name eg "F2" or "ZETA"

    # a dictionary mapping parameter number to parameter name
    parameters = {}

    numLSJ = np.sum(ndim)
    LSJlevels = []

    count = 0

    # loop over the blocks in the free ion file

    for jnumber in range(nJsub):

        line = f.readline()
        assert(line.strip() != '')
        (jsubspace_idx, jsubspace_size) = list(map(int, line.split()))
        assert(jsubspace_idx == jnumber+1)
        assert(jsubspace_size == ndim[jnumber])
        # get all the parameters for the
        # state_collection = [] #never used?
        for k in range(jsubspace_size):  # loop over all the XYJ states
            # make sure that they are what we are expecting
            line = f.readline().split()
            state = line[5]
            assert(state == statesby2J[2*jnumber+Jmin][k])
            if len(state) == 2:
                temp = "1 %s" % (state)
            else:
                temp = "%s %s" % (state[2], state[0:2])
            if Jmin == 0:
                LSJlevels.append('%s %2d  ' % (temp, jnumber))
            else:
                LSJlevels.append('%s %2d/2' % (temp, 2*jnumber+1))
            assert(line[6] == 'F')
            if nf < 7:
                assert(int(line[7]) == nf)
            else:
                assert(int(line[7]) == 14-nf)

        while(True):
            line = f.readline()
            print(line)
            if line == []:
                break
            if len(line.strip()) == 0:
                break
            # grrr the formating is such that the numbers aren't
            # always separated by free space
            (jnum, i, j, pnum) = list(map(int, line[0:24].split()[0:4]))
            line = line[24:].split()
            print(line)
            mat_element = float(line[0])
            param = line[1]
            if not(param in fi_mat):
                fi_mat[param] = np.zeros([numLSJ, numLSJ])
                parameters[pnum] = param
            assert(parameters[pnum] == param)
            II = i+count-1
            JJ = j+count-1
            fi_mat[param][II, JJ] = mat_element
        count = count+jsubspace_size

        # fill in the other triangle
    for key in fi_mat.keys():
        fi_mat[key] = fi_mat[key]+np.transpose(fi_mat[key]) - \
                      np.diag(np.diag(fi_mat[key]))
    # Who knows if this is the right way round?
    fi_mat['ALPHA'] = 100*fi_mat['.01ALPH']
    return (LSJlevels, fi_mat, LSterms, Uk, V)
