import pdb
import dieke
import numpy as np

L = 1
S = 0
J = 1
mJ = list(range(-J,J+1))
L0 = np.zeros((2*J+1, 2*J+1))
L1 = np.zeros((2*J+1, 2*J+1))

wignerlookup = dieke.WignerDict()

for ii in range(2*J+1):
    twiceL = 2*L
    twiceJ = 2*J
    twiceS = 0
    twicemJ = 2*mJ[ii]
    for jj in range(2*J+1):
        twiceLp = 2*L
        twiceJp = 2*J
        twiceSp = 0
        twicemJp = 2*mJ[jj]
        if twicemJ==twicemJp:
            sign = (-1)**((twiceJ-twicemJ)/2.0)
            L0[ii, jj] = sign * wignerlookup.w3j(twiceJ, 2, twiceJp, -twicemJ, 0, twicemJp) * \
                         dieke.reducedL(twiceS, twiceL, twiceJ, twiceSp, twiceLp, twiceJp)
            L1[ii, jj] = sign * wignerlookup.w3j(twiceJ, 2, twiceJp, -twicemJ, 2, twicemJp) * \
                         dieke.reducedL(twiceS, twiceL, twiceJ, twiceSp, twiceLp, twiceJp)

         
L = 0
S = 1
J = 1
mJ = list(range(-J,J+1))
S0 = np.zeros((2*J+1, 2*J+1))
S1 = np.zeros((2*J+1, 2*J+1))

for ii in range(2*J+1):
    twiceL = 2*L
    twiceJ = 2*J
    twiceS = 2*S
    twicemJ = 2*mJ[ii]
    for jj in range(2*J+1):
        twiceLp = 2*L
        twiceJp = 2*J
        twiceSp = 2*S
        twicemJp = 2*mJ[jj]
        if twicemJ==twicemJp:
            sign = (-1)**((twiceJ-twicemJ)/2.0)
            S0[ii, jj] = sign * wignerlookup.w3j(twiceJ, 2, twiceJp, -twicemJ, 0, twicemJp) * \
                         dieke.reducedS(twiceS, twiceL, twiceJ, twiceSp, twiceLp, twiceJp)
            S1[ii, jj] = sign * wignerlookup.w3j(twiceJ, 2, twiceJp, -twicemJ, 2, twicemJp) * \
                         dieke.reducedS(twiceS, twiceL, twiceJ, twiceSp, twiceLp, twiceJp)


            
Pr = dieke.RareEarthIon(2)
