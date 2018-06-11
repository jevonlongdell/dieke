import pdb
import dieke
import numpy as np
import qutip as q

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


            
#Pr = dieke.RareEarthIon(2)
import qutip as q


#consider two coupled spin L and S

Lval = 0.5
Sval = 0.5

[Lx,Ly,Lz] = map(lambda x: q.tensor(x, q.qeye(2*Sval+1)),q.jmat(Lval))
[Sx,Sy,Sz] = map(lambda x: q.tensor(q.qeye(2*Lval+1), x),q.jmat(Sval))

Jx = Lx + Sx
Jy = Ly + Sy
Jz = Lz + Sz

Lsq = Lx*Lx + Ly*Ly + Lz*Lz
Ssq = Sx*Sx + Sy*Sy + Sz*Sz
Jsq = Jx*Jx + Jy*Jy + Jz*Jz

g = q.simdiag([Lsq,Ssq,Jsq,Jz])
Jval = np.sqrt(4*g[0][2,:]+1)/2-0.5
mJval = g[0][3,:]

vects = g[1]

Lz3 = 0*Lz
Lz2 = 0*Lz
L1 = 0*Lz
for ii in range(len(vects)):
    L = Lval
    S = Sval
    J = Jval[ii]
    mJ = mJval[ii]
    for jj in range(len(vects)):
        Jp = Jval[jj]
        mJp = mJval[jj]
        matel1 = vects[ii].dag()*Lz*vects[jj]
        # Using
        sign = (-1)**(round(J-mJ))
        matel2 = sign * wignerlookup.w3j(round(2*J), 2, round(2*Jp), round(-2*mJ), 0, round(2*mJp)) * \
                                     dieke.reducedL(round(2*S), round(2*L), round(2*J), round(2*S), round(2*L), round(2*Jp))
        matel3 = sign * wignerlookup.w3j(round(2*J), 2, round(2*Jp), round(-2*mJ), 2, round(2*mJp)) * \
                                     dieke.reducedL(round(2*S), round(2*L), round(2*J), round(2*S), round(2*L), round(2*Jp))
        Lz2 = Lz2 + vects[ii]*matel2*vects[jj].dag()
        L1 = L1 + vects[ii]*matel3*vects[jj].dag()
        if np.abs(matel1[0,0]-matel2)>1e-6:
            print("S = %g,\nL = %g,\nJ = %g,\nJ' = %g,\nmJ = %g,\nmJ' = %g\n" % (S, L, J, Jp,mJ, mJp))
            print("W.E. gives %g, should be %g\n\n\n" % (matel2,matel1[0,0]) )
            raise Exception("oh shit")
