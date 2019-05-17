from __future__ import division
from scipy.special import factorial
from scipy import  floor, sqrt
from numpy import arange

def Wigner3j(j1,j2,j3,m1,m2,m3):
#======================================================================
# Wigner3j.m by David Terr, Raytheon, 6-17-04
#
# Compute the Wigner 3j symbol using the Racah formula [1]. 
#
# Usage: 
# from wigner import Wigner3j
# wigner = Wigner3j(j1,j2,j3,m1,m2,m3)
#
#  / j1 j2 j3 \
#  |          |  
#  \ m1 m2 m3 /
#
# Reference: Wigner 3j-Symbol entry of Eric Weinstein's Mathworld: 
# http://mathworld.wolfram.com/Wigner3j-Symbol.html
#======================================================================

    # Error checking
    if ( ( 2*j1 != floor(2*j1) ) | ( 2*j2 != floor(2*j2) ) | ( 2*j3 != floor(2*j3) ) | ( 2*m1 != floor(2*m1) ) | ( 2*m2 != floor(2*m2) ) | ( 2*m3 != floor(2*m3) ) ):
        raise('All arguments must be integers or half-integers.')
        return -1

    # Additional check if the sum of the second row equals zero
    if ( m1+m2+m3 != 0 ):
        #print '3j-Symbol unphysical'
        return 0

    if ( j1 - m1 != floor ( j1 - m1 ) ):
        #print '2*j1 and 2*m1 must have the same parity'
        return 0
    
    if ( j2 - m2 != floor ( j2 - m2 ) ):
        #print '2*j2 and 2*m2 must have the same parity'
        return; 0

    if ( j3 - m3 != floor ( j3 - m3 ) ):
        #print '2*j3 and 2*m3 must have the same parity'
        return 0
    
    if ( j3 > j1 + j2)  | ( j3 < abs(j1 - j2) ):
        #print 'j3 is out of bounds.'
        return 0

    if abs(m1) > j1:
        #print 'm1 is out of bounds.'
        return 0

    if abs(m2) > j2:
        #print 'm2 is out of bounds.'
        return 0 

    if abs(m3) > j3:
        #print 'm3 is out of bounds.'
        return 0

    t1 = j2 - m1 - j3
    t2 = j1 + m2 - j3
    t3 = j1 + j2 - j3
    t4 = j1 - m1
    t5 = j2 + m2

    tmin = max( 0, max( t1, t2 ) )
    tmax = min( t3, min( t4, t5 ) )
    tvec = arange(tmin, tmax+1, 1)

    wigner = 0

    for t in tvec:
        wigner += (-1)**t / ( factorial(t) * factorial(t-t1) * factorial(t-t2) * factorial(t3-t) * factorial(t4-t) * factorial(t5-t) )

    return wigner * (-1)**(j1-j2-m3) * sqrt( factorial(j1+j2-j3) * factorial(j1-j2+j3) * factorial(-j1+j2+j3) / factorial(j1+j2+j3+1) * factorial(j1+m1) * factorial(j1-m1) * factorial(j2+m2) * factorial(j2-m2) * factorial(j3+m3) * factorial(j3-m3) )


def Wigner6j(j1,j2,j3,J1,J2,J3):
#======================================================================
# Calculating the Wigner6j-Symbols using the Racah-Formula                
# Author: Ulrich Krohn                                            
# Date: 13th November 2009
#                                                                         
# Based upon Wigner3j.m from David Terr, Raytheon                         
# Reference: http://mathworld.wolfram.com/Wigner6j-Symbol.html            
#
# Usage: 
# from wigner import Wigner6j
# WignerReturn = Wigner6j(j1,j2,j3,J1,J2,J3)
#
#  / j1 j2 j3 \
# <            >  
#  \ J1 J2 J3 /
#
#======================================================================

    # Check that the js and Js are only integer or half integer
    if ( ( 2*j1 != round(2*j1) ) | ( 2*j2 != round(2*j2) ) | ( 2*j2 != round(2*j2) ) | ( 2*J1 != round(2*J1) ) | ( 2*J2 != round(2*J2) ) | ( 2*J3 != round(2*J3) ) ):
        raise('All arguments must be integers or half-integers.')
        return -1
    
# Check if the 4 triads ( (j1 j2 j3), (j1 J2 J3), (J1 j2 J3), (J1 J2 j3) ) satisfy the triangular inequalities
    if ( ( abs(j1-j2) > j3 ) | ( j1+j2 < j3 ) | ( abs(j1-J2) > J3 ) | ( j1+J2 < J3 ) | ( abs(J1-j2) > J3 ) | ( J1+j2 < J3 ) | ( abs(J1-J2) > j3 ) | ( J1+J2 < j3 ) ):
        #print '6j-Symbol is not triangular!'
        return 0
    
    # Check if the sum of the elements of each traid is an integer
    if ( ( 2*(j1+j2+j3) != round(2*(j1+j2+j3)) ) | ( 2*(j1+J2+J3) != round(2*(j1+J2+J3)) ) | ( 2*(J1+j2+J3) != round(2*(J1+j2+J3)) ) | ( 2*(J1+J2+j3) != round(2*(J1+J2+j3)) ) ):
        #print '6j-Symbol is not triangular!'
        return 0
    
    # Arguments for the factorials
    t1 = j1+j2+j3
    t2 = j1+J2+J3
    t3 = J1+j2+J3
    t4 = J1+J2+j3
    t5 = j1+j2+J1+J2
    t6 = j2+j3+J2+J3
    t7 = j1+j3+J1+J3

    # Finding summation borders
    tmin = max(0, max(t1, max(t2, max(t3,t4))))
    tmax = min(t5, min(t6,t7))
    tvec = arange(tmin,tmax+1,1)
        
    # Calculation the sum part of the 6j-Symbol
    WignerReturn = 0
    for t in tvec:
        WignerReturn += (-1)**t*factorial(t+1)/( factorial(t-t1)*factorial(t-t2)*factorial(t-t3)*factorial(t-t4)*factorial(t5-t)*factorial(t6-t)*factorial(t7-t) )

    # Calculation of the 6j-Symbol
    return WignerReturn*sqrt( TriaCoeff(j1,j2,j3)*TriaCoeff(j1,J2,J3)*TriaCoeff(J1,j2,J3)*TriaCoeff(J1,J2,j3) )


def TriaCoeff(a,b,c):
    # Calculating the triangle coefficient
    return factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/(factorial(a+b+c+1))
