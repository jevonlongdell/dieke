#!/usr/bin/env python
# Filename = njsymbols.py

# Copyright (C) 2013 Sebastian Horvath (sebastian.horvath@gmail.com)
# 
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
# 
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

from __future__ import division
import numpy as np
# from scipy.misc import factorial
from scipy.special import factorial



from math import fsum


def tricon_ck(a, b, c):
    r"""
    Triangular condition check; returns True if the triangular condition on the
    three integers or half-integers a, b and c is satisfied.
    """
    return(a + b >= c and c >= np.abs(a - b))


def wigner_3j(j1, j2, j3, m1, m2, m3):
    r"""
    Calculate the Wigner 3j symbol, give in terms of the Clebsch-Gordon
    coefficents as

    .. math::

        \begin{pmatrix}
        j_1 & j_2 & j_3 \\
        m_1 & m_2 & m_3
        \end{pmatrix} = \frac{(-1)^{(j_1 - j_2 - m_3)}}{\sqrt{2j_3 +1}} 
        \langle j_1 m_1 j_2 m_2 | j_3 - m_3 \rangle.

    Parameters
    ----------
    j1 : integer or half-integer 
        The value of `j_1`.
    j2 : integer or half-integer 
        The value of `j_2`.
    j3 : integer or half-integer 
        The value of `j_3`.
    m1 : integer or half-integer 
        The value of `m_1`.
    m2 : integer or half-integer 
        The value of `m_2`.
    m3 : integer or half-integer 
        The value of `m_3`.

    Returns
    -------
    result : float
        The numerical value of the 3j symbol.

    Notes
    -----
    Uses the recursive algorithm for Clebsch-Gordon coefficients from page 44,
    Edmonds - Angular Momentum in Quantum Mechanics.
    """

    if tricon_ck(j1, j2, j3) and (m1 + m2 + m3 == 0):
        phase = (-1)**(j1 - j2 - m3)
        pre_fact = np.sqrt((factorial(j1 + j2 - j3) * factorial(j1 - m1) * 
            factorial(j2 - m2) * factorial(j3 + m3) * factorial(j3 -
            m3))/(factorial(j1 + j2 + j3 + 1) * factorial(j3 + j1 - j2) *
            factorial(j3 + j2 - j1) * factorial(j1 + m1) * factorial(j2 + m2)))
        xmin = max(0, j3 - j2 - m1)
        xmax = min(j3 + m3, j1 - m1)
        xlist = np.arange(xmin, xmax + 1)
        #        try:
        result = phase * pre_fact * fsum([((-1)**(j1 - m1 - x) * (factorial(j1 + m1 + x) * factorial(j3 + j2 - m1 - x)) /(factorial(x) * factorial(j3 + m3 - x) * factorial(j1 - m1 - x) * factorial(j2 - j3 + m1 + x)))  for x in xlist])
#        except:
#            print(xlist, j1, j2,j3, m1,m2,m3)
#            import pdb; pdb.set_trace()
    else:
        result = 0

    
    return(result)


def wigner_6j(a, b, c, d, e, f):
    r"""
    Calculate the Wigner 6j symbol

    .. math::

        \begin{Bmatrix}
        j_1 & j_2 & j_3 \\
        l_1 & l_2 & l_3
        \end{Bmatrix}.

    Parameters
    ----------
    a : integer or half-integer
        The value of `j_1`.
    b : integer or half-integer
        The value of `j_2`.
    c : integer or half-integer
        The value of `j_3`.
    d : integer or half-integer
        The value of `l_1`.
    e : integer or half-integer
        The value of `l_2`.
    f : integer or half-integer
        The value of `l_3`.

    Returns
    -------
    result : float
        The numerical value of the 6j symbol.

    Notes
    -----
    Uses the algorithm on page 99 of Edmonds - Angular Momentum in Quantum
    Mechanics.
    """

    def triad(a, b, c):
        """
        Evaluate the triangular portion of the 6j symbol formula.
        """
        r = np.sqrt(factorial(a + b - c) * factorial(a - b + c) * factorial(b +
            c - a)/factorial(a + b + c + 1))
        return(r)

    if tricon_ck(a, b, c) and tricon_ck(d, b, f) and tricon_ck(d, e, c) and \
            tricon_ck(a, e, f):
        pre_fact = triad(a, b, c) * triad(a, e, f) * triad(d, b, f) * triad(d, e, c)
        xmin = max(a + b + c, a + e + f, d + b + f, d + e + c)
        xmax = min(a + b + d + e, b + c + e + f, c + a + f + d)
        xlist = np.arange(xmin, xmax + 1)
        result = pre_fact * fsum([((-1)**x * factorial(x + 1)/(factorial(x -
            a - b -c) * factorial(x - a - e - f) * factorial(x - d - b - f) *
            factorial(x - d - e - c) * factorial(a + b + d + e - x) *
            factorial(b + c + e + f - x) * factorial(c + a + f + d - x))) for x
            in xlist])
    else:
        result = 0

    return(result)


def wigner_9j(a, b, c, d, e, f, g, h, i):
    r"""
    Calculate the Wigner 9j symbol

    .. math::

        \begin{Bmatrix}
        j_{11} & j_{12} & j_{13} \\
        j_{21} & j_{22} & j_{23} \\
        j_{31} & j_{32} & j_{33}
        \end{Bmatrix}.

    Parameters
    ----------
    a : integer or half-integer
        The value of `j_{11}`.
    b : integer or half-integer
        The value of `j_{12}`.
    c : integer or half-integer
        The value of `j_{13}`.
    d : integer or half-integer
        The value of `j_{21}`.
    e : integer or half-integer
        The value of `j_{22}`.
    f : integer or half-integer
        The value of `j_{23}`.
    g : integer or half-integer
        The value of `j_{31}`.
    h : integer or half-integer
        The value of `j_{32}`.
    i : integer or half-integer
        The value of `j_{33}`.

    Returns
    -------
    result : float
        The numerical value of the 9j symbol.

    Notes
    -----
    Uses the definition on page 101 in terms of 6j symbols of Edmonds - Angular
    Momentum in Quantum Mechanics.
    """
    if tricon_ck(a, d, g) and tricon_ck(h, i, g) and tricon_ck(b, e, h) and \
            tricon_ck(d, e, f) and tricon_ck(c, f, i) and tricon_ck(c, a, b):
        xmax = min(a + i, h + d, b + f)
        xmin = max(abs(a - i), abs(h - d), abs(b - f))
        xlist = np.arange(xmin, xmax + 1)
        result = fsum([((complex(-1)**(2 * x)) * (2 * x + 1) * wigner_6j(a, d, g, h,
            i, x) * wigner_6j(b, e, h, d, x, f) * wigner_6j(c, f, i, x, a, b))
            for x in xlist])
    else:
        result = 0
    return(result)
