==============================
Format of Crosswhite Data Files
==============================

Dieke gets many of it's matrix elements from Hannah Crosswhite's data files see "data" directory for the data files as well as some ancient Fortran code used to read them.

Reduced Tensor Operators
------------------------

In files fXnm.dat where X=2,3,..,7

I will illustrate with f2nm.dat (Pr)

- First line/block

  For Pr this is

  ::
   
  "   2   7   7   2   1   3   1   3   1   2"

  The first number is the number of f-electrons, or for n>7, the number of holes. In the example this is 2.

  The second number is the number of different LS terms. 

  The third number is the number of "J subspaces", these consist of a set of LSJ levels each with a fixed J.
  
  Then follows a list of numbers one for each of the J subspaces, with the dimensions of that subspace.

-  Second line/block
   
   This contains the state labels for each of the LS terms

   For Pr this is

   ::
   
   " 3P  3F  3H  1S  1D  1G  1I"

   For n=4 (Pm) this is

   ::
   
      5S  5D  5F  5G  5I  3P1 3P2 3P3 3D1 3D2 3F1 3F2 3F3 3F4 3G1 3G2 3G3 3H1 3H2 3H3
      3H4 3I1 3I2 3K1 3K2 3L  3M  1S1 1S2 1D1 1D2 1D3 1D4 1F  1G1 1G2 1G3 1G4 1H1 1H2
      1I1 1I2 1I3 1K  1L1 1L2 1N

      The third character is an index to distinguish between different terms  with the same L and S.

- Third line/block

  For Pr this is

  ::
   
  " 3 1 3 3 3 5 1 0 1 2 1 4 1 6"

  The thrid line contains two integers for each LS state

     * The first begin the multiplicity (2S+1)
     * The second being the L value

  The values can run together when these are greater than 9 as there are only two characters per number in the format.

- Fourth line/block

  In the fourth line/block the states are listed by possible J values

  for Pr this block is

  ::
   
  3P          1S
  3P
  3P          3F          1D
  3F
  3F          3H          1G
  3H
  3H          1I

  The first line lists the LS multiplets that can have a J value =0
  The second LS mulitplets that  have a J value =1
  The thrid  J=2 etc

  For the Kramers ion J vals this starts at 1/2 rather than 0

- Fifth line/block

  Block contains a whole bunch of lines with 8 numbers the first two are the state index (i,j) both i and j are in the range [1,numLSstates]
  which refer to which LS state the matrix elements are between

  Next the three numbers which are matrix elements of the Uk (U(2),U(4),U(6)) and the three elements which are the matrix elements of the V1k (V(12),V(14),V16). Specifically the first column is U(2)



Free ion
--------
In files fXnm.dat where X=2,3,..,12

The data files consist of a number of blocks separated by a blank line, each block corresponds to a "J subspace".

The first line of a block consists of two number
- The index of the J subspace
- The size of the J subspace (number of levels) 

The next lines list the states in the J subspace (one per line)
- The state is given like "3P" and is in the 6th col

The remaining lines list give matrix elemnts the lines are of the form:
(index of j subspace) (bra state) (ket state) (parameter index) (matrix element) (parameter name)


