# z2inv

  This program calculates the parity of each band at the TR kpts
    and determines the Z2 index using the parity method
    using PRB 76, 045302 (2007)
  Requirement: 
    1. SCF calculation with SOC
    2. NSCF calculation with SOC using 8 TR kpts explicitly
    3. Both above MUST use wf_collect=.true.
    4. cd into the .save directory and run this program with a
       single parameter of the number of filled bands
  Possible caveats:
    1. No checks (assuming you follow the previous 4 steps strictly)
    2. For USPP & PAWs, one indeed needs <psi_i|P S|psi_i> where S
       is the projection S=1+Q_ij |beta_i><beta_j|, but since this is
       nothing but a parity calculation, so I did not bother.

  Author: Chao Cao @ Hangzhou Normal University (2015) (ccao1981@gmail.com)

