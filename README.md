## This program is based on IMPES method to solve the variation of pressure and water saturation with core length in one-dimensional two-phase water flooding process.
### 1) Model assumptions:
  1. Oil and water are not miscible with each other;
  2. Isothermal seepage, in line with Darcy's law of seepage;
  3. One-dimensional flow, regardless of gravity;
  4. Fluid and rock are incompressible, regardless of capillary force pcow=0;
  5. The porosity and permeability in the whole reservoir remain constant and isotropic;
### 2) Initial condition:
  The core is saturated with oil and bound water, and water is injected at the left end with a constant flow. 
  The water injection is a stable displacement, and the injection and output are Qv.
### 3) File description：
  1. IMPES_1D.m is the main program；
  2. solve_tridiagonal_matrix.m is a function for solving tridiagonal matrix equations using the catch-up method;
  3. Report.pdf is a Report based on questions and results.
  4. theory.pdf is the principle used for the entire programming, and also contains the derivation of the necessary formulas 
