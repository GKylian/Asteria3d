#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

/* Defines all the parameters set pre-compiling
<->
---> MHD: Is the gas MHD (B != 0) ?
---> CYLINDRICAL: Are we using cylindrical coordinates ?
---> HLL/HLLC/HLLD: What Riemann are we using ? (HLLD default for MHD, HLLC default for HD)
---> TRACING: Do we perform characteristic tracing in the reconstruction ?
---> NGHOST: The number of ghost cells added to N in each direction
---> NVAL: The numbe of variables in u (MHD? Iso?)
---> NWAVES: The number of characteristic waves (and the number of eigenvectors/eigenvalues)
<->
*/

#define MHD
//#define CYLINDRICAL
//#define ISO
//#define HLLD				/* Use the HLLD Riemann Solver ? */
//#define HLLC				/* Use the HLLC Riemann Solver ? */
#define HLL					/* Use the HLL Riemann Solver ? */
#define TRACING				/* Do the characteristic tracing ? */
#define PRIMBOUNDS          /* Apply boundary conditions to primitive variables */
#ifdef MHD
  #define NGHOST 3
#else
  #define NGHOST 2
#endif // MHD

#define NGHOST 3

#ifdef MHD
  #define NVAL 8
  #ifdef ISO
    #define NWAVE 6
  #else
    #define NWAVE 7
  #endif // ISO
#else
  #define NVAL 5
  #ifdef ISO
    #define NWAVE 4
  #else
    #define NWAVE 5
  #endif // ISO
#endif // MHD

enum class bounds
{
    PERIODIC,
    OUTFLOW,
    REFLECTING,
    PRESSURE    /* For RT instability, extrapolate pressure with grav. potential to ghost cells to improve conservation of energy */
};

enum class exports
{
    PRIM,
    CONS,
    Rho,
    Vvec,
    VX,
    VY,
    VZ,
    Mvec,
    MX,
    MY,
    MZ,
    E,
    P,
    Bvec,
    Bx,
    By,
    Bz,
    DIVB
};

std::string exp_names[18] = {
    "Primitive_variables",
    "Conserved_variable",
    "Density",
    "Velocity_vector",
    "Vx",
    "Vy",
    "Vz",
    "Momentum_vector",
    "Mx",
    "My",
    "Mz",
    "Energy",
    "Pressure",
    "MagnField_vector",
    "Bx",
    "By",
    "Bz",
    "B_Divergence"
};