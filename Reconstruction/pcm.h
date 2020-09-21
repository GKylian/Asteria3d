#pragma once
#include "../Defs/Arrays.h"
#include "eigen.h"
#include "../prob.h"

/* Computes the left and right interfaces of the (i, j, k) cell using PLM reconstruction
   (with characteristic variables characteritic tracing)                                 */


bool PCM(Arrays *u, int i, int j, int k) {
    int iR = 2*i; int iL = 2*i+1; int jR = 2*j; int jL = 2*j+1; int kR = 2*k; int kL = 2*k+1; //interface indices of the right (i-1/2) and left (i+1/2) interfaces (in the cell i, j, k)
    

    if (u->Nx > 1) {
        for (int n = 0; n < NVAL; n++)
        {
            u->ix(n, iR, j, k) = u->uC(n, i, j, k); u->ix(n, iL, j, k) = u->uC(n, i, j, k);
        }
#ifdef MHD
        u->ix(5, iL, j, k) = u->uC(5, i+1, j, k);
#endif // MHD
    }

    if (u->Ny > 1) {
        for (int n = 0; n < NVAL; n++)
        {
            u->iy(n, i, jR, k) = u->uC(n, i, j, k); u->iy(n, i, jL, k) = u->uC(n, i, j, k);
        }
#ifdef MHD
        u->iy(6, i, jL, k) = u->uC(6, i, j+1, k);
#endif // MHD
    }

    if (u->Nz > 1) {
        for (int n = 0; n < NVAL; n++)
        {
            u->iz(n, i, j, kR) = u->uC(n, i, j, k); u->iz(n, i, j, kL) = u->uC(n, i, j, k);
        }
#ifdef MHD
        u->iz(7, i, j, kL) = u->uC(7, i, j, k+1);
#endif // MHD
    }

    return true;
}