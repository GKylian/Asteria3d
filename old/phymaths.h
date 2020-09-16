#pragma once
#include <iostream>
#include <algorithm>

///<summary>
///Returns -1 if negative, 0 if null and 1 if positive
///</summary>
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

///<summary>
///Returns -1 if negative, and 1 if null or positive
///</summary>
template <typename T> int Sgn(T val) {
    return (val<T(0)? -1 : 1);
}


long double minmod(long double x, long double y) { // --ALL OK--
    return 0.5*(sgn(x) + sgn(y)) * std::min(fabsl(x), fabsl(y));
}

long double minmod(long double w, long double x, long double y, long double z) {
    return 0.125 * (sgn(w) + sgn(x)) * fabsl((sgn(w) + sgn(y))*(sgn(w) + sgn(z))) * std::min({ w, x, y, z });
}


/// <summary>
/// Compute the minmod limiter from r(i)
/// </summary>
/// <param name="r">r(i)</param>
/// <returns></returns>
long double lim_minmod(long double r) { // --ALL OK--
    return fmaxl(0, fminl(1, r));
}

long double lim_monocentral(long double r) {
    return fmaxl(0, fminl(fminl(2*r, 0.5*(1 + r)), 2));
}

long double lim_ospre(long double r) {
    return 1.5*(r*r+r)/(r*r+r+1);
}


bool in(long double a, long double b, long double x) {
    //Has to be between the minimum and the maximum of the two
    return(x >= fminl(a, b) && x <= fmaxl(a, b));
}

bool null(long double x) {
    return (fabsl(x) < 1e-16);
}

///
//void swap(long double *arr, int i1, int i2) {
//    long double a1 = arr[i1];
//    arr[i1] = arr[i2];
//    arr[i2]
//}