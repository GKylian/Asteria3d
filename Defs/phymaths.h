#pragma once

//#define SQ(x) x*x
//#define sgn(x) 1.0*(0.0 < x) - 1.0*(x < 0.0)
#define Sgn(x) 1.0*(0.0 <= x) - 1.0*(x < 0.0)

long double sgn(long double x) {
	return (long double)((0.0 < x) - (x < 0.0));
}

long double minmod(long double x, long double y) {
	return sgn(x)*fminl(fabsl(x), fabsl(y)) * (x*y > 0.0);
}

bool null(long double x) {
	return (fabsl(x) <= 1e-12);
}

long double SQ(long double x) {
	return x*x;
}

