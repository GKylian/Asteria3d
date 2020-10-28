#pragma once

//#define SQ(x) x*x
//#define sgn(x) 1.0*(0.0 < x) - 1.0*(x < 0.0)
#define Sgn(x) 1.0*(0.0 <= x) - 1.0*(x < 0.0)


ld sgn(ld x) {
	return (ld)((0.0 < x) - (x < 0.0));
}

ld minmod(ld x, ld y) { 
	return sgn(x)*fminl(fabsl(x), fabsl(y)) * (x*y > 0.0);
}

bool null(ld x) {
	return (fabsl(x) <= 1e-14);
}

ld SQ(ld x) {
	return x*x;
}

void pArray(ld *a, int size) {
	cout << "( " << a[0];
	for (int i = 1; i < size; i++)
		cout << ", " << a[i];
	cout << ")";
}