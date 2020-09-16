#pragma once

#define SQ(x) x*x
#define sgn(x) (0 < x) - (x < 0)
#define Sgn(x) (0 <= x) - (x < 0)
#define null(x) (x <= 1e-14)*(x >= -1e-14)