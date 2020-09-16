#pragma once

#include "Definitions.h"

/*
<-------------------- Functions -------------------->
<->
---> bool allBounds(Arrays *u)
	 <-> Applies all the boundary conditions (BCs) by calling the corresponding functions.
	 - Arrays *u: the arrays containing the conserved values we're applying the BCs to.
<->
<->
---> bool bounds_xin(Arrays *u)
	 <-> Apply the BC at the inner boundary of x (x0).
<->
---> bool bounds_xout(Arrays *u)
	 <-> Apply the BC at the outer boundary of x (xn).
<->
*/
bool allBounds(Arrays *u);
bool bounds_xin(Arrays *u);
bool bounds_xout(Arrays *u);


enum bounds
{
	PERIODIC,
	OUTFLOW,
	REFLECTING,
	PRESSURE    /* For RT instability, extrapolate pressure with grav. potential to ghost cells to improve conservation of energy */
};

bool allBounds(Arrays *u) {


	return true;
}

bool bounds_xin(Arrays *u) {

}
