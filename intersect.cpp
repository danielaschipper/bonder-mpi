#include "intersect.h"
#include <iostream>

int inline GetIntersection( double fDst1, double fDst2, vec3 P1, vec3 P2, vec3 *Hit) {
	if ( (fDst1 * fDst2) >= 0.0f) return 0;
	if ( fDst1 == fDst2) return 0; 
	*Hit = P1 + (P2-P1) * ( -fDst1/(fDst2-fDst1) );
	return 1;
}

int inline InBox(vec3* hit, vec3 B1, vec3 B2, const int Axis) {
	vec3 Hit = *hit;
	if ( Axis==1 && Hit.z > B1.z && Hit.z < B2.z && Hit.y > B1.y && Hit.y < B2.y) return 1;
	if ( Axis==2 && Hit.z > B1.z && Hit.z < B2.z && Hit.x > B1.x && Hit.x < B2.x) return 1;
	if ( Axis==3 && Hit.x > B1.x && Hit.x < B2.x && Hit.y > B1.y && Hit.y < B2.y) return 1;
	return 0;
}

// returns true if line (L1, L2) intersects with the box (B1, B2)
// returns intersection point in Hit

int CheckLineBox( vec3 B1, vec3 B2, vec3 L1, vec3 L2)
{
	vec3 hit(0,0,0);
	vec3 *Hit = &hit;
	if (L2.x < B1.x && L1.x < B1.x) return false;
	if (L2.x > B2.x && L1.x > B2.x) return false;
	if (L2.y < B1.y && L1.y < B1.y) return false;
	if (L2.y > B2.y && L1.y > B2.y) return false;
	if (L2.z < B1.z && L1.z < B1.z) return false;
	if (L2.z > B2.z && L1.z > B2.z) return false;
	if (L1.x > B1.x && L1.x < B2.x && L1.y > B1.y && L1.y < B2.y &&L1.z > B1.z && L1.z < B2.z) return true;
	if (               (GetIntersection( L1.x-B1.x, L2.x-B1.x, L1, L2, Hit) && InBox( Hit, B1, B2, 1 ))
			|| (GetIntersection( L1.y-B1.y, L2.y-B1.y, L1, L2, Hit) && InBox( Hit, B1, B2, 2 )) 
			|| (GetIntersection( L1.z-B1.z, L2.z-B1.z, L1, L2, Hit) && InBox( Hit, B1, B2, 3 )) 
			|| (GetIntersection( L1.x-B2.x, L2.x-B2.x, L1, L2, Hit) && InBox( Hit, B1, B2, 1 )) 
			|| (GetIntersection( L1.y-B2.y, L2.y-B2.y, L1, L2, Hit) && InBox( Hit, B1, B2, 2 )) 
			|| (GetIntersection( L1.z-B2.z, L2.z-B2.z, L1, L2, Hit) && InBox( Hit, B1, B2, 3 )))
		return true;

	return false;
}

/*
int CheckLineBox(vec3 B1, vec3 B2, vec3 L1, vec3 L2)
{
	std::cout << B2.x << " " << B2.y << " " << B2.z << std::endl;
	vec3 dv = (L2 - L1) * 0.0001;
	for (int i = 0; i < 10000; i++)
	{
		vec3 loc = L1 + (dv * i);
		if ((B1.x < loc.x) && (loc.x < B2.x) && (B1.y < loc.y) && (loc.y < B2.y) &&( B1.z < loc.z) && (loc.z < B2.z))
			return 1;
	}
	return 0;
}
*/
