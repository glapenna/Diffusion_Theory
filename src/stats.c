
/************************************************************************
*                            AMBER                                     **
*                                                                      **
* Copyright (c) 1986, 1991 Regents of the University of California     **
*                       All Rights Reserved.                           ** 
*                                                                      **
*  This software provided pursuant to a license agreement containing   **
*  restrictions on its disclosure, duplication, and use. This software **
*  contains confidential and proprietary information, and may not be   **
*  extracted or distributed, in whole or in part, for any purpose      **
*  whatsoever, without the express written permission of the authors.  **
*  This notice, and the associated author list, must be attached to    **
*  all copies, or extracts, of this software. Any additional           **
*  restrictions set forth in the license agreement also apply to this  **
*  software.                                                           **
*************************************************************************/
#include "head.h"

/*
 * STATISTICS
 */

void
initstat(statstruct *stat)
{
	stat->ct = 0;
	stat->weight = 0;
	stat->avg = 0.0;
	stat->avg2 = 0.0;
	stat->max = -9999999.9;
	stat->min = 9999999.9;
}

void
initCstat(Cstatstruct *stat)
{
	stat->ct = 0;
	stat->x = 0.0;
	stat->y = 0.0;
}

void
updatestat(statstruct *stat, _REAL val)
{
	stat->avg += val;
	stat->avg2 += val*val;
	if (!stat->ct++) {
		stat->max = stat->min = val;
	} else {
		if (val > stat->max)
			stat->max = val;
		if (val < stat->min)
			stat->min = val;
	}
}

void
updateCstat(Cstatstruct *stat, _REAL val)
{
	stat->ct++;
	stat->x += cos(val);
	stat->y += sin(val);
}

void
printstat(statstruct *stat, _REAL n, int convert)
{
	_REAL	avg;

	if (!stat->ct) {
		printf("     no data accumulated\n");
		return;
	}
	if (n == 0.0)
		n = (_REAL) stat->ct;

	if (n == 1.0) {
		if (convert)
			printf("single %e\n", TODEG * stat->avg);
		else
			printf("single %e\n", stat->avg);
		return;
	}
        if(stat->weight) n=1.;
	avg = stat->avg/n;
	if (convert)
		printf("avg    %e  dev  %e  max %e  min %e\n", 
				TODEG * avg, 
				TODEG * sqrt(stat->avg2/n - avg*avg),
				TODEG * stat->max, 
				TODEG * stat->min);
	else
		printf("avg    %e  dev  %e  max %e  min %e\n", 
				avg, sqrt(stat->avg2/n - avg*avg),
				stat->max, stat->min);
}

void
printCstat(Cstatstruct *stat, _REAL n)
{
	_REAL	theta, r, s;

	if (!stat->ct) {
		printf("     no data accumulated\n");
		return;
	}
	if (n == 0.0)
		n = (_REAL) stat->ct;
	stat->x /= n;
	stat->y /= n;

	theta = atan(stat->y / stat->x);
	if (stat->x < 0.0)
		theta += _PI;

	theta *= TODEG;
	if (theta > 180.0) 
		theta  -= (TODEG * 2 * _PI);

	if (n == 1.0) {
		printf("single %8.3f\n", theta);
		return;
	}
	r = sqrt(stat->x * stat->x  +  stat->y * stat->y);
	s = sqrt(2 * (1-r));

	printf("theta  %8.3f  s    %8.3f\n", theta, TODEG * s);
}
