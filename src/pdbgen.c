
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

/***********************************************************************
 							PDBGEN()
************************************************************************/

/*
 *  pdbATOM() - write an ATOM record to a pdb file
 */

pdbATOM(fp, atnum, resnum, crd, prm)
FILE		*fp;
_REAL		*crd;
parmstruct	*prm;
{
	char	atnam[5];
	char	*AtomNames = prm->AtomNames;
	char	*ResNames = prm->ResNames;

	/*
	 *  right-justify atom name
	 */
	if (AtomNames[4*atnum+3] == ' ')
		sprintf(atnam, " %.3s", &AtomNames[4*atnum]);
	else
		sprintf(atnam, "%.4s", &AtomNames[4*atnum]);

	fprintf(fp, "ATOM %6d %-.4s %-.4s %4d    %8.3f%8.3f%8.3f\n",
		atnum+1, atnam, &ResNames[4*resnum], resnum+1,
		crd[0], crd[1], crd[2]);
}

/*
 * pdbgen() - write out pdb file
 */

void
pdbgen(crd, prm, file, natoms)
_REAL		*crd;
parmstruct	*prm;
FILE		*file;
{
	int	j, k;
	_REAL	*ixyz = crd;
	char	atnam[5];
	char	*AtomNames = prm->AtomNames;
	char	*ResNames = prm->ResNames;

	for (j=0; j<prm->Nres; j++) {
#ifdef DBG
		printf("pdbgen res %d Nres %d\n", j+1, prm->Nres);
#endif
		for (k=prm->Ipres[j]-1; k<prm->Ipres[j+1]-1; k++) { /* -1 ? */
#ifdef DBG
			printf("pdbgen atom %d\n", k);
#endif
			if (k >= natoms)
				return;
			if (AtomNames[4*k+3] /*whatever*/ == ' ')
				sprintf(atnam, " %.3s", &AtomNames[4*k]);
			else
				sprintf(atnam, "%.4s", &AtomNames[4*k]);
			fprintf(file,
		    "ATOM %6d %-.4s %-.4s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
				k+1, atnam, &ResNames[4*j], j+1,
				ixyz[0], ixyz[1], ixyz[2],0.,0.);
/*
	  fprintf(file, "%8.3f%8.3f%8.3f\n",ixyz[0], ixyz[1], ixyz[2]);
*/
			ixyz += 3;
		}
	}
printf("pdbgen done\n");
}
