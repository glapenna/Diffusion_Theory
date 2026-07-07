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
%************************************************************************/
#include "head.h"

/*
 *  pointax.c - get points, axes for vector-style measurement
 */

// extern char		tok[];
// char		tok[];
extern idstruct		*firstser;
extern parmstruct	*defprm;
 
idstruct	*getstruct(), *startid();
setstruct	*getsetaddr(), *getstatset();
 


getpointsaxes(int npoints, int naxes, int type, char *name)
{
	idstruct	*objptr, *idptr;
	char		errmsg[TOKLEN];
	char		id[TOKLEN];
	int		i, atom, grptype, nitems, points = 0, axes = 0;
	void		Backbone();
 
#ifdef DBG
	printf("$$ getting items\n");
#endif
	objptr = startid (type, id);
 
	/*
	 *  see if BACKBONE
	 */
 
	gettok();
	if (!strcmp(tok, "BACKBONE")) {
		if (type != TORSION)
			inerr(name, ": BACKBONE not an option");
		objptr->id_type = BACKBONE;
 
		Backbone(objptr);
		return;
	}
	for (i=0; i<4; i++)
		objptr->idunion.work.input[i] = NULL;
 
	if (type == ANGLE)
		initstat(&objptr->idunion.work.resunion.scalar.stat);
	else if (type == TORSION)
		initCstat(&objptr->idunion.work.resunion.scalar.cstat);
	objptr->idunion.work.nitems = 0;
	for (nitems=0; nitems<(npoints > naxes ? npoints : naxes); nitems++)
		objptr->idunion.work.intype[nitems] = 0;
	if (naxes)
		sprintf(errmsg, "%s %.*s: need %d points or %d axes",
				name,TOKLEN-40,id,npoints, naxes);
	else
		sprintf(errmsg, "%s %.*s: need %d points",
				name,TOKLEN-40,id,npoints);
 
	/*
	 * get elements; points (atom or group center/cmass) or axes
	 */
 
#ifdef DBG
	printf("$$ getting items\n");
#endif
	for (nitems=0;;) {
 
#ifdef DBG
		printf("$$ item %d\n", nitems);
#endif
		if (nitems > npoints)
			inerr(errmsg,"");
 
		if (tok[0] == ';') {
#ifdef DBG
			printf("$$ ; %s: %d-%d %d-%d\n", 
				name, npoints, points, naxes, axes);
#endif
			if (points  != npoints  &&  (axes != naxes || !naxes)) 
				inerr(errmsg,"");
			break;
		} else if (tokint(&atom)) {
#ifdef DBG
			printf("$$atom # %d\n", atom);
#endif
			if (axes)
				inerr(errmsg,"");
			objptr->idunion.work.input[nitems] = 
						(char *) setatom(atom);
			objptr->idunion.work.intype[nitems] = ATOM;
			points++;
			nitems++;
			continue;
		} else if (!strcmp(tok, "ATOM")) {
#ifdef DBG
			printf("$$ATOM # \n");
#endif
			if (axes)
				inerr(errmsg,"");
			gettok();
			if (!tokint(&atom))
				inerr("Expected atom # -",name);
#ifdef DBG
			printf("$$ATOM # %d\n", atom);
#endif
			objptr->idunion.work.input[nitems] = 
							(char *) setatom(atom);
			objptr->idunion.work.intype[nitems] = ATOM;
			points++;
			nitems++;
			continue;
		} else if ((idptr = getstruct(tok)) == NULL) {
 
#ifdef DBG
			printf("$$is it an atomname\n");
#endif
			/*
			 * last chance: see if it's an atom name
			 */
 
			objptr->idunion.work.input[nitems] = 
							(char *) getsetatom();
			if (objptr->idunion.work.input[nitems] == NULL)
				inerr(errmsg,"");	/* unexpected tok */
			objptr->idunion.work.intype[nitems] = ATOM;
			points++;
			nitems++;
			continue;
		}
		switch (idptr->id_type) {
		    case GROUP:
#ifdef DBG
			printf("$$group \n");
#endif
			if (axes)
				inerr(errmsg,"");
			points++;
			grptype = grouptype(idptr);
			objptr->idunion.work.intype[nitems] = grptype;

			/*
			 * set binding comes w/ group
			 */
			if (grptype == GROUPCENTER)
				objptr->idunion.work.input[nitems] =
					(char *)idptr->idunion.group.geomcenter;
			else if (grptype == GROUPCOM)
				objptr->idunion.work.input[nitems] =
					(char *)idptr->idunion.group.masscenter;
			else
				inerr(errmsg, "expected group center");
#ifdef DBG
			printf("$$group got type: %d\n",
				objptr->idunion.work.intype[nitems]);
#endif
			nitems++;
			break;
		    case AXIS:
		    case PLANE:
#ifdef DBG
			printf("$$axis/plane\n");
#endif
			if (!naxes)
				inerr("No axes allowed for now -",name);
			if (points)
				inerr(errmsg,"");
			axes++;
			objptr->idunion.work.intype[nitems] = idptr->id_type;
 
			/*
			 *  set binding comes w/ axis or plane
			 */

			objptr->idunion.work.input[nitems] = (char *)
				idptr->idunion.work.resunion.axis.vector;
			nitems++;
			gettok();
			break;
		    default:
			inerr("Inappropriate id type","");
		}
	}
	objptr->idunion.work.nitems = nitems;
}
 
/**************************************************************************
								DIST()
***************************************************************************/
 
dist()
{
	idstruct	*objptr, *idptr;
	char		*errmsg = "DIST: 2 items required";
	char		id[TOKLEN];
	int		atom, grptype, nitems = 0;
 
	objptr = startid (DIST, id);
 
	gettok();
	for (;;) {
		if (nitems > 2)
			inerr("DIST: only 2 items","");
		if (tok[0] == ';') {
			if (nitems != 2)
				inerr(errmsg, "");
			break;
		} else if (tokint(&atom)) {
			objptr->idunion.work.input[nitems] =
							(char *) setatom(atom);
			objptr->idunion.work.intype[nitems] = ATOM;
			nitems++;
			continue;
		} else if (!strcmp(tok, "ATOM")) {
			gettok();
			if (!tokint(&atom)) 
				inerr("DIST ATOM: expected atom #","");
			objptr->idunion.work.input[nitems] =
							(char *) setatom(atom);
			objptr->idunion.work.intype[nitems] = ATOM;
			nitems++;
			continue;
		} else if ((idptr = getstruct(tok)) == NULL) {
			/*
			 * last chance: see if it's an atom name
			 */
 
			objptr->idunion.work.input[nitems] = 
							(char *) getsetatom();
			if (objptr->idunion.work.input[nitems] == NULL)
				inerr("DIST: unexpected token","");
			objptr->idunion.work.intype[nitems] = ATOM;
			nitems++;
			continue;
		}
		switch (idptr->id_type) {
		    case GROUP:
#ifdef DBG
			printf("$$group \n");
#endif
			grptype = grouptype(idptr);
			objptr->idunion.work.intype[nitems] = grptype;
			/*
			 * set binding comes w/ group
			 */
			if (grptype == GROUPCENTER)
				objptr->idunion.work.input[nitems] =
					(char *)idptr->idunion.group.geomcenter;
			else if (grptype == GROUPCOM)
				objptr->idunion.work.input[nitems] =
					(char *)idptr->idunion.group.masscenter;
			else
				inerr(errmsg, "expected group center");
#ifdef DBG
			printf("$$group got type\n");
#endif
			nitems++;
			break;
		    case AXIS:
		    case PLANE:
			inerr("axis or plane not supported for dist","");
#ifdef DBG
			printf("$$axis/plane\n");
#endif
			objptr->idunion.work.intype[nitems] = idptr->id_type;
 
			/*
			 * set binding comes w/ axis or plane
			 */
			objptr->idunion.work.input[nitems] = (char *)
				idptr->idunion.work.resunion.axis.vector;
			nitems++;
			gettok();
			break;
		    default:
			inerr("DIST: unexpected id type -", 
				" need GROUP, AXIS or PLANE");
		}
		
	}
	initstat(&objptr->idunion.work.resunion.scalar.stat);
}
 
 
 
 
/**********************************************************************
							BACKBONE()
***********************************************************************/
 
/*
 *  Backbone() - set up backbone torsion
 */
 
void
Backbone(idstruct *idptr)
{
	idstruct	*serptr;
	parmstruct	*prm;
	char		*ResName, *AtomTree;
	int		*atomarray, *resarray;
	_REAL		**addrarray;
	_REAL		*crd;
	int		atnum, natom = 0, resnum, res1 = 1, res2 = 0, start = 0;
	int		maxarr, iarr = 0, i, consecatom = 0;
	int		iarr_first = 0, iarr_last = 0; 
	int		ndih = 0;
	Cstatstruct	*stat;
 
	/*
	 *  get optional residue #s
	 */
 
	gettok();
 
	if (tokint(&resnum)) {
 
		/*
		 *  a range of residues specd
		 */
 
		res1 = resnum;
		if (res1 < 1)
			inerr("BACKBONE: 1st residue number must be > 0","");
		gettok();
		if (tokint(&resnum)) {
			res2 = resnum;
			gettok();
		} 
	}
 
	/*
	 *  determine stream
	 */
 
	if (tok[0] == ';')
		serptr = firstser;
	else {
		if ((serptr = getstruct(tok)) != NULL) {
			if (serptr->id_type != STREAM)
				inerr("BACKBONE: got non-stream id","");
			gettok();
			if (tok[0] != ';')
				inerr("BACKBONE: expected ';'","");
		} else
			inerr("BACKBONE: unexpected token","");
	}
	crd = serptr->idunion.ser.set->crd;
 
	/*
	 *  assign parm, default res #s & check res #s
	 */
 
	if ((prm = serptr->idunion.ser.set->prm) == NULL)
		inerr("BACKBONE: no parm file defined","");
	idptr->idunion.work.prm = prm;
 
	if (!res2)
		res2 = serptr->idunion.ser.set->prm->Nres;
	else if (res2 > serptr->idunion.ser.set->prm->Nres)
		inerr("BACKBONE: last residue beyond parm Nres","");
 
	if (!(res1 <= res2))
		inerr("BACKBONE: 1st residue number must be < last", "");
 
	/*
	 *  figure an upper bound on the input, allocate space,
	 *	 and build the input array of backbone atom #s
	 */
 
	maxarr = res2 - res1 + 1;
	AtomTree = &prm->AtomTree[4 * (prm->Ipres[res1-1] - 1)];
	for (atnum=prm->Ipres[res1-1]; atnum<prm->Ipres[res2]; atnum++) {
		if (!strncmp(AtomTree, "M   ", 4)) 
			maxarr++;
		AtomTree += 4;
	}
 
	atomarray = (int *) get(maxarr * sizeof(int));
	resarray = (int *) get(maxarr * sizeof(int));
	addrarray = (_REAL * *) get(maxarr * sizeof(_REAL *));
	idptr->idunion.work.input[0] = (char *) atomarray;
	idptr->idunion.work.input[1] = (char *) resarray;
	idptr->idunion.work.input[2] = (char *) addrarray;
	atomarray[0] = 0;
 
	resnum = res1;
#ifdef DBG
	ResName = &prm->ResNames[4 * (res1 - 1)];
	printf("residue %d %.4s: ", resnum, ResName);
#endif
	AtomTree = &prm->AtomTree[4 * (prm->Ipres[res1-1] - 1)];
	for (atnum=prm->Ipres[res1-1]; atnum<prm->Ipres[res2]; atnum++,
							       AtomTree += 4) {
		if (atnum == prm->Ipres[resnum]) {

			/* 
			 *  1st atom of next res 
			 */

			resnum++;
			start++;
#ifdef DBG
			ResName += 4;
			printf("\nresidue %d %.4s: ", resnum, ResName);
#endif
		}

		/*
		 *  place list of backbone atoms in atomarray, leaving 1st 
		 *	number free for indicating whether res is joined to prev
		 */
 
#ifdef DBG
		printf("%.4s ", AtomTree);
#endif
		if (!strncmp(AtomTree, "M   ", 4)) {	

			/* 
			 *  main chain atom 
			 */

			if (start) {

				/*
				 *  1st atom in (residue > 1st)
				 * 	if 1st main chain atom in res not bonded 
				 *	to last res, break the chain
				 */

				start = 0;
				if (!IsBonded(atomarray[iarr_last], 
								atnum, prm)) {
#ifdef DBG
				    printf("no bond w/ prev res\n");
#endif
				    if (consecatom < 4) {

					/*
					 *  not enough 4 torsion - back up
					 */

					printf("< 4 main atoms in chain\n");
					iarr = iarr_first;
					iarr_last = iarr;
					continue;
				    }
				    addrarray[iarr] = NULL;
				    if (IsBonded(atomarray[iarr_last], 
					  atomarray[iarr_first], prm)) {
					printf("cycle, atoms %d->%d\n",
				    		atomarray[iarr_last], 
					  	atomarray[iarr_first]);
					atomarray[iarr] = -2;
					resarray[iarr] = iarr_first;
					ndih += consecatom;
				    } else {
					atomarray[iarr] = -1;
					resarray[iarr] = 0;
					ndih += consecatom - 3;
				    }
				    iarr++;
				    iarr_first = iarr;
				    consecatom = 0;
				}
			} 
			iarr_last = iarr;
			atomarray[iarr] = atnum;
			resarray[iarr] = resnum;
			addrarray[iarr] = &crd[3 * (atnum-1)];
			iarr++; natom++; consecatom++;
		}
	}
	if (!natom) {
		printf("No main chain atoms \n");
	}

	/*
	 *  see if last residue is part of a cycle
	 */

	if (IsBonded(atomarray[iarr_first], atomarray[iarr_last], prm)){
		printf("cycle at end, atoms %d->%d\n",
			atomarray[iarr_last], atomarray[iarr_first]);
		atomarray[iarr] = -2;
		resarray[iarr] = iarr_first;
		if (consecatom > 3)
			ndih += consecatom;
	} else {
		atomarray[iarr] = -1;
		resarray[iarr] = 0;
		if (consecatom > 3)
			ndih += consecatom - 3;
	}
	iarr++;
	printf("\n%d backbone atoms, %d dihedrals\n", natom, ndih);
	idptr->idunion.work.nitems = iarr;
	idptr->idunion.work.flag = ndih;
	idptr->idunion.work.resunion.list.list = 
				(char *) get(ndih * sizeof(_REAL));
	stat = (Cstatstruct *) get(ndih * sizeof(Cstatstruct));
	idptr->idunion.work.resunion.list.stat = (char *) stat;
	for (i=0; i<ndih; i++)
		initCstat(stat++);
#ifdef DBG
	for (i=0; i<iarr; i++) {
		if (!(i%7))
			printf("\n");
		switch (atomarray[i]) {
			case -1:
				printf("  |  ");
				break;
			case -2:
				printf(" <-%d:%d  ", resarray[resarray[i]], 
					atomarray[resarray[i]]);
				break;
			default:
				printf("%2d:%2d %.4s ",
					resarray[i], atomarray[i],
					&prm->AtomNames[4*(atomarray[i]-1)]);
				break;
		}
	}
	printf("\n");
#endif
}
 
/**********************************************************************
							ISBONDED()
***********************************************************************/
 
/*
 *  IsBonded() - is the pair bonded?
 */
 
int
IsBonded(int at1, int at2, parmstruct *prm)
{
	int	*BondAt1 = prm->BondAt1;
	int	*BondAt2 = prm->BondAt2;
	int	Nbona = prm->Nbona;
	int	i;
/*int *BondHAt1 = prm->BondHAt1, *BondHAt2 = prm->BondHAt2;
printf("bonded? %d %d\n", at1, at2);
*/
	/*
	 *  convert at1, at2 to the xyz index form  used in BondAt arrays
	 */
 
	at1 = 3 * (at1 - 1);
	at2 = 3 * (at2 - 1);
 
	/*
	 *  search starts from the beginning each time..
	 */
 
	for (i=0; i<Nbona; i++) {
		if ( (*BondAt1 == at1  &&  *BondAt2 == at2)
		   || (*BondAt1 == at2  &&  *BondAt2 == at1) )
			return(1);
		BondAt1++;
		BondAt2++;
	}
/*
	for (i=0; i<prm->Nbonh; i++) {
		printf("?? %d %d\n", *BondHAt1, *BondHAt2);
		BondHAt1++; BondHAt2++;
	}
*/
	return(0);
}

