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
#include <ctype.h>
#include <math.h>
 
/*
 *  rms.c - measure rms &/or position set for min rms
 */

/* #define DBG */
// extern  char            tok[];
// char            tok[];
extern  idstruct        *idlist, *firstser;
extern	int		it;
 
extern idstruct 	*startid(), *getstruct();
extern char		*get();
extern setstruct	*getstatset();
extern _REAL		grp_total_mass();

/***********************************************************************
 							RMS()
************************************************************************/
 
/*
 * rms() - set up rms idstruct
 */
 
/* TODO: make sure group, set1, set2 are defined on same parm file */
/* TODO: make sure com not calcd for group + this */
 
void
rms()
{
	char		id[TOKLEN];
	int		i, setsize, fit = 0, grptype;
	idstruct	*rmsptr, *grptr, *setptr;

	/*
	 * get rms record in idlist & start filling in
	 *	refcrd is the reference set, 
	 *	crd2 is the stream set to be compared.
	 * for RMS FIT:
	 *	refcrdctr is reference set w/ group center of mass at origin,
	 *	crd2 will be positioned into crd2a for minimum rms of group
	 */
 
	rmsptr = startid(RMS, id);
	rmsptr->idunion.rms.fit = 0;
	rmsptr->idunion.rms.it0 = 0;
	rmsptr->idunion.rms.refstatic = 0;
	rmsptr->idunion.rms.natoms = 0;
	rmsptr->idunion.rms.rms = -123.123;
	initstat(&rmsptr->idunion.rms.stat);
	rmsptr->idunion.rms.grp = NULL;
	rmsptr->idunion.rms.refcrd = NULL;
	rmsptr->idunion.rms.refcrdctr = NULL;
	rmsptr->idunion.rms.crd2 = NULL;
	rmsptr->idunion.rms.crd2a = NULL;
 
	/*
	 *  see if FIT
	 */

	gettok();
	if (!strcmp(tok, "FIT")) {
		rmsptr->idunion.rms.fit = 1;
		fit = 1;
		gettok();
	}

	/*
	 *  look up group id
	 */
 
	if ((grptr = getstruct(tok)) == NULL)
		inerr("RMS: expected group id","");
	if (grptr->id_type != GROUP)
		inerr("RMS: expected group id","");
	rmsptr->idunion.rms.grp = grptr;
	// if geometry of group is requested:
	grptype = grouptype(grptr);

	if (fit) {

		if (grptr->idunion.group.natoms < 3)
			inerr("RMS: attempting fit of group w/ < 3 atoms","");

		/*
		 *  make sure group's center of mass is calcd
		 */
 
		setupgrpcom(grptr);
	} else {

		/*
		 *  just calc group total_mass if not already set up for com
		 */

		if (!(grptr->idunion.group.flag & COMFLAG))
			grptr->idunion.group.total_mass =
				grp_total_mass(grptr->idunion.group.natoms,
					grptr->idunion.group.Atomlist,
					grptr->idunion.group.set->prm->Masses);

	}
	if (!(grptr->idunion.group.total_mass > 0.0))
		inerr("RMS: group total mass is 0 (?)","");
 
	/*
	 *  get stream or rms id; default is to compare first stream 
	 *	to its 1st set
	 */
 
	// if geometry of group is requested, token is already got in grouptype:
	// gettok();
	if (tok[0] == ';') {
		printf("RMS: comparing default stream to its initial set\n");
		if (grptr->idunion.group.set != firstser->idunion.ser.set)
			inerr("RMS: group must be defined on stream to be ",
				fit ? "positioned" : "compared");
		setsize =  firstser->idunion.ser.set->max;
		rmsptr->idunion.rms.natoms =  setsize;
		rmsptr->idunion.rms.crd2 = firstser->idunion.ser.set->crd;
		setsize *= 3;
 
		/*
		 *  make space for 1st set and note to copy it at 1st iteration
		 */
 
		if (fit) {
			rmsptr->idunion.rms.refcrdctr = (_REAL *)
						get(setsize * sizeof(_REAL));
		} else {
			rmsptr->idunion.rms.refcrd = (_REAL *)
						get(setsize * sizeof(_REAL));
		}
		rmsptr->idunion.rms.it0++;

		if (fit) {

			/*
			 *  make space for set to be positioned into
			 */

			rmsptr->idunion.rms.crd2a = (_REAL *)
						get(setsize * sizeof(_REAL));
#ifdef DBG
			printf("got crd2a %d  size %d\n", 
				rmsptr->idunion.rms.crd2a, setsize);
#endif
		}
		return;
	}
 
	/*
	 *  (1st) stream (to be positioned/compared) is specified: look it up
	 */
 
	if ((setptr = getstruct(tok)) == NULL)
		inerr("RMS: expected id of stream/rmsfit to be ", 
					fit ? "positioned" : "compared");
	/*
	 *  1st must be STREAM or RMS
	 */
	if (setptr->id_type == STREAM) {
printf("1st is stream\n");
		if (grptr->idunion.group.set != setptr->idunion.ser.set)
			inerr("RMS: group must be defined on stream to be ",
					fit ? "positioned" : "compared");
		setsize =  setptr->idunion.ser.set->max;
		rmsptr->idunion.rms.natoms = setsize;
		rmsptr->idunion.rms.crd2 = setptr->idunion.ser.set->crd;
	} else if (setptr->id_type == RMS) {

printf("1st is rms\n");
		/*
		 *  must be RMS FIT: its refcrdctr will be compared 
		 *	to its crd2a
		 */

		if (setptr == rmsptr)
			inerr("RMS: self reference not allowed", "");
		if (fit)
			inerr("RMS: can't fit another RMS: ",
						setptr->id_name);
		if (!setptr->idunion.rms.fit)
			inerr("RMS: other RMS must be RMS FIT: ", 
						setptr->id_name);
		if (setptr->idunion.rms.grp->idunion.group.set->prm != 
				rmsptr->idunion.rms.grp->idunion.group.set->prm)
			inerr("RMS: other RMS's group must use same PARM","");
		rmsptr->idunion.rms.refcrd = setptr->idunion.rms.refcrdctr;
		rmsptr->idunion.rms.crd2 = setptr->idunion.rms.crd2a;
		rmsptr->idunion.rms.it0 = setptr->idunion.rms.it0;
		rmsptr->idunion.rms.natoms = setptr->idunion.rms.natoms;
		grptr->idunion.group.total_mass =
			grp_total_mass(grptr->idunion.group.natoms,
					grptr->idunion.group.Atomlist,
					grptr->idunion.group.set->prm->Masses);
		gettok();
		if (tok[0] != ';')
			inerr("Secondary RMS of FIT: expected\n",
				" 'RMS rmsid2  grpid  rmsid1'\n");
		return;
 	} else {
		inerr("RMS: expected id of stream/rmsfit to be ",
					fit ? "positioned" : "compared");
	}

	/*
	 *  get reference set
	 */
 
	gettok();
	if (tok[0] == ';') {
 
		/*
		 *  reference set defaults to 1st set of stream
		 */
 
		printf("RMS: %s to initial set of stream\n",
					fit ? "positioning" : "comparing");
		setsize *= 3;
		if (fit) {
			rmsptr->idunion.rms.refcrdctr = (_REAL *)
						get(setsize * sizeof(_REAL));
			rmsptr->idunion.rms.crd2a = (_REAL *)
						get(setsize * sizeof(_REAL));
		} else
			rmsptr->idunion.rms.refcrd = (_REAL *)
						get(setsize * sizeof(_REAL));
		rmsptr->idunion.rms.it0++;
		return;
	}
 
	/*
	 *  explicit reference set should be idstruct for a coord set
	 */
 
	if ((setptr = getstruct(tok)) == NULL)
		inerr("RMS: expected reference set id or ';' - got ", tok);
	switch (setptr->id_type) {
		case STREAM:

			/*
			 *  every time a fit occurs, the current set
			 *	in the stream will be centered into refcrdctr
			 */

			if (setptr->idunion.ser.set->max != setsize)
				inerr("RMS: stream and reference set",
						" are different sizes");
			rmsptr->idunion.rms.refcrd = 
						setptr->idunion.ser.set->crd;
			if (fit)
				rmsptr->idunion.rms.refcrdctr = (_REAL *)
						get(3 * setsize *sizeof(_REAL));
			break;
		case STATIC: 
		{
			_REAL	*refcrd, *refcrdctr, refcenter[3];
 

			if (setptr->idunion.stat.natoms < 
						grptr->idunion.group.maxatom)
				inerr("RMS: static reference set has fewer",
						" max atom in group");
			refcrd = (getstatset(setptr))->crd;
			rmsptr->idunion.rms.refcrd = refcrd;
			if (fit) {

				/*
				 *  calculate reference center of mass once and
				 *	center coords
				 */
 
				refcrdctr = (_REAL *) 
					get(3 * setsize * sizeof(_REAL));
				rmsptr->idunion.rms.refcrdctr = refcrdctr;
				grpcenter(refcrd,
					grptr->idunion.group.set->prm->Masses,
					grptr->idunion.group.AtomXlist,
					grptr->idunion.group.Atomlist,
					grptr->idunion.group.natoms,
					COMFLAG, (_REAL *) NULL, refcenter,
					grptr->idunion.group.total_mass, 
					(_REAL *) NULL, (_REAL *) NULL,
					grptr->idunion.group.q);
				refcrd--; refcrdctr--;
				for (i=0; i<setsize; i++) {
					*++refcrdctr = *++refcrd - refcenter[0];
					*++refcrdctr = *++refcrd - refcenter[1];
					*++refcrdctr = *++refcrd - refcenter[2];
				}
				rmsptr->idunion.rms.refstatic++;
			}
			break;
		}
		default:
			inerr("RMS: need stream or static id",
					 " for reference set");
	}
	if (fit) {
		setsize *= 3;
		rmsptr->idunion.rms.crd2a = 
				(_REAL *) get(setsize * sizeof(_REAL));
	}
	gettok();
	if (tok[0] != ';')
		inerr("RMS: expected ';'","");
}
 
 
/************************************************************************
							FITRMS()
************************************************************************/
 
/*
 *  fitrms() - compare or position coords for min rms
 */
 
void
fitrms(idstruct *idptr)
{
	_REAL 		rms_val = 0.0;
	_REAL 		*refcrd = idptr->idunion.rms.refcrd;
	_REAL 		*refcrdctr = idptr->idunion.rms.refcrdctr;
	_REAL 		*crd2 = idptr->idunion.rms.crd2;
	_REAL 		*crd2a = idptr->idunion.rms.crd2a;
	_REAL		*refcenter = idptr->idunion.rms.refcenter;
	int		i, natoms = idptr->idunion.rms.natoms;
	int		fit = idptr->idunion.rms.fit;
	idstruct	*group = idptr->idunion.rms.grp;
			/* group is defined on stream to be positioned */
	_REAL		x =  group->idunion.group.masscenter[0];
	_REAL		y =  group->idunion.group.masscenter[1];
	_REAL		z =  group->idunion.group.masscenter[2];

#ifdef DBG
	printf("rms() %s fit %d it %d  natoms %d  refcrd %d/%d crd2 %d/%d\n", 
		idptr->id_name, fit, it, natoms, idptr->idunion.rms.refcrd,
		idptr->idunion.rms.refcrdctr, idptr->idunion.rms.crd2,
		idptr->idunion.rms.crd2a);
	printf("\tcom %f %f %f\n", x, y, z);
#endif
	if (idptr->idunion.rms.it0 && it == 0) {

 		idptr->idunion.rms.rms = 0.0;

		/*
		 *  initialize by copying stream's 1st set for reference
		 */
 
#ifdef DBG
		printf("refcrdctr %d  crd2 %d\n", refcrdctr, crd2);
#endif
		if (fit) {

			/*
			 *  shift center of mass to origin 
			 */

			refcrdctr--; crd2--;
			for (i=0; i<natoms; i++) {
				*++refcrdctr = *++crd2 - x;
				*++refcrdctr = *++crd2 - y;
				*++refcrdctr = *++crd2 - z;
			}

			/*  
			 *  note that reference set is static now
			 */

			idptr->idunion.rms.refstatic++;

			/*
			 *  copy shifted version to output version in case 
			 *	a COORD wants to print it
			 */

			if ((char *)memcpy((char *)crd2a, (char *)idptr->idunion.rms.refcrdctr,
				    3*natoms*sizeof(_REAL)) != (char *) crd2a) {
				perror("memcpy");
				exit(1);
			}
		} else {

			/*
			 *  just copy 1st set to refcrd
			 */

			if ((char *)memcpy((char *)refcrd, (char *)crd2, 
				3*natoms*sizeof(_REAL)) != (char *)refcrd) {
				perror("memcpy");
				exit(1);
			}
		}
#ifdef DBG
		printf("copied.. end-of-set refcrd %d refcrdctr %d  crd2 %d\n", 
					refcrd, refcrdctr, crd2);
#endif
		return;
	}
 
	/*
	 *  if just measuring, do it and return
	 */

	rms_val = 0.0;
	if (!fit) {
		int 	k, grpnatoms = group->idunion.group.natoms;
		int	*Xlist = group->idunion.group.AtomXlist;
		int	*list = group->idunion.group.Atomlist;
		_REAL	*Masses = group->idunion.group.set->prm->Masses;
		_REAL	distancesq();
#ifdef DEBUG
		printf("just comparing, (not fitting)..\n");
#endif
		for (k=0; k<grpnatoms; k++) {
			i = Xlist[k];
			rms_val += distancesq(&refcrd[i], &crd2[i]) 
				* Masses[list[k]];
#ifdef DBG
			printf("dist^2 %f * mass %f .. accum  %f\n", 
				distancesq(&refcrd[i], &crd2[i]) ,
				Masses[list[k]], rms_val);
#endif
		}
		rms_val = sqrt(rms_val / group->idunion.group.total_mass);
		/* printf("rms = %f\n", rms_val); */
		idptr->idunion.rms.rms = rms_val;
		updatestat(&idptr->idunion.rms.stat, rms_val);
		return;
	}

	/*  FITTING: */

	/*
	 *  if reference isn't static, calculate its center of mass
	 *	and copy it, positioning it at origin
	 */
 
	if (!idptr->idunion.rms.refstatic) {
#ifdef DBG
		printf("centering ref copy..\n");
#endif
		grpcenter(idptr->idunion.rms.refcrd, 
				group->idunion.group.set->prm->Masses,
				group->idunion.group.AtomXlist,
				group->idunion.group.Atomlist,
				group->idunion.group.natoms,
				COMFLAG,  (_REAL *) NULL, refcenter,
				group->idunion.group.total_mass, 
				(_REAL *) NULL, (_REAL *) NULL,
				group->idunion.group.q);
		refcrdctr--; refcrd--;
		for (i=0; i<natoms; i++) {
			*++refcrdctr = *++refcrd - refcenter[0];
			*++refcrdctr = *++refcrd - refcenter[1];
			*++refcrdctr = *++refcrd - refcenter[2];
		}
	}
 
	/*
	 *  copy set to be positioned, translating center of mass to origin
	 */
 
	crd2a--; crd2--;
	for (i=0; i<natoms; i++) {
		*++crd2a = *++crd2 - x;
		*++crd2a = *++crd2 - y;
		*++crd2a = *++crd2 - z;
	}
 
	/*
	 *  do rotations for least square fit
	 */
 
	rotfit(group->idunion.group.natoms,
		group->idunion.group.Atomlist,
		group->idunion.group.AtomXlist,
		natoms,
		idptr->idunion.rms.crd2a,
		idptr->idunion.rms.refcrdctr,
		group->idunion.group.set->prm->Masses,
		group->idunion.group.total_mass, &rms_val);
	/* printf("rms = %f\n", rms_val); */
	idptr->idunion.rms.rms = rms_val;
	updatestat(&idptr->idunion.rms.stat, rms_val);
}
 

int 	jacobi3();
_REAL	dot();

/***********************************************************************
 							ROTFIT()
************************************************************************/
 
/*
 *  rotfit() - rotate one set of coordinates to minimize mass-weighted
 *	rms deviation between it and a target set for a group of atoms,
 *	where both sets have centers of mass placed at the origin.
 *
 *	This code was adapted from Eric Fauman's implementation of
 *	the Kabsch method described in:
 *
 *		W. Kabsch, Acta Cryst. A32, 922-923 (1976)
 *		W. Kabsch, Acta Cryst. A34, 827-828 (1978)
 */

int 
rotfit(int grpnatoms, int *list, int *Xlist, int natoms, 
       _REAL *Mov, _REAL *Ref, _REAL *Masses, _REAL total_mass, _REAL *rms_val)
{
	int 	i, j, k, i3, k3;
	_REAL 	Rmatrix[9], Rtranspose[9], rtr[9], Emat[9], b[9], U[9];
	_REAL	*Evec[3], Eigenvalue[3], cp[3];
	_REAL	*Xref, *Xmov, x, y, z;
	_REAL 	wt, sig3, ratmsq = 0.0;

#ifdef DBG
	printf("rotfit grpnatoms %d, lists %d %d, natoms %d,\n",
		grpnatoms, list, Xlist, natoms);
	printf(" Mov/Rev %d %d, Masses %d, total_mass %f \n",
		Mov, Ref, Masses, total_mass);
#endif
	/*
	 *  zero "matrices"
	 */

	for (i=0; i<9; i++) {
		Rmatrix[i] = 0;
		rtr[i] = 0;
		U[i] = 0;
	}

	/* now calculate the R matrix:	             */
	/* r(i,j) = sum(over n) y(n,i) * x(n,j)      */
	/* where x(n) and y(n) are two vector sets   */

	for (k=0; k<grpnatoms; k++) {

		/* 
		 *  get pointers to the coords of the k'th group atom 
		 *	in each set: Xmov & Xref address the X coord of the atom
		 *  also note the mass
		 */

		Xmov = &Mov[Xlist[k]];
		Xref = &Ref[Xlist[k]];
		wt   = Masses[list[k]];

		/*
		 *  coord sets are mass centered at origin
		 */

#ifdef DBG2
		printf("wt %e  mov %e %e %e\nref %e %e %e\n",
		wt, Xmov[0], Xmov[1], Xmov[2], Xref[0], Xref[1], Xref[2]);
#endif
		/*
		 *  accumulate the mass-weighted sum of the squares
		 *	of the distances from center of mass of Move 
		 *	and Reference atom
		 */
		ratmsq += wt *
		        ( Xmov[0] * Xmov[0] 	/*  x^2  */
			+ Xmov[1] * Xmov[1] 	/*  y^2  */
			+ Xmov[2] * Xmov[2]	/*  z^2  */

			+ Xref[0] * Xref[0]	/*  x^2  */
			+ Xref[1] * Xref[1]	/*  y^2  */
			+ Xref[2] * Xref[2]);	/*  z^2  */

		/*
		 *  add this atom's mass-weighted Move & Reference
		 *	cross products to the Rmatrix
		 */
		for (i=i3=0; i<3; i++,i3+=3)
			for (j=0; j<3; j++)
				Rmatrix[i3 + j] += wt * Xmov[i] * Xref[j];
	}
	ratmsq *= 0.5;
#ifdef DBG2
	printf("E0 = %e\n",ratmsq);
#endif
	/*
	 *  build the transpose of the Rmatrix
	 */
	for (i=i3=0; i<3; i++,i3+=3)
		for (j=0; j<3; j++)
			Rtranspose[i3 + j] = Rmatrix[j * 3 + i];

	/*
	 *  make transpose X R
	 */
	for (i=i3=0; i<3; i++,i3+=3) 
		for (j=0; j<3; j++) {
			for (k = 0; k<3; k++)
				rtr[i3 + j] += Rtranspose[k*3 + i]
						   * Rmatrix[j*3 + k];
		}
		

#ifdef DBG2
	printf("Rmatrix:\n");
	parr(Rmatrix);
	printf("Rtranspose:\n");
	parr(Rtranspose);
	printf("R(transpose)R:\n");
	parr(rtr);
#endif
	/*
	 *  diagonalize and sort by eigenvector
	 */
	if (!diagEsort(rtr, Emat, Evec, Eigenvalue))
		return(0);

	/* 
	 *  now set Evec(2) = ( Evec(0) cross Evec(1) ) 
	 */

	cross(Evec[2], Evec[0], Evec[1]);

#ifdef DBG2
	printf("e crossed\n%e   %e   %e\n%e   %e   %e\n%e   %e   %e\n",
		Evec[0][0], Evec[0][1], Evec[0][2],
		Evec[1][0], Evec[1][1], Evec[1][2],
		Evec[2][0], Evec[2][1], Evec[2][2]);
#endif

	for (i=i3=0; i<3; i++,i3+=3) {
		for (j=0; j<3; j++) {
			b[i3 + j] = dot(Evec[i], &Rtranspose[3*j]);
		}
		normalize(&b[i3]);
	}

	/*
	 *  calc 3rd b vector & sig3
	 */

	cross(cp, &b[0], &b[3]);
	if (dot(cp, &b[6]) < 0.0)
		sig3 = -1.0;
	else 
		sig3 = 1.0;
	b[6] = cp[0]; b[7] = cp[1]; b[8] = cp[2];

#ifdef DBG2
	printf("b crossed ..sig3 %e  3rd b\n", sig3);
	parr(b);
#endif

	/*
	 *  calc rotation matrix U
	 */

	for (k=k3=0; k<3; k++,k3+=3)
		for (i=i3=0;i<3; i++,i3+=3)
			for (j=0; j<3; j++) {
				U[i3 + j] += Evec[k][j] * b[k3 + i];
#ifdef DBG2
				printf("u[%d][%d]  +=  %e * %e = %e\n", 
					i,j, e[k][j],b[k3 + i],U[i3+j]);
#endif
			}
#ifdef DBG2
	printf("matrix U:\n");
	parr(U);
#endif

	/*
	 *  rotate coords to fit
	 */

	Xmov = Mov;
        for (k=0; k<natoms; k++) {
		x = Xmov[0] * U[0]  +  Xmov[1] * U[1]  +  Xmov[2] * U[2];
		y = Xmov[0] * U[3]  +  Xmov[1] * U[4]  +  Xmov[2] * U[5];
		z = Xmov[0] * U[6]  +  Xmov[1] * U[7]  +  Xmov[2] * U[8];
		Xmov[0] = x;
		Xmov[1] = y;
		Xmov[2] = z;
		Xmov += 3;
        }

	/*
	 *  calc least squares rms
	 */

	*rms_val = ratmsq - (_REAL)sqrt(fabs(Eigenvalue[0])) 
		      - (_REAL)sqrt(fabs(Eigenvalue[1]))
		      - sig3 * (_REAL)sqrt(fabs(Eigenvalue[2]));
	*rms_val = (_REAL)sqrt( (2.0 * *rms_val) / total_mass);
	return(1);
}

