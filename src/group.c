 
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
 *  group.c - form lists of atom numbers per set theoretic operations:
 * 	union, intersection and negation over the universe of atoms, 1..max
 */

int	debug=0;
// extern 	char	tok[];
// char	tok[];
 
extern setstruct	*getsetaddr();
extern idstruct		*startid(), *getstruct();
extern char		*get();
 
static parmstruct 	*prm;
static int		max;
 
static idstruct		*grptr;/* global to enable detection of circular defn */
static char		id[TOKLEN];
 
static int		left = 0, leaf = 0;	/* parity maintained else exit */
 
/*   forward declarations   */

void 		setupgrpcom();
 
extern void	binop(char a,char *b,char *c);
extern void	resarr(char *a,int b,char c);
extern void     atdash(), resdash();
extern char	*newarr(), *resname(), *resnum(), *solute();
 
/***********************************************************************
                                                        GROUP()
************************************************************************/
      
/*
 * group() - process group declaration - result is a list of atoms 
 */
 
void
group()
{
	char		*grparr, *grpexpr();
	int		*list1, *list2;
	int		i, j, k, nat, maxgrp;
	setstruct	*set;
        _REAL           fric_Rouse,fric_Zimm;
 
	grptr = startid(GROUP, id);
	
	/*
	 *  flag is for calculating geom center, center of mass, radius of gyr
	 */
 
	grptr->idunion.group.flag = 0;
	grptr->idunion.group.total_mass = 0;
 
	/*
	 * perform set binding
	 */
 
	set = getsetaddr();
	grptr->idunion.group.set = set;
 	prm = set->prm;
	max = set->max;
 
	/*
	 * evaluate the group expression to a binary char array 
	 *	& check termination status 
	 */
 
	grparr = grpexpr();
#ifdef DBG
	printf("group %d %s\n", left, tok);
#endif
	if (left)
		inerr("GROUP: unmatched '('","");
	/*
	 * check punctuation 
	 */
 
	gettok();
#ifdef DBG
	printf("group %d %s\n", left, tok);
#endif
	if (tok[0] == ')')
		inerr("GROUP: unmatched ')'","");
	else if (tok[0] != ';')
		inerr("GROUP: unexpected token","");
 
	/*
	 * count resulting atoms in group 
	 */
 
	if (debug)
		printf("%.*s\n", max, grparr);
 
	nat = 0;
	for (i=0; i<max; i++) 
		if (grparr[i] == 'T') 
			nat++;
	if (!nat)
		inerr("Null group defined, id=", id);
	grptr->idunion.group.natoms = nat;
 
	/* 
	 *  allocate space & write lists - one for atom #s and one for
	 *	atom*3 indexes
	 */
 
	list1 = (int *) get(sizeof(int)*nat);
	list2 = (int *) get(sizeof(int)*nat);
	grptr->idunion.group.Atomlist = list1;
	grptr->idunion.group.AtomXlist = list2;
	j = 0;
	for (i=0; i<max; i++)
		if (grparr[i] == 'T') {
			list1[j] = i;
			list2[j++] = 3 * i;
			maxgrp = i;
		}
		
	grptr->idunion.group.maxatom = maxgrp;
 
	printf("** Group %s: %d atom%c", id, nat, (nat > 1) ? 's' : ' ');
	if (debug) {
		for (i=0; i<nat; i++) {
			if (!(i%20))
				printf("\n   ");
			printf("%d ", list1[i]+1);
		}
		if (prm == NULL) {
			printf ("No parm file, can't print atom names\n");
		} else {
			for (i=0; i<nat; i++) {
				if (!(i%10))
					printf("\n");
				printf("%.4s", &prm->AtomNames[4*list1[i]]);
			}
		}
	}
	printf("\n");

/*
       grptr->idunion.group.vdw = (_REAL *)get2(nat,sizeof(_REAL));
       fric_Rouse=0.;
       fric_Zimm=0.;
       for(i=0;i<nat;i++) {
        _REAL A,C;
        j=grptr->idunion.group.Atomlist[i];
        if(prm->fric!=NULL) {
         fric_Zimm+=prm->fric[j]*prm->fric[j];
         fric_Rouse+=prm->fric[j];
        }
        grptr->idunion.group.vdw[i] = prm->ljAB[j];
*/
/*
        k=prm->Ntypes*prm->atomtype[j]+prm->atomtype[j];
        A = prm->ljAB[2*k];
        C = prm->ljAB[2*k+1];
        if(A!=0.&&C!=0) {
         A *= 1/C;
         C = pow(A,0.1666666666667);
         C *= prm->sig0r;
         C /= 2;
         if(C<=0.7) {
          C=1.0;
         }
         grptr->idunion.group.vdw[i] = C;
        } else {
         C=1.0;
         grptr->idunion.group.vdw[i] = C;
        }
*/
/*
        printf("vdw(%d)= %f\n",j,grptr->idunion.group.vdw[i]);
       }
       printf("group %s Rouse friction coefficient= %f\n",grptr->id_name,fric_Rouse);
       printf("group %s Zimm  friction coefficient= %f\n",grptr->id_name,l_sqrt(fric_Zimm));
*/
}
 
/***********************************************************************
                                                        GRPEXPR()
************************************************************************/
 
/* 
 * groupexpr() - do recursive 'set OP set' evaluation, producing a binary
 * 	chararray of size max having 'T' if an atom is in group, else '.' 
 */
 
char *
grpexpr()
{
	char	*arr, *neg(), *grpatm(), *grpres(), *grpgrp();
	char	op;

#ifdef DBG
	printf("grpexpr left %d  %s\n", left, tok);
#endif
	arr = NULL;
 
	if (tok[0] == '(') {

		/*
		 *  new subexpression: evaluate 1st part
		 */

		left++;
		gettok();
		arr = grpexpr();
		if (leaf) {

			/*
			 *  the leaf expr (ATOM etc) saw the ')'
			 */
#ifdef DBG
	printf("return leaf\n");
#endif
			leaf = 0;
			return(arr);
		}

		/*
		 *  now evaluate arr binop arr binop arr ...
		 */

		for (;;) {

			/*
			 *  see whether binary operator or done 
			 */
			gettok();
#ifdef DBG
			printf("part 2  %s\n", tok);
#endif
			switch (tok[0]) {
				case '&':
				case '|':
					op = tok[0];
					gettok();	/* for grpexpr() */
					binop(op, arr, grpexpr());
					break;
				case ')':

					/*
					 *  end of level
					 */
#ifdef DBG
	printf("end of level - %d\n", left);
#endif
					if (--left < 0)
					    inerr("GROUP: unmatched '('","");
					return(arr);
				case ';':
					inerr("GROUP: unexpected end of ",
							"expression");
					break;
				default:
					inerr("GROUP: unexpected token", "");
			}
		}
	} else if (tok[0] == '!') {
		return(neg());
	} else if (!strcmp(tok, "ATOM")) {
		leaf = 1;
#ifdef DBG
	printf("return grpatm\n");
#endif
		return(grpatm());
	} else if (!strcmp(tok, "RES")) {
		leaf = 1;
		return(grpres());
	} else if (!strcmp(tok, "GROUP")) {
		leaf = 1;
		return(grpgrp());
	} else if (!strcmp(tok, "SOLUTE")) {
		leaf = 1;
		return(solute());
	} else
		inerr("GROUP: unexpected token","");
#ifdef lint
	return((char *) NULL);
#endif
}
 
/***********************************************************************
                                                        NEG()
************************************************************************/
 
/* 
 * neg() - get a group array & return its negation 
 */
 
char *
neg()
{
 	char	*arr;
 	int 	i;
 
	gettok();
	arr = grpexpr();
 
	for (i=0; i<max; i++)
		if(arr[i] == 'T')
			arr[i] = '.';
		else
			arr[i] = 'T';
 
	return(arr);
}
 
/***********************************************************************
                                                        BINOP()
************************************************************************/
 
/* 
 * binop() - apply '&' or '|' to arrays a & b, free b, return a 
 */
 
void
binop(char op, char *a, char *b)
{
	int	i;
 
	if (debug)
		printf("%.*s\n%c\n%.*s\n", max, a, op, b);
 
	switch (op) {
		case '&':
			for(i=0; i<max; i++) 
				if (a[i] != 'T'  ||  b[i] != 'T')
					a[i] = '.';
			break;
		case '|':
			for(i=0; i<max; i++) 
				if (a[i] == 'T'  ||  b[i] == 'T')
					a[i] = 'T';
			break;
		default:
			printf("Internal error: unknown op: %c\n", op);	
			exit(1);
	}
	free (b);
}
 
/***********************************************************************
                                                        GRPATM()
************************************************************************/
 
/* 
 * grpatm() - produce group array from atom list 
 */
 
char *
grpatm()
{
	char		*atchar(), *atnums();
#ifdef DBG
	printf("grpatm\n");
#endif
	gettok();
	if (!strcmp(tok, "NAME")) {
		if (prm == NULL)
			inerr("GROUP: Need parm file for atom names","");
		return(atchar(prm->AtomNames));
/*
	} else if (!strcmp(tok, "TYPE")) {
		if (prm == NULL)
			inerr("GROUP: Need parm file for atom types","");
		return(atchar(prm->AtomSym));
*/
	} else
		return(atnums());
}
 
/***********************************************************************
                                                        ATCHAR()
************************************************************************/
 
/*
 * atchar() - get atoms matching string
 */
 
char *
atchar(char *charr)
{
	int	any = 0, i;
	char 	*arr = newarr();
 
	for (;;) {
		gettok();
		if (tok[0] == ')') {
			if (--left < 0)
				inerr("GROUP: unmatched ')'","");
			if (!any)
				inerr("GROUP ATOM NAME/TYPE: no name given","");
			return(arr);
		} else {
			/* assume it's an atom or atm type name */
			if (strlen(tok) > 4)
				inerr("GROUP ATOM NAME/TYPE: string too long",
									"");
			any = 0;
			/* TODO: use char ptr inc'd by 4; #define 4 */
			for (i=0; i<max; i++)
				if (samename(&charr[4*i])) {
					arr[i] = 'T';
					any++;
				}
			if (!any)
				inerr("GROUP ATOM NAME/TYPE: no atoms w/ name",
									"");
		}
	}
}
 
/***********************************************************************
                                                        SAMENAME()
************************************************************************/
 
/*
 * samename() - compare tok to a part of chararray, assumes tok <= 4 chars
 */
 
int
samename(char *s)
{
	int 	i;
 
	if (debug > 1)
		printf("/%s/%.4s/\n", tok, s);
 
	for (i=0; i<4; i++) {
		if (tok[i] == '\0'  &&  s[i] == ' ')	/* end of strings */
			return(1);
		else if (tok[i] == '?')			/* wild card */
			continue;
		else if (tok[i] != s[i])		/* mismatch char  */
			return(0);
	}
	return(1);				  /* both strings 4 char  */
}
 
/***********************************************************************
                                                        ATNUMS()
************************************************************************/
 
/*
 * atnums() - read numlist & return char binary array
 */
 
char *
atnums()
{
	int	at1 = 0, at2 = 0, dash = 0, num;
	char 	*arr = newarr();
 
	for (;;) {
#ifdef DBG
		printf("atnums - %d  %s\n", left, tok);
#endif
		if (tok[0] == ')') {
			if (--left < 0) {
				inerr("GROUP: unmatched ')'","");
			}
			if (dash) {
				atdash(arr, at1, at2);
			} else if (at1) {
				arr[at1-1] = 'T';
			} else
				inerr("Unexpected ')' in GROUP ATOM","");
#ifdef DBG
			printf("atnums done - %d  %s\n", left, tok);
#endif
			return(arr);
		} else if (tok[0] == ',') {
			if (dash) {
				atdash(arr, at1, at2);
				dash = 0;
			} else if (at1) {
				arr[at1-1] = 'T';
			} else
				inerr("Unexpected ',' in GROUP ATOM","");
			at1 = 0;
			at2 = 0;
		} else if (tok[0] == '-') {
			if (at2)
				inerr("Unexpected '-' in GROUP ATOM","");
			dash++;
			if (!at1)
				at1 = 1;
		} else {
			/* it'd better be a number */
			if (!sscanf(tok, "%d", &num))
				inerr("GROUP ATOM: unexpected token","");
            		if (num < 1)
     				inerr("atom number < 1 in GROUP ATOM","");
     			else if (num > max)
				inerr("atom number > parm in GROUP ATOM","");
			if (!at1) {
				at1 = num;
				at2 = 0;
			} else {
				if (!dash)
				    inerr("expected comma in GROUP ATOM","");
				at2 = num;
			}
		}
		gettok();
	}
}
 
/***********************************************************************
                                                        GRPRES()
************************************************************************/
 
/*
 * grpres() - figure out atoms in residue(s)
 */
 
char *
grpres()
{
	gettok();
	if (!strcmp(tok, "NAME")) {
		return(resname());
	} else
		return(resnum());
}
 
/***********************************************************************
                                                        RESNAME()
************************************************************************/
 
/*
 * resname() - get all atoms in residue(s) w/ given name(s)
 */
 
char *
resname()
{
	int	any = 0, i;
	char	*arr = newarr();
 
	if (prm == NULL)
		inerr("GROUP: Need parm file to interpret residues", "");
	for (;;) {
		gettok();
		if (tok[0] == ')') {
			if (--left < 0)
				inerr("GROUP: unmatched ')'","");
			if (!any)
				inerr("GROUP RES NAME: no name given","");
			return(arr);
		} else {
			/* assume it's a res name */
			if (strlen(tok) > 4)
				inerr("GROUP RES NAME: string too long", "");
			any = 0;
			/* TODO: usr ptr += 4 */
			for (i=0; i<prm->Nres; i++)
				if (samename(&prm->ResNames[i*4])) {
					any++;
					resarr(arr, i+1, 'T');
				}
			if (!any)
				inerr("GROUP RES NAME: no such residue","");
		}
	}
}
 
/***********************************************************************
                                                        RESNUM()
************************************************************************/
 
/*
 * resnum() - get all atoms in residue(s) w/ given numlist
 */
 
char *
resnum()
{
	int	res1 = 0, res2 = 0, dash = 0, num;
	char	*arr = newarr();
 
	if (prm == NULL)
		inerr("GROUP: Need parm file to interpret residues", "");
	for (;;) {
		if (tok[0] == ')') {
			if (--left < 0)
				inerr("GROUP: unmatched ')'","");
			if (dash) {
				resdash(arr, res1, res2);
			} else if (res1) {
				resarr(arr, res1, 'T');
			} else
				inerr("Unexpected ')' in GROUP RES ","");
			return(arr);
		} else if (tok[0] == ',') {
			if (dash) {
				resdash(arr, res1, res2);
				dash = 0;
			} else if (res1) {
				resarr(arr, res1, 'T');
			} else
				inerr("Unexpected ',' in GROUP RES","");
			res1 = 0;
			res2 = 0;
		} else if (tok[0] == '-') {
			if (res2)
				inerr("Unexpected '-' in GROUP RES","");
			dash++;
			if (!res1)
				res1 = 1;
		} else {
			/* it'd better be a number */
			if (!tokint(&num))
				inerr("GROUP RES: unexpected token","");
            		if (num < 1)
     				inerr("res number < 1 in GROUP RES","");
     			else if (num > prm->Nres)
				inerr("res number > parm in GROUP RES","");
			if (!res1) {
				res1 = num;
				res2 = 0;
			} else {
				if (!dash)
				    inerr("expected comma in GROUP RES","");
				res2 = num;
			}
		}
		gettok();
	}
}
 
/***********************************************************************
                                                        ATDASH()
************************************************************************/
 
/* 
 * atdash() - fill in range of atoms in array 
 */
 
void
atdash (char *arr, int at1, int at2)
{
	int	i;
 
	/* maybe range is all atoms */
	if (!at1) at1 = 1;
	if (!at2) at2 = prm->Natom;
	if (at1 > at2)
		inerr("GROUP ATOM: range is backward","");
	at2++;	
	for (i=at1; i<at2; i++)
		arr[i-1] = 'T';
}
 
/***********************************************************************
                                                        RESDASH()
************************************************************************/
 
/* 
 * resdash() - fill in range of atoms in array 
 */
 
void
resdash (char *arr, int res1, int res2)
{
	int	i;
 
	/* maybe range is all residues */
	if (!res1) res1 = 1;
	if (!res2) res2 = prm->Nres;
	if (res1 > res2)
		inerr("GROUP RES: range is backward","");
	res2++;	
	for (i=res1; i<=res2; i++)
		resarr(arr, i, 'T');
}
 
/***********************************************************************
                                                        RESARR()
************************************************************************/
 
/*
 * resarr() - set atoms in residue to given value in binchararray
 */
 
void
resarr(char *arr, int res, char val)
{
	int	i;
 
	for (i=prm->Ipres[res-1]-1; i<prm->Ipres[res]-1; i++)
		arr[i] = val;
}
 
/***********************************************************************
                                                        GRPGRP()
************************************************************************/
 
/* 
 *  grpgrp() - look up group id(s) & OR to produce logical array (union)
 */
 
char *
grpgrp()
{
	idstruct	*idptr;
	int		*list;
	int		i, one = 0;
	char		*arr = newarr();
 
	for (;;) {
		gettok();
		if (tok[0] == ')') {
			if (!one)
				inerr("GROUP GROUP: no id given","");
			if (--left < 0)
				inerr("GROUP: unmatched ')'","");
			return(arr);
		}
		one++;
		/*
	 	 * assume it's an id, look it up, then OR list into array
		 */
 
		if ((idptr = getstruct(tok)) == NULL)
			inerr("GROUP GROUP: no such group id:", tok);
		if (idptr->id_type != GROUP) 
			inerr("GROUP GROUP: non group id given","");
		if (grptr == idptr)
			inerr ("GROUP GROUP: self referencing group:", id);
 
		list = idptr->idunion.group.Atomlist;
		for (i=0; i<idptr->idunion.group.natoms; i++)
			arr[list[i]] = 'T';
	}
}
 
/***********************************************************************
                                                        SOLUTE()
************************************************************************/
 
/* 
 *  solute() - return array of solute atoms
 */
 
char *
solute()
{
	char 	*arr = newarr();
 	int 	i, natom = prm->Natom;
	int	solv = firstWAT(prm);
 
	gettok();
	if (tok[0] != ')')
		inerr("SOLUTE: expected ')'", "");
	if (--left < 0)
		inerr("GROUP: unmatched ')'", "");

	if (!solv) {

		/*
		 *  no water, so whole thing is solute
		 */

		solv = natom;
	}

	for (i=0; i<solv; i++)
		arr[i] = 'T';
	for (; i<natom; i++)
		arr[i] = ' ';
	return(arr);
}
 
/***********************************************************************
                                                        NEWARR()
************************************************************************/
 
/* 
 * newarr() - malloc an array & initialize it 
 */
 
char *
newarr()
{
	char	*arr;
	int	i;
 
	arr = get(max * sizeof(char));
	for (i=0; i<max; i++)
		arr[i] = '.';
	return(arr);
}
 
/***********************************************************************
 							GROUPTYPE()
************************************************************************/
 
/*
 * grouptype() - get token & see if a group type is specd; geometric center 
 *	is default
 */
 
int
grouptype(idstruct *grpt)
{
	int	type;
	void	setupgrpcom(), setupradgyr(),setupuz(),setupux();
 
	gettok();
	if (tok[0] == '%') {
		gettok();
		if(!strcmp(tok, "center")) {
			type = GROUPCENTER;
			grpt->idunion.group.flag |= GEOMFLAG;
		} else if (!strcmp(tok, "cmass")) {
			type = GROUPCOM;
			if (!(grpt->idunion.group.flag & COMFLAG)) {
 				setupgrpcom(grpt);
			}
		} else if (!strcmp(tok, "radgyr")) {
			type = GROUPGYR;
			printf ("sono qui.\n"); fflush(stdout);
			if (!(grpt->idunion.group.flag & GYRFLAG))
 				setupradgyr(grpt);
		} else
			inerr("Unknown group option","");
		gettok();
	} else {
		grpt->idunion.group.flag |= GEOMFLAG;
		type = GROUPCENTER;
	}
	return(type);
}

_REAL
grp_total_mass(int natoms, int *list, _REAL *Masses)
{
	_REAL	total_mass = 0.0;
	int	i;

	for (i=0; i<natoms; i++) {
		total_mass += Masses[list[i]];
	}
	return(total_mass);
}

void
setupgrpcom(idstruct *grpt)
{
	/*
	 * check for parm for set & set flag to cause com calc
	 */
 
	if (grpt->idunion.group.set->prm == NULL)
		inerr("Group center of mass requires parm file to",
				"be defined for stream");

	grpt->idunion.group.flag |= COMFLAG;
	grpt->idunion.group.total_mass = 
		grp_total_mass(grpt->idunion.group.natoms, 
			       grpt->idunion.group.Atomlist,
			       grpt->idunion.group.set->prm->Masses);
}


void
setupradgyr(idstruct *grpt)
{
	void	setupgrpcom();

	/*
	 * check for parm for set & set flag to cause com calc
	 */
 
	if (grpt->idunion.group.set->prm == NULL)
		inerr("Group radius of gyration requires parm file to",
				"be defined for stream");
	grpt->idunion.group.flag |= GYRFLAG;
	if (grpt->idunion.group.flag & COMFLAG)
		return;
	setupgrpcom(grpt);
}

/***********************************************************************
                                                        DOGRP()
************************************************************************/
 
/*
 *  dogrp() - calculate geometrical center, center of mass of group,
 *	radius of gyration as specd by flag
 */

void
dogrp(idstruct *grpt)
{
	void	grpcenter();

	if (grpt->idunion.group.flag) {
		grpcenter(grpt->idunion.group.set->crd, 
			  grpt->idunion.group.set->prm->Masses, 
			  grpt->idunion.group.AtomXlist, 
			  grpt->idunion.group.Atomlist, 
			  grpt->idunion.group.natoms, 
			  grpt->idunion.group.flag, 
			  grpt->idunion.group.geomcenter, 
			  grpt->idunion.group.masscenter, 
			  grpt->idunion.group.total_mass, 
			  &grpt->idunion.group.momin, 
			  &grpt->idunion.group.radgyr,
			  grpt->idunion.group.q);
#ifdef DEBUG
		printf("masscenter %s %f %f %f\n", grpt->id_name,
			grpt->idunion.group.masscenter[0], 
			grpt->idunion.group.masscenter[1], 
			grpt->idunion.group.masscenter[2]);
#endif
	}
}

/***********************************************************************
                                                        GRPCENTER()
************************************************************************/
 
/*
 *  grpcenter() - calculate geometrical center, center of mass of group,
 *	radius of gyration as specd by flag
 */
 
void
grpcenter(_REAL *crd, _REAL *Masses, int *Xlist, 
          int *list, int natoms, int flag, _REAL *geomcenter, 
          _REAL *masscenter, _REAL total_mass, _REAL *momin, 
          _REAL *radgyr, _REAL *q)
{
	int	i, j;
	_REAL	x, y, z, mass;
	_REAL	comX, comY, comZ, dX, dY, dZ;
        _REAL **Inert,**evec,*eval,*scp,ux[3],uy[3],uz[3],ut[3],
              ang[3],beh,alefp,alefm,
              qt[4],qtmp[4],utmp[4],
              ux0[]={0.,1.,0.,0.},
              uy0[]={0.,0.,1.,0.},
              uz0[]={0.,0.,0.,1.};
        int j1,j2,j3;
#ifdef DEBUG
	printf("grpcenter flags=%o\n", flag);
#endif

	/*
	 *  geometrical center
	 */
 
	x = 0; y = 0; z = 0; 
 
	if (flag & GEOMFLAG) {
		for (i=0; i<natoms; i++) {
			/* j is index to x crd of ith atom # */
			j = Xlist[i];
			x += crd[j++];
			y += crd[j++];
			z += crd[j];
		}
		geomcenter[0] = x / (_REAL) natoms;
		geomcenter[1] = y / (_REAL) natoms;
		geomcenter[2] = z / (_REAL) natoms;
	}
 
	/*
	 *  center of mass 
	 */
 
	if (!(flag & COMFLAG)) {
		return;
	}
 
	x = 0; y = 0; z = 0;
 
	for (i=0; i<natoms; i++) {
		j = Xlist[i];	/* get index to x crd of ith atom # */
		mass = Masses[list[i]];
		x += crd[j++] * mass;
		y += crd[j++] * mass;
		z += crd[j] * mass;
	}
	masscenter[0] = x / total_mass;
	masscenter[1] = y / total_mass;
	masscenter[2] = z / total_mass;

	/*
	 *  radius of gyration
	 */

	if ( !(flag & GYRFLAG) )
		return;

	comX = masscenter[0];
	comY = masscenter[1];
	comZ = masscenter[2];
	
	/*
	 *  calculate moment of inertia = Sum:atoms (distance_to_com^2)
	 */

        Inert=matrix(1,3,1,3);
        evec=matrix(1,3,1,3);
        eval=vector(1,3);
        scp=vector(1,3);
	*momin = 0.0;
        for(i=1;i<=3;i++) {
         for(j=i;j<=3;j++) {
          Inert[i][j]=0.;
         }
        }
	for (i=0; i<natoms; i++) {
		j = Xlist[i];	/* get index to x crd of ith atom # */
		dX = crd[j++] - comX;
		dY = crd[j++] - comY;
		dZ = crd[j] - comZ;
                Inert[1][1] += Masses[list[i]] * dX*dX;
                Inert[1][2] += Masses[list[i]] * dX*dY;
                Inert[1][3] += Masses[list[i]] * dX*dZ;
                Inert[2][2] += Masses[list[i]] * dY*dY;
                Inert[2][3] += Masses[list[i]] * dY*dZ;
                Inert[3][3] += Masses[list[i]] * dZ*dZ;
		dX *= dX;
		dY *= dY;
		dZ *= dZ;
		*momin += Masses[list[i]] * (dX + dY + dZ);
	}
	*radgyr = (_REAL) sqrt((double)(*momin / total_mass));

	printf("\tMoment of Inertia  = %8.2g g (A^2) / mol\n",*momin);
	printf("\tRadius of gyration = %8.2f Angstroms\n",*radgyr);

        Inert[2][1]=Inert[1][2];
        Inert[3][1]=Inert[1][3];
        Inert[3][2]=Inert[2][3];
        jacobi(Inert,3,eval,evec);
	j1 = Xlist[0];
	j2 = Xlist[1];
	j3 = Xlist[2];
        uz[0]=crd[j3]-crd[j1];
        uz[1]=crd[j3+1]-crd[j1+1];
        uz[2]=crd[j3+2]-crd[j1+2];
        ut[0]=crd[j1]-crd[j2];
        ut[1]=crd[j1+1]-crd[j2+1];
        ut[2]=crd[j1+2]-crd[j2+2];
        scp[1]=1./l_sqrt(uz[0]*uz[0]+uz[1]*uz[1]+uz[2]*uz[2]);
        for(i=0;i<3;i++) uz[i]*=scp[1];
        scp[1]=1./l_sqrt(ut[0]*ut[0]+ut[1]*ut[1]+ut[2]*ut[2]);
        for(i=0;i<3;i++) ut[i]*=scp[1];
        euler(uz,ut,ux,uy,ang);

/*
        beh = 0.5*ang[1];
        alefp = 0.5*(ang[0]+ang[2]);
        alefm = 0.5*(ang[0]-ang[2]);
        q[0] =  cos(beh)*cos(alefp);
        q[1] =  sin(beh)*cos(alefm);
        q[2] =  sin(beh)*sin(alefm);
        q[3] =  cos(beh)*sin(alefp);
        qt[0]= q[0];
        qt[1]=-q[1];
        qt[2]=-q[2];
        qt[3]=-q[3];
	j2 = Xlist[1];
        printf("MODEL  %5d\n",1);
        printf("ATOM  %5d  C0  MOL %5d    %8.3f%8.3f%8.3f\n",1,1,crd[j2],crd[j2+1],crd[j2+2]);
        proquat(qt,uz0,qtmp);
        proquat(qtmp,q,utmp);
        utmp[1]+=crd[j2]; utmp[2]+=crd[j2+1]; utmp[3]+=crd[j2+2];
        printf("ATOM  %5d  Cz  MOL %5d    %8.3f%8.3f%8.3f\n",2,1,utmp[1],utmp[2],utmp[3]);
        proquat(qt,ux0,qtmp);
        proquat(qtmp,q,utmp);
        utmp[1]+=crd[j2]; utmp[2]+=crd[j2+1]; utmp[3]+=crd[j2+2];
        printf("ATOM  %5d  Cx  MOL %5d    %8.3f%8.3f%8.3f\n",3,1,utmp[1],utmp[2],utmp[3]);
        proquat(qt,uy0,qtmp);
        proquat(qtmp,q,utmp);
        utmp[1]+=crd[j2]; utmp[2]+=crd[j2+1]; utmp[3]+=crd[j2+2];
        printf("ATOM  %5d  Cy  MOL %5d    %8.3f%8.3f%8.3f\n",4,1,utmp[1],utmp[2],utmp[3]);
        printf("ENDMDL\n");
*/

/*
        scp[1]=ux[0]*evec[1][1]+ux[1]*evec[2][1]+ux[2]*evec[3][1];
        scp[2]=ux[0]*evec[1][2]+ux[1]*evec[2][2]+ux[2]*evec[3][2];
        scp[3]=ux[0]*evec[1][3]+ux[1]*evec[2][3]+ux[2]*evec[3][3];
        for(i=1;i<=3;i++) scp[i]*=scp[i];
        eigsrt(scp,evec,3);

        scp[2]=uz[0]*evec[1][2]+uz[1]*evec[2][2]+uz[2]*evec[3][2];
        scp[3]=uz[0]*evec[1][3]+uz[1]*evec[2][3]+uz[2]*evec[3][3];
        if(scp[2]<scp[3]) {
         ut[0]=evec[1][2];
         ut[1]=evec[2][2];
         ut[2]=evec[3][2];
         evec[1][2]=evec[1][3];
         evec[2][2]=evec[2][3];
         evec[3][2]=evec[3][3];
         evec[1][3]=ut[0];
         evec[2][3]=ut[1];
         evec[3][3]=ut[2];
        }
        uz[0]= evec[1][1];
        uz[1]= evec[2][1];
        uz[2]= evec[3][1];
        ut[0]=-evec[1][3];
        ut[1]=-evec[2][3];
        ut[2]=-evec[3][3];
        ux[0]= evec[1][2];
        ux[1]= evec[2][2];
        ux[2]= evec[3][2];

        ut[0]= uy[0];
        ut[1]= uy[1];
        ut[2]= uy[2];
        uz[0]= evec[1][1];
        uz[1]= evec[2][1];
        uz[2]= evec[3][1];
*/
        ut[0]=-uz[0];
        ut[1]=-uz[1];
        ut[2]=-uz[2];
        uz[0]= ux[0];
        uz[1]= ux[1];
        uz[2]= ux[2];

/*
*/
        euler(uz,ut,ux,uy,ang);
        beh = 0.5*ang[1];
        alefp = 0.5*(ang[0]+ang[2]);
        alefm = 0.5*(ang[0]-ang[2]);
        q[0] =  cos(beh)*cos(alefp);
        q[1] =  sin(beh)*cos(alefm);
        q[2] =  sin(beh)*sin(alefm);
        q[3] =  cos(beh)*sin(alefp);
/*
        qt[0]= q[0];
        qt[1]=-q[1];
        qt[2]=-q[2];
        qt[3]=-q[3];

	j2 = Xlist[1];
        printf("MODEL  %5d\n",1);
        printf("ATOM  %5d  C0  MOL %5d    %8.3f%8.3f%8.3f\n",1,1,crd[j2],crd[j2+1],crd[j2+2]);
        proquat(qt,uz0,qtmp);
        proquat(qtmp,q,utmp);
        utmp[1]+=crd[j2]; utmp[2]+=crd[j2+1]; utmp[3]+=crd[j2+2];
        printf("ATOM  %5d  Cz  MOL %5d    %8.3f%8.3f%8.3f\n",2,1,utmp[1],utmp[2],utmp[3]);
        proquat(qt,ux0,qtmp);
        proquat(qtmp,q,utmp);
        utmp[1]+=crd[j2]; utmp[2]+=crd[j2+1]; utmp[3]+=crd[j2+2];
        printf("ATOM  %5d  Cx  MOL %5d    %8.3f%8.3f%8.3f\n",3,1,utmp[1],utmp[2],utmp[3]);
        proquat(qt,uy0,qtmp);
        proquat(qtmp,q,utmp);
        utmp[1]+=crd[j2]; utmp[2]+=crd[j2+1]; utmp[3]+=crd[j2+2];
        printf("ATOM  %5d  Cy  MOL %5d    %8.3f%8.3f%8.3f\n",4,1,utmp[1],utmp[2],utmp[3]);
        printf("ENDMDL\n");
*/

}
