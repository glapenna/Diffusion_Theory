/************************************************************************
*			    AMBER				     	**
*								      	**
* Copyright (c) 1986, 1991 Regents of the University of California     	**
*		       All Rights Reserved.			   	** 
*								      	**
*  This software provided pursuant to a license agreement containing   	**
*  restrictions on its disclosure, duplication, and use. This software 	**
*  contains confidential and proprietary information, and may not be   	**
*  extracted or distributed, in whole or in part, for any purpose      	**
*  whatsoever, without the express written permission of the authors.  	**
*  This notice, and the associated author list, must be attached to    	**
*  all copies, or extracts, of this software. Any additional	   	**
*  restrictions set forth in the license agreement also apply to this  	**
*  software.							   	**
*************************************************************************/

#define DBG
/* #define DBG2 */
#define NUCLEIC_DEF
#include "head.h"
#include <math.h>

// extern char	tok[];
// char	tok[];
// void *get2(int,int);
extern idstruct	*firstser;

extern idstruct	*startid(), *getstruct();
extern char	*get();
extern _REAL	*setatom(), *getsetatom();
extern _REAL	size(), dot();

static char	id[TOKLEN];

static void	initpuck(), nucpuck(), pointpuck(), setsugar(), singlepucker();
static int	allnucs(), resnums(), resnames(), resdash();

/**************************************************************************
				PUCKER ANALYSIS
 **************************************************************************/

/*
 *  Routines for setting up & calculating ring puckering coordinates.
 *
 *  	Algorithm from "A General Definition of Ring Puckering Coordinates"
 *  	D.Cremer, J.A.Pople (JACS 96:6 pp1354-1358, 1975).
 *
 *	The result of this algorithm is modified to conform (approximately)
 *	to the Sundaralingham convention by adding 90 degrees if the NUCLEIC
 *	options are used. See Altona and Sundaralingam (1972), 
 *	(JACS 94 pp8205-8212), or p. 20 of Saenger's "Principles of Nucleic 
 *	Acid Structure", Springer-Verlag, 1983.
 *
 *  Pucker calc (singlepucker()) implemented by Andre Gorin with corrections
 *	by Alexey Zenkin and Bill Ross.
 */

void initpucker()
{
	idstruct	*pkptr;
        int exists;
 
	/*
	 *  PUCKER id: start & initialize pucker structure
	 */

	pkptr = startid(PUCKER, id);

	pkptr->id_subtype = 0;
	pkptr->idunion.puck.npuckers = 0;
	pkptr->idunion.puck.npoints = 0;
	pkptr->idunion.puck.points = NULL;
	pkptr->idunion.puck.results = NULL;
	pkptr->idunion.puck.stats = NULL;
	pkptr->idunion.puck.cases = 0;
	pkptr->idunion.puck.nboxes = 1;
	pkptr->idunion.puck.min = 0.;
	pkptr->idunion.puck.max = 360.;
	
        gettok();
        if (iskeywd())
                inerr("PUCKER: expected file name","");
        strcpy(pkptr->idunion.puck.fname, tok);

        /*
         *  open the file(s)
         */

        if (fileprob(pkptr->idunion.puck.fname, 0, &exists, "PUCKER"))
                exit(1);
        if ((pkptr->idunion.puck.file =
                                  fopen(pkptr->idunion.puck.fname, "w")) == NULL) {
                perror(pkptr->idunion.puck.fname);
                exit(1);
        }
        gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/**************************************************************************
								PUCKER()
 **************************************************************************/
/*
 *  pucker() - set up pucker calc
 *
 *	PUCKER id NUCLEIC|points ;
 */

void
pucker()
{
	idstruct	*pkptr;
	int		i, j, npoints, nresults;
	statstruct	*stats;
 
	/*
	 *  PUCKER id: get & initialize pucker structure
	 */

        gettok();
        if ((pkptr = getstruct(tok)) == NULL)
                inerr("PUCKER: no such id: ", tok);
        if (pkptr->id_type != PUCKER)
                inerr("PUCKER: wrong type of id: ", tok);

	/*
	 *  NUCLEIC|points
	 */

	gettok();
	if (!strcmp(tok, "NINT")) {
         gettok();
         tokint(&pkptr->idunion.puck.nboxes);
         gettok();
         tokREAL(&pkptr->idunion.puck.min);
         gettok();
         tokREAL(&pkptr->idunion.puck.max);
         gettok();
        }
	if (!strcmp(tok, "NUCLEIC")) {
		nucpuck(pkptr);
	} else if (tokint(&npoints)) {
		pointpuck(pkptr, npoints);
	} else
		inerr("PUCKER: expected 'NUCLEIC' or number of points", "");
	

	if (pkptr->idunion.puck.npoints < 4)
		inerr("PUCKER: number of points must be >= 4","");
	printf("** Pucker: %d pucker%c of %d points %s\n", 
			pkptr->idunion.puck.npuckers,
			((pkptr->idunion.puck.npuckers > 1) ? 's' : ' '),
			pkptr->idunion.puck.npoints,
			((pkptr->idunion.puck.npuckers > 1) ? "each" : ""));

	/*
	 *  get scratch array to copy points into to orient them
	 */

	pkptr->idunion.puck.scratch = (_REAL *) 
				get(3*pkptr->idunion.puck.npoints*sizeof(_REAL));
        pkptr->idunion.puck.scale=
         (pkptr->idunion.puck.max- pkptr->idunion.puck.min)/
          pkptr->idunion.puck.nboxes;
        pkptr->idunion.puck.boxes=
         (int *)get2(pkptr->idunion.puck.nboxes,sizeof(int));
        for(i=0;i<pkptr->idunion.puck.nboxes;i++) {
         pkptr->idunion.puck.boxes[i]=0;
        }

	/*
	 *  get space for results (output) part
	 */

	nresults = pkptr->idunion.puck.npuckers * 
			(pkptr->idunion.puck.npoints - 3);
	pkptr->idunion.puck.results = (_REAL *) get(nresults * sizeof(_REAL));
	stats = (statstruct *) get(nresults * sizeof(statstruct));
	pkptr->idunion.puck.stats = stats;

	/*
	 *  initialize stats: alternate stat structs are cast to Cstatstruct
	 *	(a smaller structure than statstruct, so it is safe)
	 */

	for (j=0; j<nresults/2; j++) {
		initstat(stats);
		stats++;
		initCstat((Cstatstruct *) stats);
		stats++;
	}
	if (nresults % 2)
		initstat(stats);
}


/**************************************************************************
								NUCPUCK()
**************************************************************************/
/*
 *  nucpuck() - set up calc of a nucleic acid pucker
 *
 *	PUCKER id NUCLEIC [stream] [residue_names|residue_numbers] ;
 */

static void
nucpuck(idstruct *pkptr)
{
	idstruct	*idptr;
	setstruct	*set;
	parmstruct	*prm;
	_REAL		*crd;
	_REAL		**points;
	char		*restmp;		/* boolean scratch array */
	int		i, dum, total = 0;

	pkptr->id_subtype = NUCLEIC;

	/*
	 *  determine parm & crd according to stream
	 */

	gettok();
	if ((idptr = getstruct(tok)) != NULL) {

		/*
		 *  [stream]
		 */

		if (idptr->id_type != STREAM)
			inerr("PUCKER NUCLEIC: expected stream id","");
		set = idptr->idunion.ser.set;
		gettok();
	} else {

		/*
		 *  default stream
		 */

		set = firstser->idunion.ser.set;
	}
	crd = set->crd;
	prm = set->prm;
	if (set->max != prm->Natom)
		inerr("PUCKER NUCLEIC: owing to programmer laziness, ",
				"parm atoms must match set");

	/*
	 *  prepare scratch array for keeping track of residues
	 */

	restmp = get(prm->Nres);
	for (i=0; i<prm->Nres; i++)
		restmp[i] = ' ';

	/*
	 *  see which residues are involved
	 */

	if (tok[0] == ';') {

		/*
		 *  default: all nucleic acid sugars
		 */

		total = allnucs(restmp, prm);
	} else if (tok[0] == '-'  ||  tokint(&dum)) {

		/*
		 *  residue numbers: -a,b,c-d,e,f-g,-h
		 */

		total = resnums(restmp, prm);
	} else {

		/*
		 *  better be residue names
		 */
	
		total = resnames(restmp, prm);
	}
	if (!total)
		inerr("PUCKER NUCLEIC: no residues defined","");

	/*
	 *  allocate space and set up arrays using 'total' and 'restmp'
	 */

	pkptr->idunion.puck.npoints = 5;	/* 5 sugar atoms in nucleic  */
	pkptr->idunion.puck.npuckers = total;
	points = (_REAL **) get(5 * total * sizeof(_REAL *));
	pkptr->idunion.puck.points = points;
	pkptr->idunion.puck.results = (_REAL *) get(2 * total * sizeof(_REAL));

	for (i=0; i<prm->Nres; i++) {
		if (restmp[i] == 'T') {
			setsugar(points, i, prm, crd);
			points += 5;
		}
	}
	free(restmp);
}	

/**************************************************************************
								ALLNUCS()
**************************************************************************/
/*
 *  allnucs() - build boolean array of all nucleic acid residues
 *
 *	PUCKER NUCLEIC;
 */

static int
allnucs(char *arr, parmstruct *prm)
{
	int	i, j, total = 0;

/*
        { char rt[ATOMTOKLEN];
	 for (i=0; i<prm->Nres; i++) {
          sprintf(rt,"%4s",&prm->ResNames[i*ATOMTOKLEN]);
          printf("%s\n",rt);
         }
        }
*/
	for (i=0; i<prm->Nres; i++) {
		for (j=0; j<N_NUCLEIC; j++) {
			if (!strncmp(NucleicNames[j], &prm->ResNames[i*4], 4)) {
				arr[i] = 'T';
				total++;
				break;
			}
		}
	}
	return(total);
}

/**************************************************************************
								RESNUMS()
**************************************************************************/
/*
 *  resnums() - build boolean array of residues by number
 *
 *	PUCKER NUCLEIC residue_numbers
 */

static int
resnums(char *arr, parmstruct *prm)
{
	int	num, dash = 0, res1 = 0, res2 = 0, total = 0;

	for (;;) {
		if (tok[0] == ';') {
			if (dash) {
				total += resdash(arr, prm, res1, prm->Nres);
			} else if (res1) {
				arr[res1-1] = 'T';
				total++;
			}
			break;
		} else if (tok[0] == ',') {
			if (dash) {
				total += resdash(arr, prm, res1, res2);
				dash = 0;
			} else if (res1) {
				arr[res1-1] = 'T';
				total++;
			} else
				inerr("Unexpected ',' in ",
					"PUCKER NUCLEIC residues");
			res1 = 0;
			res2 = 0;
		} else if (tok[0] == '-') {
			if (res2)
				inerr("Unexpected '-' in",
					"PUCKER NUCLEIC residues");
			dash++;
			if (!res1)
				res1 = 1;
		} else {
			/* it'd better be a number */
			if (!tokint(&num))
				inerr("PUCKER NUCLEIC residues: ",
						"unexpected token");
	    		if (num < 1)
     				inerr("res number < 1 in ",
					"PUCKER NUCLEIC residues");
     			else if (num > prm->Nres)
				inerr("res number > parm in ",
					"PUCKER NUCLEIC residues");
			if (!res1) {
				res1 = num;
				res2 = 0;
			} else {
				if (!dash)
				    inerr("expected comma in ",
					"PUCKER NUCLEIC residues");
				res2 = num;
			}
		}
		gettok();
	}
	return(total);
}
 
static int
resdash(char *arr, parmstruct *prm, int res1, int res2)
{
	int		i;

	/* maybe range is all residues */
	if (!res1) res1 = 1;
	if (!res2) res2 = prm->Nres;
	if (res1 > res2)
		inerr("PUCKER NUCLEIC residue_numbers: range is backward","");
	res2++;	
	for (i=res1; i<res2; i++)
		arr[i-1] = 'T';
	return(res2 - res1 -1);
}
 
/**************************************************************************
								RESNAMES()
**************************************************************************/
/*
 *  resnames() - build boolean array of residues w/ names in input
 *
 *	PUCKER NUCLEIC residue_names
 */

 static int
resnames(char *arr, parmstruct *prm)
{
	int	i, k, len, found, total = 0;

	for (;  tok[0] != ';'; ) {

		len = strlen(tok);
		if (len > 4)
			inerr("PUCKER NUCLEIC resname: residue name too long: ", 
									tok);
		/*
		 *  pad resname w/ blanks if necc
		 */

		for(k=len; k<4; k++)
			tok[k] = ' ';

		/*
		 *  tok[] is current res name: collect residues
		 */

		found = 0;
		for (i=0; i<prm->Nres; i++) {
#ifdef DBG2
			printf("(%.4s)(%.4s)\n", tok, &prm->ResNames[i*4]);
#endif
			if (!strncmp(tok, &prm->ResNames[i*4], 4)) {
				found++;
				arr[i] = 'T';
			}
		}
		if (!found)
			inerr("PUCKER NUCLEIC: no residues called: ",tok);
		total += found;
		gettok();
	}
	return(total);
}

/**************************************************************************
								SETSUGAR()
**************************************************************************/
/*
 *  setsugar() - given residue, set up pointers to sugar atoms
 *	(final step of NUCLEIC option)
 */

static void
setsugar(_REAL **points, int res, parmstruct *prm, _REAL *crd)
{
	char	*atomnames = &prm->AtomNames[4*(prm->Ipres[res]-1)];
	_REAL	*atom1 = &crd[3*(prm->Ipres[res]-1)];
	int	i, natres = prm->Ipres[res+1] - prm->Ipres[res];
	char	restring[80];
#ifdef DBG2
	_REAL	**p3, **p2 = points;
	printf("res %d  1st atom %d\n", res+1, prm->Ipres[res]);
	printf(" (%.4s)(%.4s)\n", &atomnames[4*6],  &atomnames[4*3]);
#endif
	if (natres < 10)
		inerr("PUCKER NUCLEIC: residue too small to have a sugar","");

	sprintf(restring, "in residue # %d (%.4s)", res, &prm->ResNames[res*4]);

	/*
	 *  set points to sugar atoms: for atoms that have relatively
	 *	fixed order in the standard all.in/uni.in dbases, first 
	 *	check standard locations then do exhaustive search
	 */

	if (!strncmp(&atomnames[4*6], "O4' ", 4)) {
		/* all atom */
		*points = atom1 + (3 * 6);
	} else if (!strncmp(&atomnames[4*3], "O4' ", 4)) {
		/* united atom */
		*points = atom1 + (3 * 3);
	} else {
		for (i=0; i<natres; i++) {
			if (!strncmp(&atomnames[4*i], "O4' ", 4)) {
				*points = atom1 + (3 * i);
				break;
			}
		}
		if (i == natres)
			inerr("PUCKER NUCLEIC: can't find O4'", restring);
	}
	points++;
	if (!strncmp(&atomnames[4*5], "C1' ", 4)) {
		/* all atom */
		*points = atom1 + (3 * 5);
	} else if (!strncmp(&atomnames[4*4], "C1' ", 4)) {
		/* united atom */
		*points = atom1 + (3 * 4);
	} else {
		for (i=0; i<natres; i++) {
			if (!strncmp(&atomnames[4*i], "C1' ", 4)) {
				*points = atom1 + (3 * i);
				break;
			}
		}
		if (i == natres)
			inerr("PUCKER NUCLEIC: can't find C1'","");
	}
	points++;

	/*
	 *  C2' & C3' locs depend on (all-atom/united, DNA/RNA, residue), so
	 *	just do exhaustive search
	 */

	for (i=natres-1; i>-1; i--) {
		if (!strncmp(&atomnames[4*i], "C2' ", 4)) {
			*points = atom1 + (3 * i);
			break;
		}
	}
	if (i == -1)
		inerr("PUCKER NUCLEIC: can't find C2'","");
	points++;
	for (i=natres-1; i>-1; i--) {
		if (!strncmp(&atomnames[4*i], "C3' ", 4)) {
			*points = atom1 + (3 * i);
			break;
		}
	}
	if (i == -1)
		inerr("PUCKER NUCLEIC: can't find C3'","");
	points++;
	if (!strncmp(&atomnames[4*4], "C4' ", 4)) {
		/* all atom */
		*points = atom1 + (3 * 4);
	} else if (!strncmp(&atomnames[4*2], "C4' ", 4)) {
		/* united atom */
		*points = atom1 + (3 * 2);
	} else {
		for (i=0; i<natres; i++) {
			if (!strncmp(&atomnames[4*i], "C4' ", 4)) {
				*points = atom1 + (3 * i);
				break;
			}
		}
		if (i == natres)
			inerr("PUCKER NUCLEIC: can't find C4'","");
	}
#ifdef DBG2
	points = p2;
	for (i=0; i<5; i++) {
		printf(" atom %d   loc %ld->%ld\n", i, points, *points);
		points++;
	}
#endif
}

/**************************************************************************
								POINTPUCK()
**************************************************************************/
/*
 *  pointpuck() - PUCKER id npoints point1 point2 ...: single pucker
 */

static void
pointpuck(idstruct *pkptr, int npoints)
{
	idstruct	*idptr;
	int		i, atom, grptype;

	/*
	 *  allocate space
	 */

	pkptr->idunion.puck.npuckers = 1;
	pkptr->idunion.puck.npoints = npoints;
	pkptr->idunion.puck.points = (_REAL **) get(npoints * sizeof(_REAL *));
	pkptr->idunion.puck.results = (_REAL *)
			get(sizeof(_REAL) * (npoints-3));

	gettok();
	for (i=0; i<npoints; i++) {
		if (tokint(&atom)) {
#ifdef DBG2
			printf("$$atom # %d\n", atom);
#endif
			pkptr->idunion.puck.points[i] = setatom(atom);
			gettok();
		} else if ((idptr = getstruct(tok)) != NULL) {
#ifdef DBG2
			printf("$$group\n");
#endif
			if (idptr->id_type != GROUP)
				inerr("PUCKER points:",
					"expected atom or group id");
			
			grptype = grouptype(idptr);

			/*
			 * set binding comes w/ group
			 */

			if (grptype == GROUPCENTER)
				pkptr->idunion.puck.points[i] = 
					idptr->idunion.group.geomcenter;
			else if (grptype == GROUPCOM)
				pkptr->idunion.puck.points[i] = 
					idptr->idunion.group.masscenter;
			else
				inerr("PUCKER points", "expected group center");
			gettok();
 		} else {
#ifdef DBG2
			printf("$$is it an atomname\n");
#endif
			/*
			 * last chance: see if it's an atom name
			 */
 
			pkptr->idunion.puck.points[i] = getsetatom();
			if (pkptr->idunion.puck.points[i] == NULL)
				inerr("PUCKER points: not a point??", "");
		}
	}
	if (tok[0] != ';')
		inerr("PUCKER points: expected ';'", "");
}

/**************************************************************************
							CALCPUCKER()
**************************************************************************/
/*
 *  calcpucker() - calc pucker(s) on coords
 */

void
calcpucker(idstruct *idptr)
{
	int		npuck = idptr->idunion.puck.npuckers;
	int		npoints = idptr->idunion.puck.npoints;
	_REAL		**points = idptr->idunion.puck.points;
	_REAL		*scratch= idptr->idunion.puck.scratch;
	_REAL		*results = idptr->idunion.puck.results;
	statstruct	*stats = idptr->idunion.puck.stats;
	int		i, j, k, nucleic = 0,ibox;
#ifdef DBG
	printf("calc pucker %s\n", idptr->id_name);
#endif
        if (idptr->id_subtype == NUCLEIC)
		nucleic++;

        idptr->idunion.puck.cases++;

	for (i=0; i<npuck; i++) {

		/*
		 *  copy points into scratch
		 */

		for (j=0,k=0; j<npoints; j++) {
#ifdef DBG2
			printf(" atom %d loc %ld->%ld\n", j, points, *points);
			fflush(stdout);
#endif
			scratch[k++] = (*points)[0];
			scratch[k++] = (*points)[1];
			scratch[k++] = (*points)[2];
			points++;
		}

		/*
		 *  calc pucker
		 */

		singlepucker(npoints, scratch, nucleic, results, stats);

/*
*/
                if(results[1] < 0.) results[1] += 360.;
                if(results[1] > 360.) results[1] -= 360.;
/*
                if(cos(TORAD*results[1]) < 0.) results[1] += 180.;
                if(results[1] > 180.) results[1] -= 360.;
*/
 printf("result%d: %f %f\n",i,results[0],results[1]);
/* printf("theta%s: %e\n",idptr->id_name,results[0]/results[2]); */

                if(results[1]>=idptr->idunion.puck.min && 
                   results[1]<=idptr->idunion.puck.max) {
                 ibox=(int)((results[1] -idptr->idunion.puck.min)/
                              idptr->idunion.puck.scale);
                 idptr->idunion.puck.boxes[ibox]++;
                }

		/*
		 *  prepare for next result
		 */

		results += npoints - 3;
		stats += npoints - 3;
	}
}

/**************************************************************************
							SINGLEPUCKER()
**************************************************************************/
/*
 *  singlepucker() - calc pucker for a single ring
 *
 *	Input:
 *		npoints		- number of points in ring
 *		arr		- pointer to disposable array of points
 *		nucleic		- flag to add 90 degrees to angle 
 *				  (Sundaralingham convention)
 *	Output:
 *		results		- pointer to calculated puckering coordinates:
 *				  results[0] - the first amplitude
 *				  results[1] - the first phase angle
 *				  results[2] - the second amplitude
 *				  results[3] - the second phase angle
 *					...[npoints-3]
 *		stats		- stat/Cstat/stat/Cstat.. matching results
 */

static void
singlepucker(int npoints, _REAL *arr, int nucleic, 
             _REAL *results, statstruct *stats)
{
	int    		i, j, m;
	_REAL		*rr, rc[3], r1[3], r2[3], norm[3];
	_REAL		sum1, sum2;
	double		A;
	_REAL		sinA, cosA;
	_REAL		sq2n, qfi, ampl, phi, tmp_fl;

#ifdef DBG2
	_REAL distance();
	for (i=0; i<3*npoints; i+=3) 
		printf("%f  %f  %f\n", arr[i], arr[i+1], arr[i+2]);
	for (i=0; i<3*(npoints-1); i+=3) 
		printf("d %f\n", distance(&arr[i], &arr[i+3]));
	printf("d %f\n", distance(&arr[3*(npoints-1)], &arr[0]));
#endif

	for (i=0; i<3; i++) {
		rc[i] = 0;
		r1[i] = 0;
		r2[i] = 0;
	}
	sq2n = sqrt(2.0 / (double) npoints);

	/* 
	 *  Calculate center of geometry & move center to origin
	 */

	for(i=0; i<3*npoints; i+=3) 
     		dvector(rc, rc, arr+i) ;
	for(i=0; i<3; i++)
		rc[i] = -rc[i] / npoints ;

	for(i=0; i<3*npoints; i+=3)
		dvector(arr+i, arr+i, rc);

	/* 
	 *  Calculate R', R"
	 */

	rr = arr;
	for(i=0; i<npoints; i++) {
		A = (TWOPI * i) / npoints;
		sinA = sin(A);
		cosA = cos(A);

		r1[0] += *rr * sinA;
		r2[0] += *rr * cosA;
		rr++;
		r1[1] += *rr * sinA;
		r2[1] += *rr * cosA;
		rr++;
		r1[2] += *rr * sinA;
		r2[2] += *rr * cosA;
		rr++;
	}

	/* 
	 *  Calculate the unit vector normal to the mean plane 
	 */

	cross(norm, r1, r2);
	tmp_fl = size(norm);
	for (i = 0; i < 3; i++)
		norm[i] /= tmp_fl;

	/* 
	 *  Calculate new z coords (effectively rotating normal to z axis)
	 */

	for(rr=arr,i=2; i<3*npoints+2; i+=3,rr+=3)
		arr[i] = dot(rr, norm) ;
	
	/* 
	 *  Calculate the pucker coordinate q , fi 
	 */

	for (m=2; m<=(npoints-1)/2; m++) {
		sum1 = 0.0;
		sum2 = 0.0;
		for(i=0,j=2; i<npoints; i++,j+=3) {
			A = (TWOPI * m * i) / npoints;
			sum1 += arr[j] * cos(A);
			sum2 -= arr[j] * sin(A);
		} /* make q(m)*cos fi(m) and q(m)*sin fi(m) */
		qfi = sqrt(sum1 * sum1 + sum2 * sum2);

		/*
		 *  amplitude
		 */

		ampl = sq2n * qfi;
		results[0 + (m - 2) * 2] = ampl;
		updatestat(stats, ampl);
		stats++;

		/*
		 *  angle
		 */

		phi = asin(sum2 / qfi);
		if (sum1 < 0.0)
			phi = _PI - phi;
		else if (phi < 0.0)
			phi = TWOPI + phi;
		if (nucleic) {

			/*
			 *  adjust to Sundaralingham convention
			 */
			
			phi += HALFPI;
		}
		results[1 + (m - 2) * 2] = phi * TODEG;
		updateCstat((Cstatstruct *) stats, phi);
		stats++;
#ifdef DBG
		printf("ampl %f  angle %f\n", ampl, phi * TODEG);
#endif
	}

	/* 
	 *  If the number of atoms is even, calculate last amplitude 
	 */

	if (!(npoints & 1)) {
		sum1 = 0.0;
		for (i=0; i<npoints; i++) {
			if (i & 1)	/* odd atoms */
				sum1 -= arr[2 + i * 3];
			else		/* even atoms */
				sum1 += arr[2 + i * 3];
		}
		ampl = sum1 / sqrt(npoints);
		results[npoints - 4] = ampl;
		updatestat(stats, ampl);
#ifdef DBG
		printf("ampl %f\n", ampl);
#endif
	}
}

void
writepucker(idstruct *dptr)
{
        int     i, nboxes = dptr->idunion.puck.nboxes;
        int     *boxes = dptr->idunion.puck.boxes;
        int     states = dptr->idunion.puck.cases;
        FILE    *file = dptr->idunion.puck.file;
        _REAL   scale = dptr->idunion.puck.scale, min=dptr->idunion.puck.min;

        _REAL   ri,distrib,fac;

        fac = 1./ ((_REAL)states);

        for (i=0; i<nboxes; i++) {
          ri = min+(_REAL)i*scale+0.5*scale;
          distrib = (_REAL)boxes[i]*fac;
          fprintf(file, "%e %e\n", ri,distrib);
        }
        fclose(file);
}
