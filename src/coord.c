#include "head.h"
#include <ctype.h>
 
/*  coord.c - routines for output of coordinate sets */

// extern  char            tok[];
// char            tok[];
extern  idstruct        *firstser;
extern	int		it;
 
extern idstruct		*startid(), *getstruct();
extern char		*get();
extern FILE		*newpath();
extern selectstruct	*selectx();
extern int		tokint();
 
/*  forward declarations */

/***********************************************************************
 							INITCOORD()
************************************************************************/
 
/*
 * initcoord() - set up coord idstruct
 */
 
void
initcoord()
{
	idstruct	*crdptr;
	char		id[TOKLEN];
	int		gotfmt = 0;
	int		append = 0;
 
	/*
	 *  get crd record in idlist
	 */
 
	crdptr = startid(CRD, id);

	/*
	 *  initialize crd record
	 */

	gettok();
	if (iskeywd())
		inerr("Expected file name", "");
	strcpy(crdptr->idunion.crd.name, tok);
	crdptr->idunion.crd.fmt = CRD;
	crdptr->idunion.crd.blank = 0;
 
	crdptr->idunion.crd.mod = 0;
	crdptr->idunion.crd.select = NULL;
	crdptr->idunion.crd.selmod = 0;
	crdptr->idunion.crd.modsel = 0;
	crdptr->idunion.crd.option = 0;
	crdptr->idunion.crd.nsel = 0;
	crdptr->idunion.crd.conv = 0;
 
	crdptr->idunion.crd.set = NULL;
 
	/*
	 *  get options
	 */
 
	for (;;) {
		gettok();
		if (tok[0] == ';') {
			break;
		} else if (!strcmp(tok, "APPEND")) {
			append = 1;
		} else if (!strcmp(tok, "BLANK")) {
			crdptr->idunion.crd.blank = 1;
		} else if (!strcmp(tok, "RST")) {
			gotfmt++;
			crdptr->idunion.crd.fmt = RST;
		} else if (!strcmp(tok, "REVERSE_FLOAT")) {
			crdptr->idunion.crd.conv = 1;
		} else if (!strcmp(tok, "RIG")) {
			gotfmt++;
			crdptr->idunion.crd.fmt = RIG;
		} else if (!strcmp(tok, "PDB")) {
			gotfmt++;
			crdptr->idunion.crd.fmt = PDB;
		} else if (!strcmp(tok, "MOL")) {
			gotfmt++;
			crdptr->idunion.crd.fmt = MOL;
		} else if (!strcmp(tok, "MOVIE")) {
			gotfmt++;
			append = 1;
			crdptr->idunion.crd.fmt = MOVIE;
		} else if (!strcmp(tok, "XYZ")) {
			gotfmt++;
			crdptr->idunion.crd.fmt = XYZ;
		} else if (!strcmp(tok, "AMBER")) {
			gotfmt++;
			crdptr->idunion.crd.fmt = AMBER;
		} else if (!strcmp(tok, "STDF")) {
			gotfmt++;
			crdptr->idunion.crd.fmt = STDF;
		} else if (!strcmp(tok, "ORAC4_DUMP")) {
			gotfmt++;
			crdptr->idunion.crd.fmt = ORAC4_DUMP;
		} else if (!strcmp(tok, "STDD")) {
			gotfmt++;
			crdptr->idunion.crd.fmt = STDD;
		} else {
			inerr("COORD: unexpected token","");
		}
	}
	if (gotfmt > 1) 
		inerr("COORD: more than 1 format given","");
 
/*
	if ((crdptr->idunion.crd.blank  ||  append)  &&
					(crdptr->idunion.crd.fmt != CRD ||
					 crdptr->idunion.crd.fmt != MOVIE)) 
		inerr("COORD: BLANK or APPEND given w/ non-CRD format","");
*/
 
	if (crdptr->idunion.crd.fmt == AMBER || 
            crdptr->idunion.crd.fmt == XYZ ||
            crdptr->idunion.crd.fmt == MOVIE ||
            crdptr->idunion.crd.fmt == ORAC4_DUMP ||
	    crdptr->idunion.crd.fmt == STDF)                 {
 		int	exists;

		/*
		 * start file now  TODO: wait till all input read ok
	 	 */
 
		if (fileprob(crdptr->idunion.crd.name, append, &exists, 
								"COORD"))
			exit(1);
 
		if ((crdptr->idunion.crd.file = 
				fopen(crdptr->idunion.crd.name, "a+")) == NULL){
			perror(crdptr->idunion.crd.name);
			exit(1);
		}
/*
*/
		if (!exists)
		 if (crdptr->idunion.crd.fmt == AMBER) 
			fprintf(crdptr->idunion.crd.file, "MDANAL dump\n");
/* print header of car file */
		if (crdptr->idunion.crd.fmt == MOVIE) {
                 fprintf(crdptr->idunion.crd.file,"!BIOSYM archive 1\nPBC=OFF\n");
                 fflush(crdptr->idunion.crd.file);
                }


	}
}
 
/***********************************************************************
                                                        SETUPCOORD()
************************************************************************/
 
/*
 * setupcoord() - look up coord id in idqueue and fill in struct
 */
 
void
setupcoord()
{
	idstruct        *crdptr, *idptr;
	int		i,k,Idum;
        extern parmstruct *prm;
 
	/*
	 * get the coord struct
	 */
 
	gettok();
	if ((crdptr = getstruct(tok)) == NULL)
		inerr("OUTPUT COORD: no such id:", tok);
	if (crdptr->id_type != CRD)
		inerr("OUTPUT COORD: id isn't coord:", tok);
 
	if (crdptr->id_setup)
		inerr("OUTPUT COORD: id multiply referenced in OUTPUT","");
 
	/*
	 *  see if explicit crd set reference
	 */
 
	gettok();
	if ((idptr = getstruct(tok)) == NULL) {
		idptr = firstser;
		printf("COORD %s: using default stream %s\n", 
				crdptr->id_name, idptr->id_name);
	} else {
		gettok();
	}
 
	/* TODO: fix up RMS; make checkatoms not prm dep for rms' sake */
	/*
	 *  set up the set pointer 
	 */
 
	if (idptr->id_type == STREAM) {
		crdptr->idunion.crd.set = idptr->idunion.ser.set;
                crdptr->idunion.crd.prm = idptr->idunion.ser.prm;
                if(idptr->idunion.ser.iweight)
                {
                 crdptr->idunion.crd.weight = idptr->idunion.ser.w->idunion.serw.w;
                } else {
                 crdptr->idunion.crd.weight = 1.;
                }

	} else if (idptr->id_type == RMS) {
                crdptr->idunion.crd.set = firstser->idunion.ser.set;
                crdptr->idunion.crd.prm = firstser->idunion.ser.prm;

		/*
		 *  must be an RMS FIT id
		 */

		if (!idptr->idunion.rms.fit)
			inerr("COORD: RMS must be a FIT type", "");

		/*
		 *  make a setstruct for now TODO: make rms use its own setstruc
		 */

		crdptr->idunion.crd.set = (setstruct *) get(sizeof(setstruct));
		crdptr->idunion.crd.set->crd = idptr->idunion.rms.crd2a;
#ifdef DBG
		printf("crdptr %d  set %d  crd2a %d\n", 
		  crdptr, crdptr->idunion.crd.set, idptr->idunion.rms.crd2a);
#endif
		crdptr->idunion.crd.set->max = idptr->idunion.rms.natoms;
		crdptr->idunion.crd.set->prm = 
			idptr->idunion.rms.grp->idunion.group.set->prm;

	} else
		inerr("OUTPUT COORD set: id must be for STREAM", "");

	/*
	 *  set default # of atoms to print (all in set)
	 */
        crdptr->idunion.crd.natoms = crdptr->idunion.crd.set->max;

 
	if (crdptr->idunion.crd.fmt == RIG) {
         crdptr->idunion.crd.setr = (setstructr *) get(sizeof(setstructr));
         crdptr->idunion.crd.setr->crd = (_REAL *)
          get(7*(crdptr->idunion.crd.prm->Nmol)*sizeof(_REAL));
        }

	 /*  make sure set has parms if output is PDB  TODO: allow names
	 *					to come from input PDB
	 */
 
	if ((crdptr->idunion.crd.fmt == PDB || 
              crdptr->idunion.crd.fmt == MOVIE ||
              crdptr->idunion.crd.fmt == MOL

            )
				&&  crdptr->idunion.crd.prm == NULL)
	    inerr("OUTPUT COORD: set needs a parm file to output PDB or MOVIE", "");
 
	/*
	 * get options
	 */
 
/*
	printf("COORD %s: getting options\n", crdptr->id_name);
*/
	for (;;) {
		if (tok[0] == ';') {
/*
			printf("COORD %s: done\n", crdptr->id_name);
*/
			crdptr->id_setup++;
                        break;
/*
			return;
*/
		} else if (!strcmp(tok, "ATOM")) {
			crdptr->idunion.crd.option = ATOM;
                        if (idptr->id_type == STREAM) {
                                crdptr->idunion.crd.natoms =
                                        checkatoms(idptr->idunion.ser.set->prm,

                                                idptr->idunion.ser.set->max);
                                if (crdptr->idunion.crd.natoms
                                        > idptr->idunion.ser.set->max)
                                        inerr("ATOM: fewer atoms in STREAM ",
                                                "than requested for output");
                        } else {
                                /* TODO: make sense of this (RMS FIT) */
                                crdptr->idunion.crd.natoms =
                                        checkatoms(crdptr->idunion.crd.set->prm
,
                                                idptr->idunion.rms.natoms);
                                if (crdptr->idunion.crd.natoms
                                        > idptr->idunion.rms.natoms)
                                        inerr("ATOM: fewer atoms in STREAM ",
                                                "than requested for output");
                        }
/*
			printf("COORD %s: ATOM\n", crdptr->id_name);
*/
			continue;
		} else if (!strcmp(tok, "SELECT")) {
			crdptr->idunion.crd.select = selectx();
		} else if (!strcmp(tok, "MOD")) {
			mod(&crdptr->idunion.crd.mod, 
				crdptr->idunion.crd.select);
		} else 
			inerr("COORD: unexpected token","");
		gettok();
	}

}
 
 
/************************************************************************
							WRITECRD()
************************************************************************/
 
/*
 * writecrd() - if appropriate set, write coords to file or add for average
 */
 
void
writecrd(idstruct *idptr, idstruct *lastpass)
{
	_REAL 		*ixyz;
	int		i, j, k, nxyz, box,nitems;
	char		path[120];
	FILE		*file;
	int		xmod, nsel, code, item, selmod;
        float           boxc_s[3];
	selectstruct	*sel;
	char    atnam[5];
 
	/*
	 * get some long-name values in short-name vars for convenience
	 */
 
	xmod = idptr->idunion.crd.mod;
	sel = idptr->idunion.crd.select;
	ixyz = idptr->idunion.crd.set->crd;
 	nxyz = 3 * idptr->idunion.crd.natoms + 1;


/*
	if (lastpass) {
		printf("COORD %s: last pass\n", idptr->id_name);
			if (idptr->idunion.crd.box)
				 nxyz += 9;
			ixyz--;
*/
/* turns off INH2O */
/*
			idptr->idunion.crd.option = 0; 
			goto writeout;
		return;
	}
*/

	/*
	 * see if there's a reason for not writing/averaging this set
	 */
 
	if (sel != NULL) {
		nsel = idptr->idunion.crd.nsel;
		code = sel[nsel].code;
		item = sel[nsel].item;
		selmod = idptr->idunion.crd.selmod;
	}
 
	if (xmod > 0) {
 
		/*
		 * mod was defined first
		 */
 
  		/* if (it % xmod) */
  		if ((it+1) % xmod)
			return;	
 
		if (sel != NULL) {
 
			/*
			 * apply select
			 */
 
			switch (code) {
				case SEL:
					if (item != selmod)
						return;
					idptr->idunion.crd.nsel++;
					break;
				case SELTILL:
					if (selmod == item) 
						idptr->idunion.crd.nsel++;
					break;
				case SELTILLEND:
					break;
				case SELDONE:
					return;
			}
			idptr->idunion.crd.selmod++;
 
		}
	} else if (sel != NULL) {
	
		/*
		 * select was defined first
		 */
 
		switch (code) {
			case SEL:
				if (item != it+1)
					return;
				idptr->idunion.crd.nsel++;
				break;
			case SELTILL:
				if (it+1 == item) 
					idptr->idunion.crd.nsel++;
				break;
			case SELTILLEND:
				break;
			case SELDONE:
				return;
		}
		if (xmod) {
 
			/*
			 * apply mod
			 */
 
			if (idptr->idunion.crd.modsel++ % xmod)
				return;
		}
		if (idptr->idunion.crd.weight>=-10000.) {
 
			/*
			 * apply weight
			 */
			/* return; */
		}
	}
 
	/*
	 * set passed all tests: print it
	 */
 
writeout:
	/*
	 * write appropriate format
	 */
 
	switch (idptr->idunion.crd.fmt) {
	case AMBER:
		file = idptr->idunion.crd.file;
 
		for (i=1; i<nxyz; i++) {
			fprintf(file, "%8.3f", *ixyz);
			ixyz++;
			if (!(i%10))
			// if (!(i%3))
				fprintf(file, "\n");
/*
                        else if(i==nxyz-1)
				fprintf(file, "\n");
*/
		}
 
		if (--i%10)
			fprintf(file, "\n");
                 for(i=0;i<3;i++) {
                  boxc_s[i]= (float)idptr->idunion.crd.set->boxc[i];
                 // boxc_s[i]=100.;
                 }
 		 fprintf(file, "%8.3f%8.3f%8.3f\n",
                  boxc_s[0],boxc_s[1],boxc_s[2]);
/*
			fprintf(file, "%8.3f%8.3f%8.3f\n",
				idptr->idunion.crd.set->crd[box],
				idptr->idunion.crd.set->crd[box+1],
				idptr->idunion.crd.set->crd[box+2]);
*/
		if (idptr->idunion.crd.blank)
			fprintf(file, "\n");
		break;
	case XYZ:
		file = idptr->idunion.crd.file;
		fprintf(file, "%d\n", idptr->idunion.crd.natoms);
		fprintf(file, "Dump of XYZ by mdan\n");
		for (i=0; i<idptr->idunion.crd.natoms; i++) 
		{
 		 if (idptr->idunion.crd.prm->AtomNames[4*i+3] == ' ')
		  sprintf(atnam, " %.3s", &idptr->idunion.crd.prm->AtomNames[4*i]);
 		 else 
		  sprintf(atnam, " %.4s", &idptr->idunion.crd.prm->AtomNames[4*i]);
 		 boxc_s[0]=*ixyz++; boxc_s[1]=*ixyz++; boxc_s[2]=*ixyz++;
		 fprintf(file, "%-.4s ",atnam); 
		 fprintf(file, "%10.5f %10.5f %10.5f\n", 
		  boxc_s[0], boxc_s[1], boxc_s[2]);
		}
                break;
	case STDF:
		file = idptr->idunion.crd.file;
                nitems=idptr->idunion.crd.natoms;
                if(idptr->idunion.crd.conv)
                 reverse4ByteOrder(1,(char *)&nitems);
                fwrite(&nitems,sizeof(int),1,file);
                fflush(file);
                for(i=0;i<3;i++) {
                 boxc_s[i]=
                  (float)idptr->idunion.crd.set->boxc[i];
                }
                if(idptr->idunion.crd.conv)
                 reverse4ByteOrder(3,(char *)boxc_s);
                fwrite(boxc_s,3*sizeof(float),1,file);
                fflush(file);
                for(i=0;i<3*idptr->idunion.crd.natoms;i++) {
                 idptr->idunion.crd.set->crd_s[i]=
                  (float)idptr->idunion.crd.set->crd[i];
                }
                if(idptr->idunion.crd.conv)
                 reverse4ByteOrder(3*idptr->idunion.crd.natoms,
                  (char *)idptr->idunion.crd.set->crd_s);
                fwrite(idptr->idunion.crd.set->crd_s,
                 3*idptr->idunion.crd.natoms*sizeof(float),
                  1,file);
                fflush(file);
                break;

	case ORAC4_DUMP:
		file = idptr->idunion.crd.file;
                nitems=idptr->idunion.crd.natoms;
                { float dt,co_s[9];
                  int byteso,bytes,ibox=0;
                dt=(float)it;

                byteso=sizeof(float)+ibox*sizeof(float)+3*nitems*sizeof(float);
                bytes=byteso;
                if(idptr->idunion.crd.conv)
                 reverse4ByteOrder(1,(char *)&bytes);
                fwrite(&bytes, sizeof(int),1,file);
                fflush(file);

                if(idptr->idunion.crd.conv)
                 reverse4ByteOrder(1,(char *)&dt);
                fwrite(&dt,
                 sizeof(float),1,file);

                /* if(ibox) */
                {
                for(i=1,k=0;i<=3;i++) {
                 for(j=1;j<=3;j++,k++) {
                  co_s[k]=
                   (float)idptr->idunion.crd.set->cell_to_cart[i][j];
                 }
                }
                if(idptr->idunion.crd.conv)
                 reverse4ByteOrder(9,(char *)co_s);
                 fwrite(co_s,9*sizeof(float),1,file);
                }

                for(i=0;i<3*idptr->idunion.crd.natoms;i++) {
                 idptr->idunion.crd.set->crd_s[i]=
                  (float)idptr->idunion.crd.set->crd[i];
                }
                if(idptr->idunion.crd.conv)
                 reverse4ByteOrder(idptr->idunion.crd.natoms,
                  (char *)idptr->idunion.crd.set->crd_s);
                fwrite(idptr->idunion.crd.set->crd_s,
                 3*idptr->idunion.crd.natoms*sizeof(float),
                  1,file);

                fwrite(&bytes, sizeof(int),1,file);
                }
		break;

	case STDD:
		if ((file = newpath(path, idptr->idunion.crd.name)) == NULL)
			break;
                fwrite(&idptr->idunion.crd.natoms,
                 sizeof(int),1,file);
                fwrite(idptr->idunion.crd.set->boxc,3*sizeof(_REAL),1,file);
                fwrite(idptr->idunion.crd.set->crd,
                 3*idptr->idunion.crd.natoms*sizeof(_REAL),
                  1,file);
                fflush(file);
		break;
	case RST: 
		if ((file = newpath(path, idptr->idunion.crd.name)) == NULL)
			break;
		fprintf(file, "CARNAL temporary title\n%4d\n",
					idptr->idunion.crd.natoms);
		for (i=1; i<nxyz; i++) {
			fprintf(file, "%12.7f", *ixyz);
			ixyz++;
			if (!(i%6))
				fprintf(file, "\n");
		}
		if (i%6)
			fprintf(file, "\n");
		fclose(file);
		break;
	case PDB:
		if ((file = newpath(path, idptr->idunion.crd.name)) == NULL)
			break;
		pdbgen(ixyz, idptr->idunion.crd.prm, file,
		 idptr->idunion.crd.natoms,idptr->idunion.crd.fmt);
		fclose(file);
		break;
	case MOVIE:
/*
		if ((file = newpath(path, idptr->idunion.crd.name)) == NULL)
			break;
		file = idptr->idunion.crd.file;
*/
		pdbgen(ixyz, idptr->idunion.crd.prm, idptr->idunion.crd.file,
		 idptr->idunion.crd.natoms,idptr->idunion.crd.fmt);
/*
		fclose(file);
*/
		break;
	default:
		printf("Programming error: crd.fmt %d\n", 
						idptr->idunion.crd.fmt);
		exit(1);
	}
}
 
/************************************************************************
							NEWPATH()
************************************************************************/
 
/*
 * newpath() - try to create a file using iteration as extension
 */
 
FILE *
newpath(char *path, char *base)
{
	FILE	*file;
	int	exists;
 
	sprintf(path, "%s.%d", base, it+1);
	if (fileprob(path, 0, &exists, "COORD")) {
		/*
		 * another chance:
		 */
		sprintf(path, "%s.a%d", base, it);
		if (fileprob(path, 0, &exists, "COORD")) {
			printf( "COORD: file %s.%d exists: not written\n",
								path, it);
			return(NULL);
		}
	}
	printf("starting new coord output file: %s\n", path);
	if ((file = fopen(path, "w")) == NULL) {
		perror(path);
		return(NULL);
	}
	return(file);
}
 
/************************************************************************
							CLOSECRD()
************************************************************************/
 
/*
 * closecrd() - close file if necc
 */
 
void
closecrd(idstruct *crdpt)
{
	if (crdpt->idunion.crd.fmt == CRD) 
		fclose(crdpt->idunion.crd.file);
}
