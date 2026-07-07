#include "head.h"

extern int	errno;

/***********************************
	global vars
************************************/

/****	default crd set for measurements   ****/

idstruct	*firstser = NULL;

/***********************************
	external vars
************************************/

/****	util.c	****/
// extern char	tok[];
// char	tok[];

/****	doit.c	****/
extern int	it;

/****   misc routines	****/

extern idstruct	*startid(), *getstruct();
extern FILE	*genopen();
// extern char	*get2(),*get();
// void *get2();
char	*get();
extern void reverse4ByteOrder();

FILE *ftmp;
_REAL lnZ,lnZf=2.504062664748843e+04,wmax=-5.843980e+00;

/***********************************************************************
 							SETUPSTREAM()
************************************************************************/

/*
 * setupstream() - set up to read STREAM of coordinates
 */

void setupstream()
{
        extern parmstruct *prm;
	idstruct	*serptr, *idptr;
	namestruct	*nameptr, *lastname = NULL;
	char		id[TOKLEN];
	int		i,gotfmt=0;
        void startser();

	/*
	 *  put stream record in idlist
	 */

	serptr = startid(STREAM, id);
 
	/*
	 *  if this is the 1st stream defined, it becomes the default one
	 */

	if (firstser == NULL)
		firstser = serptr;

	/*
	 *  initialize the struct
	 */

        serptr->idunion.ser.ifile = 0;
        serptr->idunion.ser.fmt = STDF;
        serptr->idunion.ser.conv = 0;
        serptr->idunion.ser.reconstr = 0;
	serptr->idunion.ser.iset = 0;
        serptr->idunion.ser.nfiles = 0;
        serptr->idunion.ser.names = NULL;
	serptr->idunion.ser.fileptr = NULL;
	serptr->idunion.ser.skip = 0;
	serptr->idunion.ser.skipa = 0;
	serptr->idunion.ser.skipb = 0;
	serptr->idunion.ser.nproc = 1;
	serptr->idunion.ser.shift = 0;
	serptr->idunion.ser.nstates = 10000000;
	serptr->idunion.ser.win = 0;
	serptr->idunion.ser.winct = 0;
	serptr->idunion.ser.ibox = 1;
        serptr->idunion.ser.scale=0;
        serptr->idunion.ser.iweight=0;
        serptr->idunion.ser.atom_selected=-1;

	serptr->idunion.ser.set = (setstruct *) get2(1,sizeof(setstruct));
 	serptr->idunion.ser.set->boxc = (_REAL *) get2(3,sizeof(_REAL));
 	serptr->idunion.ser.set->boxm1 = (_REAL *) get2(3,sizeof(_REAL));
 	serptr->idunion.ser.set->cell_to_cart = (_REAL **) matrix(1,3,1,3);
 	serptr->idunion.ser.set->cart_to_cell = (_REAL **) matrix(1,3,1,3);
	/*
	 *  get options; 
	 *	WIN allows applying {skip m, read n} forever
	 */

	for (;;) {
                gettok();
                if ((idptr = getstruct(tok)) != NULL) {
                 if (idptr->id_type != STREAMW)
                  inerr("STREAM: streamw id expected\n","");
                 serptr->idunion.ser.iweight=1;
                 serptr->idunion.ser.w=idptr;
		} else if (!strcmp(tok, "STATES")) {
                        gettok();
                        if (!tokint(&serptr->idunion.ser.nstates))
                          inerr("STREAM: expected nstates", "");
		} else if (!strcmp(tok, "RED_TO_REAL")) {
                        serptr->idunion.ser.scale=1;
		} else if (!strcmp(tok, "SKIPA")) {
                        gettok();
                        tokint(&serptr->idunion.ser.skipa);
		} else if (!strcmp(tok, "SKIPB")) {
                        gettok();
                        tokint(&serptr->idunion.ser.skipb);
		} else if (!strcmp(tok, "REAL_TO_RED")) {
                        serptr->idunion.ser.scale=2;
		} else if (!strcmp(tok, "REVERSE_FLOAT")) {
                        serptr->idunion.ser.conv=1;
		} else if (!strcmp(tok, "RECONSTR")) {
                        serptr->idunion.ser.reconstr=1;
		} else if (!strcmp(tok, "SELECT_ATOM")) {
                        gettok();
                        tokint(&serptr->idunion.ser.atom_selected);
                        serptr->idunion.ser.atom_selected-=1;
		} else if (!strcmp(tok, "BOX")) {
                        serptr->idunion.ser.ibox=0;
                        gettok();
                        if (!tokREAL(&serptr->idunion.ser.set->boxc[0]))
                          inerr("STREAM: expected boxx", "");
                        gettok();
                        if (!tokREAL(&serptr->idunion.ser.set->boxc[1]))
                          inerr("STREAM: expected boxy", "");
                        gettok();
                        if (!tokREAL(&serptr->idunion.ser.set->boxc[2]))
                          inerr("STREAM: expected boxz", "");
		} else if (!strcmp(tok, "FORM")) {
                        gettok();
		        if (!strcmp(tok, "APE")) {
                          serptr->idunion.ser.fmt = APE;
		        } else if (!strcmp(tok, "XYZ")) {
                          serptr->idunion.ser.fmt = XYZ;
		        } else if (!strcmp(tok, "AMBER")) {
                          serptr->idunion.ser.fmt = AMBER;
		        } else if (!strcmp(tok, "POS")) {
                          serptr->idunion.ser.fmt = POS;
		        } else if (!strcmp(tok, "PDB")) {
                          serptr->idunion.ser.fmt = PDB;
		        } else if (!strcmp(tok, "DISCOVER2")) {
                          serptr->idunion.ser.fmt = DISCOVER2;
		        } else if (!strcmp(tok, "DISCOVER3")) {
                          serptr->idunion.ser.fmt = DISCOVER3;
		        } else if (!strcmp(tok, "ORAC3_DUMP")) {
                          serptr->idunion.ser.fmt = ORAC3_DUMP;
                        } else if (!strcmp(tok, "DLPOLY")) {
                          serptr->idunion.ser.fmt = DLPOLY;
                        } else if (!strcmp(tok, "CHARMM")) {
                          serptr->idunion.ser.fmt = CHARMM;
                        } else if (!strcmp(tok, "AXSF")) {
                          serptr->idunion.ser.fmt = AXSF;
		        } else if (!strcmp(tok, "ORAC4_DUMP")) {
                          serptr->idunion.ser.fmt = ORAC4_DUMP;
		        } else if (!strcmp(tok, "STDD")) {
                          serptr->idunion.ser.fmt = STDD;
		        } else if (!strcmp(tok, "CBMC")) {
                          serptr->idunion.ser.fmt = CBMC;
		        } else if (!strcmp(tok, "MUCA")) {
                          serptr->idunion.ser.fmt = MUCA;
                          ftmp=fopen("w.out","w");
		        } else if (!strcmp(tok, "NWCHEM")) {
                          serptr->idunion.ser.fmt = NWCHEM;
		        } else if (!strcmp(tok, "NWCHEM_RST")) {
                          serptr->idunion.ser.fmt = NWCHEM_RST;
                        } else {
                         inerr("STREAM: fmt unknown", "");
                         exit(1);
                        }
                        gotfmt++;
		} else if (!strcmp(tok, "NPROC")) {
                        gettok();
                        if (!tokint(&serptr->idunion.ser.nproc))
                          inerr("STREAM: expected nproc", "");
		} else if (!strcmp(tok, "WIN")) {

			/*
			 *  window mode = {skip m, read n} forever
			 */

			if (serptr->idunion.ser.skip)
			  inerr("STREAM: 2 WIN statements?", "");
			gettok();
			if (!tokint(&serptr->idunion.ser.skip)) 
			  inerr("STREAM/WIN: expected # of sets to skip", "");
			if (serptr->idunion.ser.skip < 0)
			  inerr("STREAM/WIN: skip # needs to be >= 0","");
			serptr->idunion.ser.skip = -serptr->idunion.ser.skip;
			serptr->idunion.ser.winct = serptr->idunion.ser.skip;
			
			gettok();
			if (!tokint(&serptr->idunion.ser.win)) 
			  inerr("STREAM/WIN: expected # of sets in window", "");
			if (serptr->idunion.ser.win < 1)
			  inerr("STREAM/WIN: window # needs to be > 0","");
		} else
			break;
	}
        if (gotfmt > 1)
                inerr("STREAM: more than 1 format given","");

	serptr->idunion.ser.set->dir[0] = 0.;
	serptr->idunion.ser.set->dir[1] = 0.;
	serptr->idunion.ser.set->dir[2] = 1.;
        serptr->idunion.ser.prm = prm;
        serptr->idunion.ser.set->prm = prm;
        serptr->idunion.ser.set->max = prm->Nsit;
        serptr->idunion.ser.set->crd_s = (float *) 
               get2(3 * serptr->idunion.ser.prm->Nsit,
                sizeof(float));

        if(serptr->idunion.ser.fmt==DISCOVER2 ||
           serptr->idunion.ser.fmt==DISCOVER3) {
         serptr->idunion.ser.set->crd = (_REAL *) 
                    get2(3 * (serptr->idunion.ser.prm->Nmov) ,sizeof(_REAL));
 	 serptr->idunion.ser.set->Enerd = (double *) 
         get2(prm->Nener,sizeof(double));
         serptr->idunion.ser.set->crd_s = (float *) 
                    get2(3 * serptr->idunion.ser.prm->Nmov , sizeof(float));
        } else {
 	 serptr->idunion.ser.set->crd = (_REAL *)
		    get2(3 * (serptr->idunion.ser.prm->Nsit) , sizeof(_REAL));
        }

        serptr->idunion.ser.shift =
         (3*serptr->idunion.ser.prm->Nsit / serptr->idunion.ser.nproc);
 
	/*
	 * get file names, starting w/ current token
	 */

	for (;;) {

		/*
		 * check for done 
		 */

		if (tok[0] == ';') {
			if (!serptr->idunion.ser.nfiles)
				inerr("No files given in STREAM stmt, id=", id);
			break;
		} 

		/*
		 *  assume token is file name: check for kwd, increment count 
		 * 	& add file name to list
		 */

		if (iskeywd())
			inerr("STREAM: expected file name, got keyword: ", tok);
		serptr->idunion.ser.nfiles++;
		nameptr = (namestruct *) get2(1,sizeof(namestruct));
		/* allow for .Z extension to name */
		nameptr->name = get(strlen(tok) + 3);
		strcpy(nameptr->name, tok);
		if (lastname == NULL)
			serptr->idunion.ser.names = nameptr;
		else
			lastname->next = nameptr;
		lastname = nameptr;
		nameptr->next = NULL;

		gettok();
	}

	/*
	 * check files: can we open them & determine format?
	 */

	for (i=1; i<serptr->idunion.ser.nfiles; i++) {
		startser(serptr, i);
		genclose(serptr->idunion.ser.fileptr, serptr->idunion.ser.name);
	}

	/*
	 *  leave the 1st one started  
	 *	TODO: leave this to a second pass over the idlist
	 */

	startser(serptr, 0);

	serptr->id_setup++;
}


/***********************************************************************
 							STARTSER()
************************************************************************/

/*
 *  startser() - open the i'th file in an input stream and read past any header
 */

void startser(idstruct *serptr, int i)
{
	namestruct 	*nameptr;
	int		j = 0;
        char title[TOKLEN];

	/*
	 *  find i'th filename
	 */

	for (nameptr=serptr->idunion.ser.names; 
	     nameptr!=NULL; 
	     nameptr=nameptr->next) 
		if (j++ == i)
			break;

	if (nameptr == NULL) {
	  printf("Internal error: stream file name %d out of range\n", i);
	  exit(1);
	}

	/*
	 *  if something is already open, close it
	 */

/*
	if (serptr->idunion.ser.fileptr != NULL) 
		genclose(serptr->idunion.ser.fileptr,serptr->idunion.ser.name);
 printf("ho chiuso la seconda volta, in startser.\n");
*/

	/*
	 *  now open the file
	 */

	serptr->idunion.ser.name = nameptr->name;
	if ((serptr->idunion.ser.fileptr = genopen(nameptr->name, "stream"))
								 == NULL)
		exit(1);

/* read only binary files without skipping any header */
/*
*/
	if (serptr->idunion.ser.fmt == AMBER) { 
	 preadln(serptr->idunion.ser.fileptr, "", title);
	 printf("%s\n", title);
        }

if (serptr->idunion.ser.fmt == APE) {
/*
*/
   int Idum[8]; _REAL Xdum[5];
   for (i=0;i<8;i++) {
    fread(Idum,sizeof(int),1,serptr->idunion.ser.fileptr);
   }
   for (i=0;i<5;i++) {
    fread(Xdum,sizeof(float),1,serptr->idunion.ser.fileptr);
   }
   printf("Xdum: %e \n",Xdum[0]);
}
/*
   if (serptr->idunion.ser.fmt == DLPOLY) {
    preadln(serptr->idunion.ser.fileptr, "", title);
    fscanf(serptr->idunion.ser.fileptr,"%*d %*d ");
   }
*/
   else if (serptr->idunion.ser.fmt == CHARMM) {
    int bytes,ntitle,icntrl[100];
    int Cdum;

    fread(&bytes,sizeof(int),1,serptr->idunion.ser.fileptr);
    if(serptr->idunion.ser.conv) reverse4ByteOrder(1,(char *)&bytes);
    bytes/=4;
    fread(icntrl,bytes*sizeof(int),1,serptr->idunion.ser.fileptr);
    if(serptr->idunion.ser.conv) reverse4ByteOrder(bytes,(char *)icntrl);
    printf("CHARMM: nfile= %d\n",icntrl[0]);
    printf("CHARMM: npriv= %d\n",icntrl[1]);
    printf("CHARMM: nsavc= %d\n",icntrl[2]);
    printf("CHARMM: nstep= %d\n",icntrl[3]);
    fread(&bytes,sizeof(int),1,serptr->idunion.ser.fileptr);

    fread(&bytes,sizeof(int),1,serptr->idunion.ser.fileptr);
    if(serptr->idunion.ser.conv) reverse4ByteOrder(1,(char *)&bytes);
    bytes/=4;
    fread(icntrl,bytes*sizeof(int),1,serptr->idunion.ser.fileptr);
    fread(&bytes,sizeof(int),1,serptr->idunion.ser.fileptr);

    fread(&bytes,sizeof(int),1,serptr->idunion.ser.fileptr);
    fread(&bytes,sizeof(int),1,serptr->idunion.ser.fileptr);
    if(serptr->idunion.ser.conv) reverse4ByteOrder(1,(char *)&bytes);
    printf("CHARMM: natom= %d\n",bytes);
    fread(&bytes,sizeof(int),1,serptr->idunion.ser.fileptr);

   }
   else if (serptr->idunion.ser.fmt == AXSF) { 
    fscanf(serptr->idunion.ser.fileptr,"%*s %*d ");
    fscanf(serptr->idunion.ser.fileptr,"%*s ");
   }
   else if (serptr->idunion.ser.fmt == NWCHEM) { 
    { int natomt,nbondt,nbondwt,nwatt,it,jt;
     fscanf(serptr->idunion.ser.fileptr,"%*s ");
     fscanf(serptr->idunion.ser.fileptr,"%*d %d %d %d %d ",
        &natomt,&nbondt,&nbondwt,&nwatt);
//     printf("natomt= %d nbondt= %d nbondwt= %d nwatt= %d\n",
//        natomt,nbondt,nbondwt,nwatt);
     preadln(serptr->idunion.ser.fileptr, "", title);
     preadln(serptr->idunion.ser.fileptr, "", title);
     preadln(serptr->idunion.ser.fileptr, "", title);
     for(i=0;i<natomt;i++)
     {
      preadln(serptr->idunion.ser.fileptr, "", title);
     }
     for(i=0;i<nbondt+nbondwt;i++)
     {
      fscanf(serptr->idunion.ser.fileptr,"%d %d ",&it,&jt);
     }
    }
   }
   else if (serptr->idunion.ser.fmt == NWCHEM_RST) { 
    { int natomt,nwatt,it;
     preadln(serptr->idunion.ser.fileptr, "", title);
     preadln(serptr->idunion.ser.fileptr, "", title);
     preadln(serptr->idunion.ser.fileptr, "", title);
     fscanf(serptr->idunion.ser.fileptr,"%*s %*s %d %*s",&it);
     for(i=0;i<it;i++)
     {
      preadln(serptr->idunion.ser.fileptr, "", title);
     }
    }

   }

	serptr->idunion.ser.fileset = 0;
}

/***********************************************************************
							NEXTSER()
************************************************************************/

/*
 * nextser() - read in next set of coords in stream
 */

int
nextser(idstruct *serptr)
{
	FILE	*fileptr = serptr->idunion.ser.fileptr;
	_REAL	*r = serptr->idunion.ser.set->crd;
	float	*r_s = serptr->idunion.ser.set->crd_s;
	_REAL	*boxc = serptr->idunion.ser.set->boxc;
	_REAL	**cell_to_cart = serptr->idunion.ser.set->cell_to_cart;
	_REAL	**cart_to_cell = serptr->idunion.ser.set->cart_to_cell;
	double	*Enerd = serptr->idunion.ser.set->Enerd;
        double Etot=serptr->idunion.ser.set->Etot;
        double Epot=serptr->idunion.ser.set->Epot;
        double Ekin=serptr->idunion.ser.set->Ekin;
        double Press=serptr->idunion.ser.set->Press;
	int  	ibox = serptr->idunion.ser.ibox;
	int Nsit = serptr->idunion.ser.prm->Nsit;
	_REAL sig0r = serptr->idunion.ser.prm->sig0r;
	int nproc = serptr->idunion.ser.nproc;
	int shift = serptr->idunion.ser.shift;
	int fmt = serptr->idunion.ser.fmt;
        int i1,i2,i3,i4,i5,i6,nitems=0,ifin=0;
        int i,j,bytes,ip,kp,itest,nstep;
	_REAL Xdum[7],fstepd;
        double Ddum[1],co[9];
        float Fdum[1],boxc_s[3],co_s[9],fstep;
        int Idum[1];
        void startser(),rdhead();
        extern parmstruct *prm;
        char Cdum[1];

	printf("nextser it %d fileset %d\n", 
		it, serptr->idunion.ser.fileset);

	/*
	 *  try reading set / open next file if necc;
	 *	keep going as req'd if WINdow mode has been specified
	 */

	for (;;) {

		/*
		 *  try to get another set from current file
		 */

		serptr->idunion.ser.fileset++;

      i1=0;i2=0;i3=0;i4=0;i5=0;i6=0;
      switch (fmt)
      {
       case STDD:
/* read C binary file */
        nitems = Nsit;
        i1=1;
        if (Nsit != 0) {
         i1=fread(&nitems,sizeof(int),1,fileptr);
         if (nitems!=Nsit) {
          printf("Nsit in stream is different from the one in input\n");
          printf("Nsit: %d nitems %d\n",Nsit,nitems);
         } 
	} else i1=0;
        if (i1!=1) {
         printf("something wrong reading nitems\n");
        }
        if (nitems != 0 && ibox) {
         i3=fread(boxc,sizeof(_REAL),3,fileptr);
/*
*/
         printf("box:\n");
         for (i=0; i<3; i++) printf("%e ",boxc[i]); printf("\n");
	} else i3=0;
        if (i3!=3 && ibox) {
         printf("something wrong reading box\n");
        }
        if (nitems != 0) {
         i5+=fread(r,sizeof(_REAL),3*Nsit,fileptr);
	} else i5=0;
        if (i5!=3*Nsit) {
         printf("something wrong reading r: %d\n",i5);
        }
        itest=1+3*Nsit;
        if(ibox) itest+=3;
        ifin=i1+i3+i4+i5;
        break;

       case APE:
/* read C binary file for many processors in APE style */
        i1=0,i3=0,i5=0;
        nitems = Nsit;
        i1=1;
/*
 printf("shift,nproc: %d %d\n",shift,nproc);
 printf("i: %d %d %d %d %d\n",shift,nproc);
*/
        if (nitems != 0) {
         for(kp=0;kp<nproc;kp++) {
          i3+=fread(&boxc[0],sizeof(_REAL),1,fileptr);
          i3+=fread(&boxc[1],sizeof(_REAL),1,fileptr);
          i3+=fread(&boxc[2],sizeof(_REAL),1,fileptr);
         }
         if(serptr->idunion.ser.conv) reverse4ByteOrder(3,(char *)boxc);
/*
*/
         printf("box:\n");
          for (i=0; i<3; i++) printf("%e ",boxc[i]); printf("\n");
	} else i3=0;
        if (i3!=3*nproc) {
         printf("something wrong reading box\n");
        }
        if (nitems != 0) {
         for(ip=0;ip<shift;ip++) {
          for(kp=0;kp<nproc;kp++) {
           i5+=fread(&r[ip+kp*shift],sizeof(_REAL),1,fileptr);
          }
         }
	} else i5=0;
        if (i5!=3*Nsit) {
         printf("something wrong reading r: %d\n",i5);
        }
        if(serptr->idunion.ser.conv) reverse4ByteOrder(3*Nsit,(char *)r);
        ifin=i3+i5;
        itest=3*nproc+3*Nsit;
        break;

       case STDF:
/* read C binary file, single precision then convert into d */
        nitems = Nsit;
        i1=1;
        if (Nsit != 0) {
         i1=fread(&nitems,sizeof(int),1,fileptr);
         if(serptr->idunion.ser.conv) reverse4ByteOrder(1,(char *)&nitems);
         if (nitems!=Nsit) {
          printf("Nsit in stream is different from the one in input\n");
          printf("Nsit: %d nitems %d\n",Nsit,nitems);
         } 
	} else i1=0;
        if (i1!=1) {
         printf("something wrong reading nitems\n");
        }
        if (nitems != 0 && ibox) {
         i3=fread(boxc_s,sizeof(float),3,fileptr);
         if(serptr->idunion.ser.conv) reverse4ByteOrder(3,(char *)boxc_s);
         for (i=0; i<3; i++) boxc[i]=(_REAL)boxc_s[i];
/*
*/
         printf("box_s: ");
         for (i=0; i<3; i++) 
          printf("%e ",boxc[i]); printf("\n");
	} else i3=0;
        if (i3!=3 && ibox) {
         printf("something wrong reading box\n");
        }
        if (nitems != 0) {
         bytes=sizeof(float)*3*serptr->idunion.ser.skipb;
         fseek(fileptr,bytes,1);
         i5+=fread(r_s,sizeof(float),3*Nsit,fileptr);
         if(serptr->idunion.ser.conv) reverse4ByteOrder(3*Nsit,(char *)r_s);
         bytes=sizeof(float)*3*serptr->idunion.ser.skipa;
         fseek(fileptr,bytes,1);
         for (i=0; i<3*Nsit; i++) r[i]=(_REAL)r_s[i];
         printf("r0= %f %f %f\n",r[0],r[1],r[2]);
         printf("rf= %f %f %f\n",r[3*Nsit-3],r[3*Nsit-2],r[3*Nsit-1]);
	} else i5=0;
        if (i5!=3*Nsit) {
         printf("something wrong reading r_s: %d\n",i5);
        }

        itest=1+3*Nsit;
        if(ibox) itest+=3;
        ifin=i1+i3+i4+i5;
        break;

       case ORAC3_DUMP:
/* read ORAC3 dump binary file, fortran single precision then convert into d */
        nitems = Nsit;
        i1=1;
        fread(&bytes,sizeof(int),1,fileptr);
        if (Nsit != 0) {
         i1=fread(&nstep,sizeof(int),1,fileptr);
         if(serptr->idunion.ser.conv) reverse4ByteOrder(1,(char )nstep);
         printf("nstep= %d\n",nstep);
	} else i1=0;
        if (i1!=1) {
         printf("something wrong reading nstep\n");
        }
        i5+=fread(r_s,sizeof(float),3*Nsit,fileptr);
        bytes=sizeof(float)*3*serptr->idunion.ser.skipa;
        fseek(fileptr,bytes,1);
        if(serptr->idunion.ser.conv) reverse4ByteOrder(3*Nsit,(char *)r_s);
        for (i=0; i<3*Nsit; i++) r[i]=(_REAL)r_s[i];
        printf("r0= %f %f %f\n",r[0],r[1],r[2]);
        if (i5!=3*Nsit) {
         printf("something wrong reading r_s: %d\n",i5);
        }
        if (ibox) {
/*
         i3=fread(co_s,sizeof(float),9,fileptr);
         if(serptr->idunion.ser.conv) reverse4ByteOrder(9,(char *)co_s);
         boxc[0]=(_REAL)co_s[0];
         boxc[1]=(_REAL)co_s[4];
         boxc[2]=(_REAL)co_s[8];
*/
         i3=fread(co,sizeof(_REAL),9,fileptr);
         boxc[0]=co[0];
         boxc[1]=co[4];
         boxc[2]=co[8];
/*
*/
         printf("box: ");
         for (i=0; i<3; i++) 
          printf("%e ",boxc[i]); printf("\n");
	} else i3=0;
        if (i3!=9 && ibox) {
         printf("something wrong reading box\n");
        }
        fread(&bytes,sizeof(int),1,fileptr);
        for (i=0; i<Nsit; i++) {
         *r++*=boxc[0];
         *r++*=boxc[1];
         *r++*=boxc[2];
        }
        r-=3*Nsit;
        for (i=0; i<3; i++) boxc[i]*=2.;
        itest=1+3*Nsit;
        if(ibox) itest+=9;
        ifin=i1+i3+i5;
        break;

       case ORAC4_DUMP:
/* read ORAC4 dump binary file, fortran single precision then convert into d */
/* It differs from ORAC3 only in the order of r/box */
        nitems = Nsit;
        i1=1;
        fread(&bytes,sizeof(int),1,fileptr);
        if (Nsit != 0) {
         i1=fread(&fstep,sizeof(float),1,fileptr);
         if(serptr->idunion.ser.conv) reverse4ByteOrder(1,(char *)&fstep);
         printf("fstep= %f\n",fstep);
	} else i1=0;
        if (i1!=1) {
         printf("something wrong reading fstep\n");
        }
        if (ibox) {
         i3=fread(co_s,sizeof(float),9,fileptr);
         if(serptr->idunion.ser.conv) reverse4ByteOrder(9,(char *)co_s);
         boxc[0]=(_REAL)co_s[0];
         boxc[1]=(_REAL)co_s[4];
         boxc[2]=(_REAL)co_s[8];
         cell_to_cart[1][1]=(_REAL)co_s[0];
         cell_to_cart[2][1]=(_REAL)co_s[1];
         cell_to_cart[3][1]=(_REAL)co_s[2];
         cell_to_cart[1][2]=(_REAL)co_s[3];
         cell_to_cart[2][2]=(_REAL)co_s[4];
         cell_to_cart[3][2]=(_REAL)co_s[5];
         cell_to_cart[1][3]=(_REAL)co_s[6];
         cell_to_cart[2][3]=(_REAL)co_s[7];
         cell_to_cart[3][3]=(_REAL)co_s[8];
/*
*/
         printf("box: ");
         for (i=0; i<3; i++) 
          printf("%e ",boxc[i]); printf("\n");
	} 
         else 
        {
         cell_to_cart[1][1]=boxc[0];
         cell_to_cart[2][1]=0.;
         cell_to_cart[3][1]=0.;
         cell_to_cart[1][2]=0.;
         cell_to_cart[2][2]=boxc[1];
         cell_to_cart[3][2]=0.;
         cell_to_cart[1][3]=0.;
         cell_to_cart[2][3]=0.;
         cell_to_cart[3][3]=boxc[2];
         i3=0;
        }
        if (i3!=9 && ibox) {
         printf("something wrong reading box\n");
        }
        bytes=sizeof(float)*3*serptr->idunion.ser.skipb;
        fseek(fileptr,bytes,1);
        i5+=fread(r_s,sizeof(float),3*Nsit,fileptr);
        bytes=sizeof(float)*3*serptr->idunion.ser.skipa;
        fseek(fileptr,bytes,1);
        if(serptr->idunion.ser.conv) reverse4ByteOrder(3*Nsit,(char *)r_s);
        for (i=0; i<3*Nsit; i++) r[i]=(_REAL)r_s[i];
/*
*/
        printf("r0= %f %f %f\n",r[0],r[1],r[2]);
        printf("rf= %f %f %f\n",r[3*Nsit-3],r[3*Nsit-2],r[3*Nsit-1]);
        if (i5!=3*Nsit) {
         printf("something wrong reading r_s: %d\n",i5);
        }
        fread(&bytes,sizeof(int),1,fileptr);
/*
        for (i=0; i<Nsit; i++) {
         *r++*=boxc[0];
         *r++*=boxc[1];
         *r++*=boxc[2];
        }
        r-=3*Nsit;
*/
        if(i3==9) for (i=0; i<3; i++) boxc[i]+=boxc[i];
        itest=1+3*Nsit;
        if(ibox) itest+=9;
        ifin=i1+i3+i5;
        break;

       case DLPOLY:
/* read DLPOLY formatted file */
        /* fscanf(fileptr,"%*s %*d "); */
        nitems = Nsit;
        i1=1;
/*
        if (Nsit != 0) {
         i1=fscanf(fileptr,"%d ",&nitems);
*/
/*
         if (nitems!=Nsit) {
          printf("Nsit in stream is different from the one in input\n");
          printf("Nsit: %d nitems %d\n",Nsit,nitems);
         } 
*/
/*
	} else i1=0;
        if (i1!=1) {
         printf("something wrong reading nitems\n");
        }
*/

        /*fscanf(fileptr,"%*d %*d %*lf");*/
        if (nitems != 0 && ibox) {
         i3+=fscanf(fileptr,"%lf %lf %lf ",&cell_to_cart[1][1],&cell_to_cart[2][1],&cell_to_cart[3][1]);
         i3+=fscanf(fileptr,"%lf %lf %lf ",&cell_to_cart[1][2],&cell_to_cart[2][2],&cell_to_cart[3][2]);
         i3+=fscanf(fileptr,"%lf %lf %lf ",&cell_to_cart[1][3],&cell_to_cart[2][3],&cell_to_cart[3][3]);
        boxc[0]=cell_to_cart[1][1];
        boxc[1]=cell_to_cart[2][2];
        boxc[2]=cell_to_cart[3][3];
         printf("box: ");
        for (i=0; i<3; i++) printf("%f ",boxc[i]); printf("\n");
	} else i3=0;
        if (i3!=9 && ibox) {
         printf("something wrong reading box\n");
        }

        printf("Nsit=%d\n",Nsit);
        for (i=0; i<serptr->idunion.ser.skipb; i++)
         fscanf(fileptr, "%*lf %*lf %*lf ");
        for (i=0,i5=0; i<3*Nsit; i++)
         i5+=fscanf(fileptr, "%lf ", r++);
        r-=3*Nsit;
        for (i=0; i<serptr->idunion.ser.skipa; i++)
         fscanf(fileptr, "%*lf %*lf %*lf ");
        ifin=i1+i3+i5;
        itest=1+9+3*Nsit;
        break;

       case NWCHEM:
        { _REAL timet,tempt,ht;
          char title[TOKLEN];
/* read NWCHEM formatted file */
        nitems = Nsit;
        i1=1;
        fscanf(fileptr,"%*s ");
	preadln(fileptr, "", title);
	sscanf(title,"%lf12.6",&timet);
	printf("timet= %e\n",timet);

        if (nitems != 0 && ibox) {
         i3+=fscanf(fileptr,"%lf %lf %lf ",
	  &cell_to_cart[1][1],&cell_to_cart[2][1],&cell_to_cart[3][1]);
         i3+=fscanf(fileptr,"%lf %lf %lf ",
	  &cell_to_cart[1][2],&cell_to_cart[2][2],&cell_to_cart[3][2]);
         i3+=fscanf(fileptr,"%lf %lf %lf ",
	  &cell_to_cart[1][3],&cell_to_cart[2][3],&cell_to_cart[3][3]);
        boxc[0]=10.*cell_to_cart[1][1];
        boxc[1]=10.*cell_to_cart[2][2];
        boxc[2]=10.*cell_to_cart[3][3];
         printf("box: ");
        for (i=0; i<3; i++) printf("%f ",boxc[i]); printf("\n");
	} else i3=0;
        if (i3!=9 && ibox) {
         printf("something wrong reading box\n");
        }

        printf("Nsit=%d\n",Nsit);
        fscanf(fileptr, "%*s %*d %*d %*d ");
        for (i=0; i<serptr->idunion.ser.skipb; i++)
         fscanf(fileptr, "%*lf %*lf %*lf ");
        for (i=0,i5=0; i<3*Nsit; i++)
         i5+=fscanf(fileptr, "%lf ", r++);
        r-=3*Nsit;
        for (i=0; i<serptr->idunion.ser.skipa; i++)
         fscanf(fileptr, "%*lf %*lf %*lf ");
        for (i=0; i<3*Nsit; i++) *r++ *= 10.;
        r-=3*Nsit;
        ifin=i1+i3+i5;
        itest=1+9+3*Nsit;
        break;
	}

       case NWCHEM_RST:
        { _REAL tempt,tempw,prest;
          char title[TOKLEN];
	  int nsolvt,at_per_solvt,at_per_solut,k;
/* read NWCHEM_RST formatted file */
     preadln(fileptr, "", title);
     printf("stream: NWCHEM_RST: title= %s\n",title);
        i1=fscanf(fileptr,"%d %*d ",&k);
     printf("stream: NWCHEM: k= %d\n",k);
        i3=0;
        if (ibox) {
         i3+=fscanf(fileptr,"%lf %lf %lf ",
	  &cell_to_cart[1][1],&cell_to_cart[2][1],&cell_to_cart[3][1]);
         i3+=fscanf(fileptr,"%lf %lf %lf ",
	  &cell_to_cart[1][2],&cell_to_cart[2][2],&cell_to_cart[3][2]);
         i3+=fscanf(fileptr,"%lf %lf %lf ",
	  &cell_to_cart[1][3],&cell_to_cart[2][3],&cell_to_cart[3][3]);
        boxc[0]=10.*cell_to_cart[1][1];
        boxc[1]=10.*cell_to_cart[2][2];
        boxc[2]=10.*cell_to_cart[3][3];
         printf("box: ");
        for (i=0; i<3; i++) printf("%f ",boxc[i]); printf("\n");
	} else i3=0;
        if (i3!=9 && ibox) {
         printf("something wrong reading box\n");
        }
        fscanf(fileptr, "%le ",&prest);
        fscanf(fileptr, "%lf %lf %*lf ",&tempt,&tempw);
        fscanf(fileptr, "%d %d %*d %d %*d %*d %*d %*d %*d ",
	 &nsolvt,&at_per_solvt,&at_per_solut);

  printf("stream: NWCHEM_RST: at_per_solvt= %d nsolvt= %d\n",
    at_per_solvt,nsolvt);
  printf("stream: NWCHEM_RST: tempt= %e prest= %e\n",tempt,prest);
  if( (at_per_solut+nsolvt*at_per_solvt)!=prm->Nsit)
  {
   printf("stream: NWCHEM_RST: lettura sbagliata.\n");
   at_per_solut=0;
   at_per_solvt=0;
   nsolvt=0;
  }
	i5=0;
	r+=3*at_per_solut;
        for (i=0; i<nsolvt; i++)
        {
         for (k=0; k<at_per_solvt; k++)
         {
          fscanf(fileptr, "%lf %lf %lf ", &r[0],&r[1],&r[2]);
	  r+=3;
          fscanf(fileptr, "%*lf %*lf %*lf ");
	  i5+=3;
         }
         fscanf(fileptr, "%*d %*lf %*lf %*lf ");
        }
        r-=3*(at_per_solut+at_per_solvt*nsolvt);
        for (i=0; i<at_per_solut; i++)
        {
         fscanf(fileptr, "%*d %lf %lf %lf ", &r[0],&r[1],&r[2]);
          fscanf(fileptr, "%*lf %*lf %*lf ");
	 r+=3;
	 i5+=3;
        }
        r-=3*at_per_solut;
        for (i=0; i<3*Nsit; i++) *r++ *= 10.;
        r-=3*(at_per_solut+at_per_solvt*nsolvt);
        ifin=i1+i3+i5;
	fflush(stdout);
        itest=1+9+3*prm->Nsit;
        break;
	}

       case CHARMM:
/* read CHARMM unformatted file */
        i1=1;
        if (ibox) {
         fread(&bytes,sizeof(int),1,fileptr);
/* single precision, never works for CHARMM */
/*
         i3=fread(co_s,sizeof(float),6,fileptr);
         if(serptr->idunion.ser.conv) reverse4ByteOrder(6,(char *)co_s);
         fread(&bytes,sizeof(int),1,fileptr);
         boxc[0]=(_REAL)co_s[0];
         boxc[1]=(_REAL)co_s[2];
         boxc[2]=(_REAL)co_s[5];
         co[0]=(_REAL)co_s[0];
         co[3]=(_REAL)co_s[2];
         co[5]=(_REAL)co_s[5];
*/
/* double precision */
         i3=fread(co,sizeof(double),6,fileptr);
         fread(&bytes,sizeof(int),1,fileptr);
         printf("CHARMM: cart_to_cell: ");
         for (i=0; i<6; i++) 
           printf("%e ",co[i]);
         printf("\n");
         boxc[0]=co[0];
         boxc[1]=co[2];
         boxc[2]=co[5];
/* only orthorhombic */
         cell_to_cart[1][1]=co[0];
         cell_to_cart[2][1]=0.;
         cell_to_cart[3][1]=0.;
         cell_to_cart[1][2]=0.;
         cell_to_cart[2][2]=co[2];
         cell_to_cart[3][2]=0.;
         cell_to_cart[1][3]=0.;
         cell_to_cart[2][3]=0.;
         cell_to_cart[3][3]=co[5];
/*
*/
         printf("CHARMM: cell_to_cart:\n");
         for (i=1; i<=3; i++) 
          for (j=1; j<=3; j++) 
           printf("%e ",cell_to_cart[i][j]);
         printf("\n");
	} else i3=0;
        if (i3!=6 && ibox) {
         printf("something wrong reading box\n");
        }

        fread(&bytes,sizeof(int),1,fileptr);
        bytes=sizeof(float)*serptr->idunion.ser.skipb;
        fseek(fileptr,bytes,1);
        for(i=0;i<Nsit;i++) {
         i5+=fread(r_s,sizeof(float),1,fileptr);
         r_s+=3;
        }
        bytes=sizeof(float)*serptr->idunion.ser.skipa;
        fseek(fileptr,bytes,1);
        fread(&bytes,sizeof(int),1,fileptr);
        r_s-=3*Nsit;
        r_s++;
        fread(&bytes,sizeof(int),1,fileptr);
        bytes=sizeof(float)*serptr->idunion.ser.skipb;
        fseek(fileptr,bytes,1);
        for(i=0;i<Nsit;i++) {
         i5+=fread(r_s,sizeof(float),1,fileptr);
         r_s+=3;
        }
        bytes=sizeof(float)*serptr->idunion.ser.skipa;
        fseek(fileptr,bytes,1);
        fread(&bytes,sizeof(int),1,fileptr);
        r_s-=3*Nsit;
        r_s++;
        fread(&bytes,sizeof(int),1,fileptr);
        bytes=sizeof(float)*serptr->idunion.ser.skipb;
        fseek(fileptr,bytes,1);
        for(i=0;i<Nsit;i++) {
         i5+=fread(r_s,sizeof(float),1,fileptr);
         r_s+=3;
        }
        bytes=sizeof(float)*serptr->idunion.ser.skipa;
        fseek(fileptr,bytes,1);
        fread(&bytes,sizeof(int),1,fileptr);
        r_s-=3*Nsit;
        r_s-=2;

        if(serptr->idunion.ser.conv) reverse4ByteOrder(3*Nsit,(char *)r_s);
        for (i=0; i<3*Nsit; i++) r[i]=(_REAL)r_s[i];
        printf("r0= %f %f %f\n",r[0],r[1],r[2]);
        if (i5!=3*Nsit) {
         printf("something wrong reading r_s: %d\n",i5);
        }
/*
        for (i=0; i<Nsit; i++) {
         *r++*=boxc[0];
         *r++*=boxc[1];
         *r++*=boxc[2];
        }
        r-=3*Nsit;
*/
        /* if(i3==6) for (i=0; i<3; i++) boxc[i]+=boxc[i]; */
        itest=3*Nsit;
        if(ibox) itest+=6;
        ifin=i3+i5;
        break;

       case PDB:
        nitems=Nsit;
        for (i=0; i<nitems; i++) {
         if (scanfmt("%lf %lf %lf ", fileptr, r)) break;
         r += 3;
        }
        r-=3*i;
        printf("r0= %f %f %f\n",r[0],r[1],r[2]);
        ifin=3*i;
        itest=3*Nsit;
        break;

       case POS:
/* read POS formatted file generated by cp.x */
        nitems=Nsit;
        fscanf(fileptr, "%*d %*f ");
        for (i=0; i<serptr->idunion.ser.skipb; i++)
         fscanf(fileptr, "%*e %*e %*e ");
        for (i=0,i5=0; i<3*nitems; i++)
                i5+=fscanf(fileptr, "%e ", r_s++);
        for (i=0; i<serptr->idunion.ser.skipa; i++)
         fscanf(fileptr, "%*e %*e %*e ");
/*
        while(getc(fileptr) != 10);
*/
        r_s-=3*nitems;
        if(ibox)
        {
         printf("stream: Box not present in POS files.\n");
         exit(1);
        }

        for (i=0; i<3*Nsit; i++) r[i]=(_REAL)r_s[i];
        printf("r0= %f %f %f\n",r[0],r[1],r[2]);
        ifin=i5;
        itest=3*Nsit;
        break;

       case AXSF:
/* read axsf formatted file (animated xcrysden, produced by QE/cppp) */
/* PRIMVEC a,b,c in cartesian coords. (Angstrom) */
        fscanf(fileptr, "%*s %*d ");
        for (i=0,i3=0; i<9; i++)
         i3+=fscanf(fileptr, "%lf ", &co[i]);
        if(i3!=9)
         printf("AXSF: something wrong reading cell components...\n");
/* calculation of cell_to_cart (see wikipedia: Fractional coordinates) */
        { _REAL cos_alpha, cos_beta, cos_gamma, sin_gamma, a,b,c,v;
         a=l_sqrt(co[0]*co[0]+co[1]*co[1]+co[2]*co[2]);
         b=l_sqrt(co[3]*co[3]+co[4]*co[4]+co[5]*co[5]);
         c=l_sqrt(co[6]*co[6]+co[7]*co[7]+co[8]*co[8]);
         cos_alpha = (co[3]*co[6]+co[4]*co[7]+co[5]*co[8])/(b*c);
         cos_beta  = (co[0]*co[6]+co[1]*co[7]+co[2]*co[8])/(a*c);
         cos_gamma = (co[0]*co[3]+co[1]*co[4]+co[2]*co[5])/(a*b);
         sin_gamma = l_sqrt(1.-cos_gamma*cos_gamma);
         v = l_sqrt(1.-cos_gamma*cos_gamma 
                      -cos_alpha*cos_alpha
                      -cos_beta*cos_beta
                      +2.*cos_alpha*cos_beta*cos_gamma );
         cell_to_cart[1][1] = a;
         cell_to_cart[1][2] = 0.;
         cell_to_cart[1][3] = 0.;
         cell_to_cart[2][1] = b*cos_gamma;
         cell_to_cart[2][2] = b*sin_gamma;
         cell_to_cart[2][3] = 0.;
         cell_to_cart[3][1] = c*cos_beta;
         cell_to_cart[3][2] = c*(cos_alpha - cos_beta*cos_gamma)/sin_gamma;
         cell_to_cart[3][3] = c*v/sin_gamma;
        }
/*
        printf ("AXSF: cell_to_cart:\n");
        for(i=1; i<=3; i++)
        {
         for(j=1; j<=3; j++)
         {
          printf ("%d %d %f\n",i,j,cell_to_cart[i][j]);
         }
        }
*/
/* CONVVEC , skip */
        fscanf(fileptr, "%*s %*d ");
        for (i=0; i<9; i++)
         fscanf(fileptr, "%*f ");
/* PRIMCOORD x,y,z cartesian (Angstrom), atoms are in the primitive cell */
        fscanf(fileptr, "%*s %*d ");
        fscanf(fileptr, "%d %*d ",&nitems);
        printf ("AXSF: nitems= %d\n",nitems);
        if(nitems != Nsit) 
        {
         printf ("AXSF: nitems= %d, different from Nsit= %d\n",nitems,Nsit);
         exit(1);
        }
        for (i=0,i5=0; i<nitems; i++)
        {
         i5+=fscanf(fileptr, "%*d %f %f %f ", &r_s[0],&r_s[1],&r_s[2]);
         r_s+=3;
        }
        r_s-=3*nitems;
        for (i=0; i<3*Nsit; i++) r[i]=(_REAL)r_s[i];
        ifin=i3+i5;
        itest=9+3*Nsit;
        break;

       case XYZ:
/* read XYZ formatted file */
        fscanf(fileptr, "%d ",&nitems);
     printf("stream: nitems= %d\n",nitems);
        if(nitems != Nsit) 
        {
         printf ("XYZ: nitems= %d, different from Nsit= %d\n",nitems,Nsit);
         exit(1);
        }
	preadln(fileptr, "", Cdum);
        for (i=0; i<serptr->idunion.ser.skipb; i++)
//         fscanf(fileptr, "%*s %*f %*f %*f %*f %*f %*f %*f ");
         fscanf(fileptr, "%*s %*f %*f %*f %*d %*d %*d ");
        for (i=0,i5=0; i<nitems; i++)
	{
                fscanf(fileptr, "%*s ");
                i5+=fscanf(fileptr, "%f ", r_s++);
                i5+=fscanf(fileptr, "%f ", r_s++);
                i5+=fscanf(fileptr, "%f ", r_s++);
//                fscanf(fileptr, "%*f %*f %*f ");
                fscanf(fileptr, "%*d %*d %*d ");
	}
        for (i=0; i<serptr->idunion.ser.skipa; i++)
//         fscanf(fileptr, "%*s %*f %*f %*f %*f %*f %*f %*f ");
         fscanf(fileptr, "%*s %*f %*f %*f %*d %*d %*d ");
        r_s-=3*nitems;
        for (i=0; i<3*Nsit; i++) r[i]=(_REAL)r_s[i];
        printf("r(0)= %f %f %f\n",r[0],r[1],r[2]);
        printf("r(1)= %f %f %f\n",r[3],r[4],r[5]);
        printf("r(N-1)= %f %f %f\n",r[3*Nsit-3],r[3*Nsit-2],r[3*Nsit-1]);
        ifin=i5;
        itest=3*Nsit;
        break;

       case AMBER:
/* read AMBER formatted trajectory file */
        nitems=Nsit;
/*
        for (i=0,i3=0; i<3; i++)
                i3+=fscanf(fileptr, "%f ", &boxc_s[i]);
        for (i=0; i<3; i++) boxc[i]=(_REAL)boxc_s[i];
        i3=3;
*/
        for (i=0; i<serptr->idunion.ser.skipb; i++)
         fscanf(fileptr, "%*f %*f %*f ");
        for (i=0,i5=0; i<3*nitems; i++)
                i5+=fscanf(fileptr, "%f ", r_s++);
        for (i=0; i<serptr->idunion.ser.skipa; i++)
         fscanf(fileptr, "%*f %*f %*f ");
/*
        while(getc(fileptr) != 10);
*/
        r_s-=3*nitems;
        if(ibox)
        {
        for (i=0,i3=0; i<3; i++)
                i3+=fscanf(fileptr, "%f ", &boxc_s[i]);
        for (i=0; i<3; i++) boxc[i]=(_REAL)boxc_s[i];
        }

/*
        {
         FILE *ftemp;
         int kk;
         ftemp=fopen("tmp.out","a");
         for (i=0,kk=0;i<3*nitems;i++) {
          fprintf(ftemp,"%8.3f",r_s[i]);
          kk++;
          if(kk%10==0) {fprintf(ftemp,"\n"); kk=0;}
         }
         fprintf(ftemp,"\n");
        }
*/

        for (i=0; i<3*Nsit; i++) r[i]=(_REAL)r_s[i];
        printf("r0= %f %f %f\n",r[0],r[1],r[2]);
        ifin=i3+i5;
        itest=3*Nsit;
        if(ibox) itest+=3;
        break;

       case DISCOVER2:
       case DISCOVER3:
  printf("reading stream...\n");
/* read discover .his f77 unformatted file */
        if(serptr->idunion.ser.fileset==1) {

         rdhead(serptr,&i3,&itest);

        } else {

/*      printf("leggo le coordinate...\n"); */
        nitems=Nsit;
        i3+=fread(&bytes,sizeof(int),1,fileptr);
        i3+=fread(Idum,sizeof(int),1,fileptr);
        i3+=fread(&bytes,sizeof(int),1,fileptr);
        i3+=fread(&bytes,sizeof(int),1,fileptr);
        i3+=fread(&Etot,sizeof(double),1,fileptr);
        i3+=fread(&Epot,sizeof(double),1,fileptr);
        i3+=fread(&Ekin,sizeof(double),1,fileptr);
        i3+=fread(Enerd,prm->Nener*sizeof(double),1,fileptr);
        for(i=0;i<prm->Nmolt;i++) i3+=fread(Ddum,sizeof(double),1,fileptr);
        for(i=0;i<prm->Nmolt*prm->Nener;i++) i3+=fread(Ddum,sizeof(double),1,fileptr);
        for(i=0;i<prm->Nmolt;i++) i3+=fread(Ddum,sizeof(double),1,fileptr);
        for(i=0;i<prm->Nmolt;i++) i3+=fread(Ddum,sizeof(double),1,fileptr);
        for(i=0;i<prm->Nmolt;i++) i3+=fread(Ddum,sizeof(double),1,fileptr);
        for(i=0;i<prm->Nmolt;i++) i3+=fread(Ddum,sizeof(double),1,fileptr);
        i3+=fread(&Press,sizeof(double),1,fileptr);
        i3+=fread(Ddum,sizeof(double),1,fileptr);
/* printf("En: %e %e %e %e\n",Etot,Epot,Ekin,Press); */
        for(i=0;i<9;i++) i3+=fread(Ddum,sizeof(double),1,fileptr);
        for(i=0;i<9;i++) i3+=fread(Ddum,sizeof(double),1,fileptr);
        for(i=0;i<9;i++) i3+=fread(Ddum,sizeof(double),1,fileptr);
        for(i=0;i<9;i++) i3+=fread(Ddum,sizeof(double),1,fileptr);
        for(i=0;i<9;i++) i3+=fread(Ddum,sizeof(double),1,fileptr);
        for(i=0;i<9;i++) i3+=fread(Ddum,sizeof(double),1,fileptr);
        i3+=fread(&bytes,sizeof(int),1,fileptr);

        i3+=fread(&bytes,sizeof(int),1,fileptr);
        if(ibox) {
         i3+=fread(boxc,3*sizeof(double),1,fileptr);
        } else {
         i3+=fread(Ddum,3*sizeof(double),1,fileptr);
        }
        for(i=0;i<3;i++) i3+=fread(Ddum,sizeof(double),1,fileptr);
        for(i=0;i<9;i++) i3+=fread(Ddum,sizeof(double),1,fileptr);
        i3+=fread(&bytes,sizeof(int),1,fileptr);

        i3+=fread(&bytes,sizeof(int),1,fileptr);
        i3+=fread(r_s,3*prm->Nmov*sizeof(float),1,fileptr);
        if(serptr->idunion.ser.conv) reverse4ByteOrder(3*prm->Nmov,(char *)r_s);
/*
         printf("stream:r_s: %e %e %e\n",r_s[0],r_s[1],r_s[2]); 
         printf("stream:r_s: %e %e %e\n",r_s[3],r_s[4],r_s[5]); 
         printf("stream:r_s: %e %e %e\n",r_s[6],r_s[7],r_s[8]); 
*/

        i3+=fread(&bytes,sizeof(int),1,fileptr);
        i3+=fread(&bytes,sizeof(int),1,fileptr);
        for(i=0;i<3*prm->Nmov;i++) i3+=fread(Fdum,sizeof(float),1,fileptr);
        i3+=fread(&bytes,sizeof(int),1,fileptr);
        for(i=0;i<3*prm->Nmov;i++) *r++ =(_REAL)(*r_s++);
        r-=3*prm->Nmov; r_s-=3*prm->Nmov;
        itest=85+5*prm->Nmolt+prm->Nmolt*prm->Nener+3*prm->Nmov;
        serptr->idunion.ser.set->Etot=Etot;
        serptr->idunion.ser.set->Epot=Epot;
        serptr->idunion.ser.set->Ekin=Ekin;
        serptr->idunion.ser.set->Press=Press;

        } /* endif fileset */

        ifin=i3;
        break;

       default:
        printf("STREAM: fmt not found\n");
        exit(1);
        break;

       case CBMC:
/* read FORTRAN binary file generated by Frenkel_Smit case study 19 (CBMC) */
        nitems = Nsit;
        i1=1;
        fread(&bytes,sizeof(int),1,fileptr);
        if (Nsit != 0) {
         i1=fread(&fstepd,sizeof(_REAL),1,fileptr);
         printf("fstep= %f\n",fstepd);
	} else i1=0;
        if (i1!=1) {
         printf("something wrong reading fstep\n");
        }
        fread(&bytes,sizeof(int),1,fileptr);
        if (ibox) {
         fread(&bytes,sizeof(int),1,fileptr);
         i3=fread(&boxc[0],sizeof(_REAL),1,fileptr);
         fread(&bytes,sizeof(int),1,fileptr);
         boxc[1]=boxc[0];
         boxc[2]=boxc[0];
/*
*/
         printf("box: ");
         for (i=0; i<3; i++) 
          printf("%e ",boxc[i]); printf("\n");
	} else i3=0;
        if (i3!=1 && ibox) {
         printf("something wrong reading box\n");
        }
        fread(&bytes,sizeof(int),1,fileptr);
        bytes=sizeof(_REAL)*3*serptr->idunion.ser.skipb;
        fseek(fileptr,bytes,1);
        i5+=fread(r,sizeof(_REAL),3*Nsit,fileptr);
        bytes=sizeof(_REAL)*3*serptr->idunion.ser.skipa;
        fseek(fileptr,bytes,1);
        printf("r0= %f %f %f\n",r[0],r[1],r[2]);
        if (i5!=3*Nsit) {
         printf("something wrong reading r_s: %d\n",i5);
        }
        fread(&bytes,sizeof(int),1,fileptr);
        itest=1+3*Nsit;
        if(ibox) itest+=1;
        ifin=i1+i3+i5;
        break;

       case MUCA:
/* read ASCII file generated by Yuji with Epot on top */
        if (Nsit != 0) {
         i1=fscanf(fileptr, "%d %le ",&nitems,&fstepd);
         printf("Natoms= %d, Epot= %f\n",nitems,fstepd);
	} else i1=0;
        if (i1!=2) {
         printf("something wrong reading Natoms and Epot.\n");
        }
        if (ibox) {
         i2=0;
         i2=fscanf(fileptr,"%lf ",&boxc[0]);
         if(i2!=1)
          i2=fscanf(fileptr,"%le ",&boxc[0]);
         if(i2==1) i3++;
         i2=0;
         i2=fscanf(fileptr,"%lf ",&boxc[1]);
         if(i2!=1)
          i2=fscanf(fileptr,"%le ",&boxc[1]);
         if(i2==1) i3++;
         i2=0;
         i2=fscanf(fileptr,"%lf ",&boxc[2]);
         if(i2!=1)
          i2=fscanf(fileptr,"%le ",&boxc[2]);
         if(i2==1) i3++;
/*
*/
         printf("box: ");
         for (i=0; i<3; i++) 
          printf("%e ",boxc[i]); printf("\n");
	} else i3=0;
        if (i3!=3 && ibox) {
         printf("something wrong reading box\n");
        }
        for (i=0; i<serptr->idunion.ser.skipb; i++)
         fscanf(fileptr, "%*s %*s %*s ");
        for (i=0,i5=0; i<3*nitems; i++)
        {
         i2=0;
         i2=fscanf(fileptr,"%lf ",r);
         if(i2!=1)
          i2=fscanf(fileptr,"%le ",r);
         if(i2==1) 
         {
          i5++; r++;
         }
        }
        r-=3*nitems;
        for (i=0; i<serptr->idunion.ser.skipa; i++)
         fscanf(fileptr, "%*s %*s %*s ");
        printf("r0= %f %f %f\n",r[0],r[1],r[2]);
        printf("rf= %f %f %f\n",r[3*nitems-3],r[3*nitems-2],r[3*nitems-1]);
        if (i5!=3*Nsit) {
         printf("something wrong reading r: %d\n",i5);
        }
        itest=2+3*Nsit;
        if(ibox) itest+=3;
        ifin=i1+i3+i5;
        break;

      }

                if (ifin==itest && 
                    it<serptr->idunion.ser.nstates) {

			int	wct;

/*
#include "decane.c"
*/
/*
  if(serptr->idunion.ser.atom_selected>0)
  {
   _REAL rt1[3],rt2[3],rt3[3],rtmp[3];
   void genquatch1();
   int atom_selected=serptr->idunion.ser.atom_selected;

   rt1[0]=r[3*atom_selected-3];
   rt1[1]=r[3*atom_selected-2];
   rt1[2]=r[3*atom_selected-1];
   rt2[0]=r[3*atom_selected+0];
   rt2[1]=r[3*atom_selected+1];
   rt2[2]=r[3*atom_selected+2];
   rt3[0]=r[3*atom_selected+3];
   rt3[1]=r[3*atom_selected+4];
   rt3[2]=r[3*atom_selected+5];
   genquatch1(rt1,rt2,rt3,rtmp);
   r[3*1+0]=rtmp[0];
   r[3*1+1]=rtmp[1];
   r[3*1+2]=rtmp[2];
  }
*/
/* GE (or MUCA) reweighting */
/*
 {
  _REAL E,a1,b1,a2,b2,S,w,
        beta=1.677, E1=-13057.5, E2=-8299.5,
        cc[]={0.,-6.372308982174522,-1.3695557599801341E-03,
              -1.0941858160981546E-07,1.1887055237102052E-13,2.0415605344344340E-16,
              -7.2741363552754789E-20,-1.0782577331525304E-23,-5.4814527404273566E-28,
              -9.9502453486683496E-33};
  int kk,nmax=9;

  E=fstepd;
  a1=0.; b1=0.;
  for(kk=0;kk<=nmax;kk++)
  {
   a1+=cc[kk]*pow(E1,(_REAL)kk);
   if(kk>0)
    b1+=(_REAL)kk*cc[kk]*pow(E1,(_REAL)(kk-1));
  }
  a2=0.; b2=0.;
  for(kk=0;kk<=nmax;kk++)
  {
   a2+=cc[kk]*pow(E2,(_REAL)kk);
   if(kk>0)
    b2+=(_REAL)kk*cc[kk]*pow(E2,(_REAL)(kk-1));
  }

  if(E<=E1)
  {
   S=a1+b1*(E-E1);
  }
  else if(E>=E2)
  {
   S=a2+b2*(E-E2);
  }
  else
  {
   S=0;
   for(kk=0;kk<=nmax;kk++)
   {
    S+=cc[kk]*pow(E,(_REAL)kk);
   }
  }
  S-=beta*E;
  S-=lnZf;
  S-=wmax;
  if(it==0) lnZ=S;
  else
  {
   a1=lnZ;
   a2=S;
   b1=MINIMUM(a1,a2);
   b2=MAXIMUM(a1,a2);
   w=log(1.+exp(b1-b2));
   w+=b2;
   lnZ=w;
  }
  printf("WW= %e %e %23.15e\n",E,S,lnZ);
  w=exp(S);
  fwrite(&w,sizeof(_REAL),1,ftmp);
 }
*/

		        printf("coord set ok:\n");
                        fflush(stdout);
/*
        printf("Thermo: %e %e %e %e\n",Etot,Epot,Ekin,Press);
*/
/* do some scaling */
               if(serptr->idunion.ser.scale==1) {
                 for(i=0;i<Nsit;i++) {
                  (*r++)*=sig0r;
                  (*r++)*=sig0r;
                  (*r++)*=sig0r;
                 }
                 r -= 3*Nsit;
                 if(i3!=0) for(i=0;i<3;i++) boxc[i] *= sig0r;
               } else if(serptr->idunion.ser.scale==2) {
                 for(i=0;i<Nsit;i++) {
                  (*r++)/=sig0r;
                  (*r++)/=sig0r;
                  (*r++)/=sig0r;
                 }
                 r -= 3*Nsit;
                 if(i3!=0) for(i=0;i<3;i++) boxc[i] /= sig0r;
               }
        if(fmt!=DLPOLY && fmt!=ORAC4_DUMP && fmt!=CHARMM && fmt!=AXSF) {
         for(i=1;i<=3;i++)
          for(j=1;j<=3;j++)
           cell_to_cart[i][j]=0.;
         cell_to_cart[1][1]=1.0*boxc[0];
         cell_to_cart[2][2]=1.0*boxc[1];
         cell_to_cart[3][3]=1.0*boxc[2];
        }
         
        serptr->idunion.ser.set->boxc[0]=cell_to_cart[1][1];
        serptr->idunion.ser.set->boxc[1]=cell_to_cart[2][2];
        serptr->idunion.ser.set->boxc[2]=cell_to_cart[3][3];
        invert(cell_to_cart,3,3,cart_to_cell);
        serptr->idunion.ser.set->boxm1[0]=cart_to_cell[1][1];
        serptr->idunion.ser.set->boxm1[1]=cart_to_cell[2][2];
        serptr->idunion.ser.set->boxm1[2]=cart_to_cell[3][3];

			/*
			 *  success: done unless in WINdow mode
			 */

			if (!serptr->idunion.ser.skip)
				return(1);

			/*
			 *  this stream is in WINdow mode - consider 
			 *	whether to skip this set
			 */

			wct = ++serptr->idunion.ser.winct;
			if (wct < 0) 
				continue;	/* skip another set */
			if (wct ==  serptr->idunion.ser.win) {
				serptr->idunion.ser.winct = 
						serptr->idunion.ser.skip;
				continue;
			}
			return(1);
		} else {
			printf("incomplete coord set: %s\n",
					serptr->idunion.ser.name); 
			printf("\tstarting another file if any\n");
		}
			
		/*
		 *  exhausted file: start next file
		 */

                rewind(fileptr);
		if (++serptr->idunion.ser.ifile 
                      == serptr->idunion.ser.nfiles) {
		    printf("STREAM %s: no more files/crds (iteration %d)\n", 
						serptr->id_name, it);
	            serptr->idunion.ser.fileset = 0;
	            serptr->idunion.ser.ifile = 0;
	            serptr->idunion.ser.winct = 
	             serptr->idunion.ser.skip;
                    startser(serptr, 0);
		    return(0);
		}
                startser(serptr, serptr->idunion.ser.ifile);
		fileptr = serptr->idunion.ser.fileptr;
		r = serptr->idunion.ser.set->crd;
		boxc = serptr->idunion.ser.set->boxc;
		Enerd = serptr->idunion.ser.set->Enerd;
		Etot=serptr->idunion.ser.set->Etot;
		Epot=serptr->idunion.ser.set->Epot;
		Ekin=serptr->idunion.ser.set->Ekin;
		Press=serptr->idunion.ser.set->Press;

		printf("nextser next file (%d) started\n", 
					serptr->idunion.ser.ifile);

		/*
		 *  continue loop: read coordinates
		 */
	}

#ifdef lint
	return(0);
#endif
}


void
rdhead(idstruct *serptr,int *i3,int *ifin)
{
   FILE    *fileptr = serptr->idunion.ser.fileptr;
   _REAL   *r = serptr->idunion.ser.set->crd;
   double Etot=serptr->idunion.ser.set->Etot;
   double Epot=serptr->idunion.ser.set->Epot;
   double Ekin=serptr->idunion.ser.set->Ekin;
   double Press=serptr->idunion.ser.set->Press;
   float   *r_s = serptr->idunion.ser.set->crd_s;
   _REAL   *boxc = serptr->idunion.ser.set->boxc;
   _REAL sig0r = serptr->idunion.ser.prm->sig0r;

   int i,ir,bytes,bytest,nitems,nrecord;
   int Idum[1],sizeAtNam;
   char Cdum[1];
   double Ddum[1];
   float Fdum[1];
   extern parmstruct *prm;

   printf("rdhead:reading header...\n");
   ir=0;
/* read headers */
/* control for first frame */
   nitems=1;
   bytest=sizeof(int);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   ir+=fread(&bytes,sizeof(int),1,fileptr);

/* Discover version */
   nitems=20;
   bytest=4*sizeof(char);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(Cdum,bytest,1,fileptr);
   nitems=1;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   ir+=fread(&bytes,sizeof(int),1,fileptr);

/* title */
   nitems=20;
   bytest=4*sizeof(char);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(Cdum,bytest,1,fileptr);
   ir+=fread(&bytes,sizeof(int),1,fileptr);

/* title */
   nitems=20;
   bytest=4*sizeof(char);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(Cdum,bytest,1,fileptr);
   ir+=fread(&bytes,sizeof(int),1,fileptr);

/* Forcefield types and masses */
   nitems=1;
   bytest=sizeof(int);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
/*  printf("rdhead:Ntypes: %d\n",*Idum);  */
   nitems=prm->Ntypes;
   bytest=4*sizeof(char);
   for(i=0;i<nitems;i++) ir+=fread(Cdum,bytest,1,fileptr);
   nitems=prm->Ntypes;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   ir+=fread(&bytes,sizeof(int),1,fileptr);

/* Residue names */
   nitems=1;
   bytest=sizeof(int);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   /* printf("rdhead:Nres: %d\n",*Idum); */
   nitems=prm->Nres;
   bytest=4*sizeof(char);
   for(i=0;i<nitems;i++) {
    ir+=fread(Cdum,bytest,1,fileptr);
   }
   ir+=fread(&bytes,sizeof(int),1,fileptr);

/* index to forcefield and atom names */
   nitems=1;
   bytest=sizeof(int);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   /* printf("rdhead:Nsit: %d\n",*Idum); */
   nitems=prm->Nsit;
   bytest=sizeof(int);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   nitems=prm->Nsit;
/* sizeAtNam = 5 for Discover 2.x, = 4 for Discover 3.0 or later */
   if(serptr->idunion.ser.fmt == DISCOVER3) {
    sizeAtNam=4; 
   } else if(serptr->idunion.ser.fmt == DISCOVER2) {
    sizeAtNam=5; 
   } else {
    printf("rdhead:DISCOVER not found.\n");
    exit(1);
   }
   bytest=sizeAtNam*sizeof(char);
   for(i=0;i<nitems;i++) {
    ir+=fread(Cdum,bytest,1,fileptr);
   }
   ir+=fread(&bytes,sizeof(int),1,fileptr);

/* moveable atoms */
   nitems=1;
   bytest=sizeof(int);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   nitems=1;
   bytest=sizeof(int);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   nitems=prm->Nmov;
   bytest=sizeof(int);
   for(i=0;i<nitems;i++) {
    ir+=fread(Idum,bytest,1,fileptr);
   }
   ir+=fread(&bytes,sizeof(int),1,fileptr);

/* molecules */
   nitems=1;
   bytest=sizeof(int);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   nitems=prm->Nmolt;
   bytest=sizeof(int);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   nitems=prm->Nmolt;
   bytest=sizeof(int);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   ir+=fread(&bytes,sizeof(int),1,fileptr);

/* residues */
   nitems=1;
   bytest=sizeof(int);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   /* printf("rdhead:nrest: %d\n",*Idum); */
   nitems=2*prm->Nrest;
   bytest=sizeof(int);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   nitems=prm->Nrest;
   bytest=sizeof(int);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   ir+=fread(&bytes,sizeof(int),1,fileptr);

/* bonds */
   nitems=1;
   bytest=sizeof(int);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   /* printf("rdhead:nbond: %d\n",*Idum); */
   ir+=fread(&bytes,sizeof(int),1,fileptr);

   if(prm->Nbond!=0) {
   nitems=2*prm->Nbond;
   bytest=sizeof(int);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   }

/* cell */
   nitems=6;
   bytest=sizeof(double);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   if (serptr->idunion.ser.ibox) {
    ir+=fread(boxc,3*bytest,1,fileptr);
   } else {
    ir+=fread(Ddum,3*bytest,1,fileptr);
   }
   for(i=0;i<3;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=3*3*196;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=3*196;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=3*3*196;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=1;
   bytest=sizeof(int);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   /* printf("rdhead:nsym: %d\n",*Idum); */
   nitems=4;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
/*
   nitems=1;
   bytest=sizeof(int);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
*/
   nitems=6;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=6;
   bytest=sizeof(int);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   ir+=fread(&bytes,sizeof(int),1,fileptr);

/* steps */
   nitems=1;
   bytest=sizeof(int);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   nitems=1;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=1;
   bytest=sizeof(int);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   nitems=1;
   bytest=sizeof(int);
   for(i=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   ir+=fread(&bytes,sizeof(int),1,fileptr);

/* energies */
   nitems=1;
   bytest=sizeof(double);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<nitems;i++) ir+=fread(&Etot,bytest,1,fileptr);
   nitems=1;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(&Epot,bytest,1,fileptr);
   nitems=1;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(&Ekin,bytest,1,fileptr);
   nitems=prm->Nener;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=prm->Nmolt;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=prm->Nener*prm->Nmolt;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=prm->Nmolt;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=prm->Nmolt;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=prm->Nmolt;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=prm->Nmolt;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=1;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(&Press,bytest,1,fileptr);
/*  printf("Thermo: %e %e %e %e\n",Etot,Epot,Ekin,Press); */
   nitems=1;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   ir+=fread(&bytes,sizeof(int),1,fileptr);

/* coordinates and velocities */
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   ir+=fread(r_s,3*prm->Nmov*sizeof(float),1,fileptr);
   if(serptr->idunion.ser.conv) reverse4ByteOrder(3*prm->Nmov,(char *)r_s);
/*
         printf("stream:r_s: %e %e %e\n",r_s[0],r_s[1],r_s[2]);
         printf("stream:r_s: %e %e %e\n",r_s[3],r_s[4],r_s[5]);
         printf("stream:r_s: %e %e %e\n",r_s[6],r_s[7],r_s[8]);
         fflush(stdout);
*/

   ir+=fread(&bytes,sizeof(int),1,fileptr);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<3*prm->Nmov;i++) ir+=fread(Fdum,sizeof(float),1,fileptr);
   ir+=fread(&bytes,sizeof(int),1,fileptr);
   for(i=0;i<3*prm->Nmov;i++) *r++ =(_REAL)(*r_s++);
   r-=3*prm->Nmov; r_s-=3*prm->Nmov;
   *i3=ir;
   nrecord=16;
   if(prm->Nbond!=0) nrecord++;
   *ifin= 3*20+2*prm->Ntypes+prm->Nres+5*prm->Nsit+1
         +prm->Nmov+7*prm->Nmolt+3*prm->Nrest+2*prm->Nbond
         +6+3*9+2*3*3*196+3*196+4+2*6+prm->Nener+prm->Nmolt*prm->Nener
         +6*9 +20 +2*nrecord-1;
/* printf("i3,ifin: %d %d\n",*i3,*ifin); */
   serptr->idunion.ser.set->Etot=Etot;
   serptr->idunion.ser.set->Epot=Epot;
   serptr->idunion.ser.set->Ekin=Ekin;
   serptr->idunion.ser.set->Press=Press;

/* end reading */
}
