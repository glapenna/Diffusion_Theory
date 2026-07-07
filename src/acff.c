#include "head.h"

#define MAXJ 1000

// extern char		tok[];
// char		tok[];
// void *get2(int,int);

/**** extern ****/
extern idstruct		*startid(),*getstruct();
extern _REAL **matrix(),*vector();
extern int it;

/************************************************************************/

void
init_acff()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(ACFF, id);

	dptr->idunion.acff.type = 0;
	dptr->idunion.acff.domega = 0.;
	dptr->idunion.acff.typo = 0;
	dptr->idunion.acff.nomega = 0;
	dptr->idunion.acff.cases = 0;
	dptr->idunion.acff.casesw = 0.;
	dptr->idunion.acff.dump = 1;
	dptr->idunion.acff.printj = 0;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("ACFF: expected file name","");
	strcpy(dptr->idunion.acff.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.acff.fname, 0, &exists, "ACFF"))
		exit(1);
	if ((dptr->idunion.acff.file = 
				  fopen(dptr->idunion.acff.fname, "w")) == NULL){
		perror(dptr->idunion.acff.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************/

void
setup_acff()
{
	idstruct	*dptr,*serptr,*bptr,*modptr;
        extern parmstruct *prm;
	int i,j,k,icount,ifound,mm;
        int resnum1,resnum2,resnum3,resnum4,resnumt;
        char atnam1[2*ATOMTOKLEN], atnam2[2*ATOMTOKLEN],
             atnam3[2*ATOMTOKLEN], atnam4[2*ATOMTOKLEN],*atnamt;
        _REAL **r1t,**r2t,**r3t,**r4t,rNH;
        _REAL **tmp;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != ACFF)
                        inerr("setupacff: expected id or keyword\n","");
        } else {
                printf("setupacff: id %s not found in idlist\n", tok);
        }

        gettok();
        if ((bptr = getstruct(tok)) != NULL) {
                if (bptr->id_type != BEAD)
                        inerr("setupacff: expected id or keyword\n","");
        } else {
                printf("setupacff: id %s not found in idlist\n", tok);
        }

        gettok();
        if ((serptr = getstruct(tok)) == NULL) 
                inerr("setupacff: no such id\n","");
        if (serptr->id_type != STREAM)
                inerr("setupacff: stream id expected\n","");

        gettok();
        if ((modptr = getstruct(tok)) == NULL) 
                inerr("setupacff: no such id\n","");
        if (   modptr->id_type != MCARIG 
            && modptr->id_type != RIG
            && modptr->id_type != MCA_P2P1)
                inerr("setupacff: some model id expected\n","");

        dptr->idunion.acff.model = modptr;
        dptr->idunion.acff.prm = serptr->idunion.ser.prm;

        gettok();
        if(!strcmp(tok,"BOND")) {

        dptr->idunion.acff.typo=0;

        rNH=1.09;
        r1t = (_REAL **) get2(MAXJ,sizeof(long));
        r2t = (_REAL **) get2(MAXJ,sizeof(long));
        tmp = matrix(1,MAXJ,0,2);

        icount=0;
        ifound=0;
        for(;;) {
         gettok();
         if (tok[0] == ';') {
          break;
         } else if(!strcmp(tok,"D2mn")) {
          dptr->idunion.acff.typo=3;
         } else if(!strcmp(tok,"DOMEGA")) {
          gettok();
          tokREAL(&dptr->idunion.acff.domega);
         } else if(!strcmp(tok,"PRINTJ")) {
          dptr->idunion.acff.printj=1;
         } else if(!strcmp(tok,"FT")) {
          dptr->idunion.acff.type=1;
          gettok();
          if(!strcmp(tok,"NUMERIC")) {
           dptr->idunion.acff.type=2;
           gettok();
          }
          tokint(&dptr->idunion.acff.nomega,tok);
          for(i=0;i<dptr->idunion.acff.nomega;i++) {
           gettok();
           tokREAL(&dptr->idunion.acff.omega0[i],tok);
          }

         } else {
          icount++;
          if(icount>MAXJ) {
           printf("setup_acff: tmp space exceeded.\n");
           exit(1);
          }

          if(icount>1) {
           r1t++;
           r2t++;
          }
          sprintf(atnam1,"%s    ",tok);
          gettok();
          tokint(&resnum1,tok);
          /* printf("atnam1: %s %d\n",atnam1,resnum1); */
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           atnamt = &serptr->idunion.ser.prm->AtomNames[ATOMTOKLEN*i];
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           if (!strncmp(atnam1,atnamt,ATOMTOKLEN)&&resnum1==resnumt) {
            *r1t = serptr->idunion.ser.set->crd;
            /* printf("atom found:\n"); */
            if(atnam1[0]=='H') {
             tmp[icount][0] =  2.6753;
            } else if(atnam1[0]=='C') {
             tmp[icount][0] =  0.673;
            } else if(atnam1[0]=='A') {
             tmp[icount][0] =  0.673;
            } else if(atnam1[0]=='O') {
             tmp[icount][0] =  0.673;
            } else if(atnam1[0]=='N') {
             tmp[icount][0] = -0.2712;
            } else if(isdigit(atnam1[0])){
             tmp[icount][0] =  2.6753; /* assumed to be H */
            } else {
             printf("setup_acff: nucleus 1 not found.\n");
             exit(1);
            }
            ifound++;
            /* printf("atom1: %s %d\n",atnam1,resnum1); */
            /* printf("atomt: %s %d\n",atnamt,resnumt); */
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;

          gettok();
          sprintf(atnam2,"%s    ",tok);
          gettok();
          tokint(&resnum2,tok);
          /* printf("atnam2: %s %d\n",atnam2,resnum2); */
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           atnamt = &serptr->idunion.ser.prm->AtomNames[ATOMTOKLEN*i];
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           if (!strncmp(atnam2,atnamt,ATOMTOKLEN)&&resnum2==resnumt) {
            *r2t = serptr->idunion.ser.set->crd;
            /* printf("atom found:\n"); */
            if(atnam2[0]=='H') {
             tmp[icount][1] =  2.6753;
            } else if(atnam2[0]=='C') {
             tmp[icount][1] =  0.673;
            } else if(atnam2[0]=='A') {
             tmp[icount][1] =  0.673;
            } else if(atnam2[0]=='O') {
             tmp[icount][1] =  0.673;
            } else if(atnam2[0]=='N') {
             tmp[icount][1] = -0.2712;
            } else if(isdigit(atnam2[0])){
             tmp[icount][1] =  2.6753; /* assumed to be H */
            } else {
             printf("setup_acff: nucleus 2 not found.\n");
             exit(1);
            }
            ifound++;
            /* printf("atom2: %s %d\n",atnam2,resnum2); */
            /* printf("atomt: %s %d\n",atnamt,resnumt); */
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;

          gettok();
          tokREAL(&tmp[icount][2],tok);
         }

        }
        r1t-=icount-1;
        r2t-=icount-1;

        if(icount!=ifound/2) {
         printf("setup_acff: atom not found.\n");
         exit(1);
        }
        dptr->idunion.acff.nbond = icount;

        dptr->idunion.acff.r1 = (_REAL **)
            get2(icount,sizeof(long));
        dptr->idunion.acff.r2 = (_REAL **)
            get2(icount,sizeof(long));
        dptr->idunion.acff.gam1 = (_REAL *)
            get2(icount,sizeof(_REAL));
        dptr->idunion.acff.gam2 = (_REAL *)
            get2(icount,sizeof(_REAL));
        dptr->idunion.acff.omegaN = matrix(0,icount-1,0,199);
        dptr->idunion.acff.omegaH = matrix(0,icount-1,0,199);
        dptr->idunion.acff.gNgH = (_REAL *)
            get2(icount,sizeof(_REAL));
        dptr->idunion.acff.csa = (_REAL *)
            get2(icount,sizeof(_REAL));
        for(i=1;i<=icount;i++) {
         *dptr->idunion.acff.r1++ = *r1t++;
         *dptr->idunion.acff.r2++ = *r2t++;
         *dptr->idunion.acff.gam1++ = tmp[i][0];
         *dptr->idunion.acff.gam2++ = tmp[i][1];
/*
         *dptr->idunion.acff.gNgH =  
          (tmp[i][0]*tmp[i][1])*(6.62618/(2.*_PI))/pow(prm->sig0r,3.);
*/
         *dptr->idunion.acff.gNgH =  
          (tmp[i][0]*tmp[i][1])*(6.62618/(2.*_PI));
         *dptr->idunion.acff.gNgH *= *dptr->idunion.acff.gNgH; 
         *dptr->idunion.acff.gNgH *= 10.;
         printf("J factor(%d)= %e\n",i,*dptr->idunion.acff.gNgH);
         dptr->idunion.acff.gNgH++;
         *dptr->idunion.acff.csa++ = tmp[i][2]*tmp[i][2]/15./1000.;
        }
        r1t-=icount;
        r2t-=icount;
/*
        free(r1t);
        free(r2t);
*/
        dptr->idunion.acff.r1-=icount;
        dptr->idunion.acff.r2-=icount;
        dptr->idunion.acff.gam1-=icount;
        dptr->idunion.acff.gam2-=icount;
        dptr->idunion.acff.gNgH-=icount;
        dptr->idunion.acff.csa-=icount;

   mm=5;
   dptr->idunion.acff.Ssq = (_REAL *)
     get2(mm*dptr->idunion.acff.nbond,sizeof(_REAL));
   for(i=0;i<mm*dptr->idunion.acff.nbond;i++)
    dptr->idunion.acff.Ssq[i]=0.;

        } else if (!strcmp(tok,"DIHE")) {

        dptr->idunion.acff.typo=1;
 
        r1t = (_REAL **) get2(MAXJ,sizeof(long));
        r2t = (_REAL **) get2(MAXJ,sizeof(long));
        r3t = (_REAL **) get2(MAXJ,sizeof(long));
        r4t = (_REAL **) get2(MAXJ,sizeof(long));

        icount=0;
        ifound=0;
        for(;;) {
         gettok();
         if (tok[0] == ';') {
          break;
         } else if(!strcmp(tok,"DOMEGA")) {
          gettok();
          tokREAL(&dptr->idunion.acff.domega);
         } else {
          icount++;

          if(icount>1) {
           r1t++;
           r2t++;
           r3t++;
           r4t++;
          }
          sprintf(atnam1,"%s    ",tok);
          gettok();
          tokint(&resnum1,tok);
          /* printf("atnam1: %s %d\n",atnam1,resnum1); */
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           atnamt = &serptr->idunion.ser.prm->AtomNames[ATOMTOKLEN*i];
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           if (!strncmp(atnam1,atnamt,ATOMTOKLEN)&&resnum1==resnumt) {
            *r1t = serptr->idunion.ser.set->crd;
            /* printf("atom found:\n"); */
            ifound++;
            /* printf("atom1: %s %d\n",atnam1,resnum1); */
            /* printf("atomt: %s %d\n",atnamt,resnumt); */
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;

          gettok();
          sprintf(atnam2,"%s    ",tok);
          gettok();
          tokint(&resnum2,tok);
          /* printf("atnam2: %s %d\n",atnam2,resnum2); */
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           atnamt = &serptr->idunion.ser.prm->AtomNames[ATOMTOKLEN*i];
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           if (!strncmp(atnam2,atnamt,ATOMTOKLEN)&&resnum2==resnumt) {
            *r2t = serptr->idunion.ser.set->crd;
            /* printf("atom found:\n"); */
            ifound++;
            /* printf("atom2: %s %d\n",atnam2,resnum2); */
            /* printf("atomt: %s %d\n",atnamt,resnumt); */
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;

          gettok();
          sprintf(atnam3,"%s    ",tok);
          gettok();
          tokint(&resnum3,tok);
          /* printf("atnam3: %s %d\n",atnam3,resnum3); */
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           atnamt = &serptr->idunion.ser.prm->AtomNames[ATOMTOKLEN*i];
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           if (!strncmp(atnam3,atnamt,ATOMTOKLEN)&&resnum3==resnumt) {
            *r3t = serptr->idunion.ser.set->crd;
            /* printf("atom found:\n"); */
            ifound++;
            /* printf("atom3: %s %d\n",atnam3,resnum3); */
            /* printf("atomt: %s %d\n",atnamt,resnumt); */
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;

          gettok();
          sprintf(atnam4,"%s    ",tok);
          gettok();
          tokint(&resnum4,tok);
          /* printf("atnam4: %s %d\n",atnam4,resnum4); */
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           atnamt = &serptr->idunion.ser.prm->AtomNames[ATOMTOKLEN*i];
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           if (!strncmp(atnam4,atnamt,ATOMTOKLEN)&&resnum4==resnumt) {
            *r4t = serptr->idunion.ser.set->crd;
            /* printf("atom found:\n"); */
            ifound++;
            /* printf("atom4: %s %d\n",atnam4,resnum4); */
            /* printf("atomt: %s %d\n",atnamt,resnumt); */
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;

         }

        }
        r1t-=icount-1;
        r2t-=icount-1;
        r3t-=icount-1;
        r4t-=icount-1;

        if(icount!=ifound/4) {
         printf("setup_acff: atom not found.\n");
         exit(1);
        }
        dptr->idunion.acff.nbond = icount;

        dptr->idunion.acff.r1 = (_REAL **)
            get2(icount,sizeof(long));
        dptr->idunion.acff.r2 = (_REAL **)
            get2(icount,sizeof(long));
        dptr->idunion.acff.r3 = (_REAL **)
            get2(icount,sizeof(long));
        dptr->idunion.acff.r4 = (_REAL **)
            get2(icount,sizeof(long));
        for(i=1;i<=icount;i++) {
         *dptr->idunion.acff.r1++ = *r1t++;
         *dptr->idunion.acff.r2++ = *r2t++;
         *dptr->idunion.acff.r3++ = *r3t++;
         *dptr->idunion.acff.r4++ = *r4t++;
        }
        r1t-=icount;
        r2t-=icount;
        r3t-=icount;
        r4t-=icount;
/*
        free(r1t);
        free(r2t);
        free(r3t);
        free(r4t);
*/
        dptr->idunion.acff.r1-=icount;
        dptr->idunion.acff.r2-=icount;
        dptr->idunion.acff.r3-=icount;
        dptr->idunion.acff.r4-=icount;

        } else if (!strcmp(tok,"NOESY")) {

        dptr->idunion.acff.typo=2;

        r1t = (_REAL **) get2(MAXJ,sizeof(long));
        r2t = (_REAL **) get2(MAXJ,sizeof(long));
        tmp = matrix(1,MAXJ,0,3);

        icount=0;
        ifound=0;
        for(;;) {
         gettok();
         if (tok[0] == ';') {
          break;
         } else if(!strcmp(tok,"DOMEGA")) {
          gettok();
          tokREAL(&dptr->idunion.acff.domega);
         } else if(!strcmp(tok,"PRINTJ")) {
          dptr->idunion.acff.printj=1;
         } else if(!strcmp(tok,"FT")) {
          dptr->idunion.acff.type=1;
          gettok();
          if(!strcmp(tok,"NUMERIC")) {
           dptr->idunion.acff.type=2;
           gettok();
          }
          tokint(&dptr->idunion.acff.nomega,tok);
          for(i=0;i<dptr->idunion.acff.nomega;i++) {
           gettok();
           tokREAL(&dptr->idunion.acff.omega0[i],tok);
          }
         } else {
          icount++;
          if(icount>MAXJ) {
           printf("setup_acff: tmp space exceeded.\n");
           exit(1);
          }

          if(icount>1) {
           r1t++;
           r2t++;
          }
          sprintf(atnam1,"%s    ",tok);
          gettok();
          tokint(&resnum1,tok);
          /* printf("atnam1: %s %d\n",atnam1,resnum1); */
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           atnamt = &serptr->idunion.ser.prm->AtomNames[ATOMTOKLEN*i];
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           if (!strncmp(atnam1,atnamt,ATOMTOKLEN)&&resnum1==resnumt) {
            *r1t = serptr->idunion.ser.set->crd;
            /* printf("atom found:\n"); */
            if(atnam1[0]=='H') {
             tmp[icount][0] =  2.6753;
            } else if(atnam1[0]=='C') {
             tmp[icount][0] =  0.673;
            } else if(atnam1[0]=='A') {
             tmp[icount][0] =  0.673;
            } else if(atnam1[0]=='N') {
             tmp[icount][0] = -0.2712;
            } else if(isdigit(atnam1[0])){
             tmp[icount][0] =  2.6753; /* assumed to be H */
            } else {
             printf("setup_acff: nucleus 1 not found.\n");
             exit(1);
            }
            ifound++;
            /* printf("atom1: %s %d\n",atnam1,resnum1); */
            /* printf("atomt: %s %d\n",atnamt,resnumt); */
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;

          gettok();
          sprintf(atnam2,"%s    ",tok);
          gettok();
          tokint(&resnum2,tok);
          /* printf("atnam2: %s %d\n",atnam2,resnum2); */
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           atnamt = &serptr->idunion.ser.prm->AtomNames[ATOMTOKLEN*i];
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           if (!strncmp(atnam2,atnamt,ATOMTOKLEN)&&resnum2==resnumt) {
            *r2t = serptr->idunion.ser.set->crd;
            /* printf("atom found:\n"); */
            if(atnam2[0]=='H') {
             tmp[icount][1] =  2.6753;
            } else if(atnam2[0]=='C') {
             tmp[icount][1] =  0.673;
            } else if(atnam2[0]=='A') {
             tmp[icount][1] =  0.673;
            } else if(atnam2[0]=='O') {
             tmp[icount][1] =  0.673;
            } else if(atnam2[0]=='N') {
             tmp[icount][1] = -0.2712;
            } else if(isdigit(atnam2[0])){
             tmp[icount][1] =  2.6753; /* assumed to be H */
            } else {
             printf("setup_acff: nucleus 2 not found.\n");
             exit(1);
            }
            ifound++;
            /* printf("atom2: %s %d\n",atnam2,resnum2); */
            /* printf("atomt: %s %d\n",atnamt,resnumt); */
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;

          gettok();
          tokREAL(&tmp[icount][2],tok);
          gettok();
          tokREAL(&tmp[icount][3],tok);
         }

        }
        r1t-=icount-1;
        r2t-=icount-1;

        if(icount!=ifound/2) {
         printf("setup_acff: atom not found.\n");
         exit(1);
        }
        dptr->idunion.acff.nbond = icount;

        dptr->idunion.acff.r1 = (_REAL **)
            get2(icount,sizeof(long));
        dptr->idunion.acff.r2 = (_REAL **)
            get2(icount,sizeof(long));
        dptr->idunion.acff.gam1 = (_REAL *)
            get2(icount,sizeof(_REAL));
        dptr->idunion.acff.gam2 = (_REAL *)
            get2(icount,sizeof(_REAL));
        dptr->idunion.acff.omegaN = matrix(0,icount-1,0,199);
        dptr->idunion.acff.omegaH = matrix(0,icount-1,0,199);
        dptr->idunion.acff.gNgH = (_REAL *)
            get2(icount,sizeof(_REAL));
        for(i=1;i<=icount;i++) {
         *dptr->idunion.acff.r1++ = *r1t++;
         *dptr->idunion.acff.r2++ = *r2t++;
         *dptr->idunion.acff.gam1++ = tmp[i][0];
         *dptr->idunion.acff.gam2++ = tmp[i][1];
/*
         *dptr->idunion.acff.gNgH =  
          (tmp[i][0]*tmp[i][1])*(6.62618/(2.*_PI))/pow(prm->sig0r,3.);
*/
         *dptr->idunion.acff.gNgH =  
          (tmp[i][0]*tmp[i][1])*(6.62618/(2.*_PI));
         *dptr->idunion.acff.gNgH *= *dptr->idunion.acff.gNgH; 
         *dptr->idunion.acff.gNgH *= 10.;
         printf("J factor(%d)= %e\n",i,*dptr->idunion.acff.gNgH);
         dptr->idunion.acff.gNgH++;
        }
        r1t-=icount;
        r2t-=icount;
/*
        free(r1t);
        free(r2t);
*/
        dptr->idunion.acff.r1-=icount;
        dptr->idunion.acff.r2-=icount;
        dptr->idunion.acff.gam1-=icount;
        dptr->idunion.acff.gam2-=icount;
        dptr->idunion.acff.gNgH-=icount;

   mm=5;
   dptr->idunion.acff.Ssq = (_REAL *)
     get2(mm*dptr->idunion.acff.nbond,sizeof(_REAL));
   for(i=0;i<mm*dptr->idunion.acff.nbond;i++)
    dptr->idunion.acff.Ssq[i]=0.;


        } else if (!strcmp(tok,"STRESS")) {
         dptr->idunion.acff.typo=4;
         for(;;) {
          gettok();
          if (tok[0] == ';') {
           break;
          }
         }
        } else {
         printf("setup_acff: type of object to correlate not found.\n");
         exit(1);
        }


   if(modptr->id_type==MCARIG) {

    dptr->idunion.acff.nbase = modptr->idunion.mcarig.n1;
    if(modptr->idunion.mcarig.sel_bonds!=0)
    {
     dptr->idunion.acff.nbase += modptr->idunion.mcarig.n2;
    }

    dptr->idunion.acff.btmp =
     (_REAL *)get2(3*dptr->idunion.acff.nbond,sizeof(_REAL));
    dptr->idunion.acff.tau_p2 =
      (_REAL *)get2(dptr->idunion.acff.nbond,sizeof(_REAL));
    for(i=0;i<dptr->idunion.acff.nbond;i++)
     dptr->idunion.acff.tau_p2[i]=0.;
    mm=5;

   } else if (modptr->id_type==MCA_P2P1) {

    dptr->idunion.acff.nbase = 5;
    dptr->idunion.acff.btmp =
     (_REAL *)get2(5*dptr->idunion.acff.nbond,sizeof(_REAL));
    mm=5;

   } else if (modptr->id_type==RIG) {

    dptr->idunion.acff.nbase = 5;
    dptr->idunion.acff.btmp =
     (_REAL *)get2(5*dptr->idunion.acff.nbond,sizeof(_REAL));
    mm=1;

   }

   dptr->idunion.acff.tmp =
     vector(1,dptr->idunion.acff.nbond);
   dptr->idunion.acff.modsq =
     vector(1,2*dptr->idunion.acff.nbond);
   for(i=1;i<=dptr->idunion.acff.nbond;i++) {
    dptr->idunion.acff.modsq[i]=0.;
   }
   if(dptr->idunion.acff.typo==3) {
    dptr->idunion.acff.proj =
      matrix(1,10*(dptr->idunion.acff.nbond-1),1,dptr->idunion.acff.nbase);
    for(i=1;i<=10*(dptr->idunion.acff.nbond-1);i++)
     for(j=1;j<=dptr->idunion.acff.nbase;j++)
      dptr->idunion.acff.proj[i][j]=0.;
   } else if(dptr->idunion.acff.typo==4) {
    dptr->idunion.acff.proj =
      matrix(1,1,1,dptr->idunion.acff.nbase);
    for(j=1;j<=dptr->idunion.acff.nbase;j++)
     dptr->idunion.acff.proj[1][j]=0.;
   } else {
    dptr->idunion.acff.proj =
      matrix(1,dptr->idunion.acff.nbond,1,dptr->idunion.acff.nbase);
    for(i=1;i<=dptr->idunion.acff.nbond;i++)
     for(j=1;j<=dptr->idunion.acff.nbase;j++)
      dptr->idunion.acff.proj[i][j]=0.;
   }
/*
   dptr->idunion.acff.projc =
     Cmatrix(1,dptr->idunion.acff.nbond,1,dptr->idunion.acff.nbase);
   for(i=1;i<=dptr->idunion.acff.nbond;i++)
    for(j=1;j<=dptr->idunion.acff.nbase;j++)
     dptr->idunion.acff.projc[i][j]=Complex(0.,0.);
*/
/*
   dptr->idunion.acff.proj2 =
     matrix(1,dptr->idunion.acff.nbond,1,dptr->idunion.acff.nbase);
   for(i=1;i<=dptr->idunion.acff.nbond;i++)
    for(j=1;j<=dptr->idunion.acff.nbase;j++)
     dptr->idunion.acff.proj2[i][j]=0.;
*/
   dptr->idunion.acff.domega *= 2.*_PI;
   if(dptr->idunion.acff.domega==0.)
   {
    printf("setup_acff: WARNING: domega is zero.\n");
    // exit(1);
   }
   dptr->idunion.acff.nt = 1000;
/*
   dptr->idunion.acff.nt = 32768;
*/

   dptr->idunion.acff.dtime = _PI/
        ((_REAL)dptr->idunion.acff.nt*dptr->idunion.acff.domega);
   dptr->idunion.acff.dtime /= bptr->idunion.bead.factns;

   if(dptr->idunion.acff.type>0) {
 
    dptr->idunion.acff.relax_Sz =
     matrix(1,dptr->idunion.acff.nbond,1,dptr->idunion.acff.nomega);
    dptr->idunion.acff.relax_Sxy =
     matrix(1,dptr->idunion.acff.nbond,1,dptr->idunion.acff.nomega);
    dptr->idunion.acff.relax_Iz_Sz =
     matrix(1,dptr->idunion.acff.nbond,1,dptr->idunion.acff.nomega);
    dptr->idunion.acff.nOe_Iz_Sz =
     matrix(1,dptr->idunion.acff.nbond,1,dptr->idunion.acff.nomega);

    if(dptr->idunion.acff.type==2) {
     dptr->idunion.acff.nzero = 20000;
     dptr->idunion.acff.epsilon = 0.3;
     dptr->idunion.acff.T = dptr->idunion.acff.dtime*
         (_REAL)dptr->idunion.acff.nt;
#ifdef IMSL
     dptr->idunion.acff.workft = imsl_d_fft_real_init(dptr->idunion.acff.nt);
#endif
     dptr->idunion.acff.p2t = vector(1,dptr->idunion.acff.nt);
    }
    dptr->idunion.acff.domega *= 2.;

    dptr->idunion.acff.gamH =  2.6753;
/*
    dptr->idunion.acff.gNgH =  
     (dptr->idunion.acff.gamN*dptr->idunion.acff.gamH)*(6.62618/(2.*_PI))/
      pow(rNH,3.);
    printf("rNH= %e\n",pow(rNH,3.));
*/

    if( dptr->idunion.acff.typo==2)
    {

    for(i=0;i<dptr->idunion.acff.nomega;i++)
    {
     dptr->idunion.acff.omega0[i] *=  2.*_PI*1.e-03;
     for(j=0;j<dptr->idunion.acff.nbond;j++)
     {
      dptr->idunion.acff.omegaN[j][i] =
       dptr->idunion.acff.gam1[j]*dptr->idunion.acff.omega0[i]*
        (1.-tmp[j+1][2]*1.e-06)
        /dptr->idunion.acff.gamH  ;
      dptr->idunion.acff.omegaH[j][i] =
       dptr->idunion.acff.gam2[j]*dptr->idunion.acff.omega0[i]*
        (1.-tmp[j+1][3]*1.e-06)
        /dptr->idunion.acff.gamH  ;
     }
    }

   } else {

    for(i=0;i<dptr->idunion.acff.nomega;i++) {
     dptr->idunion.acff.omega0[i] *=  2.*_PI*1.e-03;
     for(j=0;j<dptr->idunion.acff.nbond;j++) {
      dptr->idunion.acff.omegaN[j][i] =
       dptr->idunion.acff.gam1[j]*dptr->idunion.acff.omega0[i]
        /dptr->idunion.acff.gamH  ;
      dptr->idunion.acff.omegaH[j][i] =
       dptr->idunion.acff.gam2[j]*dptr->idunion.acff.omega0[i]
        /dptr->idunion.acff.gamH  ;
     }
    }

   }

   }

   printf("setup_acff: exiting...\n");
   fflush(stdout);
}

/************************************************************************/
void
adjust_acff(idstruct *dptr)
{
   idstruct	*modptr;
   int i,j;

/*
   free_matrix(dptr->idunion.acff.proj,
     1,dptr->idunion.acff.nbond,1,dptr->idunion.acff.nbase);
*/
   modptr=dptr->idunion.acff.model;

   if(modptr->id_type==MCA_P2P1) {
    if(modptr->idunion.mca_p2p1.type==1) {
     dptr->idunion.acff.nbase = modptr->idunion.mca_p2p1.n12;
    } else if(modptr->idunion.mca_p2p1.type==2) {
     dptr->idunion.acff.nbase = modptr->idunion.mca_p2p1.n12
                               +modptr->idunion.mca_p2p1.n1234;
    } else if(modptr->idunion.mca_p2p1.type==3) {
     dptr->idunion.acff.nbase = modptr->idunion.mca_p2p1.n12
                               +modptr->idunion.mca_p2p1.n1234
                               +modptr->idunion.mca_p2p1.n123456;
    }
   } else if(modptr->id_type==MCARIG) {
    printf("adjust_acff: case no implemented yet.\n");
   }

   printf("adjust_acff: new nbase= %d\n",dptr->idunion.acff.nbase);
   printf("adjust_acff: new nbase2= %d\n",dptr->idunion.acff.nbase2);
   fflush(stdout);

   if(modptr->id_type==MCA_P2P1 && dptr->idunion.acff.typo==4)
   { 
    dptr->idunion.acff.proj =
      matrix(1,1,1,dptr->idunion.acff.nbase);
    for(j=1;j<=dptr->idunion.acff.nbase;j++)
       dptr->idunion.acff.proj[1][j]=0.;
   }
   else 
   {
    dptr->idunion.acff.proj =
      matrix(1,dptr->idunion.acff.nbond,1,dptr->idunion.acff.nbase);
    for(i=1;i<=dptr->idunion.acff.nbond;i++)
     for(j=1;j<=dptr->idunion.acff.nbase;j++)
       dptr->idunion.acff.proj[i][j]=0.;
    for(i=1;i<=dptr->idunion.acff.nbond;i++)
     dptr->idunion.acff.modsq[i]=0.;
   }

   dptr->idunion.acff.cases=0;

}


/************************************************************************/
void
do_acff2_stress(idstruct *dptr, idstruct *serptr)
{
   idstruct *modptr=dptr->idunion.acff.model;
   _REAL bt[9],stress[5];
   int i,j,j2,k,mm;
   _REAL b,c,scp;
   _REAL weight;

   if(serptr->idunion.ser.iweight)
   {
    weight=serptr->idunion.ser.w->idunion.serw.w;
   } else {
    weight=1.;
   }
   mm=5;
   b=l_sqrt(3./8.);
   c=0.5;
/*
   b=0.5;
   c=1./l_sqrt(6.);
*/

   ++dptr->idunion.acff.cases;
   dptr->idunion.acff.casesw+=weight;

// getting stress tensor from weight into bt ...
   if(serptr->idunion.ser.iweight)
   {
    bt[0]=serptr->idunion.ser.w->idunion.serw.w2[0]; // x^2
    bt[1]=serptr->idunion.ser.w->idunion.serw.w2[1]; // xy
    bt[2]=serptr->idunion.ser.w->idunion.serw.w2[2]; // xz
    bt[3]=serptr->idunion.ser.w->idunion.serw.w2[3]; // yx
    bt[4]=serptr->idunion.ser.w->idunion.serw.w2[4]; // y^2
    bt[5]=serptr->idunion.ser.w->idunion.serw.w2[5]; // yz
    bt[6]=serptr->idunion.ser.w->idunion.serw.w2[6]; // zx
    bt[7]=serptr->idunion.ser.w->idunion.serw.w2[7]; // zy
    bt[8]=serptr->idunion.ser.w->idunion.serw.w2[8]; // z^2
// from cartesian to 2nd rank spherical components:
    scp=bt[0]+bt[4]+bt[8];
    stress[0]=c*(3.*bt[8] - scp);
    stress[1]=b*(bt[0] - bt[4]);
    stress[2]=b*(bt[1] + bt[3]);
    stress[3]=b*(bt[2] + bt[6]);
    stress[4]=b*(bt[5] + bt[7]);
/*
    scp=stress[0]*stress[0];
    scp+=2.*stress[1]*stress[1];
    scp+=2.*stress[2]*stress[2];
    scp+=2.*stress[3]*stress[3];
    scp+=2.*stress[4]*stress[4];
    scp*=0.2;
    printf("WW2 %e\n",scp);
*/
   } else {
    printf("do_acff2_stress: no stress in weight.\n");
    fflush(stdout);
    exit(1);
   }

    for (j=1;j<=modptr->idunion.mca_p2p1.n12;j++) 
    {
     scp=0.;
     scp+= modptr->idunion.mca_p2p1.l12->phi[0]*stress[0];
     scp+= 2.*modptr->idunion.mca_p2p1.l12->phi[1]*stress[1];
     scp+= 2.*modptr->idunion.mca_p2p1.l12->phi[2]*stress[2];
     scp+= 2.*modptr->idunion.mca_p2p1.l12->phi[3]*stress[3];
     scp+= 2.*modptr->idunion.mca_p2p1.l12->phi[4]*stress[4];
     scp*= 0.2;
     dptr->idunion.acff.proj[1][j] += weight*scp;
     // printf("scp= %e weight= %e\n",scp,weight); fflush(stdout);
     /* dptr->idunion.acff.proj[i][j] += weight*lsq*scp; */
     modptr->idunion.mca_p2p1.l12++;
    }
    modptr->idunion.mca_p2p1.l12-=modptr->idunion.mca_p2p1.n12;

    if(modptr->idunion.mca_p2p1.type>1) 
    {

    for (j=1;j<=modptr->idunion.mca_p2p1.n1234;j++) 
    {
     j2=modptr->idunion.mca_p2p1.n12+j;
     scp=0.;
     scp+= modptr->idunion.mca_p2p1.l1234->phi[0]*stress[0];
     scp+= 2.*modptr->idunion.mca_p2p1.l1234->phi[1]*stress[1];
     scp+= 2.*modptr->idunion.mca_p2p1.l1234->phi[2]*stress[2];
     scp+= 2.*modptr->idunion.mca_p2p1.l1234->phi[3]*stress[3];
     scp+= 2.*modptr->idunion.mca_p2p1.l1234->phi[4]*stress[4];
     scp*= 0.2;
     dptr->idunion.acff.proj[1][j2] += weight*scp;
     modptr->idunion.mca_p2p1.l1234++;
    }
    modptr->idunion.mca_p2p1.l1234-=modptr->idunion.mca_p2p1.n1234;

    }

    if(modptr->idunion.mca_p2p1.type>2) 
    {

    for (j=1;j<=modptr->idunion.mca_p2p1.n123456;j++) 
    {
     j2=modptr->idunion.mca_p2p1.n12+modptr->idunion.mca_p2p1.n1234+j;
     scp=0.;
     scp+= modptr->idunion.mca_p2p1.l123456->phi[0]*stress[0];
     scp+= 2.*modptr->idunion.mca_p2p1.l123456->phi[1]*stress[1];
     scp+= 2.*modptr->idunion.mca_p2p1.l123456->phi[2]*stress[2];
     scp+= 2.*modptr->idunion.mca_p2p1.l123456->phi[3]*stress[3];
     scp+= 2.*modptr->idunion.mca_p2p1.l123456->phi[4]*stress[4];
     scp*= 0.2;
     dptr->idunion.acff.proj[1][j2] += weight*scp;
     /* dptr->idunion.acff.proj[i][j2] += weight*lsq*scp; */
     modptr->idunion.mca_p2p1.l123456++;
    }
    modptr->idunion.mca_p2p1.l123456-=modptr->idunion.mca_p2p1.n123456;

    }


}



/************************************************************************/
void
write_acff_stress(idstruct *dptr)
{
 extern parmstruct *prm;
 idstruct *modptr=dptr->idunion.acff.model;
 _REAL **proj=dptr->idunion.acff.proj;
 _REAL *eval;
 _REAL **evec,**cc;
 _REAL tau,fact,sum,norm;
 _REAL factns;
 int i,j,k,l,nf,mm,bytes;
 FILE *file=dptr->idunion.acff.file;

/* fact=1./((_REAL)dptr->idunion.acff.cases); */
 fact=1./dptr->idunion.acff.casesw;
 mm=5;

 if(dptr->idunion.acff.model->id_type==MCA_P2P1) 
 {
  nf = 1;
  mm=5;
  eval=modptr->idunion.mca_p2p1.eval;
  evec=modptr->idunion.mca_p2p1.evec;
 }
 else 
 {
  printf("write_acff_stress: model not implemented yet.\n");
  fflush(stdout);
  exit(1);
 }
 norm=1./((_REAL)mm);

 printf("write_acff_stress: number of basis functions, nbase= %d\n",dptr->idunion.acff.nbase);
 printf("write_acff_stress: nuber of functions, nf= %d\n",nf);
 printf("write_acff_stress: number of 2nd-rank components, mm= %d\n",mm);
 fflush(stdout);

 for(j=1;j<=dptr->idunion.acff.nbase;j++) {
  proj[1][j]*=fact;
 }
 bytes=sizeof(_REAL)*dptr->idunion.acff.nbase;
 fwrite(&bytes,sizeof(int),1,file);
 fwrite(&proj[1][1],dptr->idunion.acff.nbase*sizeof(_REAL),1,file);
 fwrite(&bytes,sizeof(int),1,file);
 fclose(file);
/*
 for(j=1;j<=dptr->idunion.acff.nbase;j++) {
  fprintf(file,"%e\n",proj[1][j]);
 }
*/
 return;

}
