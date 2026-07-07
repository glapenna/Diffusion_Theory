#include "head.h"

// extern char		tok[];
// char		tok[];
// void *get2(int,int);

rdiadstruct **rdiad_matrix();
_REAL **matrix();

/**** extern ****/
extern idstruct		*startid(),*getstruct();

/************************************************************************/
void
setup_bead_fric()
{
	idstruct	*dptr,*serptr;
        char id[TOKLEN];
	int i,j,k,icount,ifound;
        int resnum1,resnumt;
        _REAL z,temp,eta,eps0r,mass0r,sig0r;
        char atnam1[ATOMTOKLEN],*atnamt;

	dptr = startid(BEAD, id);

	dptr->idunion.bead.nbond = 0;
	dptr->idunion.bead.nbead = 0;
	dptr->idunion.bead.hydrod = OSEEN;
	dptr->idunion.bead.screen = 0.;
	dptr->idunion.bead.zr = -1.;
	dptr->idunion.bead.bav = 0.;
	dptr->idunion.bead.bavt = 0.;
	dptr->idunion.bead.bavt2 = 0.;
	dptr->idunion.bead.lambda = 1.;
	dptr->idunion.bead.cutsq = 0.0001;
	dptr->idunion.bead.stick_slip = 6.;
	eta = 0.;

/* get id */
        gettok();
        if ((serptr = getstruct(tok)) == NULL) 
                inerr("setup_bead_fric: no such id\n","");
        if (serptr->id_type != STREAM)
                inerr("setup_bead_fric: stream id expected\n","");

        printf("setup_bead_fric: Nsit= %d\n",serptr->idunion.ser.prm->Nsit); fflush(stdout);
        dptr->idunion.bead.zt = (_REAL *)
         get2(serptr->idunion.ser.prm->Nsit,sizeof(_REAL));

        for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
         *dptr->idunion.bead.zt++=0.;
        }
        dptr->idunion.bead.zt-=serptr->idunion.ser.prm->Nsit;
        serptr->idunion.ser.prm->fric=(_REAL *)
         get2(serptr->idunion.ser.prm->Nsit,sizeof(_REAL));
        for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
         serptr->idunion.ser.prm->fric[i]=0.;
        }

        icount=0;
        ifound=0;
        for(;;) {
         gettok();
         if (tok[0] == ';') {
          break;
         } else if(!strcmp(tok,"TEMP")) {
          gettok();
          tokREAL(&temp,tok);
         } else if(!strcmp(tok,"SOLVISC")) {
          gettok();
          tokREAL(&eta,tok);
         } else if(!strcmp(tok,"SLIP")) {
          dptr->idunion.bead.stick_slip=4.;
         } else if(!strcmp(tok,"ZR")) {
          gettok();
          tokREAL(&dptr->idunion.bead.zr,tok);
         } else if(!strcmp(tok,"OSEEN")) {
          dptr->idunion.bead.hydrod=OSEEN;
         } else if(!strcmp(tok,"SCREEN")) {
          gettok();
          tokREAL(&dptr->idunion.bead.screen,tok);
         } else if(!strcmp(tok,"RP")) {
          dptr->idunion.bead.hydrod=ROTNE_PRAGER;
         } else if(!strcmp(tok,"OTTINGER")) {
          dptr->idunion.bead.hydrod=OTTINGER;
         } else if(!strcmp(tok,"LAMBDA")) {
          gettok();
          tokREAL(&dptr->idunion.bead.lambda,tok);
         } else if(!strcmp(tok,"RCUT")) {
          gettok();
          tokREAL(&dptr->idunion.bead.cutsq,tok);
          dptr->idunion.bead.cutsq *= dptr->idunion.bead.cutsq;
         } else if(!strcmp(tok,"BAV")) {
          gettok();
          tokREAL(&dptr->idunion.bead.bav,tok);
         } else {
          icount++;
          // sprintf(atnam1,"%s    ",tok);
          // strncpy(atnam1,tok,TOKLEN);
          strcpy(atnam1,tok);

          gettok();
          tokint(&resnum1,tok);
          // printf("atnam1: %s %d\n",atnam1,resnum1);
          gettok();
          tokREAL(&z,tok);
/* find atomnames in list */
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           atnamt = &serptr->idunion.ser.prm->AtomNames[ATOMTOKLEN*i];
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           if (!strncmp(atnam1,atnamt,ATOMTOKLEN)&&resnum1==resnumt) {
/*
            { _REAL err,dev,drand48(void );
              err=0.025;
              dev=(1.-2.*drand48())*err;
              z = z*(1.+dev);
            }
*/
            dptr->idunion.bead.zt[i]=z;

            ifound++;
/*
             printf("atom found: %d %d\n",icount,ifound);
             printf("atom1: %s %d\n",atnam1,resnum1);
             printf("atomt: %s %d\n",atnamt,resnumt);
*/
           }
          }
         }
        }

        if(icount>0) {
         if(icount!=ifound) {
          printf("setup_bead_fric: atom not found.\n");
          printf("icount= %d, ifound= %d\n",icount,ifound);
          exit(1);
         }
         dptr->idunion.bead.nbead = icount;
        } else {
         printf("setup_bead_fric: you must give some beads.\n");
         exit(1);
        }
        if(eta == 0.) {
         printf("setup_bead_fric: you must give eta.\n");
         exit(1);
        }
        if(dptr->idunion.bead.zr<0.) {
         printf("setup_bead_fric: you must give zr.\n");
         exit(1);
        }

        eps0r = serptr->idunion.ser.prm->eps0r;
        mass0r = serptr->idunion.ser.prm->mass0r;
        sig0r = serptr->idunion.ser.prm->sig0r;
/*
*/
        if(serptr->idunion.ser.scale==1) {
         sig0r=1.;
        } else if(serptr->idunion.ser.scale==2) {
         if(dptr->idunion.bead.bav>0.) {
          dptr->idunion.bead.bav /= sig0r;
         }
         for(i=0; i<serptr->idunion.ser.prm->Nsit;i++) {
          dptr->idunion.bead.zt[i] /= sig0r;
         }
        }
        for(i=0,dptr->idunion.bead.zav=0.,j=0;
             i<serptr->idunion.ser.prm->Nsit;i++) {
         if(dptr->idunion.bead.zt[i]!=0.) {
          j++;
          dptr->idunion.bead.zav += dptr->idunion.bead.zt[i];
         }
        }
        dptr->idunion.bead.zav /= (_REAL)j;
/*
        for(i=0; i<serptr->idunion.ser.prm->Nsit;i++) {
         dptr->idunion.bead.zt[i] *= 2.0;
        }
*/
        for(i=0; i<serptr->idunion.ser.prm->Nsit;i++) {
         serptr->idunion.ser.prm->fric[i]=dptr->idunion.bead.zt[i];
        }

/*
        dptr->idunion.bead.zr *= dptr->idunion.bead.stick_slip*
         dptr->idunion.bead.bav/dptr->idunion.bead.zav/8.;
*/
/*
*/
        dptr->idunion.bead.zr *= dptr->idunion.bead.stick_slip/8.;

        printf("zr  = %e\n",dptr->idunion.bead.zr);
        printf("zav = %e\n",dptr->idunion.bead.zav);
        printf("bav = %e\n",dptr->idunion.bead.bav);
        printf("sig0r = %e\n",serptr->idunion.ser.prm->sig0r);
        printf("cutsq = %e\n",dptr->idunion.bead.cutsq);

/* store bead information */
        dptr->idunion.bead.z = (_REAL *)
         get2(dptr->idunion.bead.nbead,sizeof(_REAL));
        dptr->idunion.bead.r = (_REAL **) 
         get2(dptr->idunion.bead.nbead,sizeof(long));
        for(i=0,j=0; i<serptr->idunion.ser.prm->Nsit;i++) {
         if(dptr->idunion.bead.zt[i]!=0.) {
          j++;
          *dptr->idunion.bead.z++ = dptr->idunion.bead.zt[i];
          *dptr->idunion.bead.r++ = serptr->idunion.ser.set->crd;
         }
         serptr->idunion.ser.set->crd += 3;
        }
        if(j!=dptr->idunion.bead.nbead) {
         printf("bead: wrong number of beads: j= %d, and nbead= %d.\n",
                 j,dptr->idunion.bead.nbead);
         exit(1);
        }
        dptr->idunion.bead.z -= dptr->idunion.bead.nbead;
        dptr->idunion.bead.r -= dptr->idunion.bead.nbead;
        serptr->idunion.ser.set->crd -= 3*serptr->idunion.ser.prm->Nsit;

/* stick boundary */
        dptr->idunion.bead.sigma = 1.e+07*(BOLTZ*temp/AVOG)/
          (dptr->idunion.bead.stick_slip*_PI*eta*sig0r*sig0r*sig0r);
        dptr->idunion.bead.factns = 1.e+09/dptr->idunion.bead.sigma;
        printf("factns = %e\n",dptr->idunion.bead.factns);
}
/************************************************************************/

void
setup_bead_bond()
{
	idstruct	*dptr,*serptr;
	int i,j,k,icount,icount_max,ifound;
        int resnum1,resnum2,resnumt;
        char atnam1[ATOMTOKLEN], atnam2[ATOMTOKLEN],*atnamt;
        _REAL **r1t,
              **r2t,
              *z1t,
              *z2t;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != BEAD)
                        inerr("setup_bead_bond: expected id or keyword\n","");
        } else {
                printf("setup_bead_bond: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((serptr = getstruct(tok)) == NULL) 
                inerr("setup_bead_bond: no such id\n","");
        if (serptr->id_type != STREAM)
                inerr("setup_bead_bond: stream id expected\n","");

        dptr->idunion.bead.prm = serptr->idunion.ser.prm;
        icount_max=15000;
        r1t = (_REAL **) get2(icount_max,sizeof(long));
        r2t = (_REAL **) get2(icount_max,sizeof(long));
        z1t = (_REAL *) get2(icount_max,sizeof(_REAL));
        z2t = (_REAL *) get2(icount_max,sizeof(_REAL));

        icount=0;
        ifound=0;
        for(;;) {
         gettok();
         if (tok[0] == ';') {
          break;
         } else {
          icount++;
          if(icount>icount_max) {
           printf("setup_bead_bond: icount_max exceeded.\n");
           exit(1);
          }

/*
          if(icount>1) {
           r1t++;
           r2t++;
           z1t++;
           z2t++;
          }
*/
          // sprintf(atnam1,"%s    ",tok);
          strcpy(atnam1,tok);
          gettok();
          tokint(&resnum1,tok);

          atnamt = serptr->idunion.ser.prm->AtomNames;
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           // printf("atomt: %.4s %d\n",atnamt,resnumt);
           if (!strncmp(atnam1,atnamt,ATOMTOKLEN)&&resnum1==resnumt) {
            r1t[icount-1] = serptr->idunion.ser.set->crd;
            z1t[icount-1] = *dptr->idunion.bead.zt;
            /* printf("atom 1 found:\n"); */
            /* printf("atom1: %.4s %d\n",atnam1,resnum1); */
            ifound++;
           }
           atnamt += ATOMTOKLEN;
           serptr->idunion.ser.set->crd+=3;
           dptr->idunion.bead.zt++;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;
          dptr->idunion.bead.zt-=serptr->idunion.ser.prm->Nsit;

          gettok();
          // sprintf(atnam2,"%s    ",tok);
          strcpy(atnam2,tok);
          gettok();
          tokint(&resnum2,tok);

          atnamt = serptr->idunion.ser.prm->AtomNames;
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           // printf("atomt: %.4s %d\n",atnamt,resnumt);
           if (!strncmp(atnam2,atnamt,ATOMTOKLEN)&&resnum2==resnumt) {
            r2t[icount-1] = serptr->idunion.ser.set->crd;
            z2t[icount-1] = *dptr->idunion.bead.zt;
             /* printf("atom 2 found:\n"); */
             /* printf("atom2: %.4s %d\n",atnam2,resnum2); */
            ifound++;
           }
           atnamt += ATOMTOKLEN;
           serptr->idunion.ser.set->crd+=3;
           dptr->idunion.bead.zt++;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;
          dptr->idunion.bead.zt-=serptr->idunion.ser.prm->Nsit;
         }
        }
/*
        r1t-=icount-1;
        r2t-=icount-1;
        z1t-=icount-1;
        z2t-=icount-1;
*/

        if(icount!=ifound/2) {
         printf("setup_bead_bond: atom not found.\n");
         exit(1);
        }
        dptr->idunion.bead.nbond = icount;

        dptr->idunion.bead.r1 = (_REAL **)
            get2(icount,sizeof(long));
        dptr->idunion.bead.r2 = (_REAL **)
            get2(icount,sizeof(long));
        for(i=0;i<icount;i++) {
         *dptr->idunion.bead.r1++ = *r1t++;
         *dptr->idunion.bead.r2++ = *r2t++;
        }
/*
        r1t-=icount;
        r2t-=icount;
        free(r1t);
        free(r2t);
*/
        free(dptr->idunion.bead.zt);
        dptr->idunion.bead.r1-=icount;
        dptr->idunion.bead.r2-=icount;

        dptr->idunion.bead.l = (bondstruct *)
          get2(dptr->idunion.bead.nbond,sizeof(bondstruct));
/*
        dptr->idunion.bead.g = (rdiadstruct *)
          get2(dptr->idunion.bead.nbond,sizeof(rdiadstruct));
*/

/* stores friction coefficients */
   for (i=0;i<dptr->idunion.bead.nbond;i++) {
    dptr->idunion.bead.l->zeta1 = *z1t;
    dptr->idunion.bead.l->zeta2 = *z2t;
    dptr->idunion.bead.l++; z1t++; z2t++;
   }
   dptr->idunion.bead.l -= dptr->idunion.bead.nbond;
/*
   z1t -= dptr->idunion.bead.nbond;
   z2t -= dptr->idunion.bead.nbond;
   free(z1t);
   free(z2t);
*/

/*
   dptr->idunion.bead.F1 =
     rdiad_matrix(1,dptr->idunion.bead.nbond,1,dptr->idunion.bead.nbond);
   dptr->idunion.bead.M1 =
     rdiad_matrix(1,dptr->idunion.bead.nbond,1,dptr->idunion.bead.nbond); 
*/

}
