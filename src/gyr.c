#include "head.h"

// extern char		tok[];
// char		tok[];

/**** extern ****/
extern idstruct		*startid(),*getstruct();

/************************************************************************/

void
setupgyr()
{
	idstruct	*dptr,*serptr;
        extern parmstruct *prm;
	char		id[TOKLEN];
	int		exists;
	int i,j,k,icount,ifound;
        int resnum1,resnumt;
        char atnam1[2*ATOMTOKLEN],*atnamt;

	dptr = startid(GYR, id);

	dptr->idunion.work3.type = 0;
	dptr->idunion.work3.cases = 0;

/* get id */
        gettok();
        if ((serptr = getstruct(tok)) == NULL) 
                inerr("setupgyr: no such id\n","");
        if (serptr->id_type != STREAM)
                inerr("setupgyr: stream id expected\n","");

        dptr->idunion.work3.type = 1;
        dptr->idunion.work3.prm = serptr->idunion.ser.prm;

        dptr->idunion.work3.r1 = (_REAL **)
            get(dptr->idunion.work3.prm->Nsit*sizeof(long));
        if(dptr->idunion.work3.r1==NULL) {
         printf("setupgyr: error allocating r1.\n");
         exit(1);
        }
        *dptr->idunion.work3.r1 = (_REAL *)serptr->idunion.ser.set->crd;
        for(i=0;i<dptr->idunion.work3.prm->Nsit;i++) {
         *dptr->idunion.work3.r1++=serptr->idunion.ser.set->crd;
         serptr->idunion.ser.set->crd+=3;
        }
        dptr->idunion.work3.r1-=dptr->idunion.work3.prm->Nsit;
        serptr->idunion.ser.set->crd-=3*dptr->idunion.work3.prm->Nsit;

        icount=0;
        ifound=0;
        for(;;) {
         gettok();
         if (tok[0] == ';') {
          break;
         } else if (!strcmp(tok, "ATOMS")) {
           gettok();
           if (!tokint(&dptr->idunion.work3.nlist))
            inerr("GYR: expected integer # of atoms", "");
           icount=-1;
           *dptr->idunion.work3.r1 = serptr->idunion.ser.set->crd;
         } else {
          icount++;

          if(icount>1) {
           dptr->idunion.work3.r1++;
          }
          sprintf(atnam1,"%s    ",tok);
          gettok();
          tokint(&resnum1,tok);
          printf("atnam1: %s %d\n",atnam1,resnum1);
/* find atomnames in list */
          for(i=0;i<prm->Nsit;i++) {
           atnamt = &prm->AtomNames[ATOMTOKLEN*i];
           resnumt=prm->ResNum[i];
           if (!strncmp(atnam1,atnamt,ATOMTOKLEN)&&resnum1==resnumt) {
            *dptr->idunion.work3.r1 = serptr->idunion.ser.set->crd;
            printf("atom found:\n");
            ifound++;
            printf("atom1: %s %d\n",atnam1,resnum1);
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*prm->Nsit;
         }
        }

        if(icount==0) {
         dptr->idunion.work3.nlist = dptr->idunion.work3.prm->Nsit;
         icount = dptr->idunion.work3.prm->Nsit;
        } else if(icount>0) {
         if(icount!=ifound) {
          printf("gyr: atom not found.\n");
          exit(1);
         }
         dptr->idunion.work3.nlist = icount;
         dptr->idunion.work3.r1-=icount-1;
        }
}

/************************************************************************/

void
gyr(idstruct *objptr,idstruct *serptr)
{
   parmstruct   *prm = serptr->idunion.ser.prm;
   int	nsit=prm->Nsit;

   _REAL   *box = serptr->idunion.ser.set->boxc;

   _REAL **r1 = objptr->idunion.work3.r1;

   int i,j,k;
   _REAL *rt1;
   _REAL x1,y1,z1,xc,yc,zc,dx,dy,dz,dst;

/*
   printf("gyr:nlist: %d\n",objptr->idunion.work3.nlist); fflush(stdout);
*/

   ++objptr->idunion.work3.cases;
   objptr->idunion.work3.gyr=0.;
   xc=0.; yc=0.; zc=0.;

   for (i=0;i<objptr->idunion.work3.nlist;i++,r1++) {
    rt1=*r1;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
    xc += x1; yc += y1; zc += z1;
   }
   xc /= objptr->idunion.work3.nlist;
   yc /= objptr->idunion.work3.nlist;
   zc /= objptr->idunion.work3.nlist;
   r1 -= objptr->idunion.work3.nlist;

   for (i=0;i<objptr->idunion.work3.nlist;i++,r1++) {
    rt1=*r1;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
/*
   printf("gyr:r1: %e %e %e\n",x1,y1,z1);
*/
   fflush(stdout);

    dx=x1-xc; dy=y1-yc; dz=z1-zc;
    dst = dx*dx+dy*dy+dz*dz;
    objptr->idunion.work3.gyr += dst;
   }
   objptr->idunion.work3.gyr /= (objptr->idunion.work3.nlist);
   objptr->idunion.work3.gyr = l_sqrt(objptr->idunion.work3.gyr);
   r1 -= objptr->idunion.work3.nlist;
}
