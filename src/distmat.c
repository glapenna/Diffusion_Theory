#include "head.h"

// extern char		tok[];
// char		tok[];

/**** extern ****/
extern idstruct		*startid(),*getstruct();

/************************************************************************/

void
initdistmat()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(DISTMAT, id);

	dptr->idunion.work3.type = 0;
	dptr->idunion.work3.cases = 0;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("DISTMAT: expected file name","");
	strcpy(dptr->idunion.work3.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.work3.fname, 0, &exists, "DISTMAT"))
		exit(1);
	if ((dptr->idunion.work3.file = 
				  fopen(dptr->idunion.work3.fname, "w")) == NULL){
		perror(dptr->idunion.work3.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************
							setupdistmat()
*************************************************************************/

void
setupdistmat()
{
	idstruct	*dptr,*serptr;
        extern parmstruct *prm;
	int i,j,k,icount,ifound;
        int resnum1,resnumt;
        char atnam1[ATOMTOKLEN],atnamt[ATOMTOKLEN];

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != DISTMAT)
                        inerr("setupdistmat: expected id or keyword\n","");
        } else {
                printf("setupdistmat: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((serptr = getstruct(tok)) == NULL) 
                inerr("setupdistmat: no such id\n","");
        if (serptr->id_type != STREAM)
                inerr("setupdistmat: stream id expected\n","");

        dptr->idunion.work3.type = 1;
        dptr->idunion.work3.prm = serptr->idunion.ser.prm;

        dptr->idunion.work3.r1 = (_REAL **)
            get(dptr->idunion.work3.prm->Nmov*sizeof(int));
        if(dptr->idunion.work3.r1==NULL) {
         printf("setupdistmat: error allocating r1.\n");
         exit(1);
        }
        *dptr->idunion.work3.r1 = (_REAL *)serptr->idunion.ser.set->crd;
        for(i=0;i<dptr->idunion.work3.prm->Nmov;i++) {
         *dptr->idunion.work3.r1++=serptr->idunion.ser.set->crd;
         serptr->idunion.ser.set->crd+=3;
        }
        dptr->idunion.work3.r1-=dptr->idunion.work3.prm->Nmov;
        serptr->idunion.ser.set->crd-=3*dptr->idunion.work3.prm->Nmov;

        icount=0;
        ifound=0;
        for(;;) {
         gettok();
         if (tok[0] == ';') {
          break;
         } else if (!strcmp(tok, "ATOMS")) {
           gettok();
           if (!tokint(&dptr->idunion.work3.nlist))
            inerr("DISTMAT: expected integer # of atoms", "");
           icount=-1;
           *dptr->idunion.work3.r1 = serptr->idunion.ser.set->crd;
         } else {
          icount++;

          if(icount>1) {
           dptr->idunion.work3.r1++;
          }
          for (i=0; i<ATOMTOKLEN; i++) {
                atnam1[i] = (char )NULL;
          }
          strcat(atnam1,tok);
          for (k=strlen(tok)+1;k<ATOMTOKLEN;k++)
            strcat(atnam1," ");
          gettok();
          tokint(&resnum1,tok);
          printf("atnam1: %s %d\n",atnam1,resnum1);
/* find atomnames in list */
          for(i=0;i<prm->Nsit;i++) {
           for (j=0; j<ATOMTOKLEN;j++)
            atnamt[j]= prm->AtomNames[ATOMTOKLEN*i+j];
           resnumt=prm->ResNum[i];
           if (!strncmp(atnam1,atnamt,ATOMTOKLEN)&&resnum1==resnumt) {
            *dptr->idunion.work3.r1 = serptr->idunion.ser.set->crd;
            printf("atom found:\n");
            ifound++;
            printf("atom1: %s %d\n",atnam1,resnum1);
            printf("atomt: %s %d\n",atnamt,resnumt);
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*prm->Nsit;
         }
        }

        if(icount==0) {
         dptr->idunion.work3.nlist = dptr->idunion.work3.prm->Nmov;
         icount = dptr->idunion.work3.prm->Nmov;
        } else if(icount>0) {
         if(icount!=ifound) {
          printf("distmat: atom not found.\n");
          exit(1);
         }
         dptr->idunion.work3.nlist = icount;
         dptr->idunion.work3.r1-=icount-1;
        }
        dptr->idunion.work3.udim=(dptr->idunion.work.nlist)*
                                 (dptr->idunion.work.nlist-1)/2;
        dptr->idunion.work3.umat = (_REAL *)
          get(dptr->idunion.work3.udim*sizeof(_REAL));
        for(i=0;i<dptr->idunion.work3.udim;i++) {
         *dptr->idunion.work3.umat++ =0.;
        }
        dptr->idunion.work3.umat-=dptr->idunion.work3.udim;
}


/************************************************************************/

void
distribdistmat(idstruct *objptr,idstruct *serptr)
{
   parmstruct   *prm = serptr->idunion.ser.prm;
   int	nsit=prm->Nsit;

   _REAL   *box = serptr->idunion.ser.set->boxc;

   _REAL **r1 = objptr->idunion.work3.r1;
   _REAL **r2;

   int i,j,k;
   _REAL *rt1,*rt2;
   _REAL x1,y1,z1,x2,y2,z2,dx,dy,dz,dst;

/*
   printf("distribdistmat:nlist: %d\n",objptr->idunion.work3.nlist);
*/
   ++objptr->idunion.work3.cases;
   objptr->idunion.work3.gyr=0.;

   for (i=0;i<objptr->idunion.work3.nlist-1;i++,r1++) {
    rt1=*r1;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];

    for (j=i+1,r2=r1+1;j<objptr->idunion.work3.nlist;j++,r2++) {
/*
   printf("distribdistmat:i,j: %d %d\n",i,j);
*/
     rt2=*r2;
     x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
/*
   printf("distribdistmat:r2: %e %e %e\n",x2,y2,z2);
   printf("distribdistmat:r1: %e %e %e\n",x1,y1,z1);
   fflush(stdout);
*/
     dx=x1-x2;
     dy=y1-y2;
     dz=z1-z2;
     dst = dx*dx+dy*dy+dz*dz;
     objptr->idunion.work3.gyr += dst;
     dst = l_sqrt(dst);
     dst=1./dst;
     *objptr->idunion.work3.umat++ += dst;
    }
   }
   objptr->idunion.work3.umat-=(objptr->idunion.work3.nlist)*
                                (objptr->idunion.work3.nlist-1)/2;
   objptr->idunion.work3.gyr /= (objptr->idunion.work3.nlist)*
                                (objptr->idunion.work3.nlist);
   objptr->idunion.work3.gyr = l_sqrt(objptr->idunion.work3.gyr);
}

/************************************************************************/

void
writedistmat(idstruct *dptr)
{
	int	i,j;
	int	states = dptr->idunion.work3.cases;
	FILE 	*file = dptr->idunion.work3.file;

	_REAL	rm1,fac;

	fac = 1./ ((_REAL)states);
        printf("number of states: %d\n",states);

        for (i=0; i<dptr->idunion.work3.nlist-1; i++) {
         for (j=i+1; j<dptr->idunion.work3.nlist; j++) {
           rm1 = (*dptr->idunion.work3.umat++)*fac;
           fprintf(file, "%3d%4d%18.10e\n", i+1,j+1,rm1);
         }
        }
        fclose(file);
}
