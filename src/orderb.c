#include "head.h"

// extern char		tok[];
// char		tok[];

/**** extern ****/
// extern idstruct		*startid(),*getstruct();
// extern _REAL		*getsetatom();

// FILE *ftmp;

/************************************************************************/

void
initorderb()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(ORDERB, id);

	dptr->idunion.work2.type = 0;
	dptr->idunion.work2.nboxesa = 0;
	dptr->idunion.work2.nboxesb = 0;
	dptr->idunion.work2.cases = 0;
	dptr->idunion.work2.scalea = 0.;
	dptr->idunion.work2.scaleb = 0.;
	dptr->idunion.work2.mina = 0.;
	dptr->idunion.work2.minb = 1.;
	dptr->idunion.work2.maxa = TWOPI;
	dptr->idunion.work2.maxb = -1.;
        dptr->idunion.work2.resunion.scalar.value = 0.;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("ORDERB: expected file name","");
	strcpy(dptr->idunion.work2.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.work2.fname, 0, &exists, "ORDERB"))
		exit(1);
	if ((dptr->idunion.work2.file = 
				  fopen(dptr->idunion.work2.fname, "w")) == NULL){
		perror(dptr->idunion.work2.fname);
		exit(1);
	}
	gettok();
//        ftmp=fopen("wign.time","w");

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************
							setuporderb()
*************************************************************************/

void
setuporderb()
{
	idstruct	*dptr,*serptr;
        extern parmstruct *prm;
	int i,j,k;
        char atnam1[2*ATOMTOKLEN], atnam2[2*ATOMTOKLEN], atnam3[2*ATOMTOKLEN],
             atnamt[2*ATOMTOKLEN];
        int atom1,atom2,atom3;

        for (i=0; i<2*ATOMTOKLEN; i++) {
                atnam1[i] = (char )NULL;
                atnam2[i] = (char )NULL;
                atnam3[i] = (char )NULL;
        }

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != ORDERB)
                        inerr("setuporderb: expected id or keyword\n","");
        } else {
                printf("setuporder: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((serptr = getstruct(tok)) == NULL) 
                inerr("setuporderb: no such id\n","");
        if (serptr->id_type != STREAM)
                inerr("setuporderb: stream id expected\n","");

        dptr->idunion.work2.type = 0;

               for(;;) {
	        gettok();
        	if (tok[0] == ';') {
                 break;
                } else if (!strcmp(tok, "NINTA")) {
/*
                if (!strcmp(tok, "NINTA")) {
*/
                        gettok();
                        if (!tokint(&dptr->idunion.work2.nboxesa))
                         inerr("ORDERB: expected integer # of boxes", "");
                } else if (!strcmp(tok, "MINA")) {
                        gettok();
                        if (!tokREAL(&dptr->idunion.work2.mina))
                         inerr("ORDERB: expected mina", "");
                } else if (!strcmp(tok, "MAXA")) {
                        gettok();
                        if (!tokREAL(&dptr->idunion.work2.maxa))
                         inerr("ORDERB: expected maxa", "");
                } else if (!strcmp(tok, "NINTB")) {
                        gettok();
                        if (!tokint(&dptr->idunion.work2.nboxesb))
                         inerr("ORDERB: expected integer # of boxes", "");
                } else if (!strcmp(tok, "MINB")) {
                        gettok();
                        if (!tokREAL(&dptr->idunion.work2.minb))
                         inerr("ORDERB: expected minb", "");
                } else if (!strcmp(tok, "MAXB")) {
                        gettok();
                        if (!tokREAL(&dptr->idunion.work2.maxb))
                         inerr("ORDERB: expected maxb", "");

		} else if (!strcmp(tok, "ATOMNAME")) {
                  dptr->idunion.work2.type = 1;
                        gettok();
                        sprintf(atnam1,"%s    ",tok); 
		} else if (!strcmp(tok, "ATOMNAMES")) {
                  dptr->idunion.work2.type = 3;
                        gettok();
                        dptr->idunion.work2.r1 = getsetatom();
                        dptr->idunion.work2.r2 = getsetatom();
                        dptr->idunion.work2.r3 = getsetatom();
                        break;
/*
                        sprintf(atnam1,"%s    ",tok); 
                        gettok();
                        sprintf(atnam2,"%s    ",tok); 
                        gettok();
                        sprintf(atnam3,"%s    ",tok); 
*/
		} else if (!strcmp(tok, "ATOMNUM")) {
                  dptr->idunion.work2.type = 2;
                        gettok(); tokint(&atom1);
		} else if (!strcmp(tok, "ATOMNUMS")) {
                  dptr->idunion.work2.type = 4;
                        gettok(); tokint(&atom1);
                        gettok(); tokint(&atom2);
                        gettok(); tokint(&atom3);
		}
/*
                else inerr("setuporderb: option not found", "");
*/
	       }

        dptr->idunion.work2.jump=3*prm->Natom;
        dptr->idunion.work2.prm = prm;
/* printf("natom: %d\n",prm->Natom);*/

/* allocate and initialize */
         dptr->idunion.work2.mina *= TORAD;
         dptr->idunion.work2.maxa *= TORAD;
         dptr->idunion.work2.scalea =
          ( dptr->idunion.work2.maxa - dptr->idunion.work2.mina)
                                / dptr->idunion.work2.nboxesa;
/*
         dptr->idunion.work2.minb *= TORAD;
         dptr->idunion.work2.maxb *= TORAD;
*/

         dptr->idunion.work2.scaleb =
          ( dptr->idunion.work2.maxb - dptr->idunion.work2.minb)
                                / dptr->idunion.work2.nboxesb;
         dptr->idunion.work2.nboxes = dptr->idunion.work2.nboxesa*
                                      dptr->idunion.work2.nboxesb;

		dptr->idunion.work2.boxes = 
                 (_REAL *) get(dptr->idunion.work2.nboxes*sizeof(_REAL));
		for (i=0; i<dptr->idunion.work2.nboxes; i++) {
		 dptr->idunion.work2.boxes[i] = 0.;
                }

       if(dptr->idunion.work2.type==1) {
        for (i=0; i< prm->Natom; i++) {
         serptr->idunion.ser.set->crd+=3;
         for (j=0; j<ATOMTOKLEN;j++) atnamt[j]= prm->AtomNames[ATOMTOKLEN*i+j];
         if (!strncmp(atnam1,atnamt,ATOMTOKLEN)) {
          serptr->idunion.ser.set->crd-=3;
          dptr->idunion.work2.r1 = serptr->idunion.ser.set->crd;
          serptr->idunion.ser.set->crd+=3;
         }
        }
        serptr->idunion.ser.set->crd-=3*prm->Natom;
       } else if(dptr->idunion.work2.type==3) {
        for (i=0; i< prm->Natom; i++) {
          serptr->idunion.ser.set->crd+=3;
         for (j=0; j<ATOMTOKLEN;j++) atnamt[j]= prm->AtomNames[ATOMTOKLEN*i+j];
         if (!strncmp(atnam1,atnamt,ATOMTOKLEN)) {
          serptr->idunion.ser.set->crd-=3;
          dptr->idunion.work2.r1 = serptr->idunion.ser.set->crd;
          serptr->idunion.ser.set->crd+=3;
         }
         if (!strncmp(atnam2,atnamt,ATOMTOKLEN)) {
          serptr->idunion.ser.set->crd-=3;
          dptr->idunion.work2.r2 = serptr->idunion.ser.set->crd;
          serptr->idunion.ser.set->crd+=3;
         }
         if (!strncmp(atnam3,atnamt,ATOMTOKLEN)) {
          serptr->idunion.ser.set->crd-=3;
          dptr->idunion.work2.r3 = serptr->idunion.ser.set->crd;
          serptr->idunion.ser.set->crd+=3;
         }
        }
        serptr->idunion.ser.set->crd-=3*prm->Natom;
/*
       } else if(dptr->idunion.work2.type==2) {
        for (i=0; i< prm->Natom; i++) {
         serptr->idunion.ser.set->crd+=3;
         if (i+1==atom1) {
          serptr->idunion.ser.set->crd-=3;
          dptr->idunion.work2.r1 = serptr->idunion.ser.set->crd;
          serptr->idunion.ser.set->crd+=3;
         }
        }
        serptr->idunion.ser.set->crd-=3*prm->Natom;
*/
       } else if(dptr->idunion.work2.type==4) {
        for (i=0; i< prm->Natom; i++) {
         serptr->idunion.ser.set->crd+=3;
         if (i+1==atom1) {
          serptr->idunion.ser.set->crd-=3;
          dptr->idunion.work2.r1 = serptr->idunion.ser.set->crd;
          serptr->idunion.ser.set->crd+=3;
         }
         if (i+1==atom2) {
          serptr->idunion.ser.set->crd-=3;
          dptr->idunion.work2.r2 = serptr->idunion.ser.set->crd;
          serptr->idunion.ser.set->crd+=3;
         }
         if (i+1==atom3) {
          serptr->idunion.ser.set->crd-=3;
          dptr->idunion.work2.r3 = serptr->idunion.ser.set->crd;
          serptr->idunion.ser.set->crd+=3;
         }
        }
        serptr->idunion.ser.set->crd-=3*prm->Natom;
       }

       dptr->idunion.work2.jump = 3*(prm->Natom);
}


/************************************************************************/

void
distriborderb(idstruct *objptr,idstruct *serptr)
{
   parmstruct   *prm = serptr->idunion.ser.prm;
   int	nmol=prm->Nmol;
   _REAL phase();
   void genquatch2(),genquatch1();

   _REAL   *r = serptr->idunion.ser.set->crd;
   _REAL   *box = serptr->idunion.ser.set->boxc;
   _REAL   *dir = serptr->idunion.ser.set->dir;

   _REAL *r1 = objptr->idunion.work2.r1;
   _REAL *r2 = objptr->idunion.work2.r2;
   _REAL *r3 = objptr->idunion.work2.r3;
   int	jump = objptr->idunion.work2.jump;
   int     nboxesa = objptr->idunion.work2.nboxesa;
   int     nboxesb = objptr->idunion.work2.nboxesb;
   _REAL   scalea = objptr->idunion.work2.scalea, 
           mina = objptr->idunion.work2.mina,
           scaleb = objptr->idunion.work2.scaleb, 
           minb = objptr->idunion.work2.minb;
   _REAL   *boxes = objptr->idunion.work2.boxes;

   int im,k,ial,ibe,ibox,iat1,iat2;
   _REAL *rt1,*rt2,*rt3,dir1[3],dst,cosa,cosb,cosc,cosm,cosmt,cosv[64],cosvar;
   _REAL boxm1[3],q1[4],q2[4],qt[4],qtemp[4],u0[4],u1[4],beta,alpha,Scd;

   for(k=0;k<3;k++) boxm1[k]=1./box[k];
   for(k=0;k<4;k++) u0[k]=0.; u0[3]=1.;
   rt1=r1; rt2=r2; rt3=r3;

/*
*/
   initstat(&objptr->idunion.work.resunion.scalar.stat);
  
   for (im=0;im<nmol;im++) {

    ++objptr->idunion.work2.cases;

/*
    dst = l_sqrt((rt3[2]-rt1[2])*(rt3[2]-rt1[2])+(rt3[1]-rt1[1])*
                 (rt3[1]-rt1[1])+(rt3[0]-rt1[0])*(rt3[0]-rt1[0]));
    cosa = (rt3[0]-rt1[0])/dst; 
    cosb = (rt3[1]-rt1[1])/dst; 
    cosc = (rt3[2]-rt1[2])/dst; 
    cosm = cosa*dir[0]+cosb*dir[1]+cosc*dir[2];
    Scd = 0.5*(3*cosm*cosm-1);
*/

/*
    iat1=627;
    iat2=1881;
    iat1--;
    iat2--;
    dir[0]=r1[3*iat2]-r1[3*iat1];
    dir[1]=r1[3*iat2+1]-r1[3*iat1+1];
    dir[2]=r1[3*iat2+2]-r1[3*iat1+2];
    dst = 1./l_sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
    dir[0]*=dst;
    dir[1]*=dst;
    dir[2]*=dst;
*/
    genquatch2(rt1,rt2,rt3,dir,q1,q2);
    qt[0]=q1[0]; for(k=1;k<4;k++) qt[k]=-q1[k];
    proquat(qt,u0,qtemp);
    proquat(qtemp,q1,u1);
/*
    fprintf(ftmp,"C %e %e %e\n",rt1[0],rt1[1],rt1[2]);
    fprintf(ftmp,"N %e %e %e\n",rt2[0],rt2[1],rt2[2]);
    fprintf(ftmp,"N %e %e %e\n",rt3[0],rt3[1],rt3[2]);
    fprintf(ftmp,"H %e %e %e\n",rt1[0]+u1[1],rt1[1]+u1[2],rt1[2]+u1[3]);
    fprintf(ftmp,"%e %e %e\n",u1[1],u1[2],u1[3]);
*/
/*
    genquatch1(rt1,rt2,rt3,q1);
    u1[1]=q1[0];
    u1[2]=q1[1];
    u1[3]=q1[2];
*/

/*
    { _REAL x,y,z,a=0.7071067814,c=0.5,b=0.6123724357,
            D_1_00,D_1_10r,D_1_10i,D_2_00,D_2_10r,D_2_10i,D_2_20r,D_2_20i;
     x=u1[1];
     y=u1[2];
     z=u1[3];
     D_1_00  = z;
     D_1_10r = a*x;
     D_1_10i = a*y;
     D_2_00  = c*(3.*z*z-1.);
     D_2_10r = b*(x*x-y*y);
     D_2_10i = b*(x*y+y*x);
     D_2_20r = b*(x*z+z*x);
     D_2_20i = b*(y*z+z*y);
     fprintf(ftmp,"%e %e %e %e %e %e %e %e\n",
      D_1_00,D_1_10r,D_1_10i,D_2_00,D_2_10r,D_2_10i,D_2_20r,D_2_20i);
     fwrite(&D_1_00 ,sizeof(_REAL),1,ftmp);
     fwrite(&D_1_10r,sizeof(_REAL),1,ftmp);
     fwrite(&D_1_10i,sizeof(_REAL),1,ftmp);
     fwrite(&D_2_00 ,sizeof(_REAL),1,ftmp);
     fwrite(&D_2_10r,sizeof(_REAL),1,ftmp);
     fwrite(&D_2_10i,sizeof(_REAL),1,ftmp);
     fwrite(&D_2_20r,sizeof(_REAL),1,ftmp);
     fwrite(&D_2_20i,sizeof(_REAL),1,ftmp);
    }
*/

    Scd = q1[0]*q1[0]-q1[1]*q1[1]-q1[2]*q1[2]+q1[3]*q1[3];
    Scd = 0.5*(3.*Scd*Scd -1.);
    updatestat(&objptr->idunion.work.resunion.scalar.stat, Scd);
    Scd = q2[0]*q2[0]-q2[1]*q2[1]-q2[2]*q2[2]+q2[3]*q2[3];
    Scd = 0.5*(3.*Scd*Scd -1.);
    updatestat(&objptr->idunion.work.resunion.scalar.stat, Scd);
    alpha = phase(u1[1],u1[2]);
/*
    beta = acos(u1[3]);
*/
    beta = u1[3];
//printf("distriborderb: beta= %f\n",beta); fflush(stdout);

    ial = (int)((alpha-mina)/scalea);
    ibe = (int)((beta-minb)/scaleb);
    ibox = nboxesa*ibe+ial;
    boxes[ibox]+=1.;
    rt1+= jump; rt2+=jump; rt3+=jump;

   } /* enddo im */
/*
   cosmt = cosmt/(float)nmol;
   cosvar = 0.0;
   for (im=0;im<nmol;im++) {
    cosvar = cosvar + (cosmt - cosv[im])*(cosmt - cosv[im]);
   }
   cosvar = (l_sqrt(cosvar)/(float)(nmol-1));
   printf("cosm%f %f \n",cosmt,cosvar);
*/
/*
   objptr->idunion.work.resunion.scalar.stat.avg /= (_REAL)(nmol);
   objptr->idunion.work.resunion.scalar.stat.avg /= 2.;
*/
/*   printf("Scd(%s): %e\n",
    objptr->id_name,objptr->idunion.work.resunion.scalar.stat.avg); 
*/
}

/************************************************************************/

void
writeorderb(idstruct *dptr)
{
	int	i,j,ibox;
        int     nboxesa = dptr->idunion.work2.nboxesa,
                nboxesb = dptr->idunion.work2.nboxesb;
	_REAL   *boxes = dptr->idunion.work2.boxes;
	int	states = dptr->idunion.work2.cases;
	FILE 	*file = dptr->idunion.work2.file;
	_REAL	scalea = dptr->idunion.work2.scalea,
                mina=dptr->idunion.work2.mina,
	        scaleb = dptr->idunion.work2.scaleb,
                minb=dptr->idunion.work2.minb;

	_REAL	alpha,beta,distrib,fac;

	fac = 1./ ((_REAL)states);

        for (i=0,ibox=0; i<nboxesb; i++) {
         beta = minb+(_REAL)i*scaleb;
/*
         beta *= TODEG;
*/
         for (j=0; j<nboxesa; j++,ibox++) {
           alpha = TODEG*(mina+(_REAL)j*scalea);
           distrib = boxes[ibox]*fac;
           fprintf(file, "%e %e %e\n", beta,alpha,distrib);
         }
         fprintf(file, "\n");
        }
        fclose(file);
//        fclose(ftmp);
}
