#include "head.h"

// extern char		tok[];
// char		tok[];
// void *get2(int,int);

/**** extern ****/
extern idstruct		*startid(),*getstruct();

/************************************************************************/

void
init_mca_dihe()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(MCA_DIHE, id);

	dptr->idunion.mca_dihe.cases = 0;
	dptr->idunion.mca_dihe.casesw = 0.;
	dptr->idunion.mca_dihe.dump = 1;
	dptr->idunion.mca_dihe.order = 1;
	dptr->idunion.mca_dihe.rigid = 0;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("MCA_DIHE: expected file name","");
	strcpy(dptr->idunion.mca_dihe.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.mca_dihe.fname, 0, &exists, "MCA_DIHE"))
		exit(1);
	if ((dptr->idunion.mca_dihe.file = 
				  fopen(dptr->idunion.mca_dihe.fname, "w")) == NULL){
		perror(dptr->idunion.mca_dihe.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************/

void
setup_bead_dihe()
{
	idstruct	*dptr,*serptr,*bptr;
	int i,j,k,icount,icount_max,ifound;
        int resnum1,resnum2,resnum3,resnum4,resnumt;
        char atnam1[2*ATOMTOKLEN], atnam2[2*ATOMTOKLEN],
             atnam3[2*ATOMTOKLEN], atnam4[2*ATOMTOKLEN],*atnamt;
        _REAL **r1t, **r2t, **r3t, **r4t,
              *z1t, *z2t, *z3t, *z4t;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != MCA_DIHE)
                        inerr("setup_bead_dihe: expected mca_dihe id or keyword\n","");
        } else {
                printf("setup_bead_dihe: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((bptr = getstruct(tok)) == NULL) 
                inerr("setup_bead_dihe: no such id\n","");
        if (bptr->id_type != BEAD)
                inerr("setup_bead_dihe: bead id expected\n","");
        gettok();
        if ((serptr = getstruct(tok)) == NULL) 
                inerr("setup_bead_dihe: no such id\n","");
        if (serptr->id_type != STREAM)
                inerr("setup_bead_dihe: stream id expected\n","");

        dptr->idunion.mca_dihe.bead = bptr;
        bptr->idunion.bead.prm = serptr->idunion.ser.prm;
        icount_max=10000;
        r1t = (_REAL **) get2(icount_max,sizeof(long));
        r2t = (_REAL **) get2(icount_max,sizeof(long));
        r3t = (_REAL **) get2(icount_max,sizeof(long));
        r4t = (_REAL **) get2(icount_max,sizeof(long));
        z1t = (_REAL *) get2(icount_max,sizeof(long));
        z2t = (_REAL *) get2(icount_max,sizeof(long));
        z3t = (_REAL *) get2(icount_max,sizeof(long));
        z4t = (_REAL *) get2(icount_max,sizeof(long));

        icount=0;
        ifound=0;
        for(;;) {
         gettok();
         if (tok[0] == ';') {
          break;
         } else if(!strcmp(tok,"ORDER")) {
           gettok();
           tokint(&dptr->idunion.mca_dihe.order,tok);
         } else if(!strcmp(tok,"RIGID")) {
	   dptr->idunion.mca_dihe.rigid = 1;
         } else {
          icount++;
          if(icount>icount_max) {
           printf("setup_bead_dihe: icount_max exceeded.\n");
           exit(1);
          }

          if(icount>1) {
           r1t++; r2t++; r3t++; r4t++;
           z1t++; z2t++; z3t++; z4t++;
          }
          sprintf(atnam1,"%s    ",tok);
          gettok();
          tokint(&resnum1,tok);

          atnamt = serptr->idunion.ser.prm->AtomNames;
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           /* printf("atomt: %.4s %d\n",atnamt,resnumt); */
           if (!strncmp(atnam1,atnamt,ATOMTOKLEN)&&resnum1==resnumt) {
            *r1t = serptr->idunion.ser.set->crd;
            *z1t = *bptr->idunion.bead.zt;
            /* printf("atom 1 found:\n"); */
            /* printf("atom1: %.4s %d\n",atnam1,resnum1); */
            ifound++;
           }
           atnamt += ATOMTOKLEN;
           serptr->idunion.ser.set->crd+=3;
           bptr->idunion.bead.zt++;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;
          bptr->idunion.bead.zt-=serptr->idunion.ser.prm->Nsit;

          gettok();
          sprintf(atnam2,"%s    ",tok);
          gettok();
          tokint(&resnum2,tok);

          atnamt = serptr->idunion.ser.prm->AtomNames;
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           /* printf("atomt: %.4s %d\n",atnamt,resnumt); */
           if (!strncmp(atnam2,atnamt,ATOMTOKLEN)&&resnum2==resnumt) {
            *r2t = serptr->idunion.ser.set->crd;
            *z2t = *bptr->idunion.bead.zt;
             /* printf("atom 2 found:\n"); */
             /* printf("atom2: %.4s %d\n",atnam2,resnum2); */
            ifound++;
           }
           atnamt += ATOMTOKLEN;
           serptr->idunion.ser.set->crd+=3;
           bptr->idunion.bead.zt++;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;
          bptr->idunion.bead.zt-=serptr->idunion.ser.prm->Nsit;

          gettok();
          sprintf(atnam3,"%s    ",tok);
          gettok();
          tokint(&resnum3,tok);

          atnamt = serptr->idunion.ser.prm->AtomNames;
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           /* printf("atomt: %.4s %d\n",atnamt,resnumt); */
           if (!strncmp(atnam3,atnamt,ATOMTOKLEN)&&resnum3==resnumt) {
            *r3t = serptr->idunion.ser.set->crd;
            *z3t = *bptr->idunion.bead.zt;
             /* printf("atom 2 found:\n"); */
             /* printf("atom2: %.4s %d\n",atnam2,resnum2); */
            ifound++;
           }
           atnamt += ATOMTOKLEN;
           serptr->idunion.ser.set->crd+=3;
           bptr->idunion.bead.zt++;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;
          bptr->idunion.bead.zt-=serptr->idunion.ser.prm->Nsit;

          gettok();
          sprintf(atnam4,"%s    ",tok);
          gettok();
          tokint(&resnum4,tok);

          atnamt = serptr->idunion.ser.prm->AtomNames;
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           /* printf("atomt: %.4s %d\n",atnamt,resnumt); */
           if (!strncmp(atnam4,atnamt,ATOMTOKLEN)&&resnum4==resnumt) {
            *r4t = serptr->idunion.ser.set->crd;
            *z4t = *bptr->idunion.bead.zt;
             /* printf("atom 2 found:\n"); */
             /* printf("atom2: %.4s %d\n",atnam2,resnum2); */
            ifound++;
           }
           atnamt += ATOMTOKLEN;
           serptr->idunion.ser.set->crd+=3;
           bptr->idunion.bead.zt++;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;
          bptr->idunion.bead.zt-=serptr->idunion.ser.prm->Nsit;

         }

        }
        r1t-=icount-1; r2t-=icount-1; r3t-=icount-1; r4t-=icount-1;
        z1t-=icount-1; z2t-=icount-1; z3t-=icount-1; z4t-=icount-1;

        if(icount!=ifound/4) {
         printf("setup_bead_dihe: atom not found.\n");
         exit(1);
        }
        dptr->idunion.mca_dihe.ndihe = icount;
        dptr->idunion.mca_dihe.nbase =
         1+dptr->idunion.mca_dihe.ndihe*dptr->idunion.mca_dihe.order;
        dptr->idunion.mca_dihe.dim =
         (dptr->idunion.mca_dihe.nbase)*(dptr->idunion.mca_dihe.nbase+1)/2;
        printf("setup_bead_dihe: nbase: %d\n",dptr->idunion.mca_dihe.nbase);
        printf("setup_bead_dihe: dim:   %d\n",dptr->idunion.mca_dihe.dim);

        dptr->idunion.mca_dihe.fc = 
		matrix(1,dptr->idunion.mca_dihe.nbase-1,
                       1,dptr->idunion.mca_dihe.nbase-1);

        dptr->idunion.mca_dihe.s = (_REAL *)
          get2(dptr->idunion.mca_dihe.dim,sizeof(_REAL));
        dptr->idunion.mca_dihe.f = (_REAL *)
          get2(dptr->idunion.mca_dihe.dim,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_dihe.dim;i++) {
         *dptr->idunion.mca_dihe.s++ =0.;
         *dptr->idunion.mca_dihe.f++ =0.;
        }
        dptr->idunion.mca_dihe.s -=dptr->idunion.mca_dihe.dim;
        dptr->idunion.mca_dihe.f -=dptr->idunion.mca_dihe.dim;

        dptr->idunion.mca_dihe.lba = (dihestruct *)
          get2(dptr->idunion.mca_dihe.ndihe,sizeof(dihestruct));

/* stores friction coefficients and ids */
   for (i=0;i<dptr->idunion.mca_dihe.ndihe;i++) {
    dptr->idunion.mca_dihe.lba->id[0] = (int)(*r1t);
    dptr->idunion.mca_dihe.lba->id[1] = (int)(*r2t);
    dptr->idunion.mca_dihe.lba->id[2] = (int)(*r3t);
    dptr->idunion.mca_dihe.lba->id[3] = (int)(*r4t);
    dptr->idunion.mca_dihe.lba->r[0] = *r1t++;
    dptr->idunion.mca_dihe.lba->r[1] = *r2t++;
    dptr->idunion.mca_dihe.lba->r[2] = *r3t++;
    dptr->idunion.mca_dihe.lba->r[3] = *r4t++;
    dptr->idunion.mca_dihe.lba->zeta[0] = *z1t++;
    dptr->idunion.mca_dihe.lba->zeta[1] = *z2t++;
    dptr->idunion.mca_dihe.lba->zeta[2] = *z3t++;
    dptr->idunion.mca_dihe.lba->zeta[3] = *z4t++;
    dptr->idunion.mca_dihe.lba++; 
   }
   dptr->idunion.mca_dihe.lba -= dptr->idunion.mca_dihe.ndihe;
   r1t -= dptr->idunion.mca_dihe.ndihe;
   r2t -= dptr->idunion.mca_dihe.ndihe;
   r3t -= dptr->idunion.mca_dihe.ndihe;
   r4t -= dptr->idunion.mca_dihe.ndihe;
   z1t -= dptr->idunion.mca_dihe.ndihe;
   z2t -= dptr->idunion.mca_dihe.ndihe;
   z3t -= dptr->idunion.mca_dihe.ndihe;
   z4t -= dptr->idunion.mca_dihe.ndihe;
   free(r1t); free(r2t); free(r3t); free(r4t);
   free(z1t); free(z2t); free(z3t); free(z4t);

}


/************************************************************************/
void do_dihe_basis(idstruct *objptr,idstruct *serptr)
{
  idstruct *bptr;
  _REAL *r = serptr->idunion.ser.set->crd;
  _REAL *box = serptr->idunion.ser.set->boxc;

  int ndihe;
  dihestruct *lb,*lb1,*lb2;

  int i,j,k,order,overlap,i1,j1,j1_0;
  _REAL
       x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,dx,dy,dz,boxm1[3],
       x21,y21,z21,x32,y32,z32,x43,y43,z43,x41,y41,z41,
       Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,rA,rB,rC,
       dcosdAx,dcosdAy,dcosdAz,dcosdBx,dcosdBy,dcosdBz,
       cosphi,cph1,cph2,cph3,cph4,cph5,cph6,
       fc1x,fc1y,fc1z,fc2x,fc2y,fc2z,fc4x,fc4y,fc4z,
       E[21],dEdc[21],weight,scp,f;
   _REAL *r1, *r2, *r3, *r4, 
         zi,zj,zr,cutsq;

   ndihe=objptr->idunion.mca_dihe.ndihe;
   bptr=objptr->idunion.mca_dihe.bead;
   zr=bptr->idunion.bead.zr;
   cutsq=bptr->idunion.bead.cutsq;
   order=objptr->idunion.mca_dihe.order;
   for(k=0;k<3;k++) {
    boxm1[k]=1./box[k];
   }
   if(serptr->idunion.ser.iweight)
   {
    weight=serptr->idunion.ser.w->idunion.serw.w;
   } else {
    weight=1.;
   }
   ++objptr->idunion.mca_dihe.cases;
   objptr->idunion.mca_dihe.casesw+=weight;

   lb=objptr->idunion.mca_dihe.lba;

   for(i=0;i<ndihe;i++,lb++) 
   {
    r1=lb->r[0];
    r2=lb->r[1];
    r3=lb->r[2];
    r4=lb->r[3];
    x1=r1[0]; y1=r1[1]; z1=r1[2];
    x2=r2[0]; y2=r2[1]; z2=r2[2];
    x3=r3[0]; y3=r3[1]; z3=r3[2];
    x4=r4[0]; y4=r4[1]; z4=r4[2];
    lb->x[0]=x1; lb->y[0]=y1; lb->z[0]=z1;
    lb->x[1]=x2; lb->y[1]=y2; lb->z[1]=z2;
    lb->x[2]=x3; lb->y[2]=y3; lb->z[2]=z3;
    lb->x[3]=x4; lb->y[3]=y4; lb->z[3]=z4;
/*
   printf("idihe(%d)= %d %d %d %d %d\n",ii,ib,jb,kb,lb,idihe);
*/

    dx = x2-x1; dy = y2-y1; dz = z2-z1;
#ifdef PBC
     dx -= box[0]*rint(dx*boxm1[0]);
     dy -= box[1]*rint(dy*boxm1[1]);
     dz -= box[2]*rint(dz*boxm1[2]);
#endif
    x21 = dx; y21 = dy; z21 = dz;
    dx = x3-x2; dy = y3-y2; dz = z3-z2;
#ifdef PBC
     dx -= box[0]*rint(dx*boxm1[0]);
     dy -= box[1]*rint(dy*boxm1[1]);
     dz -= box[2]*rint(dz*boxm1[2]);
#endif
    x32 = dx; y32 = dy; z32 = dz;
    dx = x4-x3; dy = y4-y3; dz = z4-z3;
#ifdef PBC
     dx -= box[0]*rint(dx*boxm1[0]);
     dy -= box[1]*rint(dy*boxm1[1]);
     dz -= box[2]*rint(dz*boxm1[2]);
#endif
    x43 = dx; y43 = dy; z43 = dz;

      Ax = y21*z32 - z21*y32;
      Ay = z21*x32 - x21*z32;
      Az = x21*y32 - y21*x32;
      rA = l_sqrt(Ax*Ax+Ay*Ay+Az*Az);
      Bx = y32*z43 - z32*y43;
      By = z32*x43 - x32*z43;
      Bz = x32*y43 - y32*x43;
      rB = l_sqrt(Bx*Bx+By*By+Bz*Bz);
      cosphi = (Ax*Bx+Ay*By+Az*Bz)/(rA*rB);
      rB = 1./rB;
      Bx *= rB;
      By *= rB;
      Bz *= rB;
      if (cosphi>1.0) cosphi = 1.0;
      else if (cosphi < -1.0) cosphi = -1.0;
/* printf("cos_phi(%d)= %e\n",i,cosphi); */
      rA = 1/rA;
      Ax *= rA;
      Ay *= rA;
      Az *= rA;
      dcosdAx = rA*(Bx - cosphi*Ax);
      dcosdAy = rA*(By - cosphi*Ay);
      dcosdAz = rA*(Bz - cosphi*Az);
      dcosdBx = rB*(Ax - cosphi*Bx);
      dcosdBy = rB*(Ay - cosphi*By);
      dcosdBz = rB*(Az - cosphi*Bz);

      E[0]=1.; dEdc[0]=0.;
      for(k=1;k<=order;k++)
      {
       E[k]=E[k-1]*cosphi;
       dEdc[k]=(_REAL)(k)*E[k-1];
       // printf("dE(%d,%d)= %e\n",i+1,k,dEdc[k]);
      }

      fc1x = (              z32*dcosdAy -y32*dcosdAz);
      fc1y = (-z32*dcosdAx              +x32*dcosdAz);
      fc1z = ( y32*dcosdAx -x32*dcosdAy             );
      fc2x = (             -z21*dcosdAy +y21*dcosdAz              
                           +z43*dcosdBy -y43*dcosdBz);
      fc2y = ( z21*dcosdAx              -x21*dcosdAz
              -z43*dcosdBx              +x43*dcosdBz);
      fc2z = (-y21*dcosdAx +x21*dcosdAy 
              +y43*dcosdBx -x43*dcosdBy              );
      fc4x = (              z32*dcosdBy -y32*dcosdBz);
      fc4y = (-z32*dcosdBx              +x32*dcosdBz);
      fc4z = ( y32*dcosdBx -x32*dcosdBy);
      lb->fx[0] = fc1x;
      lb->fy[0] = fc1y;
      lb->fz[0] = fc1z;
      lb->fx[1] = (fc2x-fc1x);
      lb->fy[1] = (fc2y-fc1y);
      lb->fz[1] = (fc2z-fc1z);
      lb->fx[2] = (-fc4x-fc2x);
      lb->fy[2] = (-fc4y-fc2y);
      lb->fz[2] = (-fc4z-fc2z);
      lb->fx[3] = -fc4x;
      lb->fy[3] = -fc4y;
      lb->fz[3] = -fc4z;
/*
   printf("f1(%d) = %e %e %e\n",i,lb->fx[0],lb->fy[0],lb->fz[0]);
   printf("f2(%d) = %e %e %e\n",i,lb->fx[1],lb->fy[1],lb->fz[1]);
   printf("f3(%d) = %e %e %e\n",i,lb->fx[2],lb->fy[2],lb->fz[2]);
   printf("f4(%d) = %e %e %e\n",i,lb->fx[3],lb->fy[3],lb->fz[3]);
*/
      for(k=0;k<=order;k++)
      {
       lb->base[k] = E[k];
       lb->dbase_dc[k] = dEdc[k];
      }
   } /* enddo i and lb */
   lb-=ndihe;

   lb=objptr->idunion.mca_dihe.lba;
   for(i=1,lb1=lb;i<=ndihe;i++,lb1++)
   {
    for(j=i,lb2=lb1;j<=ndihe;j++,lb2++)
    {
#include "fdihe.c"
// printf("f= %d %d %e\n",i,j,f);
     objptr->idunion.mca_dihe.fc[i][j] = f;
     objptr->idunion.mca_dihe.fc[j][i] = f;
    }
   }

   k=1;
   *objptr->idunion.mca_dihe.s++ += weight;
   *objptr->idunion.mca_dihe.f++ += 0.;
   k++;

   lb1=objptr->idunion.mca_dihe.lba;
   for(i=1;i<=ndihe;i++,lb1++)
   {
    for(i1=1;i1<=order;i1++,k++)
    {
     scp=lb1->base[i1];
     *objptr->idunion.mca_dihe.s++ += weight*scp;
     *objptr->idunion.mca_dihe.f++ += 0.;
    }
   }

   lb1=objptr->idunion.mca_dihe.lba;
   for(i=1;i<=ndihe;i++,lb1++)
   {
    for(i1=1;i1<=order;i1++)
    {
     for(j=i,lb2=lb1;j<=ndihe;j++,lb2++)
     {
      if(i==j)
      {
       j1_0=i1;
      }
      else
      {
       j1_0=1;
      }
      for(j1=j1_0;j1<=order;j1++,k++)
      {
       scp=lb1->base[i1]*lb2->base[j1];
/*
 printf("i,i1,j,j1= %d %d %d %d\n",i,i1,j,j1);
 printf("fc= %e\n",objptr->idunion.mca_dihe.fc[i][j]);
 printf("dbase_dc= %e %e\n",lb1->dbase_dc[i1],lb2->dbase_dc[j1]);
*/
       f = 
        lb1->dbase_dc[i1]*lb2->dbase_dc[j1]*objptr->idunion.mca_dihe.fc[i][j];
       *objptr->idunion.mca_dihe.s++ += weight*scp;
       *objptr->idunion.mca_dihe.f++ += weight*f;
// printf("f(i,i1,j,j1)= %d %d %d %d %e\n",i,i1,j,j1,f);
      }
     }
    }
   }
   printf("last k= %d, dim= %d\n",k-1,objptr->idunion.mca_dihe.dim);
   objptr->idunion.mca_dihe.s-=objptr->idunion.mca_dihe.dim;
   objptr->idunion.mca_dihe.f-=objptr->idunion.mca_dihe.dim;

}

/************************************************************************/
void
write_mca_dihe(idstruct *dptr)
{
       int i,j,k,l,m,i2,j2,i3,j3,nf,nbasep1;
       FILE *file = dptr->idunion.mcarig.file,*dump;
       idstruct *bptr=dptr->idunion.mcarig.bead;
       void solve_mateq();

       nbasep1=dptr->idunion.mca_dihe.nbase;
       dptr->idunion.mca_dihe.eval = (_REAL *)
        get2(nbasep1,sizeof(_REAL));
       dptr->idunion.mca_dihe.evec = matrix(
         1,nbasep1,1,nbasep1);

       solve_mateq(dptr->idunion.mca_dihe.s,
                   dptr->idunion.mca_dihe.f,
                   nbasep1,
                   dptr->idunion.mca_dihe.dim,
                   dptr->idunion.mca_dihe.cases,
                   dptr->idunion.mca_dihe.casesw,
                   dptr->idunion.mca_dihe.evec,
                   dptr->idunion.mca_dihe.eval,
                   dptr->idunion.mca_dihe.rigid,
                   bptr->idunion.bead.factns,
                   dptr->idunion.mca_dihe.file,
                   &dptr->idunion.mca_dihe.nf);
        printf("write_mca_dihe: nf= %d\n",dptr->idunion.mca_dihe.nf);

}
