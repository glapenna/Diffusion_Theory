#include "head.h"

// extern char		tok[];
// char		tok[];
// void *get2(int,int);

/**** extern ****/
extern idstruct		*startid(),*getstruct();

rotstruct **rot_matrix();
_REAL S_switch[1000];
void setup_switch(void);

/************************************************************************/

void
init_mcarig()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(MCARIG, id);

	dptr->idunion.mcarig.cases = 0;
	dptr->idunion.mcarig.casesw = 0.;
	dptr->idunion.mcarig.dump = 1;
	dptr->idunion.mcarig.rigid = 0;
	dptr->idunion.mcarig.sel_modes = 3;
	dptr->idunion.mcarig.list = 0;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("MCARIG: expected file name","");
	strcpy(dptr->idunion.mcarig.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.mcarig.fname, 0, &exists, "MCARIG"))
		exit(1);
	if ((dptr->idunion.mcarig.file = 
				  fopen(dptr->idunion.mcarig.fname, "w")) == NULL){
		perror(dptr->idunion.mcarig.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************/

void
setup_mcarig_m1()
{
	idstruct	*dptr,*bptr;
	int i,j,k,offset_r,nt;
        FILE *ftmp;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != MCARIG)
                        inerr("setup_mcarig_m1: expected id or keyword\n","");
        } else {
                printf("setup_mcarig_m1: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((bptr = getstruct(tok)) == NULL) 
                inerr("setup_mcarig_m1: no such id\n","");
        if (bptr->id_type != BEAD)
                inerr("setup_mcarig_m1: bead id expected\n","");

        dptr->idunion.mcarig.type = 0;
        dptr->idunion.mcarig.sel_bonds = 0;

        gettok();
        for(;;) {
         if(tok[0]==';') {
          break;
         } else if(!strcmp(tok,"READFSI")) {
          dptr->idunion.mcarig.dump=-1;
          gettok();
         } else if(!strcmp(tok,"LOCAL")) {
          dptr->idunion.mcarig.type=1;
          gettok();
         } else if(!strcmp(tok,"RIGID")) {
          dptr->idunion.mcarig.rigid=1;
          gettok();
         } else if(!strcmp(tok,"NMODES")) {
          gettok();
          tokint(&dptr->idunion.mcarig.sel_modes,tok);
          dptr->idunion.mcarig.type=2;
          gettok();
         } else if(!strcmp(tok,"NBONDS")) {
          gettok();
          tokint(&dptr->idunion.mcarig.sel_bonds,tok);
          gettok();
         } else if(!strcmp(tok,"LIST")) {
          dptr->idunion.mcarig.list=1;
          gettok();
         } else {
          printf("setup_mcarig_m1: option not found.\n");
          exit(1);
         }
        }

        dptr->idunion.mcarig.bead = bptr;

        dptr->idunion.mcarig.n1=
         bptr->idunion.bead.nbond;
        dptr->idunion.mcarig.dim1=
         bptr->idunion.bead.nbond*(bptr->idunion.bead.nbond+1)/2;

        dptr->idunion.mcarig.s1 = (_REAL *)
          get2(dptr->idunion.mcarig.dim1,sizeof(_REAL));
        dptr->idunion.mcarig.f1 = (_REAL *)
          get2(dptr->idunion.mcarig.dim1,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mcarig.dim1;i++) {
         *dptr->idunion.mcarig.s1++ =0.;
         *dptr->idunion.mcarig.f1++ =0.;
        }
        dptr->idunion.mcarig.s1-=dptr->idunion.mcarig.dim1;
        dptr->idunion.mcarig.f1-=dptr->idunion.mcarig.dim1;


       if(dptr->idunion.mcarig.sel_bonds!=0)
       {

        if(dptr->idunion.mcarig.list==0)
        {

        bptr->idunion.bead.l->offset_l = (int *)
          get2(dptr->idunion.mcarig.n1,sizeof(int));
        bptr->idunion.bead.l->doffset = (int *)
          get2(dptr->idunion.mcarig.n1,sizeof(int));

        nt=0;
        for(i=0;i<bptr->idunion.bead.nbond;i++) {
         bptr->idunion.bead.l->offset_l[i]= dptr->idunion.mcarig.sel_bonds -1;
                            offset_r   = dptr->idunion.mcarig.sel_bonds -1;
         bptr->idunion.bead.l->doffset[i]=  offset_r
                                        +bptr->idunion.bead.l->offset_l[i] +1;

         if(i-bptr->idunion.bead.l->offset_l[i] < 0)
         {
          bptr->idunion.bead.l->offset_l[i]-=dptr->idunion.mcarig.sel_bonds-i-1;
         }
         if(i+offset_r > bptr->idunion.bead.nbond-1)
         {
                            offset_r   -=i+1
                                        +dptr->idunion.mcarig.sel_bonds-1
                                        -bptr->idunion.bead.nbond;
         }
         bptr->idunion.bead.l->doffset[i]=  offset_r
                                        +bptr->idunion.bead.l->offset_l[i] +1;
/*
         printf("offset(%d)= %d %d\n",
          i,bptr->idunion.bead.l->offset_l[i],bptr->idunion.bead.l->doffset[i]);
         fflush(stdout);
*/

         nt+=bptr->idunion.bead.l->doffset[i]*
             (bptr->idunion.bead.l->doffset[i]+1)/2;

         }

        }
        else
        {
          ftmp=fopen("j1j2j3.selected","r");
          nt=0;
          while( !feof(ftmp) ) 
          {
           fscanf(ftmp,"%*d %*d %*d %*d %*d %*le %*le ");
           nt++;
          }
          printf("setup_mcarig_m1: reading %d entries from file j1j2j3.selected\n",nt);
          dptr->idunion.mcarig.nlist=nt;
          dptr->idunion.mcarig.j1=(int *)get2(nt,sizeof(int));
          dptr->idunion.mcarig.j2=(int *)get2(nt,sizeof(int));
          dptr->idunion.mcarig.j3=(int *)get2(nt,sizeof(int));
          rewind(ftmp);
          nt=0;
          while( !feof(ftmp) ) 
          {
           fscanf(ftmp,"%*d %*d %d %d %d %*le %*le ",
              &dptr->idunion.mcarig.j1[nt],
              &dptr->idunion.mcarig.j2[nt],
              &dptr->idunion.mcarig.j3[nt]);
           dptr->idunion.mcarig.j1[nt]--;
           dptr->idunion.mcarig.j2[nt]--;
           dptr->idunion.mcarig.j3[nt]--;
           nt++;
          }
          if(nt!=dptr->idunion.mcarig.nlist)
          {
           printf("setup_mcarig_m1: error reading file j1j2j3.selected.\n");
           exit(1);
          }
          fclose(ftmp);

        }

        dptr->idunion.mcarig.n2=nt;

        dptr->idunion.mcarig.dim12=
         dptr->idunion.mcarig.n1*dptr->idunion.mcarig.n2;
        dptr->idunion.mcarig.dim2=
         dptr->idunion.mcarig.n2*(dptr->idunion.mcarig.n2+1)/2;

        dptr->idunion.mcarig.s12 = (_REAL *)
          get2(dptr->idunion.mcarig.dim12,sizeof(_REAL));
        dptr->idunion.mcarig.f12 = (_REAL *)
          get2(dptr->idunion.mcarig.dim12,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mcarig.dim12;i++) {
         *dptr->idunion.mcarig.s12++ =0.;
         *dptr->idunion.mcarig.f12++ =0.;
        }
        dptr->idunion.mcarig.s12-=dptr->idunion.mcarig.dim12;
        dptr->idunion.mcarig.f12-=dptr->idunion.mcarig.dim12;

        dptr->idunion.mcarig.s2 = (_REAL *)
          get2(dptr->idunion.mcarig.dim2,sizeof(_REAL));
        dptr->idunion.mcarig.f2 = (_REAL *)
          get2(dptr->idunion.mcarig.dim2,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mcarig.dim2;i++) {
         *dptr->idunion.mcarig.s2++ =0.;
         *dptr->idunion.mcarig.f2++ =0.;
        }
        dptr->idunion.mcarig.s2-=dptr->idunion.mcarig.dim2;
        dptr->idunion.mcarig.f2-=dptr->idunion.mcarig.dim2;

        dptr->idunion.mcarig.l123 = (mod123struct *)
          get2(dptr->idunion.mcarig.n2,sizeof(mod123struct));

       }

        dptr->idunion.mcarig.nn=dptr->idunion.mcarig.n1+dptr->idunion.mcarig.n2;

        printf("setup_mcarig: n1:    %d\n",dptr->idunion.mcarig.n1);
        printf("setup_mcarig: n2:    %d\n",dptr->idunion.mcarig.n2);
        printf("setup_mcarig: nn:    %d\n",dptr->idunion.mcarig.nn);
        printf("setup_mcarig: dim1:  %d\n",dptr->idunion.mcarig.dim1);
        printf("setup_mcarig: dim12: %d\n",dptr->idunion.mcarig.dim12);
        printf("setup_mcarig: dim2:  %d\n",dptr->idunion.mcarig.dim2);

        dptr->idunion.mcarig.eval = (_REAL *)
         get2(dptr->idunion.mcarig.nn,sizeof(_REAL));
        dptr->idunion.mcarig.evec = matrix(
              1,dptr->idunion.mcarig.nn,
              1,dptr->idunion.mcarig.nn);
/*
        dptr->idunion.mcarig.Rm1 = matrix(
              1,bptr->idunion.bead.nbead,
              1,bptr->idunion.bead.nbead);
        for(i=1;i<=bptr->idunion.bead.nbead;i++) {
         for(j=1;j<=bptr->idunion.bead.nbead;j++) {
          dptr->idunion.mcarig.Rm1[i][j] =0.;
         }
        }
*/
/*
        dptr->idunion.mcarig.Fb = rot_matrix
         (0,bptr->idunion.bead.nbond-1,0,bptr->idunion.bead.nbond-1);
*/

        /* setup_switch(); */

}

/************************************************************************/

do_first_m1(idstruct *objptr,idstruct *serptr)
{
   _REAL *box=serptr->idunion.ser.set->boxc;
   idstruct *bptr;
   int nbond, *offset_l, *doffset;
   int i,j,kk,l,i1,i2,i3,j1,j2,j3,j1_0,j2_0,j3_0,ib,jb,k,m1,m2,overlap,
       i1_0,i2_0,i3_0,i1_1,i2_1,i3_1,j1_1,j2_1,j3_1,ii,jj,jj0;
   _REAL **r1, **r2;
   _REAL *rt1,*rt2,zi,zj,zr,cutsq;
   _REAL boxm1[3],x1,y1,z1,x2,y2,z2,dx,dy,dz,scp;
   bondstruct *lb1, *lb2, *lb,*p1,*p2,*p3;
   mod123struct *l123,*l123_1,*l123_2;
   _REAL sig0r,f[3][3],weight;
   extern _REAL S_switch[1000];

   sig0r=serptr->idunion.ser.prm->sig0r;
   bptr=objptr->idunion.mcarig.bead;
   zr=bptr->idunion.bead.zr;
   cutsq=bptr->idunion.bead.cutsq;
   nbond=bptr->idunion.bead.nbond;
   offset_l=bptr->idunion.bead.l->offset_l;
   doffset=bptr->idunion.bead.l->doffset;

   for(k=0;k<3;k++) {
    boxm1[k]=1./box[k];
   }

   if(serptr->idunion.ser.iweight)
   {
    weight=serptr->idunion.ser.w->idunion.serw.w;
   } else {
    weight=1.;
   }

   ++objptr->idunion.mcarig.cases;
   objptr->idunion.mcarig.casesw+=weight;

/*
   r1=bptr->idunion.bead.r;
   r2=bptr->idunion.bead.r;
   for (i=1;i<=bptr->idunion.bead.nbead-1;i++,r1++) {
    rt1=*r1;;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
    for (j=i+1,r2=r1+1;j<=bptr->idunion.bead.nbead;j++,r2++) {
     rt2=*r2;
     x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
     dx=x1-x2;
     dy=y1-y2;
     dz=z1-z2;
     scp = dx*dx+dy*dy+dz*dz;
     scp = 1./l_sqrt(scp);
     objptr->idunion.mcarig.Rm1[i][j]+=scp;
    }
   }
   return;
*/

   lb=bptr->idunion.bead.l;
   r1=bptr->idunion.bead.r1;
   r2=bptr->idunion.bead.r2;
   for (i=0;i<bptr->idunion.bead.nbond;i++,r1++,r2++) {
    rt1=*r1; rt2=*r2;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
    x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
/*
   printf("do_first_m1:r1(%d): %e %e %e\n",i+1,sig0r*x1,sig0r*y1,sig0r*z1);
   printf("do_first_m1:r2(%d): %e %e %e\n",i+1,sig0r*x2,sig0r*y2,sig0r*z2);
   fflush(stdout);
*/

    dx=x1-x2;
    dy=y1-y2;
    dz=z1-z2;
    lb->id1 = (int)(*r1);
    lb->id2 = (int)(*r2);
    lb->x1 = x1;
    lb->x2 = x2;
    lb->y1 = y1;
    lb->y2 = y2;
    lb->z1 = z1;
    lb->z2 = z2;
#ifdef PBC
    dx -= box[0]*rint(dx*boxm1[0]);
    dy -= box[1]*rint(dy*boxm1[1]);
    dz -= box[2]*rint(dz*boxm1[2]);
#endif
    lb->x = dx;
    lb->y = dy;
    lb->z = dz;
    lb->sq = dx*dx+dy*dy+dz*dz;
    bptr->idunion.bead.bavt += lb->sq;
    bptr->idunion.bead.bavt2 += lb->sq*lb->sq;
    lb++;
   }
   r1-=bptr->idunion.bead.nbond;
   r2-=bptr->idunion.bead.nbond;

   if(bptr->idunion.bead.bav==0) return;

   if(objptr->idunion.mcarig.sel_bonds!=0)
   {

   l123=objptr->idunion.mcarig.l123;
   lb=bptr->idunion.bead.l;
   for (i1=0,p1=lb;i1<nbond;i1++,p1++) {
    for (i2=0,p2=p1-offset_l[i1];i2<doffset[i1];i2++,p2++) {
     for (i3=i2,p3=p2;i3<doffset[i1];i3++,p3++) {
       scp=p2->x*p3->x+p2->y*p3->y+p2->z*p3->z;
       l123->phi[0]=p1->x*scp;
       l123->phi[1]=p1->y*scp;
       l123->phi[2]=p1->z*scp;
       l123->dphidm1[0][0]=rdelta(0,0)*scp;
       l123->dphidm2[0][0]=p1->x*p3->x;
       l123->dphidm3[0][0]=p1->x*p2->x;
       l123->dphidm1[0][1]=rdelta(0,1)*scp;
       l123->dphidm2[0][1]=p1->x*p3->y;
       l123->dphidm3[0][1]=p1->x*p2->y;
       l123->dphidm1[0][2]=rdelta(0,2)*scp;
       l123->dphidm2[0][2]=p1->x*p3->z;
       l123->dphidm3[0][2]=p1->x*p2->z;
       l123->dphidm1[1][0]=rdelta(1,0)*scp;
       l123->dphidm2[1][0]=p1->y*p3->x;
       l123->dphidm3[1][0]=p1->y*p2->x;
       l123->dphidm1[1][1]=rdelta(1,1)*scp;
       l123->dphidm2[1][1]=p1->y*p3->y;
       l123->dphidm3[1][1]=p1->y*p2->y;
       l123->dphidm1[1][2]=rdelta(1,2)*scp;
       l123->dphidm2[1][2]=p1->y*p3->z;
       l123->dphidm3[1][2]=p1->y*p2->z;
       l123->dphidm1[2][0]=rdelta(2,0)*scp;
       l123->dphidm2[2][0]=p1->z*p3->x;
       l123->dphidm3[2][0]=p1->z*p2->x;
       l123->dphidm1[2][1]=rdelta(2,1)*scp;
       l123->dphidm2[2][1]=p1->z*p3->y;
       l123->dphidm3[2][1]=p1->z*p2->y;
       l123->dphidm1[2][2]=rdelta(2,2)*scp;
       l123->dphidm2[2][2]=p1->z*p3->z;
       l123->dphidm3[2][2]=p1->z*p2->z;
/*
 printf("i= %d,%d,%d\n",i1+1,i2+1,i3+1);
 printf("l123= %e\n",l123->dphidm1[0][0]);
 printf("l123= %e\n",l123->dphidm2[0][0]);
 printf("l123= %e\n",l123->dphidm3[0][0]);
 fflush(stdout);
*/
      l123++;
     }
    }
   }

   }
   return;

   lb=bptr->idunion.bead.l;
/* I-I matrices */
   for (ib=1,lb1=lb;ib<=nbond;ib++,lb1++) {
    for (jb=ib,lb2=lb1;jb<=nbond;jb++,lb2++) {

#include "fbb.c"
     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mcarig.f1++ += weight*scp/3.;
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[0][0]=f[0][0];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[0][1]=f[0][1];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[0][2]=f[0][2];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[1][0]=f[0][1];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[1][1]=f[1][1];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[1][2]=f[1][2];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[2][0]=f[0][2];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[2][1]=f[1][2];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[2][2]=f[2][2];

     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[0][0]=f[0][0];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[0][1]=f[0][1];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[0][2]=f[0][2];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[1][0]=f[0][1];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[1][1]=f[1][1];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[1][2]=f[1][2];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[2][0]=f[0][2];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[2][1]=f[1][2];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[2][2]=f[2][2];

/*
     zr=0.;
#include "fbb0.c"
     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mcarig.f0++ += scp/3.;
     zr=bptr->idunion.bead.zr;
*/
     scp= lb1->x*lb2->x+lb1->y*lb2->y+lb1->z*lb2->z;
     /* scp *= S_switch[jb-ib]; */
     *objptr->idunion.mcarig.s1++ += weight*scp/3.;
    }
   }
   objptr->idunion.mcarig.s1 -= objptr->idunion.mcarig.dim1;
   objptr->idunion.mcarig.f1 -= objptr->idunion.mcarig.dim1;
/*   objptr->idunion.mcarig.f0 -= objptr->idunion.mcarig.dim1; */


   if(objptr->idunion.mcarig.sel_bonds==0) return;


   lb=bptr->idunion.bead.l;
   l123=objptr->idunion.mcarig.l123;
   for (i1=0,lb1=lb,i=0,kk=0;i1<nbond;i1++,lb1++) {

   j1_0=0;
   j1_1=nbond;
   for (j1=j1_0,j=0,kk=0,l123_2=l123;j1<j1_1;j1++) {
    j2_0=j1-offset_l[j1];
    j2_1=j2_0+doffset[j1];
    for (j2=j2_0,jj=0;j2<j2_1;j2++,jj++) {
     j3_0=j2;
     j3_1=j3_0+doffset[j1]-jj;
     for (j3=j3_0;j3<j3_1;j3++,j++,l123_2++) {

/*
 printf("i1,j1,j2,j3= %d %d %d %d\n",i1+1,j1+1,j2+1,j3+1);
 fflush(stdout);
*/

     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<3;m1++) {
       {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          rdelta(m1,k)*
          objptr->idunion.mcarig.Fb[i1][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          rdelta(m1,k)*
          objptr->idunion.mcarig.Fb[i1][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          rdelta(m1,k)*
          objptr->idunion.mcarig.Fb[i1][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mcarig.f12 += weight*scp/((_REAL)3);
     objptr->idunion.mcarig.f12++;

     scp=0.;
     scp+=lb1->x*l123_2->phi[0];
     scp+=lb1->y*l123_2->phi[1];
     scp+=lb1->z*l123_2->phi[2];
     *objptr->idunion.mcarig.s12 += weight*scp/((_REAL)3);
     objptr->idunion.mcarig.s12++;

     }
    }
   }

   }
   objptr->idunion.mcarig.s12 -= objptr->idunion.mcarig.dim12;
   objptr->idunion.mcarig.f12 -= objptr->idunion.mcarig.dim12;

   l123=objptr->idunion.mcarig.l123;
   i1_0=0;
   i1_1=nbond;
   for (i1=i1_0,i=0,kk=0,l123_1=l123;i1<i1_1;i1++) {
    i2_0=i1-offset_l[i1];
    i2_1=i2_0+doffset[i1];
    for (i2=i2_0,ii=0;i2<i2_1;i2++,ii++) {
     i3_0=i2;
     i3_1=i3_0+doffset[i1]-ii;
     for (i3=i3_0;i3<i3_1;i3++,i++,l123_1++) {

     j1_0=i1;
     j2_0=i2;
     j3_0=i3;
     j1_1=i1_1;
     j2_1=i2_1;
     j3_1=i3_1;

/*
  printf("i1_0,i2_0,i3_0,j1_0,j2_0,j3_0= %d %d %d %d %d %d\n", 
   i1_0,i2_0,i3_0,j1_0,j2_0,j3_0);
  printf("i1_1,i2_1,i3_1,j1_1,j2_1,j3_1= %d %d %d %d %d %d\n", 
   i1_1,i2_1,i3_1,j1_1,j2_1,j3_1);
*/

   jj0=ii;
   for (j1=j1_0,j=i,l123_2=l123_1;j1<j1_1;j1++,j2_0=j1-offset_l[j1],
    j2_1=j2_0+doffset[j1],j3_0=j2_0,j3_1=j2_1,jj0=0)
   {
    
    for (j2=j2_0,jj=jj0;j2<j2_1;j2++,jj++,j3_0=j2,j3_1=j3_0+doffset[j1]-jj)
    {
/* printf("jj,j3_1= %d %d\n",jj,j3_1); */

     for (j3=j3_0;j3<j3_1;j3++,j++,l123_2++,kk++) 
     {
/*

  printf("i1,i2,i3,j1,j2,j3= %d %d %d %d %d %d\n", 
   i1+1,i2+1,i3+1,j1+1,j2+1,j3+1);
*/
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<3;m1++) {
       {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l123_1->dphidm1[m1][k] *
          objptr->idunion.mcarig.Fb[i1][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm1[m1][k] *
          objptr->idunion.mcarig.Fb[i1][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm1[m1][k] *
          objptr->idunion.mcarig.Fb[i1][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mcarig.Fb[i2][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mcarig.Fb[i2][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mcarig.Fb[i2][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mcarig.Fb[i3][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mcarig.Fb[i3][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mcarig.Fb[i3][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mcarig.f2 += weight*scp/((_REAL)3);
     objptr->idunion.mcarig.f2++;

     scp=0.;
     scp+=l123_1->phi[0]*l123_2->phi[0];
     scp+=l123_1->phi[1]*l123_2->phi[1];
     scp+=l123_1->phi[2]*l123_2->phi[2];
     *objptr->idunion.mcarig.s2 += weight*scp/((_REAL)3);
     objptr->idunion.mcarig.s2++;

     }
    }
   }
     }
    }
   }
   objptr->idunion.mcarig.s2 -= objptr->idunion.mcarig.dim2;
   objptr->idunion.mcarig.f2 -= objptr->idunion.mcarig.dim2;

}


/************************************************************************/

do_first_m1_list(idstruct *objptr,idstruct *serptr)
{
   _REAL *box=serptr->idunion.ser.set->boxc;
   idstruct *bptr;
   int nbond,nlist;
   int i,j,kk,l,i1,i2,i3,j1,j2,j3,j1_0,j2_0,j3_0,ib,jb,k,m1,m2,overlap,
       i1_0,i2_0,i3_0,i1_1,i2_1,i3_1,j1_1,j2_1,j3_1,ii,jj,jj0;
   _REAL **r1, **r2;
   _REAL *rt1,*rt2,zi,zj,zr,cutsq;
   _REAL boxm1[3],x1,y1,z1,x2,y2,z2,dx,dy,dz,scp;
   bondstruct *lb1, *lb2, *lb,p1,p2,p3;
   mod123struct *l123,*l123_1,*l123_2;
   _REAL sig0r,f[3][3],weight;

   sig0r=serptr->idunion.ser.prm->sig0r;
   bptr=objptr->idunion.mcarig.bead;
   zr=bptr->idunion.bead.zr;
   cutsq=bptr->idunion.bead.cutsq;
   nbond=bptr->idunion.bead.nbond;
   nlist=objptr->idunion.mcarig.nlist;

   for(k=0;k<3;k++) {
    boxm1[k]=1./box[k];
   }

/*
    printf("do_first_m1_list: nbond: %d\n",bptr->idunion.bead.nbond);
    printf("do_first_m1_list: nbond: %d\n",objptr->idunion.mcarig.nlist);
    fflush(stdout);
*/

   if(serptr->idunion.ser.iweight)
   {
    weight=serptr->idunion.ser.w->idunion.serw.w;
   } else {
    weight=1.;
   }

   ++objptr->idunion.mcarig.cases;
   objptr->idunion.mcarig.casesw+=weight;

   lb=bptr->idunion.bead.l;
   r1=bptr->idunion.bead.r1;
   r2=bptr->idunion.bead.r2;
   for (i=0;i<bptr->idunion.bead.nbond;i++,r1++,r2++) {
    rt1=*r1; rt2=*r2;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
    x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
/*
   printf("do_first_m1:r1(%d): %e %e %e\n",i+1,sig0r*x1,sig0r*y1,sig0r*z1);
   printf("do_first_m1:r2(%d): %e %e %e\n",i+1,sig0r*x2,sig0r*y2,sig0r*z2);
   fflush(stdout);
*/

    dx=x1-x2;
    dy=y1-y2;
    dz=z1-z2;
    lb->id1 = (int)(*r1);
    lb->id2 = (int)(*r2);
    lb->x1 = x1;
    lb->x2 = x2;
    lb->y1 = y1;
    lb->y2 = y2;
    lb->z1 = z1;
    lb->z2 = z2;
#ifdef PBC
    dx -= box[0]*rint(dx*boxm1[0]);
    dy -= box[1]*rint(dy*boxm1[1]);
    dz -= box[2]*rint(dz*boxm1[2]);
#endif
    lb->x = dx;
    lb->y = dy;
    lb->z = dz;
    lb->sq = dx*dx+dy*dy+dz*dz;
    bptr->idunion.bead.bavt += lb->sq;
    bptr->idunion.bead.bavt2 += lb->sq*lb->sq;
    lb++;
   }
   r1-=bptr->idunion.bead.nbond;
   r2-=bptr->idunion.bead.nbond;

   if(objptr->idunion.mcarig.sel_bonds!=0)
   {

   l123=objptr->idunion.mcarig.l123;
   lb=bptr->idunion.bead.l;
   for (i=0;i<nlist;i++) {
    i1=objptr->idunion.mcarig.j1[i];
    i2=objptr->idunion.mcarig.j2[i];
    i3=objptr->idunion.mcarig.j3[i];
    p1=bptr->idunion.bead.l[i1];
    p2=bptr->idunion.bead.l[i2];
    p3=bptr->idunion.bead.l[i3];
       scp=p2.x*p3.x+p2.y*p3.y+p2.z*p3.z;
       l123->phi[0]=p1.x*scp;
       l123->phi[1]=p1.y*scp;
       l123->phi[2]=p1.z*scp;
       l123->dphidm1[0][0]=rdelta(0,0)*scp;
       l123->dphidm2[0][0]=p1.x*p3.x;
       l123->dphidm3[0][0]=p1.x*p2.x;
       l123->dphidm1[0][1]=rdelta(0,1)*scp;
       l123->dphidm2[0][1]=p1.x*p3.y;
       l123->dphidm3[0][1]=p1.x*p2.y;
       l123->dphidm1[0][2]=rdelta(0,2)*scp;
       l123->dphidm2[0][2]=p1.x*p3.z;
       l123->dphidm3[0][2]=p1.x*p2.z;
       l123->dphidm1[1][0]=rdelta(1,0)*scp;
       l123->dphidm2[1][0]=p1.y*p3.x;
       l123->dphidm3[1][0]=p1.y*p2.x;
       l123->dphidm1[1][1]=rdelta(1,1)*scp;
       l123->dphidm2[1][1]=p1.y*p3.y;
       l123->dphidm3[1][1]=p1.y*p2.y;
       l123->dphidm1[1][2]=rdelta(1,2)*scp;
       l123->dphidm2[1][2]=p1.y*p3.z;
       l123->dphidm3[1][2]=p1.y*p2.z;
       l123->dphidm1[2][0]=rdelta(2,0)*scp;
       l123->dphidm2[2][0]=p1.z*p3.x;
       l123->dphidm3[2][0]=p1.z*p2.x;
       l123->dphidm1[2][1]=rdelta(2,1)*scp;
       l123->dphidm2[2][1]=p1.z*p3.y;
       l123->dphidm3[2][1]=p1.z*p2.y;
       l123->dphidm1[2][2]=rdelta(2,2)*scp;
       l123->dphidm2[2][2]=p1.z*p3.z;
       l123->dphidm3[2][2]=p1.z*p2.z;
/*
 printf("i= %d,%d,%d\n",i1+1,i2+1,i3+1);
 fflush(stdout);
*/
      l123++;
   }

   }

   lb=bptr->idunion.bead.l;
/* I-I matrices */
   for (ib=1,lb1=lb;ib<=nbond;ib++,lb1++) {
    for (jb=ib,lb2=lb1;jb<=nbond;jb++,lb2++) {

#include "fbb.c"
     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mcarig.f1++ += weight*scp/3.;
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[0][0]=f[0][0];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[0][1]=f[0][1];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[0][2]=f[0][2];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[1][0]=f[0][1];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[1][1]=f[1][1];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[1][2]=f[1][2];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[2][0]=f[0][2];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[2][1]=f[1][2];
     objptr->idunion.mcarig.Fb[ib-1][jb-1].f[2][2]=f[2][2];

     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[0][0]=f[0][0];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[0][1]=f[0][1];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[0][2]=f[0][2];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[1][0]=f[0][1];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[1][1]=f[1][1];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[1][2]=f[1][2];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[2][0]=f[0][2];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[2][1]=f[1][2];
     objptr->idunion.mcarig.Fb[jb-1][ib-1].f[2][2]=f[2][2];

/*
     zr=0.;
#include "fbb0.c"
     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mcarig.f0++ += scp/3.;
     zr=bptr->idunion.bead.zr;
*/
     scp= lb1->x*lb2->x+lb1->y*lb2->y+lb1->z*lb2->z;
     /* scp *= S_switch[jb-ib]; */
     *objptr->idunion.mcarig.s1++ += weight*scp/3.;

    }
   }
   objptr->idunion.mcarig.s1 -= objptr->idunion.mcarig.dim1;
   objptr->idunion.mcarig.f1 -= objptr->idunion.mcarig.dim1;
/*   objptr->idunion.mcarig.f0 -= objptr->idunion.mcarig.dim1; */


   if(objptr->idunion.mcarig.sel_bonds==0) return;


   lb=bptr->idunion.bead.l;
   l123=objptr->idunion.mcarig.l123;
   for (i1=0,lb1=lb,kk=0;i1<nbond;i1++,lb1++) {

   for (j=0,l123_2=l123;j<nlist;j++,l123_2++,kk++) {
    j1=objptr->idunion.mcarig.j1[j];
    j2=objptr->idunion.mcarig.j2[j];
    j3=objptr->idunion.mcarig.j3[j];

/*
 printf("kk,i1,j1,j2,j3= %d %d %d %d %d\n",kk,i1+1,j1+1,j2+1,j3+1);
 fflush(stdout);
*/

     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<3;m1++) {
       {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          rdelta(m1,k)*
          objptr->idunion.mcarig.Fb[i1][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          rdelta(m1,k)*
          objptr->idunion.mcarig.Fb[i1][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          rdelta(m1,k)*
          objptr->idunion.mcarig.Fb[i1][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mcarig.f12 += weight*scp/((_REAL)3);
     objptr->idunion.mcarig.f12++;

     scp=0.;
     scp+=lb1->x*l123_2->phi[0];
     scp+=lb1->y*l123_2->phi[1];
     scp+=lb1->z*l123_2->phi[2];
     *objptr->idunion.mcarig.s12 += weight*scp/((_REAL)3);
     objptr->idunion.mcarig.s12++;

   }

   }
   objptr->idunion.mcarig.s12 -= objptr->idunion.mcarig.dim12;
   objptr->idunion.mcarig.f12 -= objptr->idunion.mcarig.dim12;

   l123=objptr->idunion.mcarig.l123;
   for (i=0,kk=0,l123_1=l123;i<nlist;i++,l123_1++) {
    i1=objptr->idunion.mcarig.j1[i];
    i2=objptr->idunion.mcarig.j2[i];
    i3=objptr->idunion.mcarig.j3[i];

   for (j=i,l123_2=l123_1;j<nlist;j++,l123_2++,kk++)
   {
    
/*
  printf("i1,i2,i3,j1,j2,j3= %d %d %d %d %d %d\n", 
   i1+1,i2+1,i3+1,j1+1,j2+1,j3+1);
*/

     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<3;m1++) {
       {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l123_1->dphidm1[m1][k] *
          objptr->idunion.mcarig.Fb[i1][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm1[m1][k] *
          objptr->idunion.mcarig.Fb[i1][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm1[m1][k] *
          objptr->idunion.mcarig.Fb[i1][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mcarig.Fb[i2][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mcarig.Fb[i2][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mcarig.Fb[i2][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mcarig.Fb[i3][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mcarig.Fb[i3][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mcarig.Fb[i3][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mcarig.f2 += weight*scp/((_REAL)3);
     objptr->idunion.mcarig.f2++;

     scp=0.;
     scp+=l123_1->phi[0]*l123_2->phi[0];
     scp+=l123_1->phi[1]*l123_2->phi[1];
     scp+=l123_1->phi[2]*l123_2->phi[2];
     *objptr->idunion.mcarig.s2 += weight*scp/((_REAL)3);
     objptr->idunion.mcarig.s2++;

    }
   }
   objptr->idunion.mcarig.s2 -= objptr->idunion.mcarig.dim2;
   objptr->idunion.mcarig.f2 -= objptr->idunion.mcarig.dim2;

}


/************************************************************************/

void
write_first(idstruct *dptr)
{
	int	i,j,k,l,m,i2,j2;
	int	states = dptr->idunion.mcarig.cases;
	FILE 	*file = dptr->idunion.mcarig.file,*dump;
        idstruct *bptr=dptr->idunion.mcarig.bead;

	_REAL fact;

        int n1,n2,nn,nf,rank,mm;
        _REAL **F,**S,**Fm1,*w,**tmp,**F0,**tmp0;
        _REAL sum,sum2,tol,RKR;
        _REAL *u,*v;
#ifdef IMSL
        d_complex *vi,*wi;
#else
        _REAL *vi;
#endif

        return;
	/* fact = 1./ ((_REAL)states); */
	fact = 1./ dptr->idunion.mcarig.casesw;

        printf("write_first: number of states: %d\n",states);
/*
        for(i=1,sum=0.;i<=bptr->idunion.bead.nbead-1;i++) {
         for(j=i+1;j<=bptr->idunion.bead.nbead;j++) {
          dptr->idunion.mcarig.Rm1[i][j] *= fact;
          sum+=2.*bptr->idunion.bead.z[i-1]*bptr->idunion.bead.z[j-1]*
                  dptr->idunion.mcarig.Rm1[i][j];
         }
        }
        RKR = ((_REAL)bptr->idunion.bead.nbead)*bptr->idunion.bead.zav/
          (1. + (4./3.)*zr*sum/
                 ((_REAL)bptr->idunion.bead.nbead)/bptr->idunion.bead.zav );
        printf("RKR= %e\n",RKR);
        exit(0);
*/

        if(bptr->idunion.bead.bav==0.) {
         bptr->idunion.bead.bav =
          bptr->idunion.bead.bavt*fact/bptr->idunion.bead.nbond;
         bptr->idunion.bead.bav = bptr->idunion.bead.prm->sig0r*
          l_sqrt(bptr->idunion.bead.bav);
         printf("write_first: bav: %e\n",bptr->idunion.bead.bav);
         bptr->idunion.bead.bav =
          bptr->idunion.bead.bavt2*fact/bptr->idunion.bead.nbond;
         bptr->idunion.bead.bav =
          bptr->idunion.bead.prm->sig0r*
          bptr->idunion.bead.prm->sig0r*
          l_sqrt(bptr->idunion.bead.bav);
         printf("write_first: bav2: %e\n",bptr->idunion.bead.bav);
         exit(0);
        }

        for (i=1; i<=bptr->idunion.bead.nbond; i++) {
         for (j=i; j<=bptr->idunion.bead.nbond; j++) {
           *dptr->idunion.mcarig.s1++ *= fact;
         }
        }
        dptr->idunion.mcarig.s1 -= dptr->idunion.mcarig.dim1;

        for (i=1; i<=bptr->idunion.bead.nbond; i++) {
         for (j=i; j<=bptr->idunion.bead.nbond; j++) {
           *dptr->idunion.mcarig.f1++ *= fact;
         }
        }
        dptr->idunion.mcarig.f1 -= dptr->idunion.mcarig.dim1;

/*
        for (i=1; i<=bptr->idunion.bead.nbond; i++) {
         for (j=i; j<=bptr->idunion.bead.nbond; j++) {
           *dptr->idunion.mcarig.f0++ *= fact;
         }
        }
        dptr->idunion.mcarig.f0 -= dptr->idunion.mcarig.dim1;
*/

        if(dptr->idunion.mcarig.sel_bonds!=0)
        {

        for (i=1; i<=dptr->idunion.mcarig.n1; i++) {
         for (j=1; j<=dptr->idunion.mcarig.n2; j++) {
           *dptr->idunion.mcarig.s12++ *= fact;
         }
        }
        dptr->idunion.mcarig.s12 -= dptr->idunion.mcarig.dim12;

        for (i=1; i<=dptr->idunion.mcarig.n1; i++) {
         for (j=1; j<=dptr->idunion.mcarig.n2; j++) {
           *dptr->idunion.mcarig.f12++ *= fact;
         }
        }
        dptr->idunion.mcarig.f12 -= dptr->idunion.mcarig.dim12;

        for (i=1; i<=dptr->idunion.mcarig.n2; i++) {
         for (j=i; j<=dptr->idunion.mcarig.n2; j++) {
           *dptr->idunion.mcarig.s2++ *= fact;
         }
        }
        dptr->idunion.mcarig.s2 -= dptr->idunion.mcarig.dim2;

        for (i=1; i<=dptr->idunion.mcarig.n2; i++) {
         for (j=i; j<=dptr->idunion.mcarig.n2; j++) {
           *dptr->idunion.mcarig.f2++ *= fact;
         }
        }
        dptr->idunion.mcarig.f2 -= dptr->idunion.mcarig.dim2;

        }

        if(dptr->idunion.mcarig.dump==1) {
         dump=fopen("fsI.dump","w");
         fwrite(dptr->idunion.mcarig.f1,dptr->idunion.mcarig.dim1*sizeof(_REAL),
                1,dump);
         fwrite(dptr->idunion.mcarig.s1,dptr->idunion.mcarig.dim1*sizeof(_REAL),
                1,dump);
         fwrite(dptr->idunion.mcarig.f12,dptr->idunion.mcarig.dim12*sizeof(_REAL),
                1,dump);
         fwrite(dptr->idunion.mcarig.s12,dptr->idunion.mcarig.dim12*sizeof(_REAL),
                1,dump);
         fwrite(dptr->idunion.mcarig.f2,dptr->idunion.mcarig.dim2*sizeof(_REAL),
                1,dump);
         fwrite(dptr->idunion.mcarig.s2,dptr->idunion.mcarig.dim2*sizeof(_REAL),
                1,dump);
         fflush(dump);
         fclose(dump);
        } else if(dptr->idunion.mcarig.dump==-1) {
         dump=fopen("fsI.dump","r");
         fread(dptr->idunion.mcarig.f1,dptr->idunion.mcarig.dim1*sizeof(_REAL),
                1,dump);
         fread(dptr->idunion.mcarig.s1,dptr->idunion.mcarig.dim1*sizeof(_REAL),
                1,dump);
         fread(dptr->idunion.mcarig.f12,dptr->idunion.mcarig.dim12*sizeof(_REAL),
                1,dump);
         fread(dptr->idunion.mcarig.s12,dptr->idunion.mcarig.dim12*sizeof(_REAL),
                1,dump);
         fread(dptr->idunion.mcarig.f2,dptr->idunion.mcarig.dim2*sizeof(_REAL),
                1,dump);
         fread(dptr->idunion.mcarig.s2,dptr->idunion.mcarig.dim2*sizeof(_REAL),
                1,dump);
         fclose(dump);
        }


/* invert F */
        n1 = dptr->idunion.mcarig.n1;
        n2 = dptr->idunion.mcarig.n2;
        nn = dptr->idunion.mcarig.nn;
        /* F0= matrix(1,nn,1,nn); */
        F = matrix(1,nn,1,nn);
        S = matrix(1,nn,1,nn);

        for(i=1;i<=n1;i++) {
         for(j=i;j<=n1;j++) {
          F[i][j]=(*dptr->idunion.mcarig.f1++);
          F[j][i]=F[i][j];
/*
          F0[i][j]=(*dptr->idunion.mcarig.f0++);
          F0[j][i]=F0[i][j];
*/
          S[i][j]=*dptr->idunion.mcarig.s1++;
          S[j][i]=S[i][j];
         }
        }
        dptr->idunion.mcarig.s1 -= dptr->idunion.mcarig.dim1;
        dptr->idunion.mcarig.f1 -= dptr->idunion.mcarig.dim1;
        /* dptr->idunion.mcarig.f0 -= dptr->idunion.mcarig.dim1; */

        if(dptr->idunion.mcarig.sel_bonds!=0) {

        for(i=1;i<=n1;i++) {
         for(j=1;j<=n2;j++) {
          j2=n1+j;
          F[i][j2]=*dptr->idunion.mcarig.f12++;
          F[j2][i]=F[i][j2];
          S[i][j2]=*dptr->idunion.mcarig.s12++;
          S[j2][i]=S[i][j2];
         }
        }
        dptr->idunion.mcarig.s12 -= dptr->idunion.mcarig.dim12;
        dptr->idunion.mcarig.f12 -= dptr->idunion.mcarig.dim12;

        for(i=1;i<=n2;i++) {
         for(j=i;j<=n2;j++) {
          i2=n1+i;
          j2=n1+j;
          F[i2][j2]=*dptr->idunion.mcarig.f2++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*dptr->idunion.mcarig.s2++;
          S[j2][i2]=S[i2][j2];
         }
        }
        dptr->idunion.mcarig.s2 -= dptr->idunion.mcarig.dim2;
        dptr->idunion.mcarig.f2 -= dptr->idunion.mcarig.dim2;

        }

/*
        fprintf(file,"S,F,F0:\n");
        for(i=1;i<=nn;i++)
         for(j=1;j<=nn;j++)
          fprintf(file,"%d %d %e %e\n",i,j,S[i][j],F[i][j]);
        fflush(file);
*/
        

        for(i=1;i<=nn;i++) {
         if(S[i][i]<ZERO) {
          printf("write_first:\n");
          printf("S element is wrong: S(%d,%d) = %e\n",i,i,S[i][i]);
         }
        }

        Fm1 = matrix(1,nn,1,nn);
#ifdef IMSL
        tol = TOLSVD;
        u = vector(0,nn*nn-1);
        if(dptr->idunion.mcarig.rigid==0) {
         for(i=1,k=0;i<=nn;i++)
          for(j=1;j<=nn;j++,k++)
           u[k]=F[i][j];
        } else {
         printf("not implemented...\n");
        }
        w = imsl_d_lin_svd_gen(nn,nn,u,
              IMSL_RANK, tol, &rank,
              IMSL_INVERSE,&v,0);
        free(w);
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
          Fm1[i][j]=v[k];
         }
        }
        free(v);
#else
        if(dptr->idunion.mcarig.rigid==0) {
         rank=invert(F,nn,nn,Fm1);
        } else {
         printf("not implemented...\n");
        }
#endif
        fprintf(file,"rank of F: %d\n",rank);

        tmp = matrix(1,nn,1,nn);
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=Fm1[i][k]*S[k][j];
          }
          tmp[i][j]=sum;
         }
        }

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*tmp[k][j];
          }
          Fm1[i][j]=sum;
         }
        }

        for(i=1;i<=nn;i++) {
         if(Fm1[i][i]<=ZERO) {
          printf("write_first:\n");
          printf("SF^-1S diagonal element is wrong: F(%d,%d) = %e\n",
           i,i,Fm1[i][i]);
         }
         if(Nabs(S[i][i])>ZERO) {
          Fm1[i][i] /= S[i][i];
         }
        }

        for(i=1;i<=nn;i++) {
         fprintf(file,"tau(%d)= %e ",i,1.*Fm1[i][i]);
         fprintf(file,"%e\n",bptr->idunion.bead.factns*Fm1[i][i]);
        }
        fflush(file);

#ifdef IMSL
        tol = TOLSVD;
        u = vector(0,nn*nn-1);
        if(dptr->idunion.mcarig.rigid==0) {
         for(i=1,k=0;i<=nn;i++)
          for(j=1;j<=nn;j++,k++)
           u[k]=S[i][j];
        } else {
         for(i=1,k=0;i<=nn;i++)
          for(j=1;j<=nn;j++,k++)
           u[k]=F[i][j];
        }
        w = imsl_d_lin_svd_gen(nn,nn,u,
              IMSL_RANK, tol, &rank,
              IMSL_INVERSE,&v,0);
        free(w);
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
          Fm1[i][j]=v[k];
         }
        }
        free(v);
#else
        if(dptr->idunion.mcarig.rigid==0) {
         rank=invert(S,nn,nn,Fm1);
        } else {
         rank=invert(F,nn,nn,Fm1);
        }
#endif
        fprintf(file,"rank of S: %d\n",rank);

/*
        fprintf(file,"Sm1*S:\n");
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=Fm1[i][k]*S[k][j];
          }
          fprintf(file,"%d %d %e\n",i,j,sum);
          fflush(file);
         }
        }
*/

        /* tmp0 = matrix(1,nn,1,nn); */

        if(dptr->idunion.mcarig.rigid==0) {

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.,sum2=0.;k<=nn;k++) {
           sum+=Fm1[i][k]*F[k][j];
           /* sum2+=Fm1[i][k]*F0[k][j]; */
          }
          tmp[i][j]=sum;
          /* tmp0[i][j]=sum2; */
         }
        }

        } else {

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=Fm1[i][k]*S[k][j];
          }
          tmp[i][j]=sum;
         }
        }

        }

/* diagonalize A = S-1*F0 or A = F0^-1*S if rigidi */
        
        w=vector(1,nn);
/*
#ifdef IMSL
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
          u[k]=tmp0[i][j];
         }
        }
        wi = imsl_d_eig_gen(nn,u,IMSL_VECTORS,&vi,0);
        if(dptr->idunion.mcarig.rigid==0) {
         for(i=0;i<nn;i++) {
          w[i+1] = 1./wi[i].re;
         }
        } else {
         for(i=0;i<nn;i++) {
          w[i+1] = wi[i].re;
         }
        }
        free(wi);
#else
        vi = (_REAL *)get2(nn*nn,sizeof(_REAL));
#ifdef EISPACK
        eigvec_ev(nn,tmp0,vi,w);
        if(dptr->idunion.mcarig.rigid==0) {
         for(i=0;i<nn;i++) {
          w[i+1] = 1./w[i+1];
         }
        }
#else
        eigvec_nr(tmp0,nn,w,vi);
#endif
#endif

        for(i=1,sum=0.,fact=0.,nf=0;i<=nn;i++) {
         sum+=w[i];
         if(Nabs(w[i])>ZERO)
         {
          nf++;
          fprintf(file,"lambda0(%d)= %e\n",i,1./w[i]);
          fact+= w[i];
         }
        }
        fprintf(file,"\n"); fflush(file);

        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
#ifdef IMSL
          dptr->idunion.mcarig.evec0[i][j] = vi[k].re;
#else
          dptr->idunion.mcarig.evec0[i][j] = vi[k];
#endif
         }
        }

        eigsrt(w,dptr->idunion.mcarig.evec0,nn);
        for(i=1;i<=nn;i++) {
         dptr->idunion.mcarig.eval0[i-1] = w[i];
        }
*/
/*
        free_vector(w,1,nn);
        free(vi);
*/
/*
#ifdef IMSL
        tol = TOLSVD;
        for(i=1,k=0;i<=nn;i++)
         for(j=1;j<=nn;j++,k++)
          u[k]=dptr->idunion.mcarig.evec0[i][j];
        w = imsl_d_lin_svd_gen(nn,nn,u,
              IMSL_RANK, tol, &rank,
              IMSL_INVERSE,&v,0);
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
          Fm1[i][j]=v[k];
         }
        }
        free(v);
#else
        rank=invert(dptr->idunion.mcarig.evec0,nn,nn,Fm1);
#endif


        fprintf(file,"G matrix:\n");
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           for(l=1;l<=nn;l++) {
            sum+= Fm1[i][k]*
                  tmp[k][l]*
                  dptr->idunion.mcarig.evec0[l][j];
           }
          }
          fprintf(file,"%d %d %e\n",i,j,sum);
          if(i==j) 
          w[i]=sum;
         }
        }

        for(i=1;i<=nn;i++) {
         fprintf(file,"lambda1(%d)= %e\n",i,w[i]);
        }
        fflush(stdout);
*/

/* diagonalize A = S-1*F or A = F^-1*S if rigidi */
        
        /* w=vector(1,nn); */
#ifdef IMSL
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
          u[k]=tmp[i][j];
         }
        }
        wi = imsl_d_eig_gen(nn,u,IMSL_VECTORS,&vi,0);
        if(dptr->idunion.mcarig.rigid==0) {
         for(i=0;i<nn;i++) {
          w[i+1] = 1./wi[i].re;
         }
        } else {
         for(i=0;i<nn;i++) {
          w[i+1] = wi[i].re;
         }
        }
        free(wi);
        free_vector(u,0,nn*nn-1);
#else
        vi = (_REAL *)get2(nn*nn,sizeof(_REAL));
#ifdef EISPACK
        eigvec_ev(nn,tmp,vi,w);
        if(dptr->idunion.mcarig.rigid==0) {
         for(i=0;i<nn;i++) {
          w[i+1] = 1./w[i+1];
         }
        }
#else
        eigvec_nr(tmp,nn,w,vi);
#endif
#endif
        for(i=1,sum=0.,fact=0.,nf=0;i<=nn;i++) {
         sum+=w[i];
         if(Nabs(w[i])>ZERO)
         {
          nf++;
          fprintf(file,"lambda2(%d)= %e\n",i,1./w[i]);
          fact+= w[i];
         }
        }
        fprintf(file,"\n"); fflush(file);

        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
#ifdef IMSL
          dptr->idunion.mcarig.evec[i][j] = vi[k].re;
#else
          dptr->idunion.mcarig.evec[i][j] = vi[k];
#endif
         }
        }
        eigsrt(w,dptr->idunion.mcarig.evec,nn);
        for(i=1;i<=nn;i++) {
         dptr->idunion.mcarig.eval[i-1] = w[i];
        }
/*
        free_vector(w,1,nn);
        free(vi);
*/

/* copy evec0 to evec */
/*
        { FILE *ftmp;
         ftmp=fopen("evec.out","w");
         for(i=1;i<=bptr->idunion.bead.nbond;i++) {
          fprintf(ftmp,"%d %f %f %f %f\n",i,
           dptr->idunion.mcarig.evec[i][1],
           dptr->idunion.mcarig.evec[i][2],
           dptr->idunion.mcarig.evec[i][3],
           dptr->idunion.mcarig.evec[i][4]);
           dptr->idunion.mcarig.evec0[i][1],
           dptr->idunion.mcarig.evec0[i][2],
           dptr->idunion.mcarig.evec0[i][3]);
         }
        }
*/
/*
        for(i=1;i<=bptr->idunion.bead.nbond;i++) {
         for(j=1;j<=bptr->idunion.bead.nbond;j++) {
          dptr->idunion.mcarig.evec[i][j]=dptr->idunion.mcarig.evec0[i][j];
         }
        }
*/
/*
        fprintf(file,"evec:\n");
        for(i=1;i<=bptr->idunion.bead.nbond;i++) {
         for(j=1;j<=bptr->idunion.bead.nbond;j++) {
          fprintf(file,"%d %d %e\n",i,j,dptr->idunion.mcarig.evec[i][j]);
         }
        }
        fflush(file);
*/

        if(dptr->idunion.mcarig.type==1) {
         nf=nn;
        } if(dptr->idunion.mcarig.type==2) {
         nf=dptr->idunion.mcarig.sel_modes;
        }

        if(dptr->idunion.mcarig.rigid==1) {
         if(dptr->idunion.mcarig.type==1) {
          nf=3;
         } if(dptr->idunion.mcarig.type==2) {
          nf=dptr->idunion.mcarig.sel_modes;
         }
        }

        dptr->idunion.mcarig.nf = nf;
/*
        dptr->idunion.mcarig.dim1=
         dptr->idunion.mcarig.nf*(dptr->idunion.mcarig.nf+1)/2;
*/
/*
        dptr->idunion.mcarig.dim1=
         dptr->idunion.mcarig.n1*(dptr->idunion.mcarig.n1+1)/2;
        free(dptr->idunion.mcarig.s1);
        free(dptr->idunion.mcarig.f1);
        dptr->idunion.mcarig.s1 = (_REAL *)
          get2(dptr->idunion.mcarig.dim1,sizeof(_REAL));
        dptr->idunion.mcarig.f1 = (_REAL *)
          get2(dptr->idunion.mcarig.dim1,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mcarig.dim1;i++) {
         *dptr->idunion.mcarig.s1++ =0.;
         *dptr->idunion.mcarig.f1++ =0.;
        }
        dptr->idunion.mcarig.s1-=dptr->idunion.mcarig.dim1;
        dptr->idunion.mcarig.f1-=dptr->idunion.mcarig.dim1;
*/

        /* nn=bptr->idunion.bead.nbond; */

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mcarig.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mcarig.evec[k][i]*tmp[k][j];
          }
          Fm1[i][j]=sum;
         }
        }
/*
        fprintf(file,"Ct*S(I)*C:\n");
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          fprintf(file,"%d %d %e\n",i,j,Fm1[i][j]);
         }
        }
        fflush(file);
*/

/* normalize eigenvectors to Ct*S*C=I */
        for(j=1,k=0;j<=nn;j++) {
         if(Nabs(Fm1[j][j]) > ZERO) {
          k++;
          sum=Fm1[j][j];
         } else {
          sum=1.;
         }
         sum=l_sqrt(sum);
         for(i=1;i<=nn;i++) {
          dptr->idunion.mcarig.evec[i][j]/=sum;
         }
        }
        if(k!=nf) {
         printf("write_first: missed zero in Ct*S*C.\n");
        }
/*
#ifdef IMSL

#endif

#ifdef EISPACK
        adjust_ev(nn,Fm1,tmp,dptr->idunion.mcarig.evec);
#endif
*/

/* manually swap two eigenvectors */
/*
   for(i=1;i<=nn;i++) {
    sum=dptr->idunion.mcarig.evec[i][5];
    dptr->idunion.mcarig.evec[i][5]=dptr->idunion.mcarig.evec[i][4];
    dptr->idunion.mcarig.evec[i][4]=sum;
   }
*/

/* modify manually eigenvectors */
/*
         for(i=1;i<=nn;i++) {
          dptr->idunion.mcarig.evec[i][4]=0.;
         }
         dptr->idunion.mcarig.evec[10][4]=1.;
         for(i=1;i<=nn;i++) {
          dptr->idunion.mcarig.evec[i][5]=0.;
         }
         dptr->idunion.mcarig.evec[10][5]=1.;
         for(i=1;i<=nn;i++) {
          dptr->idunion.mcarig.evec[i][6]=0.;
         }
         dptr->idunion.mcarig.evec[11][6]=1.;
*/
/*
*/
        { FILE *fevec;
        fevec=fopen("evec.m1","w");
        fprintf(file,"writing evec in file evec.m1\n");
        for(j=1;j<=nn;j++) {
         for(i=1;i<=nn;i++) {
          fprintf(fevec,"evec(%d,%d)= %e\n",
           i,j,dptr->idunion.mcarig.evec[i][j]);
         }
        }
        fclose(fevec);
        }

/* transform F in Ct*F*C and S in Ct*S*C */
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=F[i][k]*dptr->idunion.mcarig.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mcarig.evec[k][i]*tmp[k][j];
          }
          F[i][j]=sum;
         }
        }
/*
        fprintf(file,"Ct*F1*C:\n");
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          fprintf(file,"%d %d %e\n",i,j,F[i][j]);
         }
        }
        fflush(file);
*/

/*
        for(i=1;i<=dptr->idunion.mcarig.n1;i++) {
         for(j=i;j<=dptr->idunion.mcarig.n1;j++) {
          *dptr->idunion.mcarig.f1++ = F[i][j];
         }
        }
        dptr->idunion.mcarig.f1 -= dptr->idunion.mcarig.dim1;
*/

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mcarig.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
/*
        fprintf(file,"S1*C:\n");
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          fprintf(file,"%d %d %e\n",i,j,tmp[i][j]);
         }
        }
        fflush(file);
*/

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mcarig.evec[k][i]*tmp[k][j];
          }
          S[i][j]=sum;
         }
        }

/*
        fprintf(file,"Ct*S1*C:\n");
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          fprintf(file,"%d %d %e\n",i,j,S[i][j]);
         }
        }
        fflush(file);
*/

/*
        for(i=1;i<=dptr->idunion.mcarig.n1;i++) {
         for(j=i;j<=dptr->idunion.mcarig.n1;j++) {
          *dptr->idunion.mcarig.s1++ = S[i][j];
         }
        }
        dptr->idunion.mcarig.s1 -= dptr->idunion.mcarig.dim1;
*/

        free_matrix(F,1,nn,1,nn);
        free_matrix(S,1,nn,1,nn);
        free_matrix(Fm1,1,nn,1,nn);
        free_matrix(tmp,1,nn,1,nn);
}

void setup_switch()
{
 int i;
 int iRc,iRl;
 _REAL S,aa,bb,R,R2;
 extern _REAL S_switch[1000];

 iRc=24;
 iRl=12;
 aa=1./(_REAL)(iRc-iRl);
 bb=(_REAL)(iRl)/(_REAL)(iRc-iRl);

 for(i=0;i<1000;i++)
 {
/*
  if(i>=iRc)
  {
   S=0.;
  } 
  else if (i<=iRl)
  {
   S=1.;
  }
  else
  {
   R=aa*i-bb;
   R2=R*R;
   S=R2*(2.*R-3)+1.;
  }
*/
  R=(_REAL)i/(_REAL)iRl;
  R2=R*R;
  S=exp(-R2*0.5);
  /* printf("i,S= %d %e\n",i,S); */
  S_switch[i]=S;
 }

}
