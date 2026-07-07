#include "head.h"
#include <time.h>

// extern char tok[];
// char tok[];
// void *get2(int,int);

extern idstruct *startid(),*getstruct();
rotstruct **rot_matrix();
rot1struct **rot1_matrix();

/************************************************************************/

void
init_mca_p2p1mm()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(MCA_P2P1MM, id);

	dptr->idunion.mca_p2p1mm.cases = 0;
	dptr->idunion.mca_p2p1mm.dump = 1;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("MCA_P2P1MM: expected file name","");
	strcpy(dptr->idunion.mca_p2p1mm.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.mca_p2p1mm.fname, 0, &exists, "MCA_P2P1MM"))
		exit(1);
	if ((dptr->idunion.mca_p2p1mm.file = 
			 fopen(dptr->idunion.mca_p2p1mm.fname, "w")) == NULL){
		perror(dptr->idunion.mca_p2p1mm.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************/
void
setup_mca_p2p1mm()
{
	idstruct	*dptr,*bptr,*mptr;
	int i,j,k,icount;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != MCA_P2P1MM)
                        inerr("setup_mca_p2p1mm: expected id or keyword\n","");
        } else {
                printf("setup_mca_p2p1mm: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((mptr = getstruct(tok)) == NULL) 
                inerr("setup_mca_p2p1mm: no such id\n","");
        if (mptr->id_type != MCARIG)
                inerr("setup_mca_p2p1mm: mcarig id expected\n","");
        gettok();
        if ((bptr = getstruct(tok)) == NULL) 
                inerr("setup_mca_p2p1mm: no such id\n","");
        if (bptr->id_type != BEAD)
                inerr("setup_mca_p2p1mm: bead id expected\n","");

        gettok();
        for(;;) {
         if(tok[0]==';') {
          break;
         } else if(!strcmp(tok,"READFSII")) {
          dptr->idunion.mca_p2p1mm.dump=-1;
          gettok();
         } else {
          printf("setup_mca_p2p1: option not found.\n");
          exit(1);
         }
        }

        dptr->idunion.mca_p2p1mm.mca = mptr;
        dptr->idunion.mca_p2p1mm.bead = bptr;

        dptr->idunion.mca_p2p1mm.type = 0;

        dptr->idunion.mca_p2p1mm.n0=mptr->idunion.mcarig.nf;
        dptr->idunion.mca_p2p1mm.n12=dptr->idunion.mca_p2p1mm.n0*
         (dptr->idunion.mca_p2p1mm.n0+1)/2;
        dptr->idunion.mca_p2p1mm.n1234=dptr->idunion.mca_p2p1mm.n0*
         (dptr->idunion.mca_p2p1mm.n0+1)/2
                                 *dptr->idunion.mca_p2p1mm.n0;

        dptr->idunion.mca_p2p1mm.dim1=
         dptr->idunion.mca_p2p1mm.n12*(dptr->idunion.mca_p2p1mm.n12+1)/2;
        dptr->idunion.mca_p2p1mm.dim12=
         dptr->idunion.mca_p2p1mm.n12*dptr->idunion.mca_p2p1mm.n1234;
        dptr->idunion.mca_p2p1mm.dim2=
         dptr->idunion.mca_p2p1mm.n1234*(dptr->idunion.mca_p2p1mm.n1234+1)/2;

        dptr->idunion.mca_p2p1mm.s1 = (_REAL *)
          get2(dptr->idunion.mca_p2p1mm.dim1,sizeof(_REAL));
        dptr->idunion.mca_p2p1mm.f1 = (_REAL *)
          get2(dptr->idunion.mca_p2p1mm.dim1,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p2p1mm.dim1;i++) {
         *dptr->idunion.mca_p2p1mm.s1++ =0.;
         *dptr->idunion.mca_p2p1mm.f1++ =0.;
        }
        dptr->idunion.mca_p2p1mm.s1-=dptr->idunion.mca_p2p1mm.dim1;
        dptr->idunion.mca_p2p1mm.f1-=dptr->idunion.mca_p2p1mm.dim1;

        dptr->idunion.mca_p2p1mm.s12 = (_REAL *)
          get2(dptr->idunion.mca_p2p1mm.dim12,sizeof(_REAL));
        dptr->idunion.mca_p2p1mm.f12 = (_REAL *)
          get2(dptr->idunion.mca_p2p1mm.dim12,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p2p1mm.dim12;i++) {
         *dptr->idunion.mca_p2p1mm.s12++ =0.;
         *dptr->idunion.mca_p2p1mm.f12++ =0.;
        }
        dptr->idunion.mca_p2p1mm.s12-=dptr->idunion.mca_p2p1mm.dim12;
        dptr->idunion.mca_p2p1mm.f12-=dptr->idunion.mca_p2p1mm.dim12;

        dptr->idunion.mca_p2p1mm.s2 = (_REAL *)
          get2(dptr->idunion.mca_p2p1mm.dim2,sizeof(_REAL));
        dptr->idunion.mca_p2p1mm.f2 = (_REAL *)
          get2(dptr->idunion.mca_p2p1mm.dim2,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p2p1mm.dim2;i++) {
         *dptr->idunion.mca_p2p1mm.s2++ =0.;
         *dptr->idunion.mca_p2p1mm.f2++ =0.;
        }
        dptr->idunion.mca_p2p1mm.s2-=dptr->idunion.mca_p2p1mm.dim2;
        dptr->idunion.mca_p2p1mm.f2-=dptr->idunion.mca_p2p1mm.dim2;

        dptr->idunion.mca_p2p1mm.m = (modstruct *)
          get2(dptr->idunion.mca_p2p1mm.n0,sizeof(modstruct));
        dptr->idunion.mca_p2p1mm.l12 = (mod12struct *)
          get2(dptr->idunion.mca_p2p1mm.n12,sizeof(mod12struct));
        dptr->idunion.mca_p2p1mm.l1234 = (mod1234struct *)
          get2(dptr->idunion.mca_p2p1mm.n1234,sizeof(mod1234struct));

        dptr->idunion.mca_p2p1mm.Fb = rot_matrix
         (0,bptr->idunion.bead.nbond-1,0,bptr->idunion.bead.nbond-1);
        dptr->idunion.mca_p2p1mm.Fm = rot_matrix
         (0,dptr->idunion.mca_p2p1mm.n0-1,0,dptr->idunion.mca_p2p1mm.n0-1);

        printf("setup_mca_p2p1mm: n0:    %d\n",dptr->idunion.mca_p2p1mm.n0);
        printf("setup_mca_p2p1mm: n12:   %d\n",dptr->idunion.mca_p2p1mm.n12);
        printf("setup_mca_p2p1mm: n1234: %d\n",dptr->idunion.mca_p2p1mm.n1234);
        printf("setup_mca_p2p1mm: dim1:  %d\n",dptr->idunion.mca_p2p1mm.dim1);
        printf("setup_mca_p2p1mm: dim12: %d\n",dptr->idunion.mca_p2p1mm.dim12);
        printf("setup_mca_p2p1mm: dim2:  %d\n",dptr->idunion.mca_p2p1mm.dim2);
}

/************************************************************************/
do_first_p2p1mm(idstruct *objptr,idstruct *serptr)
{
   _REAL *box=serptr->idunion.ser.set->boxc;
   idstruct *bptr,*mptr;
   int i1,i2,i3,i4,j1_0,j2_0,j3_0,j4_0,j1,j2,j3,j4,
       i,j,k,l,ib,jb,kk,im,m1,m2,nbondm,nmode,overlap;
   _REAL **r1, **r2;
   _REAL *rt1,*rt2,r,zi,zj,zr,
         boxm1[3],x1,y1,z1,x2,y2,z2,dx,dy,dz,scp;
   bondstruct *lb,*lb1,*lb2;
   modstruct *p,*p1,*p2,*p3,*p4;
   mod12struct *l12,*l12_1,*l12_2;
   mod1234struct *l1234,*l1234_1,*l1234_2;
   _REAL xm,ym,zm,f[5][5],Ct[3],csum,D[3][3],b,c;

   mptr=objptr->idunion.mca_p2p1mm.mca;
   bptr=mptr->idunion.mcarig.bead;
   zr = bptr->idunion.bead.zr;
   nmode=objptr->idunion.mca_p2p1mm.n0;
   nbondm=bptr->idunion.bead.nbond;
   for(k=0;k<3;k++) {
    boxm1[k]=1./box[k];
   }
   b=l_sqrt(3./8.);
   c=0.5;
/*
   b=0.5;
   c=1./l_sqrt(6.);
*/

   ++objptr->idunion.mca_p2p1mm.cases;

   lb=bptr->idunion.bead.l;
   r1=bptr->idunion.bead.r1;
   r2=bptr->idunion.bead.r2;
   for (i=0;i<nbondm;i++,r1++,r2++) {
    rt1=*r1; rt2=*r2;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
    x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
    dx=x2-x1; dy=y2-y1; dz=z2-z1;
    lb->id1 = (int)(*r1);
    lb->id2 = (int)(*r2);
    lb->x1 = x1; lb->y1 = y1; lb->z1 = z1;
    lb->x2 = x2; lb->y2 = y2; lb->z2 = z2;
    lb++;
   }

   p=objptr->idunion.mca_p2p1mm.m;
   for (j=0;j<nmode;j++) {
    xm=0.; ym=0.; zm=0.;
    r1=bptr->idunion.bead.r1;
    r2=bptr->idunion.bead.r2;
    for (i=0;i<nbondm;i++,r1++,r2++) {
     rt1=*r1; rt2=*r2;
     x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
     x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
     dx=x2-x1; dy=y2-y1; dz=z2-z1;
     xm += mptr->idunion.mcarig.evec[i+1][j+1]*dx;
     ym += mptr->idunion.mcarig.evec[i+1][j+1]*dy;
     zm += mptr->idunion.mcarig.evec[i+1][j+1]*dz;
    }
    p->mod=l_sqrt(xm*xm+ym*ym+zm*zm);
    p->m[0]=xm;
    p->m[1]=ym;
    p->m[2]=zm;
    p++;
   }

   p=objptr->idunion.mca_p2p1mm.m;
   l12=objptr->idunion.mca_p2p1mm.l12;
   for (i1=0,p1=p;i1<nmode;i1++,p1++) {
    for (i2=i1,p2=p1;i2<nmode;i2++,p2++) {
     scp = p1->m[0]*p2->m[0]+p1->m[1]*p2->m[1]+p1->m[2]*p2->m[2];
     l12->phi[0]=c*(3.*p1->m[2]*p2->m[2] - scp);
     l12->dphidm1[0][0]=-c*p2->m[0];
     l12->dphidm1[0][1]=-c*p2->m[1];
     l12->dphidm1[0][2]= 2.*c*p2->m[2];
     l12->dphidm2[0][0]=-c*p1->m[0];
     l12->dphidm2[0][1]=-c*p1->m[1];
     l12->dphidm2[0][2]= 2.*c*p1->m[2];
     l12->phi[1]=b*(p1->m[0]*p2->m[0] - p1->m[1]*p2->m[1]);
     l12->dphidm1[1][0]= b*p2->m[0];
     l12->dphidm1[1][1]=-b*p2->m[1];
     l12->dphidm1[1][2]= 0.;
     l12->dphidm2[1][0]= b*p1->m[0];
     l12->dphidm2[1][1]=-b*p1->m[1];
     l12->dphidm2[1][2]= 0.;
     l12->phi[2]=b*(p1->m[0]*p2->m[1] + p1->m[1]*p2->m[0]);
     l12->dphidm1[2][0]= b*p2->m[1];
     l12->dphidm1[2][1]= b*p2->m[0];
     l12->dphidm1[2][2]= 0.;
     l12->dphidm2[2][0]= b*p1->m[1];
     l12->dphidm2[2][1]= b*p1->m[0];
     l12->dphidm2[2][2]= 0.;
     l12->phi[3]=b*(p1->m[0]*p2->m[2] + p1->m[2]*p2->m[0]);
     l12->dphidm1[3][0]= b*p2->m[2];
     l12->dphidm1[3][1]= 0.;
     l12->dphidm1[3][2]= b*p2->m[0];
     l12->dphidm2[3][0]= b*p1->m[2];
     l12->dphidm2[3][1]= 0.;
     l12->dphidm2[3][2]= b*p1->m[0];
     l12->phi[4]=b*(p1->m[1]*p2->m[2] + p1->m[2]*p2->m[1]);
     l12->dphidm1[4][0]= 0.;
     l12->dphidm1[4][1]= b*p2->m[2];
     l12->dphidm1[4][2]= b*p2->m[1];
     l12->dphidm2[4][0]= 0.;
     l12->dphidm2[4][1]= b*p1->m[2];
     l12->dphidm2[4][2]= b*p1->m[1];
     l12++;
    }
   }

/*
*/
   l1234=objptr->idunion.mca_p2p1mm.l1234;
   l12=objptr->idunion.mca_p2p1mm.l12;
   p=objptr->idunion.mca_p2p1mm.m;
   for (i1=0,p1=p;i1<nmode;i1++,p1++) {
    for (i2=i1,p2=p1;i2<nmode;i2++,p2++,l12++) {
     for (i3=0,p3=p,p4=p;i3<nmode;i3++,p3++,p4++) {
       scp=p3->m[0]*p4->m[0]+p3->m[1]*p4->m[1]+p3->m[2]*p4->m[2];
       for(k=0;k<5;k++) {
        l1234->phi[k]=l12->phi[k]*scp;
        for(kk=0;kk<3;kk++) {
         l1234->dphidm1[k][kk]=l12->dphidm1[k][kk]*scp;
         l1234->dphidm2[k][kk]=l12->dphidm2[k][kk]*scp;
         l1234->dphidm3[k][kk]=l12->phi[k]*p4->m[kk];
         l1234->dphidm4[k][kk]=l12->phi[k]*p3->m[kk];
        }
       }
       l1234++;
     }
    }
   }
       
   lb=bptr->idunion.bead.l;
   for (ib=0,lb1=lb;ib<nbondm;ib++,lb1++) {
    for (jb=ib,lb2=lb1;jb<nbondm;jb++,lb2++) {
#include "fbb.c"
     for(m1=0;m1<3;m1++) {
      for(m2=m1;m2<3;m2++) {
       f[m2][m1]=f[m1][m2];
      }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_p2p1mm.Fb[ib][jb].f[m1][m2]=f[m1][m2];
       objptr->idunion.mca_p2p1mm.Fb[jb][ib].f[m2][m1]=f[m1][m2];
      }
     }
    }
   }

   for (i=0;i<nmode;i++) {
    for (j=0;j<nmode;j++) {
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       f[m1][m2]=0.;
      }
     }
     for (ib=0;ib<nbondm;ib++) {
      for (jb=0;jb<nbondm;jb++) {
       for(m1=0;m1<3;m1++) {
        for(m2=0;m2<3;m2++) {
         f[m1][m2]+=
          mptr->idunion.mcarig.evec[ib+1][i+1] *
          objptr->idunion.mca_p2p1mm.Fb[ib][jb].f[m1][m2] *
          mptr->idunion.mcarig.evec[jb+1][j+1];
        }
       }
      }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_p2p1mm.Fm[i][j].f[m1][m2]=f[m1][m2];
      }
     }
    }
   }

   l12=objptr->idunion.mca_p2p1mm.l12;
   for (i1=0,i=0,kk=0,l12_1=l12;i1<nmode;i1++) {
    for (i2=i1;i2<nmode;i2++,i++,l12_1++) {
     j1_0=i1;
     j2_0=i2;

   for (j1=j1_0,j=i,l12_2=l12_1;j1<nmode;j1++,j2_0=j1) {
    for (j2=j2_0;j2<nmode;j2++,j++,l12_2++) {

     for(m1=0;m1<5;m1++) {
      for(m2=0;m2<5;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<5;m1++) {
      {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i1][j1].f[k][l]  *
          l12_2->dphidm1[m2][l];
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i1][j2].f[k][l]  *
          l12_2->dphidm2[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i2][j1].f[k][l]  *
          l12_2->dphidm1[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i2][j2].f[k][l]  *
          l12_2->dphidm2[m2][l];
        }
       }
  
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     scp+=2.*f[3][3];
     scp+=2.*f[4][4];
     *objptr->idunion.mca_p2p1mm.f1 += scp/((_REAL)5);
     objptr->idunion.mca_p2p1mm.f1++;

     scp=0.;
     scp+=l12_1->phi[0]*l12_2->phi[0];
     scp+=2.*l12_1->phi[1]*l12_2->phi[1];
     scp+=2.*l12_1->phi[2]*l12_2->phi[2];
     scp+=2.*l12_1->phi[3]*l12_2->phi[3];
     scp+=2.*l12_1->phi[4]*l12_2->phi[4];
     *objptr->idunion.mca_p2p1mm.s1 += scp/((_REAL)5);
     objptr->idunion.mca_p2p1mm.s1++;

     kk++;
    }
   }
    }
   }
   objptr->idunion.mca_p2p1mm.s1 -= objptr->idunion.mca_p2p1mm.dim1;
   objptr->idunion.mca_p2p1mm.f1 -= objptr->idunion.mca_p2p1mm.dim1;

/*
*/
   l12=objptr->idunion.mca_p2p1mm.l12;
   l1234=objptr->idunion.mca_p2p1mm.l1234;
   for (i1=0,i=0,kk=0,l12_1=l12;i1<nmode;i1++) {
    for (i2=i1;i2<nmode;i2++,i++,l12_1++) {

   for (j1=0,j=0,l1234_2=l1234;j1<nmode;j1++) {
    for (j2=j1;j2<nmode;j2++) {
     for (j3=0;j3<nmode;j3++,l1234_2++,kk++) {
      j4=j3;

     for(m1=0;m1<5;m1++) {
      for(m2=0;m2<5;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<5;m1++) {
       {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i1][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i1][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i1][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i1][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i2][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i2][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i2][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i2][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     scp+=2.*f[3][3];
     scp+=2.*f[4][4];
     *objptr->idunion.mca_p2p1mm.f12 += scp/((_REAL)5);
     objptr->idunion.mca_p2p1mm.f12++;

     scp=0.;
     scp+=l12_1->phi[0]*l1234_2->phi[0];
     scp+=2.*l12_1->phi[1]*l1234_2->phi[1];
     scp+=2.*l12_1->phi[2]*l1234_2->phi[2];
     scp+=2.*l12_1->phi[3]*l1234_2->phi[3];
     scp+=2.*l12_1->phi[4]*l1234_2->phi[4];
     *objptr->idunion.mca_p2p1mm.s12 += scp/((_REAL)5);
     objptr->idunion.mca_p2p1mm.s12++;

     }
    }
   }
    }
   }
   objptr->idunion.mca_p2p1mm.s12 -= objptr->idunion.mca_p2p1mm.dim12;
   objptr->idunion.mca_p2p1mm.f12 -= objptr->idunion.mca_p2p1mm.dim12;

   l1234=objptr->idunion.mca_p2p1mm.l1234;
   for (i1=0,i=0,kk=0,l1234_1=l1234;i1<nmode;i1++) {
    for (i2=i1;i2<nmode;i2++) {
     for (i3=0;i3<nmode;i3++,l1234_1++,i++) {
     i4=i3;
     j1_0=i1;
     j2_0=i2;
     j3_0=i3;
     j4_0=i4;
   for (j1=j1_0,j=i,l1234_2=l1234_1;j1<nmode;j1++,j2_0=j1) {
    for (j2=j2_0;j2<nmode;j2++,j4_0=0,j3_0=0) {
     for (j3=j3_0;j3<nmode;j3++,j4_0=j3,l1234_2++,j++) {
      j4=j3;

     for(m1=0;m1<5;m1++) {
      for(m2=0;m2<5;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<5;m1++) {
       {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i1][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i1][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i1][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i1][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i2][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i2][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i2][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i2][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i3][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i3][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i3][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i3][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i4][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i4][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i4][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1mm.Fm[i4][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     scp+=2.*f[3][3];
     scp+=2.*f[4][4];
     *objptr->idunion.mca_p2p1mm.f2 += scp/((_REAL)5);
     objptr->idunion.mca_p2p1mm.f2++;

     scp=0.;
     scp+=l1234_1->phi[0]*l1234_2->phi[0];
     scp+=2.*l1234_1->phi[1]*l1234_2->phi[1];
     scp+=2.*l1234_1->phi[2]*l1234_2->phi[2];
     scp+=2.*l1234_1->phi[3]*l1234_2->phi[3];
     scp+=2.*l1234_1->phi[4]*l1234_2->phi[4];
     *objptr->idunion.mca_p2p1mm.s2 += scp/((_REAL)5);
     objptr->idunion.mca_p2p1mm.s2++;

       kk++;
     }
    }
   }
     }
    }
   }
   objptr->idunion.mca_p2p1mm.s2 -= objptr->idunion.mca_p2p1mm.dim2;
   objptr->idunion.mca_p2p1mm.f2 -= objptr->idunion.mca_p2p1mm.dim2;

}

/************************************************************************/
void
write_first_p2p1mm(idstruct *dptr)
{
	int	i,j,i2,j2,k,l,m;
	int	states = dptr->idunion.mca_p2p1mm.cases;
	FILE 	*file = dptr->idunion.mca_p2p1mm.file,*dump;
        idstruct *mptr,*bptr;

	_REAL fact,tol,factns,sum,sum2;
        int rank,nbondm;

        int n12,n1234,nn,nf;
        _REAL **S;
        _REAL **F,**Fm1,**tmp,csum,ctmp;
        _REAL *ui,*v,*w;
#ifdef IMSL
        d_complex *wi,*vi;
#else
        _REAL *vi;
#endif

        n12=dptr->idunion.mca_p2p1mm.n12;
        n1234=dptr->idunion.mca_p2p1mm.n1234;
        mptr=dptr->idunion.mca_p2p1mm.mca;
        bptr=mptr->idunion.mcarig.bead;
	fact = 1./ ((_REAL)states);
        factns = bptr->idunion.bead.factns;

        printf("write_first_p2p1: number of states: %d\n",states);
        for (i=1; i<=n12; i++) {
         for (j=i; j<=n12; j++) {
           *dptr->idunion.mca_p2p1mm.s1++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1mm.s1 -= dptr->idunion.mca_p2p1mm.dim1;
        for (i=1; i<=n12; i++) {
         for (j=i; j<=n12; j++) {
           *dptr->idunion.mca_p2p1mm.f1++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1mm.f1 -= dptr->idunion.mca_p2p1mm.dim1;
/*
*/
        for (i=1; i<=n12; i++) {
         for (j=1; j<=n1234; j++) {
           *dptr->idunion.mca_p2p1mm.s12++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1mm.s12 -= dptr->idunion.mca_p2p1mm.dim12;
        for (i=1; i<=n12; i++) {
         for (j=1; j<=n1234; j++) {
           *dptr->idunion.mca_p2p1mm.f12++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1mm.f12 -= dptr->idunion.mca_p2p1mm.dim12;
        for (i=1; i<=n1234; i++) {
         for (j=i; j<=n1234; j++) {
           *dptr->idunion.mca_p2p1mm.s2++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1mm.s2 -= dptr->idunion.mca_p2p1mm.dim2;
        for (i=1; i<=n1234; i++) {
         for (j=i; j<=n1234; j++) {
           *dptr->idunion.mca_p2p1mm.f2++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1mm.f2 -= dptr->idunion.mca_p2p1mm.dim2;

        if(dptr->idunion.mca_p2p1mm.dump==1) {
         dump=fopen("fsII.dump","w");
         fwrite(dptr->idunion.mca_p2p1mm.f1,dptr->idunion.mca_p2p1mm.dim1*sizeof(_REAL),
                1,dump);
         fwrite(dptr->idunion.mca_p2p1mm.s1,dptr->idunion.mca_p2p1mm.dim1*sizeof(_REAL),
                1,dump);
/*
*/
         fwrite(dptr->idunion.mca_p2p1mm.f12,dptr->idunion.mca_p2p1mm.dim12*sizeof(_REAL),
                1,dump);
         fwrite(dptr->idunion.mca_p2p1mm.s12,dptr->idunion.mca_p2p1mm.dim12*sizeof(_REAL),
                1,dump);
         fwrite(dptr->idunion.mca_p2p1mm.f2,dptr->idunion.mca_p2p1mm.dim2*sizeof(_REAL),
                1,dump);
         fwrite(dptr->idunion.mca_p2p1mm.s2,dptr->idunion.mca_p2p1mm.dim2*sizeof(_REAL),
                1,dump);
         fflush(dump);
         fclose(dump);
        } else if(dptr->idunion.mca_p2p1mm.dump==-1) {
         dump=fopen("fsII.dump","r");
         fread(dptr->idunion.mca_p2p1mm.f1,dptr->idunion.mca_p2p1mm.dim1*sizeof(_REAL),
                1,dump);
         fread(dptr->idunion.mca_p2p1mm.s1,dptr->idunion.mca_p2p1mm.dim1*sizeof(_REAL),
                1,dump);
/*
*/
         fread(dptr->idunion.mca_p2p1mm.f12,dptr->idunion.mca_p2p1mm.dim12*sizeof(_REAL),
                1,dump);
         fread(dptr->idunion.mca_p2p1mm.s12,dptr->idunion.mca_p2p1mm.dim12*sizeof(_REAL),
                1,dump);
         fread(dptr->idunion.mca_p2p1mm.f2,dptr->idunion.mca_p2p1mm.dim2*sizeof(_REAL),
                1,dump);
         fread(dptr->idunion.mca_p2p1mm.s2,dptr->idunion.mca_p2p1mm.dim2*sizeof(_REAL),
                1,dump);
         fclose(dump);
        }

/* invert F */
        nn = n12+n1234;
        /* nn = n12; */
        F = matrix(1,nn,1,nn);
        S = matrix(1,nn,1,nn);

        for(i=1;i<=n12;i++) {
         for(j=i;j<=n12;j++) {
          F[i][j]=*dptr->idunion.mca_p2p1mm.f1++;
          F[j][i]=F[i][j];
          S[i][j]=*dptr->idunion.mca_p2p1mm.s1++;
          S[j][i]=S[i][j];
         }
        }
        dptr->idunion.mca_p2p1mm.s1 -= dptr->idunion.mca_p2p1mm.dim1;
        dptr->idunion.mca_p2p1mm.f1 -= dptr->idunion.mca_p2p1mm.dim1;

/*
*/
        for(i=1;i<=n12;i++) {
         for(j=1;j<=n1234;j++) {
          j2=n12+j;
          F[i][j2]=*dptr->idunion.mca_p2p1mm.f12++;
          F[j2][i]=F[i][j2];
          S[i][j2]=*dptr->idunion.mca_p2p1mm.s12++;
          S[j2][i]=S[i][j2];
         }
        }
        dptr->idunion.mca_p2p1mm.s12 -= dptr->idunion.mca_p2p1mm.dim12;
        dptr->idunion.mca_p2p1mm.f12 -= dptr->idunion.mca_p2p1mm.dim12;

        for(i=1;i<=n1234;i++) {
         for(j=i;j<=n1234;j++) {
          i2=n12+i;
          j2=n12+j;
          F[i2][j2]=*dptr->idunion.mca_p2p1mm.f2++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*dptr->idunion.mca_p2p1mm.s2++;
          S[j2][i2]=S[i2][j2];
         }
        }
        dptr->idunion.mca_p2p1mm.s2 -= dptr->idunion.mca_p2p1mm.dim2;
        dptr->idunion.mca_p2p1mm.f2 -= dptr->idunion.mca_p2p1mm.dim2;

/*
*/
        fprintf(file,"S,F:\n");
        for(i=1;i<=nn;i++)
         for(j=1;j<=nn;j++)
          fprintf(file,"%d %d %e %e\n",i,j,S[i][j],F[i][j]);
        fflush(file);

        for(i=1;i<=nn;i++) {
         if(S[i][i]<ZERO) {
          printf("write_first:\n");
          printf("S element is wrong: S(%d,%d) = %e\n",i,i,S[i][i]);
         }
        }

        fflush(stdout);
        Fm1 = matrix(1,nn,1,nn);
#ifdef IMSL
        tol=TOLSVD;
        ui = (_REAL *)get2(nn*nn,sizeof(_REAL));
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
          ui[k]=F[i][j];
         }
        }
        w = imsl_d_lin_svd_gen(nn,nn,ui,
               IMSL_RANK,tol,&rank,
               IMSL_INVERSE,&v,0);
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
          Fm1[i][j]=v[k];
         }
        }
        free(v);
        free(w);
#else
        rank=invert(F,nn,nn,Fm1);
#endif
        fprintf(file,"write_first_p2: rank of F= %d\n",rank); fflush(file);

        tmp = matrix(1,nn,1,nn);

/*
        fprintf(file,"Fm1*F:\n");
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          csum=0.;
          for(k=1;k<=nn;k++) {
           csum+=Fm1[i][k]*F[k][j];
          }
          fprintf(file,"%d %d %e\n",i,j,csum);
         }
        }
        fflush(file);
*/

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          csum=0.;
          for(k=1;k<=nn;k++) {
           csum+=Fm1[i][k]*S[k][j];
          }
          tmp[i][j]=csum;
         }
        }
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          csum=0.;
          for(k=1;k<=nn;k++) {
           csum+=S[i][k]*tmp[k][j];
          }
          Fm1[i][j]=csum;
         }
        }
        for(i=1;i<=nn;i++) {
         tol=1./S[i][i];
         csum=tol*Fm1[i][i];
         fprintf(file,"tau(%d)= %e\n",i,factns*csum);
        }
        fflush(file);

        w=vector(1,nn);
#ifdef IMSL
        for(i=1,l=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,l++) {
          ui[l]=tmp[i][j];
         }
        }
        wi = imsl_d_eig_gen(nn,ui,IMSL_VECTORS,&vi,0);
        for(i=0;i<nn;i++) {
         w[i+1]=wi[i].re;
        }
        free(ui);
        free(wi);
#else
        vi=(_REAL *)get2(nn*nn,sizeof(_REAL));
#ifdef EISPACK
        eigvec_ev(nn,tmp,vi,w);
#else
        eigvec_nr(tmp,nn,w,vi);
#endif
#endif

        fprintf(file,"lambda: ");
        nf=0;
        for(i=1;i<=nn;i++) {
         fprintf(file,"%e ",1./w[i]);
         if(Nabs(w[i])>=ZERO) {
          nf++;
         }
        }
        fprintf(file,"\n");
        fflush(file);

        if(nf!=dptr->idunion.mca_p2p1mm.n12) {
         printf("write_first_p2p1: WARNING: nf=%d is different from n12=%d.\n",
                 nf,n12);
        }
        fflush(file);

        dptr->idunion.mca_p2p1mm.nf = nf;
        nf=nn;
        dptr->idunion.mca_p2p1mm.eval= (_REAL *)get2(nf,sizeof(_REAL));
        dptr->idunion.mca_p2p1mm.evec= matrix(1,nn,1,nf);

        nf=0;
        for(i=1;i<=nn;i++) {
         dptr->idunion.mca_p2p1mm.eval[i-1]=1./w[i];
        }
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
#ifdef IMSL
          dptr->idunion.mca_p2p1mm.evec[i][j]=vi[k].re;
#else
          dptr->idunion.mca_p2p1mm.evec[i][j]=vi[k];
#endif
         }
        }
/*
        free_vector(w,1,nn);
        free(vi);
*/
/*
        for(i=1,k=0;i<=nn;i++,k++) {
         if(Nabs(wi[i-1].re)>=ZERO)
         {
          nf++;
          dptr->idunion.mca_p2p1mm.eval[nf-1]=1./wi[i-1].re;
          vt=&vi[k];
          for(j=1;j<=nn;j++) {
           dptr->idunion.mca_p2p1mm.evec[j][nf]=(*vt).re;
           vt+=nn;
          }
         }
        }
*/

/*
*/
        fprintf(file,"evec:\n");
        for(j=1;j<=5;j++) {
         for(i=1;i<=nn;i++) {
          fprintf(file,"%d %d %e\n",i,j,dptr->idunion.mca_p2p1mm.evec[i][j]);
         }
        }

        nf=nn;
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mca_p2p1mm.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mca_p2p1mm.evec[k][i]*tmp[k][j];
          }
          Fm1[i][j]=sum;
         }
        }
/*
        fprintf(file,"Ct*S*C:\n");
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nf;j++) {
          fprintf(file,"%d %d %e\n",i,j,Fm1[i][j]);
         }
        }
        fflush(file);
*/

/*
        for(i=1,k=0;i<=nf;i++) {
         for(j=1;j<=nf;j++,k++) {
          ui[k]=Fm1[i][j];
         }
        }
        w = imsl_d_eig_sym(nf,ui,IMSL_VECTORS,&v,0);
        for(i=1,k=0;i<=nf;i++) {
         if(Nabs(w[i-1]) > ZERO)
         {
          w[i-1]=1./l_sqrt(w[i-1]);
          for(j=1;j<=nf;j++,k++) {
           tmp[i][j]=v[k];
          }
         } else {
          w[i-1]=1.;
          for(j=1;j<=nf;j++,k++) {
           tmp[i][j]=0.;
          }
          tmp[i][i]=1.;
         }
        }
        for(i=1,m=0;i<=nf;i++) {
         for(j=1;j<=nf;j++,m++) {
          for(k=1,sum=0.;k<=nf;k++) {
           for(l=1,sum2=0.;l<=nf;l++) {
            sum2+=tmp[k][l]*w[l-1]*tmp[j][l];
           }
           sum+=dptr->idunion.mca_p2p1mm.evec[i][k]*sum2;
          }
          v[m]=sum;
         }
        }
        for(i=1,k=0;i<=nf;i++) {
         for(j=1;j<=nf;j++,k++) {
          dptr->idunion.mca_p2p1mm.evec[i][j]=v[k];
         }
        }
*/

        for(j=1,k=0;j<=nn;j++) {
         if(Nabs(Fm1[j][j]) > ZERO) {
          k++;
          sum=Fm1[j][j];
         } else {
          sum=1.;
         }
         sum=l_sqrt(sum);
         for(i=1;i<=nn;i++) {
          dptr->idunion.mca_p2p1mm.evec[i][j]/=sum;
         }
        }
        if(k!=nf) {
         printf("write_first: missed zero in Ct*S*C.\n");
        }

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mca_p2p1mm.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mca_p2p1mm.evec[k][i]*tmp[k][j];
          }
          S[i][j]=sum;
         }
        }
/*
        fprintf(file,"Ct*S1*C:\n");
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nf;j++) {
          fprintf(file,"%d %d %e\n",i,j,S[i][j]);
         }
        }
        fflush(file);
*/

        free_matrix(F,1,nn,1,nn);
        free_matrix(S,1,nn,1,nn);
        free_matrix(tmp,1,nn,1,nn);
}
