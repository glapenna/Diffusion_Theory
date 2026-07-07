#include "head.h"

// extern char tok[];
// char tok[];
// void *get2(int,int);

extern idstruct *startid(),*getstruct();
rotstruct **rot_matrix();
rot1struct **rot1_matrix();

/************************************************************************/

void
init_mca_p1l()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(MCA_P1L, id);

	dptr->idunion.mca_p1l.cases = 0;
	dptr->idunion.mca_p1l.dump = 0;
	dptr->idunion.mca_p1l.type = 1;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("MCA_P1L: expected file name","");
	strcpy(dptr->idunion.mca_p1l.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.mca_p1l.fname, 0, &exists, "MCA_P1L"))
		exit(1);
	if ((dptr->idunion.mca_p1l.file = 
			 fopen(dptr->idunion.mca_p1l.fname, "w")) == NULL){
		perror(dptr->idunion.mca_p1l.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************/
void
setup_mca_p1l()
{
	idstruct	*dptr,*bptrm,*bptrl,*mptr;
	int i,j,k,icount;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != MCA_P1L)
                        inerr("setup_mca_p1l: expected id or keyword\n","");
        } else {
                printf("setup_mca_p1l: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((mptr = getstruct(tok)) == NULL) 
                inerr("setup_mca_p1l: no such id\n","");
        if (mptr->id_type != MCARIG)
                inerr("setup_mca_p1l: mcarig id expected\n","");
        gettok();
        if ((bptrm = getstruct(tok)) == NULL) 
                inerr("setup_mca_p1l: no such id\n","");
        if (bptrm->id_type != BEAD)
                inerr("setup_mca_p1l: beadm id expected\n","");
        gettok();
        if ((bptrl = getstruct(tok)) == NULL) 
                inerr("setup_mca_p1l: no such id\n","");
        if (bptrl->id_type != BEAD)
                inerr("setup_mca_p1l: beadl id expected\n","");

        gettok();
        for(;;) {
         if(tok[0]==';') {
          break;
         } else if(!strcmp(tok,"READFSII")) {
          dptr->idunion.mca_p1l.dump=-1;
          gettok();
         } else if(!strcmp(tok,"RM1-I")) {
          dptr->idunion.mca_p1l.type=1;
          gettok();
         } else if(!strcmp(tok,"RM1-II")) {
          dptr->idunion.mca_p1l.type=2;
          gettok();
         } else {
          printf("setup_mca_p1l: option not found.\n");
          exit(1);
         }
        }

        dptr->idunion.mca_p1l.mca = mptr;
        dptr->idunion.mca_p1l.beadm = bptrm;
        dptr->idunion.mca_p1l.beadl = bptrl;

        dptr->idunion.mca_p1l.n0m=mptr->idunion.mcarig.nf;
        dptr->idunion.mca_p1l.n0l=bptrl->idunion.bead.nbond;
        dptr->idunion.mca_p1l.n1= dptr->idunion.mca_p1l.n0m
          +dptr->idunion.mca_p1l.n0l;
        dptr->idunion.mca_p1l.n123m= dptr->idunion.mca_p1l.n0m
                                  *dptr->idunion.mca_p1l.n0m
                                  *(dptr->idunion.mca_p1l.n0m+1)/2;

        dptr->idunion.mca_p1l.dim_m_m=
         dptr->idunion.mca_p1l.n0m*(dptr->idunion.mca_p1l.n0m+1)/2;
        dptr->idunion.mca_p1l.dim_m_l=
         dptr->idunion.mca_p1l.n0m*dptr->idunion.mca_p1l.n0l;
        dptr->idunion.mca_p1l.dim_l_l=
         dptr->idunion.mca_p1l.n0l*(dptr->idunion.mca_p1l.n0l+1)/2;
        dptr->idunion.mca_p1l.dim_m_m2=
         dptr->idunion.mca_p1l.n0m*dptr->idunion.mca_p1l.n123m;
        dptr->idunion.mca_p1l.dim_l_m2=
         dptr->idunion.mca_p1l.n0l*dptr->idunion.mca_p1l.n123m;
        dptr->idunion.mca_p1l.dim_m2_m2=
         dptr->idunion.mca_p1l.n123m*(dptr->idunion.mca_p1l.n123m+1)/2;

        dptr->idunion.mca_p1l.s_m_m = (_REAL *)
          get2(dptr->idunion.mca_p1l.dim_m_m,sizeof(_REAL));
        dptr->idunion.mca_p1l.f_m_m = (_REAL *)
          get2(dptr->idunion.mca_p1l.dim_m_m,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p1l.dim_m_m;i++) {
         *dptr->idunion.mca_p1l.s_m_m++ =0.;
         *dptr->idunion.mca_p1l.f_m_m++ =0.;
        }
        dptr->idunion.mca_p1l.s_m_m-=dptr->idunion.mca_p1l.dim_m_m;
        dptr->idunion.mca_p1l.f_m_m-=dptr->idunion.mca_p1l.dim_m_m;

        dptr->idunion.mca_p1l.s_m_l = (_REAL *)
          get2(dptr->idunion.mca_p1l.dim_m_l,sizeof(_REAL));
        dptr->idunion.mca_p1l.f_m_l = (_REAL *)
          get2(dptr->idunion.mca_p1l.dim_m_l,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p1l.dim_m_l;i++) {
         *dptr->idunion.mca_p1l.s_m_l++ =0.;
         *dptr->idunion.mca_p1l.f_m_l++ =0.;
        }
        dptr->idunion.mca_p1l.s_m_l-=dptr->idunion.mca_p1l.dim_m_l;
        dptr->idunion.mca_p1l.f_m_l-=dptr->idunion.mca_p1l.dim_m_l;

        dptr->idunion.mca_p1l.s_l_l = (_REAL *)
          get2(dptr->idunion.mca_p1l.dim_l_l,sizeof(_REAL));
        dptr->idunion.mca_p1l.f_l_l = (_REAL *)
          get2(dptr->idunion.mca_p1l.dim_l_l,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p1l.dim_l_l;i++) {
         *dptr->idunion.mca_p1l.s_l_l++ =0.;
         *dptr->idunion.mca_p1l.f_l_l++ =0.;
        }
        dptr->idunion.mca_p1l.s_l_l-=dptr->idunion.mca_p1l.dim_l_l;
        dptr->idunion.mca_p1l.f_l_l-=dptr->idunion.mca_p1l.dim_l_l;

        dptr->idunion.mca_p1l.s_m_m2 = (_REAL *)
          get2(dptr->idunion.mca_p1l.dim_m_m2,sizeof(_REAL));
        dptr->idunion.mca_p1l.f_m_m2 = (_REAL *)
          get2(dptr->idunion.mca_p1l.dim_m_m2,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p1l.dim_m_m2;i++) {
         *dptr->idunion.mca_p1l.s_m_m2++ =0.;
         *dptr->idunion.mca_p1l.f_m_m2++ =0.;
        }
        dptr->idunion.mca_p1l.s_m_m2-=dptr->idunion.mca_p1l.dim_m_m2;
        dptr->idunion.mca_p1l.f_m_m2-=dptr->idunion.mca_p1l.dim_m_m2;

        dptr->idunion.mca_p1l.s_l_m2 = (_REAL *)
          get2(dptr->idunion.mca_p1l.dim_l_m2,sizeof(_REAL));
        dptr->idunion.mca_p1l.f_l_m2 = (_REAL *)
          get2(dptr->idunion.mca_p1l.dim_l_m2,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p1l.dim_l_m2;i++) {
         *dptr->idunion.mca_p1l.s_l_m2++ =0.;
         *dptr->idunion.mca_p1l.f_l_m2++ =0.;
        }
        dptr->idunion.mca_p1l.s_l_m2-=dptr->idunion.mca_p1l.dim_l_m2;
        dptr->idunion.mca_p1l.f_l_m2-=dptr->idunion.mca_p1l.dim_l_m2;

        dptr->idunion.mca_p1l.s_m2_m2 = (_REAL *)
          get2(dptr->idunion.mca_p1l.dim_m2_m2,sizeof(_REAL));
        dptr->idunion.mca_p1l.f_m2_m2 = (_REAL *)
          get2(dptr->idunion.mca_p1l.dim_m2_m2,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p1l.dim_m2_m2;i++) {
         *dptr->idunion.mca_p1l.s_m2_m2++ =0.;
         *dptr->idunion.mca_p1l.f_m2_m2++ =0.;
        }
        dptr->idunion.mca_p1l.s_m2_m2-=dptr->idunion.mca_p1l.dim_m2_m2;
        dptr->idunion.mca_p1l.f_m2_m2-=dptr->idunion.mca_p1l.dim_m2_m2;

        dptr->idunion.mca_p1l.m = (modstruct *)
          get2(dptr->idunion.mca_p1l.n0m,sizeof(modstruct));
        dptr->idunion.mca_p1l.l1 = (mod1struct *)
          get2(dptr->idunion.mca_p1l.n1,sizeof(mod1struct));
        dptr->idunion.mca_p1l.l123 = (mod123struct *)
          get2(dptr->idunion.mca_p1l.n123m,sizeof(mod123struct));

        printf("setup_mca_p1l: boundaries of data.\n");
        printf("n0m:    %d\n",dptr->idunion.mca_p1l.n0m);
        printf("n123m:  %d\n",dptr->idunion.mca_p1l.n123m);
        printf("n0l:    %d\n",dptr->idunion.mca_p1l.n0l);

        dptr->idunion.mca_p1l.Fb_m_m = rot_matrix
         (0,bptrm->idunion.bead.nbond-1,0,bptrm->idunion.bead.nbond-1);
        dptr->idunion.mca_p1l.Fb_m_l = rot_matrix
         (0,bptrm->idunion.bead.nbond-1,0,bptrl->idunion.bead.nbond-1);
        dptr->idunion.mca_p1l.Fb_l_l = rot_matrix
         (0,bptrl->idunion.bead.nbond-1,0,bptrl->idunion.bead.nbond-1);

        dptr->idunion.mca_p1l.Fm_m_m = rot_matrix
         (0,dptr->idunion.mca_p1l.n0m-1,0,dptr->idunion.mca_p1l.n0m-1);
        dptr->idunion.mca_p1l.Fm_m_l = rot_matrix
         (0,dptr->idunion.mca_p1l.n0m-1,0,bptrl->idunion.bead.nbond-1);
        dptr->idunion.mca_p1l.Fm_l_m = rot_matrix
         (0,bptrl->idunion.bead.nbond-1,0,dptr->idunion.mca_p1l.n0m-1);

        printf("dim_m_m:   %d\n",dptr->idunion.mca_p1l.dim_m_m);
        printf("dim_m_l:   %d\n",dptr->idunion.mca_p1l.dim_m_l);
        printf("dim_m_m2:  %d\n",dptr->idunion.mca_p1l.dim_m_m2);
        printf("dim_l_l:   %d\n",dptr->idunion.mca_p1l.dim_l_l);
        printf("dim_l_m2:  %d\n",dptr->idunion.mca_p1l.dim_l_m2);
        printf("dim_m2_m2: %d\n",dptr->idunion.mca_p1l.dim_m2_m2);

}

/************************************************************************/
do_first_p1l(idstruct *objptr,idstruct *serptr)
{
   _REAL *box=serptr->idunion.ser.set->boxc;
   idstruct *bptrm,*bptrl,*bptr,*mptr;
   int i1,i2,i3,i4,j1_0,j2_0,j3_0,j4_0,j1,j2,j3,j4,
       i,j,k,l,ib,jb,kk,im,m1,m2,nbondm,nbondl,nmode,overlap;
   _REAL **r1, **r2;
   _REAL *rt1,*rt2,r,zi,zj,zr,cutsq,
         boxm1[3],x1,y1,z1,x2,y2,z2,dx,dy,dz,scp;
   bondstruct *lb,*lb1,*lb2,*lb1_0,*lb2_0;
   modstruct *p,*p1,*p2,*p3,*p4;
   mod1struct *l1,*l1_1,*l1_2,*l1m_0,*l1l_0;
   mod123struct *l123,*l123_1,*l123_2;
   _REAL xm,ym,zm,f[3][3],Ct[3],csum,D[3][3],b,c;

   mptr=objptr->idunion.mca_p1l.mca;
   bptrm=objptr->idunion.mca_p1l.beadm;
   bptrl=objptr->idunion.mca_p1l.beadl;
   bptr=bptrm;
   zr=bptr->idunion.bead.zr;
   cutsq=bptr->idunion.bead.cutsq;
   nmode=objptr->idunion.mca_p1l.n0m;
   nbondm=bptrm->idunion.bead.nbond;
   nbondl=bptrl->idunion.bead.nbond;
   for(k=0;k<3;k++) {
    boxm1[k]=1./box[k];
   }
   b=1./l_sqrt(2.);

   ++objptr->idunion.mca_p1l.cases;

   lb=bptrm->idunion.bead.l;
   r1=bptrm->idunion.bead.r1;
   r2=bptrm->idunion.bead.r2;
   for (i=0;i<nbondm;i++,r1++,r2++) {
    rt1=*r1; rt2=*r2;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
    x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
    dx=x2-x1; dy=y2-y1; dz=z2-z1;
    lb->x=dx; lb->y=dy; lb->z=dz;
    lb->id1 = (int)(*r1);
    lb->id2 = (int)(*r2);
    lb->x1 = x1; lb->y1 = y1; lb->z1 = z1;
    lb->x2 = x2; lb->y2 = y2; lb->z2 = z2;
    lb++;
   }
   lb=bptrl->idunion.bead.l;
   r1=bptrl->idunion.bead.r1;
   r2=bptrl->idunion.bead.r2;
   for (i=0;i<nbondl;i++,r1++,r2++) {
    rt1=*r1; rt2=*r2;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
    x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
    dx=x2-x1; dy=y2-y1; dz=z2-z1;
    lb->x=dx; lb->y=dy; lb->z=dz;
    lb->id1 = (int)(*r1);
    lb->id2 = (int)(*r2);
    lb->x1 = x1; lb->y1 = y1; lb->z1 = z1;
    lb->x2 = x2; lb->y2 = y2; lb->z2 = z2;
    lb++;
   }

   p=objptr->idunion.mca_p1l.m;
   for (j=0;j<nmode;j++) {
    xm=0.; ym=0.; zm=0.;
    lb=bptrm->idunion.bead.l;
/*
    r1=bptrm->idunion.bead.r1;
    r2=bptrm->idunion.bead.r2;
*/
    for (i=0;i<nbondm;i++,r1++,r2++) {
/*
     rt1=*r1; rt2=*r2;
     x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
     x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
     dx=x2-x1; dy=y2-y1; dz=z2-z1;
*/
     dx=lb->x; dy=lb->y; dz=lb->z;
     xm += mptr->idunion.mcarig.evec[i+1][j+1]*dx;
     ym += mptr->idunion.mcarig.evec[i+1][j+1]*dy;
     zm += mptr->idunion.mcarig.evec[i+1][j+1]*dz;
     lb++;
    }
    p->mod=l_sqrt(xm*xm+ym*ym+zm*zm);
    p->m[0]=xm;
    p->m[1]=ym;
    p->m[2]=zm;
    p++;
   }

   p=objptr->idunion.mca_p1l.m;
   l1=objptr->idunion.mca_p1l.l1;
   l1m_0=l1;
   for (i1=0,p1=p;i1<nmode;i1++,p1++) {
     l1->phi[0]=p1->m[2];
     l1->dphidm1[0][0]=0.;
     l1->dphidm1[0][1]=0.;
     l1->dphidm1[0][2]=1.;
     l1->phi[1]=b*p1->m[0];
     l1->dphidm1[1][0]= b;
     l1->dphidm1[1][1]= 0.;
     l1->dphidm1[1][2]= 0.;
     l1->phi[2]=b*p1->m[1];
     l1->dphidm1[2][0]= 0.;
     l1->dphidm1[2][1]= b;
     l1->dphidm1[2][2]= 0.;
     l1++;
   }

   if(objptr->idunion.mca_p1l.type>1) {
   l123=objptr->idunion.mca_p1l.l123;
   p=objptr->idunion.mca_p1l.m;
   for (i1=0,k=0,p1=p;i1<nmode;i1++,p1++) {
    for (i2=0,p2=p;i2<nmode;i2++,p2++) {
     for (i3=i2,p3=p2;i3<nmode;i3++,p3++,k++) {
      scp=p2->m[0]*p3->m[0]+p2->m[1]*p3->m[1]+p2->m[2]*p3->m[2];
      for(k=0;k<3;k++) {
       l123->phi[k]=p1->m[k]*scp;
       for(kk=0;kk<3;kk++) {
        l123->dphidm1[k][kk]=rdelta(k,kk)*scp;
        l123->dphidm2[k][kk]=p1->m[k]*p3->m[kk];
        l123->dphidm3[k][kk]=p1->m[k]*p2->m[kk];
       }
      }
      l123++;
     }
    }
   }
   }  

   lb=bptrl->idunion.bead.l;
   l1l_0=l1;
   for (i1=0,lb1=lb;i1<nbondl;i1++,lb1++) {
     l1->phi[0]=lb1->z;
     l1->dphidm1[0][0]= 0.;
     l1->dphidm1[0][1]= 0.;
     l1->dphidm1[0][2]= 1.;
     l1->phi[1]=b*lb1->x;
     l1->dphidm1[1][0]= b;
     l1->dphidm1[1][1]= 0.;
     l1->dphidm1[1][2]= 0.;
     l1->phi[2]=b*lb1->y;
     l1->dphidm1[2][0]= 0.;
     l1->dphidm1[2][1]= b;
     l1->dphidm1[2][2]= 0.;
     l1++;
   }

   lb=bptrm->idunion.bead.l;
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
       objptr->idunion.mca_p1l.Fb_m_m[ib][jb].f[m1][m2]=f[m1][m2];
       objptr->idunion.mca_p1l.Fb_m_m[jb][ib].f[m2][m1]=f[m1][m2];
      }
     }
    }
   }

   lb1_0=bptrm->idunion.bead.l;
   lb2_0=bptrl->idunion.bead.l;
   for (ib=0,lb1=lb1_0;ib<nbondm;ib++,lb1++) {
    for (jb=0,lb2=lb2_0;jb<nbondl;jb++,lb2++) {
#include "fbb.c"
     for(m1=0;m1<3;m1++) {
      for(m2=m1;m2<3;m2++) {
       f[m2][m1]=f[m1][m2];
      }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_p1l.Fb_m_l[ib][jb].f[m1][m2]=f[m1][m2];
       /* objptr->idunion.mca_p1l.Fb_m_l[jb][ib].f[m2][m1]=f[m1][m2]; */
      }
     }
    }
   }

   lb=bptrl->idunion.bead.l;
   for (ib=0,lb1=lb;ib<nbondl;ib++,lb1++) {
    for (jb=ib,lb2=lb1;jb<nbondl;jb++,lb2++) {
#include "fbb.c"
     for(m1=0;m1<3;m1++) {
      for(m2=m1;m2<3;m2++) {
       f[m2][m1]=f[m1][m2];
      }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_p1l.Fb_l_l[ib][jb].f[m1][m2]=f[m1][m2];
       objptr->idunion.mca_p1l.Fb_l_l[jb][ib].f[m2][m1]=f[m1][m2];
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
          objptr->idunion.mca_p1l.Fb_m_m[ib][jb].f[m1][m2] *
          mptr->idunion.mcarig.evec[jb+1][j+1];
        }
       }
      }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_p1l.Fm_m_m[i][j].f[m1][m2]=f[m1][m2];
      }
     }
    }
   }

   for (i=0;i<nmode;i++) {
    for (j=0;j<nbondl;j++) {
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       f[m1][m2]=0.;
      }
     }
     for (ib=0;ib<nbondm;ib++) {
       for(m1=0;m1<3;m1++) {
        for(m2=0;m2<3;m2++) {
         f[m1][m2]+=
          mptr->idunion.mcarig.evec[ib+1][i+1] *
          objptr->idunion.mca_p1l.Fb_m_l[ib][j].f[m1][m2];
        }
       }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_p1l.Fm_m_l[i][j].f[m1][m2]=f[m1][m2];
      }
     }
    }
   }

   for (i=0;i<nbondl;i++) {
    for (j=0;j<nmode;j++) {
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       f[m1][m2]=0.;
      }
     }
     for (ib=0;ib<nbondm;ib++) {
       for(m1=0;m1<3;m1++) {
        for(m2=0;m2<3;m2++) {
         f[m1][m2]+=
          objptr->idunion.mca_p1l.Fb_m_l[ib][i].f[m1][m2] *
          mptr->idunion.mcarig.evec[ib+1][j+1];
        }
       }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_p1l.Fm_l_m[i][j].f[m1][m2]=f[m1][m2];
      }
     }
    }
   }

   for (i1=0,i=0,kk=0,l1_1=l1m_0;i1<nmode;i1++,l1_1++) {
    j1_0=i1;
    for (j1=j1_0,j=i,l1_2=l1_1;j1<nmode;j1++,kk++,l1_2++) {

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
          l1_1->dphidm1[m1][k] *
          objptr->idunion.mca_p1l.Fm_m_m[i1][j1].f[k][l]  *
          l1_2->dphidm1[m2][l];
        }
       }
  
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     *objptr->idunion.mca_p1l.f_m_m += scp/((_REAL)3);
     objptr->idunion.mca_p1l.f_m_m++;

     scp=0.;
     scp+=l1_1->phi[0]*l1_2->phi[0];
     scp+=2.*l1_1->phi[1]*l1_2->phi[1];
     scp+=2.*l1_1->phi[2]*l1_2->phi[2];
     *objptr->idunion.mca_p1l.s_m_m += scp/3;
     objptr->idunion.mca_p1l.s_m_m++;

    }
   }
   objptr->idunion.mca_p1l.s_m_m -= objptr->idunion.mca_p1l.dim_m_m;
   objptr->idunion.mca_p1l.f_m_m -= objptr->idunion.mca_p1l.dim_m_m;

   for (i1=0,i=0,kk=0,l1_1=l1m_0;i1<nmode;i1++,l1_1++) {
    for (j1=0,j=0,l1_2=l1l_0;j1<nbondl;j1++,l1_2++,kk++) {

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
          l1_1->dphidm1[m1][k] *
          objptr->idunion.mca_p1l.Fm_m_l[i1][j1].f[k][l]  *
          l1_2->dphidm1[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     *objptr->idunion.mca_p1l.f_m_l += scp/((_REAL)3);
     objptr->idunion.mca_p1l.f_m_l++;

     scp=0.;
     scp+=l1_1->phi[0]*l1_2->phi[0];
     scp+=2.*l1_1->phi[1]*l1_2->phi[1];
     scp+=2.*l1_1->phi[2]*l1_2->phi[2];
     *objptr->idunion.mca_p1l.s_m_l += scp/((_REAL)3);
     objptr->idunion.mca_p1l.s_m_l++;

    }
   }
   objptr->idunion.mca_p1l.s_m_l -= objptr->idunion.mca_p1l.dim_m_l;
   objptr->idunion.mca_p1l.f_m_l -= objptr->idunion.mca_p1l.dim_m_l;

   for (i1=0,i=0,kk=0,l1_1=l1l_0;i1<nbondl;i1++,l1_1++) {
    for (j1=i1,j=i,l1_2=l1_1;j1<nbondl;j1++,l1_2++,kk++) {

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
          l1_1->dphidm1[m1][k] *
          objptr->idunion.mca_p1l.Fb_l_l[i1][j1].f[k][l]  *
          l1_2->dphidm1[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     *objptr->idunion.mca_p1l.f_l_l += scp/((_REAL)3);
     objptr->idunion.mca_p1l.f_l_l++;

     scp=0.;
     scp+=l1_1->phi[0]*l1_2->phi[0];
     scp+=2.*l1_1->phi[1]*l1_2->phi[1];
     scp+=2.*l1_1->phi[2]*l1_2->phi[2];
     *objptr->idunion.mca_p1l.s_l_l += scp/((_REAL)3);
     objptr->idunion.mca_p1l.s_l_l++;

    }
   }
   objptr->idunion.mca_p1l.s_l_l -= objptr->idunion.mca_p1l.dim_l_l;
   objptr->idunion.mca_p1l.f_l_l -= objptr->idunion.mca_p1l.dim_l_l;


   if(objptr->idunion.mca_p1l.type>1) {
   p=objptr->idunion.mca_p1l.m;
   l123=objptr->idunion.mca_p1l.l123;
   for (i1=0,i=0,kk=0,p1=p;i1<nmode;i1++,p1++) {

   for (j1=0,j=0,l123_2=l123;j1<nmode;j1++) {
     for (j2=0;j2<nmode;j2++) {
      for (j3=j2;j3<nmode;j3++,j++,l123_2++,kk++) {

     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<3;m1++) {
       {
       m2=m1;
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          objptr->idunion.mca_p1l.Fm_m_m[i1][j1].f[m1][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          objptr->idunion.mca_p1l.Fm_m_m[i1][j2].f[m1][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          objptr->idunion.mca_p1l.Fm_m_m[i1][j3].f[m1][l]  *
          l123_2->dphidm3[m2][l];
        }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mca_p1l.f_m_m2 += scp/((_REAL)3);
     objptr->idunion.mca_p1l.f_m_m2++;

     scp=0.;
     scp+=p1->m[0]*l123_2->phi[0];
     scp+=p1->m[1]*l123_2->phi[1];
     scp+=p1->m[2]*l123_2->phi[2];
     *objptr->idunion.mca_p1l.s_m_m2 += scp/((_REAL)3);
     objptr->idunion.mca_p1l.s_m_m2++;

     }
    }
   }

   }
   objptr->idunion.mca_p1l.s_m_m2 -= objptr->idunion.mca_p1l.dim_m_m2;
   objptr->idunion.mca_p1l.f_m_m2 -= objptr->idunion.mca_p1l.dim_m_m2;

   p=objptr->idunion.mca_p1l.m;
   l123=objptr->idunion.mca_p1l.l123;
   for (i1=0,i=0,kk=0,l1=l1l_0;i1<nbondl;i1++,l1++) {

   for (j1=0,j=0,l123_2=l123;j1<nmode;j1++) {
     for (j2=0;j2<nmode;j2++) {
      for (j3=j2;j3<nmode;j3++,j++,l123_2++,kk++) {

     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<3;m1++) {
       {
       m2=m1;
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l1_1->dphidm1[m1][k] *
          objptr->idunion.mca_p1l.Fm_l_m[i1][j1].f[m1][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1_1->dphidm1[m1][k] *
          objptr->idunion.mca_p1l.Fm_l_m[i1][j2].f[m1][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1_1->dphidm1[m1][k] *
          objptr->idunion.mca_p1l.Fm_l_m[i1][j3].f[m1][l]  *
          l123_2->dphidm3[m2][l];
        }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mca_p1l.f_l_m2 += scp/((_REAL)3);
     objptr->idunion.mca_p1l.f_l_m2++;

     scp=0.;
     scp+=l1->phi[0]*l123_2->phi[0];
     scp+=l1->phi[1]*l123_2->phi[1];
     scp+=l1->phi[2]*l123_2->phi[2];
     *objptr->idunion.mca_p1l.s_l_m2 += scp/((_REAL)3);
     objptr->idunion.mca_p1l.s_l_m2++;

     }
    }
   }

   }
   objptr->idunion.mca_p1l.s_l_m2 -= objptr->idunion.mca_p1l.dim_l_m2;
   objptr->idunion.mca_p1l.f_l_m2 -= objptr->idunion.mca_p1l.dim_l_m2;



   l123=objptr->idunion.mca_p1l.l123;
   for (i1=0,i=0,kk=0,l123_1=l123;i1<nmode;i1++) {
    for (i2=0;i2<nmode;i2++) {
     for (i3=i2;i3<nmode;i3++,i++,l123_1++) {
     j1_0=i1;
     j2_0=i2;
     j3_0=i3;
   for (j1=j1_0,j=i,l123_2=l123_1;j1<nmode;j1++,j3_0=0,j2_0=0) {
    for (j2=j2_0;j2<nmode;j2++,j3_0=j2) {
     for (j3=j3_0;j3<nmode;j3++,j++,l123_2++,kk++) {

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
          objptr->idunion.mca_p1l.Fm_m_m[i1][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm1[m1][k] *
          objptr->idunion.mca_p1l.Fm_m_m[i1][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm1[m1][k] *
          objptr->idunion.mca_p1l.Fm_m_m[i1][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mca_p1l.Fm_m_m[i2][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mca_p1l.Fm_m_m[i2][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mca_p1l.Fm_m_m[i2][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mca_p1l.Fm_m_m[i3][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mca_p1l.Fm_m_m[i3][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mca_p1l.Fm_m_m[i3][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mca_p1l.f_m2_m2 += scp/((_REAL)3);
     objptr->idunion.mca_p1l.f_m2_m2++;

     scp=0.;
     scp+=l123_1->phi[0]*l123_2->phi[0];
     scp+=l123_1->phi[1]*l123_2->phi[1];
     scp+=l123_1->phi[2]*l123_2->phi[2];
     *objptr->idunion.mca_p1l.s_m2_m2 += scp/((_REAL)3);
     objptr->idunion.mca_p1l.s_m2_m2++;

     }
    }
   }
     }
    }
   }
   objptr->idunion.mca_p1l.s_m2_m2 -= objptr->idunion.mca_p1l.dim_m2_m2;
   objptr->idunion.mca_p1l.f_m2_m2 -= objptr->idunion.mca_p1l.dim_m2_m2;

   }

}

/************************************************************************/
void
write_first_p1l(idstruct *dptr)
{
	int	i,j,i2,j2,k,l,m;
	int	states = dptr->idunion.mca_p1l.cases;
	FILE 	*file = dptr->idunion.mca_p1l.file,*dump;
        idstruct *mptr,*bptrm;

	_REAL fact,tol,factns,sum,sum2;
        int rank,nbondm;

        int n0m,n0l,n123m,nn,nf;
        _REAL **S;
        _REAL **F,**Fm1,**tmp,csum,ctmp;
        _REAL *ui,*v,*w;
#ifdef IMSL
        d_complex *wi,*vi;
#else
        _REAL *vi;
#endif

        n0m=dptr->idunion.mca_p1l.n0m;
        n0l=dptr->idunion.mca_p1l.n0l;
        n123m=dptr->idunion.mca_p1l.n123m;
        mptr=dptr->idunion.mca_p1l.mca;
        bptrm=mptr->idunion.mcarig.bead;
	fact = 1./ ((_REAL)states);
        factns = bptrm->idunion.bead.factns;

        printf("write_first_p1l: number of states: %d\n",states);
        for (i=1,k=0; i<=n0m; i++) {
         for (j=i; j<=n0m; j++,k++) {
           *dptr->idunion.mca_p1l.s_m_m++ *= fact;
           *dptr->idunion.mca_p1l.f_m_m++ *= fact;
         }
        }
        dptr->idunion.mca_p1l.s_m_m -= dptr->idunion.mca_p1l.dim_m_m;
        dptr->idunion.mca_p1l.f_m_m -= dptr->idunion.mca_p1l.dim_m_m;

        for (i=1,k=0; i<=n0m; i++) {
         for (j=1; j<=n0l; j++,k++) {
           *dptr->idunion.mca_p1l.s_m_l++ *= fact;
           *dptr->idunion.mca_p1l.f_m_l++ *= fact;
         }
        }
        dptr->idunion.mca_p1l.s_m_l -= dptr->idunion.mca_p1l.dim_m_l;
        dptr->idunion.mca_p1l.f_m_l -= dptr->idunion.mca_p1l.dim_m_l;

        for (i=1,k=0; i<=n0l; i++) {
         for (j=i; j<=n0l; j++,k++) {
           *dptr->idunion.mca_p1l.s_l_l++ *= fact;
           *dptr->idunion.mca_p1l.f_l_l++ *= fact;
         }
        }
        dptr->idunion.mca_p1l.s_l_l -= dptr->idunion.mca_p1l.dim_l_l;
        dptr->idunion.mca_p1l.f_l_l -= dptr->idunion.mca_p1l.dim_l_l;

        if(dptr->idunion.mca_p1l.type>1) {

        for (i=1,k=0; i<=n0m; i++) {
         for (j=1; j<=n123m; j++,k++) {
           *dptr->idunion.mca_p1l.s_m_m2++ *= fact;
           *dptr->idunion.mca_p1l.f_m_m2++ *= fact;
         }
        }
        dptr->idunion.mca_p1l.s_m_m2 -= dptr->idunion.mca_p1l.dim_m_m2;
        dptr->idunion.mca_p1l.f_m_m2 -= dptr->idunion.mca_p1l.dim_m_m2;
        for (i=1,k=0; i<=n0l; i++) {
         for (j=1; j<=n123m; j++,k++) {
           *dptr->idunion.mca_p1l.s_l_m2++ *= fact;
           *dptr->idunion.mca_p1l.f_l_m2++ *= fact;
         }
        }
        dptr->idunion.mca_p1l.s_l_m2 -= dptr->idunion.mca_p1l.dim_l_m2;
        dptr->idunion.mca_p1l.f_l_m2 -= dptr->idunion.mca_p1l.dim_l_m2;
        for (i=1,k=0; i<=n123m; i++) {
         for (j=i; j<=n123m; j++,k++) {
           *dptr->idunion.mca_p1l.s_m2_m2++ *= fact;
           *dptr->idunion.mca_p1l.f_m2_m2++ *= fact;
         }
        }
        dptr->idunion.mca_p1l.s_m2_m2 -= dptr->idunion.mca_p1l.dim_m2_m2;
        dptr->idunion.mca_p1l.f_m2_m2 -= dptr->idunion.mca_p1l.dim_m2_m2;
        }

        if(dptr->idunion.mca_p1l.type==1) {
         /* nn = n0m+n0l; */
         nn = n0l;
        } else if(dptr->idunion.mca_p1l.type==2) {
         /* nn = n0m+n0l+n123m; */
         nn = n0l+n123m;
        } else {
         printf("write_first_p1l: type not found.\n");
         exit(1);
        }  
        printf("write_first_p1l: F/S dimension= %d\n",nn);

        F = matrix(1,nn,1,nn);
        S = matrix(1,nn,1,nn);
        printf("write_first_p1l: allocated F/S ...\n");

/*
        for(i=1,k=0;i<=n0m;i++) {
         for(j=i;j<=n0m;j++,k++) {
          F[i][j]=*dptr->idunion.mca_p1l.f_m_m++;
          F[j][i]=F[i][j];
          S[i][j]=*dptr->idunion.mca_p1l.s_m_m++;
          S[j][i]=S[i][j];
         }
        }
        dptr->idunion.mca_p1l.s_m_m -= dptr->idunion.mca_p1l.dim_m_m;
        dptr->idunion.mca_p1l.f_m_m -= dptr->idunion.mca_p1l.dim_m_m;

        for(i=1,k=0;i<=n0m;i++) {
         for(j=1;j<=n0l;j++,k++) {
          j2=n0m+j;
          F[i][j2]=*dptr->idunion.mca_p1l.f_m_l++;
          F[j2][i]=F[i][j2];
          S[i][j2]=*dptr->idunion.mca_p1l.s_m_l++;
          S[j2][i]=S[i][j2];
         }
        }
        dptr->idunion.mca_p1l.s_m_l -= dptr->idunion.mca_p1l.dim_m_l;
        dptr->idunion.mca_p1l.f_m_l -= dptr->idunion.mca_p1l.dim_m_l;
*/

        for(i=1,k=0;i<=n0l;i++) {
         for(j=i;j<=n0l;j++,k++) {
/*
          i2=n0m+i;
          j2=n0m+j;
*/
          i2=i;
          j2=j;
          F[i2][j2]=*dptr->idunion.mca_p1l.f_l_l++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*dptr->idunion.mca_p1l.s_l_l++;
          S[j2][i2]=S[i2][j2];
         }
        }
        dptr->idunion.mca_p1l.s_l_l -= dptr->idunion.mca_p1l.dim_l_l;
        dptr->idunion.mca_p1l.f_l_l -= dptr->idunion.mca_p1l.dim_l_l;

        if(dptr->idunion.mca_p1l.type>1) {

/*
        for(i=1,k=0;i<=n0m;i++) {
         for(j=1;j<=n123m;j++,k++) {
          i2=i;
          j2=n0m+n0l+j;
          F[i2][j2]=*dptr->idunion.mca_p1l.f_m_m2++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*dptr->idunion.mca_p1l.s_m_m2++;
          S[j2][i2]=S[i2][j2];
         }
        }
        dptr->idunion.mca_p1l.s_m_m2 -= dptr->idunion.mca_p1l.dim_m_m2;
        dptr->idunion.mca_p1l.f_m_m2 -= dptr->idunion.mca_p1l.dim_m_m2;
*/

        for(i=1,k=0;i<=n0l;i++) {
         for(j=1;j<=n123m;j++,k++) {
/*
          i2=n0m+i;
          j2=n0m+n0l+j;
*/
          i2=i;
          j2=n0l+j;
          F[i2][j2]=*dptr->idunion.mca_p1l.f_l_m2++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*dptr->idunion.mca_p1l.s_l_m2++;
          S[j2][i2]=S[i2][j2];
         }
        }
        dptr->idunion.mca_p1l.s_l_m2 -= dptr->idunion.mca_p1l.dim_l_m2;
        dptr->idunion.mca_p1l.f_l_m2 -= dptr->idunion.mca_p1l.dim_l_m2;

        for(i=1,k=0;i<=n123m;i++) {
         for(j=i;j<=n123m;j++,k++) {
/*
          i2=n0m+n0l+i;
          j2=n0m+n0l+j;
*/
          i2=n0l+i;
          j2=n0l+j;
          F[i2][j2]=*dptr->idunion.mca_p1l.f_m2_m2++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*dptr->idunion.mca_p1l.s_m2_m2++;
          S[j2][i2]=S[i2][j2];
         }
        }
        dptr->idunion.mca_p1l.s_m2_m2 -= dptr->idunion.mca_p1l.dim_m2_m2;
        dptr->idunion.mca_p1l.f_m2_m2 -= dptr->idunion.mca_p1l.dim_m2_m2;

        }

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
          ui[k]=S[i][j];
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
        rank=invert(S,nn,nn,Fm1);
#endif
        fprintf(file,"write_first_p1: rank of S= %d\n",rank); fflush(file);

        tmp = matrix(1,nn,1,nn);

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          csum=0.;
          for(k=1;k<=nn;k++) {
           csum+=Fm1[i][k]*F[k][j];
          }
          tmp[i][j]=csum;
         }
        }

        w=vector(1,nn);
#ifdef IMSL
        for(i=1,l=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,l++) {
          ui[l]=tmp[i][j];
         }
        }
        wi = imsl_d_eig_gen(nn,ui,IMSL_VECTORS,&vi,0);
        for(i=0;i<nn;i++) {
         w[i+1]=1./wi[i].re;
        }
        free(ui);
        free(wi);
#else
        vi=(_REAL *)get2(nn*nn,sizeof(_REAL));
#ifdef EISPACK
        eigvec_ev(nn,tmp,vi,w);
        for(i=0;i<nn;i++) {
         w[i+1]=1./w[i+1];
        }
#else
        eigvec_nr(tmp,nn,w,vi);
#endif
#endif

        fprintf(file,"Sm1*F eigenvalues:\n");
        nf=0;
        for(i=1;i<=nn;i++) {
         fprintf(file,"w(%d)= %e\n",i,w[i]);
         if(Nabs(w[i])>=ZERO) {
          nf++;
         }
        }
        fflush(file);

        dptr->idunion.mca_p1l.nf = nf;
        nf=nn;
        dptr->idunion.mca_p1l.eval= (_REAL *)get2(nf,sizeof(_REAL));
        dptr->idunion.mca_p1l.evec= matrix(1,nn,1,nf);

        nf=0;
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
#ifdef IMSL
          dptr->idunion.mca_p1l.evec[i][j]=vi[k].re;
#else
          dptr->idunion.mca_p1l.evec[i][j]=vi[k];
#endif
         }
        }
        eigsrt(w,dptr->idunion.mca_p1l.evec,nn);
        for(i=1;i<=nn;i++) {
         dptr->idunion.mca_p1l.eval[i-1]=1./w[i];
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
          dptr->idunion.mca_p1l.eval[nf-1]=1./wi[i-1].re;
          vt=&vi[k];
          for(j=1;j<=nn;j++) {
           dptr->idunion.mca_p1l.evec[j][nf]=(*vt).re;
           vt+=nn;
          }
         }
        }
*/

/*
*/
        fprintf(file,"evec:\n");
        for(j=1;j<=nn;j++) {
         for(i=1;i<=nn;i++) {
          fprintf(file,"%d %d %e\n",i,j,dptr->idunion.mca_p1l.evec[i][j]);
         }
        }

        nf=nn;
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mca_p1l.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mca_p1l.evec[k][i]*tmp[k][j];
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
      /* ortho(dptr->idunion.mca_p1l.eval,dptr->idunion.mca_p1l.evec,
            Fm1,nn); */

/*
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
          dptr->idunion.mca_p1l.evec[i][j]/=sum;
         }
        }
        if(k!=nf) {
         printf("write_first: missed zero in Ct*S*C.\n");
        }

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mca_p1l.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mca_p1l.evec[k][i]*tmp[k][j];
          }
          S[i][j]=sum;
         }
        }
/*
*/
        fprintf(file,"Ct*S1*C:\n");
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nf;j++) {
          fprintf(file,"%d %d %e\n",i,j,S[i][j]);
         }
        }
        fflush(file);

        free_matrix(F,1,nn,1,nn);
        free_matrix(S,1,nn,1,nn);
        free_matrix(tmp,1,nn,1,nn);
}
