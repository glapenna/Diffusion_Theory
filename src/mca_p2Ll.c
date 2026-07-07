#include "head.h"

// extern char tok[];
// char tok[];
// void *get2(int,int);

extern idstruct *startid(),*getstruct();
rotstruct **rot_matrix();

/************************************************************************/

void
init_mca_p2Ll()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(MCA_P2LL, id);

	dptr->idunion.mca_p2Ll.cases = 0;
	dptr->idunion.mca_p2Ll.dump = 1;
	dptr->idunion.mca_p2Ll.type = 1;
	dptr->idunion.mca_p2Ll.rigid = 0;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("MCA_P2LL: expected file name","");
	strcpy(dptr->idunion.mca_p2Ll.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.mca_p2Ll.fname, 0, &exists, "MCA_P2LL"))
		exit(1);
	if ((dptr->idunion.mca_p2Ll.file = 
			 fopen(dptr->idunion.mca_p2Ll.fname, "w")) == NULL){
		perror(dptr->idunion.mca_p2Ll.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************/
void
setup_mca_p2Ll()
{
	idstruct	*dptr,*mptrL,*mptrl,*bptrL,*bptrl;
	int i,j,k,icount,nn;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != MCA_P2LL)
                        inerr("setup_mca_p2Ll: expected id or keyword\n","");
        } else {
                printf("setup_mca_p2Ll: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((mptrL = getstruct(tok)) == NULL) 
                inerr("setup_mca_p2Ll: no such id\n","");
        if (mptrL->id_type != MCA_P2)
                inerr("setup_mca_p2Ll: mca_p2L id expected\n","");
        gettok();
        if ((mptrl = getstruct(tok)) == NULL) 
                inerr("setup_mca_p2Ll: no such id\n","");
        if (mptrl->id_type != MCA_P2)
                inerr("setup_mca_p2Ll: mca_p2l id expected\n","");

        gettok();
        for(;;) {
         if(tok[0]==';') {
          break;
         }
        }

        dptr->idunion.mca_p2Ll.mcaL = mptrL;
        dptr->idunion.mca_p2Ll.mcal = mptrl;

        dptr->idunion.mca_p2Ll.n1=
         mptrL->idunion.mca_p2.n1+mptrl->idunion.mca_p2.n1;
        dptr->idunion.mca_p2Ll.n2=
         mptrL->idunion.mca_p2.n2+mptrl->idunion.mca_p2.n2;

        dptr->idunion.mca_p2Ll.dim_L12_l12=
         mptrL->idunion.mca_p2.n1*mptrl->idunion.mca_p2.n1;
        dptr->idunion.mca_p2Ll.dim_L12_l1234=
         mptrL->idunion.mca_p2.n1*mptrl->idunion.mca_p2.n2;
        dptr->idunion.mca_p2Ll.dim_L1234_l12=
         mptrL->idunion.mca_p2.n2*mptrl->idunion.mca_p2.n1;
        dptr->idunion.mca_p2Ll.dim_L1234_l1234=
         mptrL->idunion.mca_p2.n2*mptrl->idunion.mca_p2.n2;

        dptr->idunion.mca_p2Ll.s_L12_l12 = (_REAL *)
          get2(dptr->idunion.mca_p2Ll.dim_L12_l12,sizeof(_REAL));
        dptr->idunion.mca_p2Ll.f_L12_l12 = (_REAL *)
          get2(dptr->idunion.mca_p2Ll.dim_L12_l12,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p2Ll.dim_L12_l12;i++) {
         *dptr->idunion.mca_p2Ll.s_L12_l12++ =0.;
         *dptr->idunion.mca_p2Ll.f_L12_l12++ =0.;
        }
        dptr->idunion.mca_p2Ll.s_L12_l12-=
         dptr->idunion.mca_p2Ll.dim_L12_l12;
        dptr->idunion.mca_p2Ll.f_L12_l12-=
         dptr->idunion.mca_p2Ll.dim_L12_l12;

        dptr->idunion.mca_p2Ll.s_L12_l1234 = (_REAL *)
          get2(dptr->idunion.mca_p2Ll.dim_L12_l1234,sizeof(_REAL));
        dptr->idunion.mca_p2Ll.f_L12_l1234 = (_REAL *)
          get2(dptr->idunion.mca_p2Ll.dim_L12_l1234,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p2Ll.dim_L12_l1234;i++) {
         *dptr->idunion.mca_p2Ll.s_L12_l1234++ =0.;
         *dptr->idunion.mca_p2Ll.f_L12_l1234++ =0.;
        }
        dptr->idunion.mca_p2Ll.s_L12_l1234-=
         dptr->idunion.mca_p2Ll.dim_L12_l1234;
        dptr->idunion.mca_p2Ll.f_L12_l1234-=
         dptr->idunion.mca_p2Ll.dim_L12_l1234;

        dptr->idunion.mca_p2Ll.s_L1234_l12 = (_REAL *)
          get2(dptr->idunion.mca_p2Ll.dim_L1234_l12,sizeof(_REAL));
        dptr->idunion.mca_p2Ll.f_L1234_l12 = (_REAL *)
          get2(dptr->idunion.mca_p2Ll.dim_L1234_l12,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p2Ll.dim_L1234_l12;i++) {
         *dptr->idunion.mca_p2Ll.s_L1234_l12++ =0.;
         *dptr->idunion.mca_p2Ll.f_L1234_l12++ =0.;
        }
        dptr->idunion.mca_p2Ll.s_L1234_l12-=
         dptr->idunion.mca_p2Ll.dim_L1234_l12;
        dptr->idunion.mca_p2Ll.f_L1234_l12-=
         dptr->idunion.mca_p2Ll.dim_L1234_l12;

        dptr->idunion.mca_p2Ll.s_L1234_l1234 = (_REAL *)
          get2(dptr->idunion.mca_p2Ll.dim_L1234_l1234,sizeof(_REAL));
        dptr->idunion.mca_p2Ll.f_L1234_l1234 = (_REAL *)
          get2(dptr->idunion.mca_p2Ll.dim_L1234_l1234,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p2Ll.dim_L1234_l1234;i++) {
         *dptr->idunion.mca_p2Ll.s_L1234_l1234++ =0.;
         *dptr->idunion.mca_p2Ll.f_L1234_l1234++ =0.;
        }
        dptr->idunion.mca_p2Ll.s_L1234_l1234-=
         dptr->idunion.mca_p2Ll.dim_L1234_l1234;
        dptr->idunion.mca_p2Ll.f_L1234_l1234-=
         dptr->idunion.mca_p2Ll.dim_L1234_l1234;

        bptrL=mptrL->idunion.mca_p2.bead;
        bptrl=mptrl->idunion.mca_p2.bead;

        dptr->idunion.mca_p2Ll.Fb_L_l = rot_matrix
         (0,bptrL->idunion.bead.nbond-1,0,bptrl->idunion.bead.nbond-1);

        nn=mptrL->idunion.mca_p2.n1+mptrl->idunion.mca_p2.n1;

        printf("setup_mca_p2Ll: boundaries of data.\n");
        printf("dim_L12_l12:      %d\n",
         dptr->idunion.mca_p2Ll.dim_L12_l12);
        printf("dim_L12_l1234:    %d\n",
         dptr->idunion.mca_p2Ll.dim_L12_l1234);
        printf("dim_L1234_l12:    %d\n",
         dptr->idunion.mca_p2Ll.dim_L1234_l12);
        printf("dim_L1234_l1234:  %d\n",
         dptr->idunion.mca_p2Ll.dim_L1234_l1234);

}

/************************************************************************/
do_first_p2Ll(idstruct *objptr,idstruct *serptr)
{
   _REAL *box=serptr->idunion.ser.set->boxc;
   idstruct *bptrL,*bptrl,*bptr,*mptrL,*mptrl;
   int i1,i1_0,i1_1,i2,i2_0,i2_1,i3,i3_0,i3_1,i4,i4_0,i4_1,
       j1_0,j2_0,j1,j2,j1_1,j2_1,j3,j3_0,j3_1,j4,j4_0,j4_1,
       i,j,k,l,ib,jb,ii,jj,jj0,kk,im,m1,m2,nbondL,nbondl,overlap;
   int *doffset_l, *doffset_L, *offset_l2_l, *doffset2_l,
       *offset_l2_L, *doffset2_L;
   _REAL **r1, **r2;
   _REAL *rt1,*rt2,r,zi,zj,zr,cutsq,
         boxm1[3],x1,y1,z1,x2,y2,z2,dx,dy,dz,scp;
   bondstruct *lb,*lL,*ll,*lb1,*lb2;
   mod12struct *l12_l,*l12_L,*l12_1,*l12_2;
   mod1234struct *l1234_l,*l1234_L,*l1234_1,*l1234_2;
   _REAL xm,ym,zm,f[5][5],Ct[3],csum,D[3][3],b,c,weight;

   mptrL=objptr->idunion.mca_p2Ll.mcaL;
   mptrl=objptr->idunion.mca_p2Ll.mcal;
   bptrL=mptrL->idunion.mca_p2.bead;
   bptrl=mptrl->idunion.mca_p2.bead;
   bptr=bptrL;
   doffset_L=bptrL->idunion.bead.l->doffset;
   doffset_l=bptrl->idunion.bead.l->doffset;
   offset_l2_L=bptrL->idunion.bead.l->offset_l2;
   doffset2_L=bptrL->idunion.bead.l->doffset2;
   offset_l2_l=bptrl->idunion.bead.l->offset_l2;
   doffset2_l=bptrl->idunion.bead.l->doffset2;
   zr=bptrL->idunion.bead.zr;
   cutsq=bptrL->idunion.bead.cutsq;
   nbondL=bptrL->idunion.bead.nbond;
   nbondl=bptrl->idunion.bead.nbond;
   if(serptr->idunion.ser.iweight)
   {
    weight=serptr->idunion.ser.w->idunion.serw.w;
   } else {
    weight=1.;
   }
   for(k=0;k<3;k++) {
    boxm1[k]=1./box[k];
   }
   b=l_sqrt(3./8.);
   c=0.5;

   ++objptr->idunion.mca_p2Ll.cases;
   objptr->idunion.mca_p2Ll.casesw+=weight;

   lL=bptrL->idunion.bead.l;
   ll=bptrl->idunion.bead.l;

   for (ib=0,lb1=lL;ib<nbondL;ib++,lb1++) {
    for (jb=0,lb2=ll;jb<nbondl;jb++,lb2++) {
#include "fbb.c"
     for(m1=0;m1<3;m1++) {
      for(m2=m1;m2<3;m2++) {
       f[m2][m1]=f[m1][m2];
      }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_p2Ll.Fb_L_l[ib][jb].f[m1][m2]=f[m1][m2];
       /* objptr->idunion.mca_p2Ll.Fb_L_l[jb][ib].f[m2][m1]=f[m1][m2]; */
      }
     }
    }
   }


   l12_L=mptrL->idunion.mca_p2.l12;
   i1_0=0;
   i1_1=nbondL;

   for (i1=i1_0,i=0,kk=0,l12_1=l12_L;i1<i1_1;i1++) {
    i2_0=i1;
    i2_1=i2_0+doffset_L[i1];
    for (i2=i2_0,ii=0;i2<i2_1;i2++,ii++,i++,l12_1++) {

     l12_l=mptrl->idunion.mca_p2.l12;
     j1_0=0;
     j1_1=nbondl;

/*
  printf("i1_0,i2_0,j1_0,j2_0= %d %d %d %d\n", 
   i1_0,i2_0,j1_0,j2_0);
  printf("i1_1,i2_1,j1_1,j2_1= %d %d %d %d\n", 
   i1_1,i2_1,j1_1,j2_1);
*/

   for (j1=j1_0,l12_2=l12_l;j1<j1_1;j1++) 
   {
    j2_0=j1;
    j2_1=j2_0+doffset_l[j1];
    
    for (j2=j2_0;j2<j2_1;j2++,l12_2++)
    {

/*
  printf("i1,i2,j1,j2= %d %d %d %d\n", 
   i1+1,i2+1,j1+1,j2+1);
  fflush(stdout);
*/

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
          objptr->idunion.mca_p2Ll.Fb_L_l[i1][j1].f[k][l]  *
          l12_2->dphidm1[m2][l];
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i1][j2].f[k][l]  *
          l12_2->dphidm2[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i2][j1].f[k][l]  *
          l12_2->dphidm1[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i2][j2].f[k][l]  *
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
     *objptr->idunion.mca_p2Ll.f_L12_l12 += weight*scp/5.;
     objptr->idunion.mca_p2Ll.f_L12_l12++;

     scp=0.;
     scp+=l12_1->phi[0]*l12_2->phi[0];
     scp+=2.*l12_1->phi[1]*l12_2->phi[1];
     scp+=2.*l12_1->phi[2]*l12_2->phi[2];
     scp+=2.*l12_1->phi[3]*l12_2->phi[3];
     scp+=2.*l12_1->phi[4]*l12_2->phi[4];
     *objptr->idunion.mca_p2Ll.s_L12_l12 += weight*scp/5.;
     objptr->idunion.mca_p2Ll.s_L12_l12++;

    }
   }
    }
   }
   objptr->idunion.mca_p2Ll.s_L12_l12 -= 
    objptr->idunion.mca_p2Ll.dim_L12_l12;
   objptr->idunion.mca_p2Ll.f_L12_l12 -= 
    objptr->idunion.mca_p2Ll.dim_L12_l12;

   if(mptrl->idunion.mca_p2.order==2)
   {

/* L12/l1234 */

   l12_L=mptrL->idunion.mca_p2.l12;
   i1_0=0;
   i1_1=nbondL;
   for (i1=i1_0,l12_1=l12_L;i1<i1_1;i1++) {
    i2_0=i1;
    i2_1=i2_0+doffset_L[i1];
    for (i2=i2_0;i2<i2_1;i2++,i++,l12_1++) {

     l1234_l=mptrl->idunion.mca_p2.l1234;
     j1_0=0;
     j1_1=nbondl;

     for (j1=j1_0,l1234_2=l1234_l;j1<j1_1;j1++) {
      j2_0=j1;
      j2_1=j2_0+doffset_l[j1];
      for (j2=j2_0;j2<j2_1;j2++) {
       j3_0=j1-offset_l2_l[j1];
       j3_1=j3_0+doffset2_l[j1];
       for (j3=j3_0;j3<j3_1;j3++) {
        j4_0=j3;
        j4_1=j3_0+doffset2_l[j1];
        for (j4=j4_0;j4<j4_1;j4++,l1234_2++) {

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
          objptr->idunion.mca_p2Ll.Fb_L_l[i1][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i1][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i1][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i1][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i2][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i2][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i2][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i2][j4].f[k][l]  *
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
     *objptr->idunion.mca_p2Ll.f_L12_l1234 += weight*scp/5.;
     objptr->idunion.mca_p2Ll.f_L12_l1234++;

     scp=0.;
     scp+=l12_1->phi[0]*l1234_2->phi[0];
     scp+=2.*l12_1->phi[1]*l1234_2->phi[1];
     scp+=2.*l12_1->phi[2]*l1234_2->phi[2];
     scp+=2.*l12_1->phi[3]*l1234_2->phi[3];
     scp+=2.*l12_1->phi[4]*l1234_2->phi[4];
     *objptr->idunion.mca_p2Ll.s_L12_l1234 += weight*scp/5.;
     objptr->idunion.mca_p2Ll.s_L12_l1234++;

        }
       }
      }
     }
    }
   }
   objptr->idunion.mca_p2Ll.s_L12_l1234 -= 
    objptr->idunion.mca_p2Ll.dim_L12_l1234;
   objptr->idunion.mca_p2Ll.f_L12_l1234 -= 
    objptr->idunion.mca_p2Ll.dim_L12_l1234;

    if(mptrL->idunion.mca_p2.order==2)
    {

/* L1234/l12 */

     l1234_L=mptrL->idunion.mca_p2.l1234;
     i1_0=0; 
     i1_1=nbondL;
     for (i1=i1_0,l1234_1=l1234_L;i1<i1_1;i1++) {
      i2_0=i1;
      i2_1=i2_0+doffset_L[i1];
      for (i2=i2_0;i2<i2_1;i2++) {
       i3_0=i1-offset_l2_L[i1];
       i3_1=i3_0+doffset2_L[i1];
       for (i3=i3_0;i3<i3_1;i3++) {
        i4_0=i3;
        i4_1=i3_0+doffset2_L[i1];
        for (i4=i4_0;i4<i4_1;i4++,l1234_1++) {

         l12_l=mptrl->idunion.mca_p2.l12;
         j1_0=0; 
         j1_1=nbondl;
         for (j1=j1_0,l12_2=l12_l;j1<j1_1;j1++) {
          j2_0=j1;
          j2_1=j2_0+doffset_l[j1];
          for (j2=j2_0;j2<j2_1;j2++,l12_2++) {

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
          objptr->idunion.mca_p2Ll.Fb_L_l[i1][j1].f[k][l]  *
          l12_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i1][j2].f[k][l]  *
          l12_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i2][j1].f[k][l]  *
          l12_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i2][j2].f[k][l]  *
          l12_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i3][j1].f[k][l]  *
          l12_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i3][j2].f[k][l]  *
          l12_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i4][j1].f[k][l]  *
          l12_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i4][j2].f[k][l]  *
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
     *objptr->idunion.mca_p2Ll.f_L1234_l12 += weight*scp/5.;
     objptr->idunion.mca_p2Ll.f_L1234_l12++;

     scp=0.;
     scp+=l1234_1->phi[0]*l12_2->phi[0];
     scp+=2.*l1234_1->phi[1]*l12_2->phi[1];
     scp+=2.*l1234_1->phi[2]*l12_2->phi[2];
     scp+=2.*l1234_1->phi[3]*l12_2->phi[3];
     scp+=2.*l1234_1->phi[4]*l12_2->phi[4];
     *objptr->idunion.mca_p2Ll.s_L1234_l12 += weight*scp/5.;
     objptr->idunion.mca_p2Ll.s_L1234_l12++;

        }
       }
      }
     }
    }
   }
   objptr->idunion.mca_p2Ll.s_L1234_l12 -= 
    objptr->idunion.mca_p2Ll.dim_L1234_l12;
   objptr->idunion.mca_p2Ll.f_L1234_l12 -= 
    objptr->idunion.mca_p2Ll.dim_L1234_l12;

/* L1234/l1234 */

     l1234_L=mptrL->idunion.mca_p2.l1234;
     i1_0=0; 
     i1_1=nbondL;
     for (i1=i1_0,l1234_1=l1234_L;i1<i1_1;i1++) {
      i2_0=i1;
      i2_1=i2_0+doffset_L[i1];
      for (i2=i2_0;i2<i2_1;i2++) {
       i3_0=i1-offset_l2_L[i1];
       i3_1=i3_0+doffset2_L[i1];
       for (i3=i3_0;i3<i3_1;i3++) {
        i4_0=i3;
        i4_1=i3_0+doffset2_L[i1];
        for (i4=i4_0;i4<i4_1;i4++,l1234_1++) {

     l1234_l=mptrl->idunion.mca_p2.l1234;
     j1_0=0;
     j1_1=nbondl;
     for (j1=j1_0,l1234_2=l1234_l;j1<j1_1;j1++) {
      j2_0=j1;
      j2_1=j2_0+doffset_l[j1];
      for (j2=j2_0;j2<j2_1;j2++) {
       j3_0=j1-offset_l2_l[j1];
       j3_1=j3_0+doffset2_l[j1];
       for (j3=j3_0;j3<j3_1;j3++) {
        j4_0=j3;
        j4_1=j3_0+doffset2_l[j1];
        for (j4=j4_0;j4<j4_1;j4++,l1234_2++) {

/*
  printf("i1_0,i2_0,j1_0,j2_0= %d %d %d %d\n", 
   i1_0,i2_0,j1_0,j2_0);
  printf("i1_1,i2_1,j1_1,j2_1= %d %d %d %d\n", 
   i1_1,i2_1,j1_1,j2_1);
*/

/*
  printf("i1,i2,i3,i4= %d %d %d %d,  ", 
   i1+1,i2+1,i3+1,i4+1);
  printf("j1,j2,j3,j4= %d %d %d %d\n", 
   j1+1,j2+1,j3+1,j4+1);
*/

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
          objptr->idunion.mca_p2Ll.Fb_L_l[i1][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i1][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i1][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i1][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i2][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i2][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i2][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i2][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i3][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i3][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i3][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i3][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i4][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i4][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i4][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2Ll.Fb_L_l[i4][j4].f[k][l]  *
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
     *objptr->idunion.mca_p2Ll.f_L1234_l1234 += weight*scp/5.;
     objptr->idunion.mca_p2Ll.f_L1234_l1234++;

     scp=0.;
     scp+=l1234_1->phi[0]*l1234_2->phi[0];
     scp+=2.*l1234_1->phi[1]*l1234_2->phi[1];
     scp+=2.*l1234_1->phi[2]*l1234_2->phi[2];
     scp+=2.*l1234_1->phi[3]*l1234_2->phi[3];
     scp+=2.*l1234_1->phi[4]*l1234_2->phi[4];
     *objptr->idunion.mca_p2Ll.s_L1234_l1234 += weight*scp/5.;
     objptr->idunion.mca_p2Ll.s_L1234_l1234++;

      }
     }
    }
   }
      }
     }
    }
   }
   objptr->idunion.mca_p2Ll.s_L1234_l1234 -= 
    objptr->idunion.mca_p2Ll.dim_L1234_l1234;
   objptr->idunion.mca_p2Ll.f_L1234_l1234 -= 
    objptr->idunion.mca_p2Ll.dim_L1234_l1234;

    }

   }

}

/************************************************************************/

void
write_p2Ll(idstruct *dptr)
{
	int	i,j,k,l,m,i2,j2;
	int	states = dptr->idunion.mca_p2Ll.cases;
	FILE 	*file = dptr->idunion.mca_p2Ll.file,*dump;
        idstruct *mptrL,*mptrl,*bptrL,*bptrl;

	_REAL fact;

        int n1_L,n2_L,n1_l,n2_l,nn,nf,rank;
        _REAL **F,**S,**Fm1,*w,**tmp;
        _REAL sum,sum2,tol;
        _REAL *u,*v;
#ifdef IMSL
        d_complex *vi,*wi;
#else
        _REAL *vi;
#endif

        fact = 1./ dptr->idunion.mca_p2Ll.casesw;
        mptrL=dptr->idunion.mca_p2Ll.mcaL;
        mptrl=dptr->idunion.mca_p2Ll.mcal;
        bptrL=mptrL->idunion.mca_p2.bead;
        bptrl=mptrl->idunion.mca_p2.bead;

        printf("write_p2Ll: number of states: %d\n",states);

        for (i=1; i<=dptr->idunion.mca_p2Ll.dim_L12_l12; i++) {
          *dptr->idunion.mca_p2Ll.s_L12_l12++ *= fact;
          *dptr->idunion.mca_p2Ll.f_L12_l12++ *= fact;
        }
        dptr->idunion.mca_p2Ll.s_L12_l12 -= 
         dptr->idunion.mca_p2Ll.dim_L12_l12;
        dptr->idunion.mca_p2Ll.f_L12_l12 -= 
         dptr->idunion.mca_p2Ll.dim_L12_l12;

        if(mptrl->idunion.mca_p2.order==2)
        {
         for (i=1; i<=dptr->idunion.mca_p2Ll.dim_L12_l1234; i++) {
          *dptr->idunion.mca_p2Ll.s_L12_l1234++ *= fact;
          *dptr->idunion.mca_p2Ll.f_L12_l1234++ *= fact;
         }
         dptr->idunion.mca_p2Ll.s_L12_l1234 -= 
          dptr->idunion.mca_p2Ll.dim_L12_l1234;
         dptr->idunion.mca_p2Ll.f_L12_l1234 -= 
          dptr->idunion.mca_p2Ll.dim_L12_l1234;

         if(mptrL->idunion.mca_p2.order==2)
         {
          for (i=1; i<=dptr->idunion.mca_p2Ll.dim_L1234_l12; i++) {
           *dptr->idunion.mca_p2Ll.s_L1234_l12++ *= fact;
           *dptr->idunion.mca_p2Ll.f_L1234_l12++ *= fact;
          }
          dptr->idunion.mca_p2Ll.s_L1234_l12 -= 
           dptr->idunion.mca_p2Ll.dim_L1234_l12;
          dptr->idunion.mca_p2Ll.f_L1234_l12 -= 
           dptr->idunion.mca_p2Ll.dim_L1234_l12;

          for (i=1; i<=dptr->idunion.mca_p2Ll.dim_L1234_l1234; i++) {
           *dptr->idunion.mca_p2Ll.s_L1234_l1234++ *= fact;
           *dptr->idunion.mca_p2Ll.f_L1234_l1234++ *= fact;
          }
          dptr->idunion.mca_p2Ll.s_L1234_l1234 -= 
           dptr->idunion.mca_p2Ll.dim_L1234_l1234;
          dptr->idunion.mca_p2Ll.f_L1234_l1234 -= 
           dptr->idunion.mca_p2Ll.dim_L1234_l1234;

         }
        }
 
        if(dptr->idunion.mca_p2Ll.dump==1) {
         dump=fopen("fsLl.dump","w");
         fwrite(dptr->idunion.mca_p2Ll.f_L12_l12,
           dptr->idunion.mca_p2Ll.dim_L12_l12*sizeof(_REAL),1,dump);
         fwrite(dptr->idunion.mca_p2Ll.s_L12_l12,
           dptr->idunion.mca_p2Ll.dim_L12_l12*sizeof(_REAL),1,dump);
         fwrite(dptr->idunion.mca_p2Ll.f_L12_l1234,
           dptr->idunion.mca_p2Ll.dim_L12_l1234*sizeof(_REAL),1,dump);
         fwrite(dptr->idunion.mca_p2Ll.s_L12_l1234,
           dptr->idunion.mca_p2Ll.dim_L12_l1234*sizeof(_REAL),1,dump);
         fwrite(dptr->idunion.mca_p2Ll.f_L1234_l12,
           dptr->idunion.mca_p2Ll.dim_L1234_l12*sizeof(_REAL),1,dump);
         fwrite(dptr->idunion.mca_p2Ll.s_L1234_l12,
           dptr->idunion.mca_p2Ll.dim_L1234_l12*sizeof(_REAL),1,dump);
         fwrite(dptr->idunion.mca_p2Ll.f_L1234_l1234,
           dptr->idunion.mca_p2Ll.dim_L1234_l1234*sizeof(_REAL),1,dump);
         fwrite(dptr->idunion.mca_p2Ll.s_L1234_l1234,
           dptr->idunion.mca_p2Ll.dim_L1234_l1234*sizeof(_REAL),1,dump);
         fflush(dump);
         fclose(dump);
        } else if(dptr->idunion.mca_p2Ll.dump==-1) {
         printf("write_p2Ll: not implemented yet.\n");
         exit(1);
        }
       /* if(dptr->idunion.mca_p2Ll.print_and_exit) return; */
         return;

/* invert F */
        n1_L = mptrL->idunion.mca_p2.n1;
        n2_L = mptrL->idunion.mca_p2.n2;
        n1_l = mptrl->idunion.mca_p2.n1;
        n2_l = mptrl->idunion.mca_p2.n2;
        nn = n1_L+n2_L+n1_l+n2_l;
        F = matrix(1,nn,1,nn);
        S = matrix(1,nn,1,nn);

        for(i=1;i<=n1_L;i++) {
         for(j=i;j<=n1_L;j++) {
          F[i][j]=*mptrL->idunion.mca_p2.f1++;
          F[j][i]=F[i][j];
          S[i][j]=*mptrL->idunion.mca_p2.s1++;
          S[j][i]=S[i][j];
         }
        }
        mptrL->idunion.mca_p2.s1 -= mptrL->idunion.mca_p2.dim1;
        mptrL->idunion.mca_p2.f1 -= mptrL->idunion.mca_p2.dim1;

        if(mptrL->idunion.mca_p2.order==2)
        {

        for(i=1;i<=n1_L;i++) {
         for(j=1;j<=n2_L;j++) {
          j2=n1_L+j;
          F[i][j2]=*mptrL->idunion.mca_p2.f12++;
          F[j2][i]=F[i][j2];
          S[i][j2]=*mptrL->idunion.mca_p2.s12++;
          S[j2][i]=S[i][j2];
         }
        }
        mptrL->idunion.mca_p2.s12 -= mptrL->idunion.mca_p2.dim12;
        mptrL->idunion.mca_p2.f12 -= mptrL->idunion.mca_p2.dim12;

        }

        for(i=1;i<=n1_L;i++) {
         for(j=1;j<=n1_l;j++) {
          j2=n1_L+n2_L+j;
          F[i][j2]=*dptr->idunion.mca_p2Ll.f_L12_l12++;
          F[j2][i]=F[i][j2];
          S[i][j2]=*dptr->idunion.mca_p2Ll.s_L12_l12++;
          S[j2][i]=S[i][j2];
         }
        }
        dptr->idunion.mca_p2Ll.s_L12_l12 -= 
         dptr->idunion.mca_p2Ll.dim_L12_l12;
        dptr->idunion.mca_p2Ll.f_L12_l12 -= 
         dptr->idunion.mca_p2Ll.dim_L12_l12;

        if(mptrl->idunion.mca_p2.order==2)
        {

        for(i=1;i<=n1_L;i++) {
         for(j=1;j<=n2_l;j++) {
          j2=n1_L+n2_L+n1_l+j;
          F[i][j2]=*dptr->idunion.mca_p2Ll.f_L12_l1234++;
          F[j2][i]=F[i][j2];
          S[i][j2]=*dptr->idunion.mca_p2Ll.s_L12_l1234++;
          S[j2][i]=S[i][j2];
         }
        }
        dptr->idunion.mca_p2Ll.s_L12_l1234 -= 
         dptr->idunion.mca_p2Ll.dim_L12_l1234;
        dptr->idunion.mca_p2Ll.f_L12_l1234 -= 
         dptr->idunion.mca_p2Ll.dim_L12_l1234;

        }

        if(mptrL->idunion.mca_p2.order==2)
        {

        for(i=1;i<=n2_L;i++) {
         for(j=i;j<=n2_L;j++) {
          i2=n1_L+i;
          j2=n1_L+j;
          F[i2][j2]=*mptrL->idunion.mca_p2.f2++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*mptrL->idunion.mca_p2.s2++;
          S[j2][i2]=S[i2][j2];
         }
        }
        mptrL->idunion.mca_p2.s2 -= mptrL->idunion.mca_p2.dim2;
        mptrL->idunion.mca_p2.f2 -= mptrL->idunion.mca_p2.dim2;

        for(i=1;i<=n2_L;i++) {
         for(j=1;j<=n1_l;j++) {
          i2=n1_L+i;
          j2=n1_L+n2_L+j;
          F[i2][j2]=*dptr->idunion.mca_p2Ll.f_L1234_l12++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*dptr->idunion.mca_p2Ll.s_L1234_l12++;
          S[j2][i2]=S[i2][j2];
         }
        }
        dptr->idunion.mca_p2Ll.s_L1234_l12 -= 
         dptr->idunion.mca_p2Ll.dim_L1234_l12;
        dptr->idunion.mca_p2Ll.f_L1234_l12 -= 
         dptr->idunion.mca_p2Ll.dim_L1234_l12;

        if(mptrl->idunion.mca_p2.order==2)
        {
        for(i=1;i<=n2_L;i++) {
         for(j=1;j<=n2_l;j++) {
          i2=n1_L+i;
          j2=n1_L+n2_L+n1_l+j;
          F[i2][j2]=*dptr->idunion.mca_p2Ll.f_L1234_l1234++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*dptr->idunion.mca_p2Ll.s_L1234_l1234++;
          S[j2][i2]=S[i2][j2];
         }
        }
        dptr->idunion.mca_p2Ll.s_L1234_l1234 -= 
         dptr->idunion.mca_p2Ll.dim_L1234_l1234;
        dptr->idunion.mca_p2Ll.f_L1234_l1234 -= 
         dptr->idunion.mca_p2Ll.dim_L1234_l1234;

        }
        }

        for(i=1;i<=n1_l;i++) {
         for(j=i;j<=n1_l;j++) {
          i2=n1_L+n2_L+i;
          j2=n1_L+n2_L+j;
          F[i2][j2]=*mptrl->idunion.mca_p2.f1++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*mptrl->idunion.mca_p2.s1++;
          S[j2][i2]=S[i2][j2];
         }
        }
        mptrl->idunion.mca_p2.s1 -= mptrl->idunion.mca_p2.dim1;
        mptrl->idunion.mca_p2.f1 -= mptrl->idunion.mca_p2.dim1;

        if(mptrl->idunion.mca_p2.order==2)
        {
        for(i=1;i<=n1_l;i++) {
         for(j=1;j<=n2_l;j++) {
          i2=n1_L+n2_L+i;
          j2=n1_L+n2_L+n1_l+j;
          F[i2][j2]=*mptrl->idunion.mca_p2.f12++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*mptrl->idunion.mca_p2.s12++;
          S[j2][i2]=S[i2][j2];
         }
        }
        mptrl->idunion.mca_p2.s12 -= mptrl->idunion.mca_p2.dim12;
        mptrl->idunion.mca_p2.f12 -= mptrl->idunion.mca_p2.dim12;

        for(i=1;i<=n2_l;i++) {
         for(j=i;j<=n2_l;j++) {
          i2=n1_L+n2_L+n1_l+i;
          j2=n1_L+n2_L+n1_l+j;
          F[i2][j2]=*mptrl->idunion.mca_p2.f2++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*mptrl->idunion.mca_p2.s2++;
          S[j2][i2]=S[i2][j2];
         }
        }
        mptrl->idunion.mca_p2.s2 -= mptrl->idunion.mca_p2.dim2;
        mptrl->idunion.mca_p2.f2 -= mptrl->idunion.mca_p2.dim2;

        }

/*
        fprintf(file,"S,F:\n");
        for(i=1;i<=nn;i++)
         for(j=1;j<=nn;j++)
          fprintf(file,"%d %d %e %e\n",i,j,S[i][j],F[i][j]);
        fflush(file);
*/
        dptr->idunion.mca_p2Ll.eval = (_REAL *)
         get2(nn,sizeof(_REAL));
        dptr->idunion.mca_p2Ll.evec = matrix( 1,nn, 1,nn);

        for(i=1;i<=nn;i++) {
         if(S[i][i]<ZERO) {
          printf("write_p2Ll:\n");
          printf("S element is wrong: S(%d,%d) = %e\n",i,i,S[i][i]);
         }
        }

        Fm1 = matrix(1,nn,1,nn);
#ifdef IMSL
        tol = TOLSVD;
        u = vector(0,nn*nn-1);
        if(dptr->idunion.mca_p2Ll.rigid==0) {
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
        if(dptr->idunion.mca_p2Ll.rigid==0) {
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
          printf("write_p1Ll:\n");
          printf("SF^-1S diagonal element is wrong: F(%d,%d) = %e\n",
           i,i,Fm1[i][i]);
         }
         if(Nabs(S[i][i])>ZERO) {
          Fm1[i][i] /= S[i][i];
         }
        }

        for(i=1;i<=nn;i++) {
         fprintf(file,"tau(%d)= %e ",i,1.*Fm1[i][i]);
         fprintf(file,"%e\n",bptrL->idunion.bead.factns*Fm1[i][i]);
        }
        fflush(file);

#ifdef IMSL
        tol = TOLSVD;
        u = vector(0,nn*nn-1);
        if(dptr->idunion.mca_p2Ll.rigid==0) {
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
        if(dptr->idunion.mca_p2Ll.rigid==0) {
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

        if(dptr->idunion.mca_p2Ll.rigid==0) {

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.,sum2=0.;k<=nn;k++) {
           sum+=Fm1[i][k]*F[k][j];
          }
          tmp[i][j]=sum;
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

/* diagonalize A = S-1*F or A = F^-1*S if rigidi */
        
        w=vector(1,nn);
#ifdef IMSL
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
          u[k]=tmp[i][j];
         }
        }
        wi = imsl_d_eig_gen(nn,u,IMSL_VECTORS,&vi,0);
        if(dptr->idunion.mca_p2Ll.rigid==0) {
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
        if(dptr->idunion.mca_p2Ll.rigid==0) {
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
          dptr->idunion.mca_p2Ll.evec[i][j] = vi[k].re;
#else
          dptr->idunion.mca_p2Ll.evec[i][j] = vi[k];
#endif
         }
        }
        eigsrt(w,dptr->idunion.mca_p2Ll.evec,nn);
        for(i=1;i<=nn;i++) {
         dptr->idunion.mca_p2Ll.eval[i-1] = 1./w[i];
        }
/*
        free_vector(w,1,nn);
        free(vi);
*/

/*
        fprintf(file,"evec:\n");
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          fprintf(file,"%d %d %e\n",i,j,dptr->idunion.mca_p2Ll.evec[i][j]);
         }
        }
        fflush(file);
*/

        if(dptr->idunion.mca_p2Ll.type==1) {
         nf=nn;
        } else {
         nf=dptr->idunion.mca_p2Ll.sel_modes;
        }
        if(dptr->idunion.mca_p2Ll.rigid==1) {
         if(dptr->idunion.mca_p2Ll.type==1) {
          nf=3;
         } else {
          nf=dptr->idunion.mca_p2Ll.sel_modes;
         }
        }

        dptr->idunion.mca_p2Ll.nf = nf;
/*
        dptr->idunion.mca_p2Ll.dim1=
         dptr->idunion.mca_p2Ll.nf*(dptr->idunion.mca_p2Ll.nf+1)/2;
*/

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mca_p2Ll.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mca_p2Ll.evec[k][i]*tmp[k][j];
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
          dptr->idunion.mca_p2Ll.evec[i][j]/=sum;
         }
        }
        if(k!=nf) {
         printf("write_p2Ll: missed zero in Ct*S*C.\n");
        }
/*
#ifdef IMSL

#endif

#ifdef EISPACK
        adjust_ev(nn,Fm1,tmp,dptr->idunion.mca_p2Ll.evec);
#endif
*/

/* transform F in Ct*F*C and S in Ct*S*C */
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=F[i][k]*dptr->idunion.mca_p2Ll.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mca_p2Ll.evec[k][i]*tmp[k][j];
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

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mca_p2Ll.evec[k][j];
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
           sum+=dptr->idunion.mca_p2Ll.evec[k][i]*tmp[k][j];
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


        free_matrix(F,1,nn,1,nn);
        free_matrix(S,1,nn,1,nn);
        free_matrix(Fm1,1,nn,1,nn);
        free_matrix(tmp,1,nn,1,nn);
}
