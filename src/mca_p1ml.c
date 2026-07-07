#include "head.h"

// extern char tok[];
// char tok[];
// void *get2(int,int);

extern idstruct *startid(),*getstruct();
rotstruct **rot_matrix();

/************************************************************************/

void
init_mca_p1ml()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(MCA_P1ML, id);

	dptr->idunion.mca_p1ml.cases = 0;
	dptr->idunion.mca_p1ml.dump = 1;
	dptr->idunion.mca_p1ml.type = 1;
	dptr->idunion.mca_p1ml.rigid = 0;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("MCA_P1ML: expected file name","");
	strcpy(dptr->idunion.mca_p1ml.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.mca_p1ml.fname, 0, &exists, "MCA_P1ML"))
		exit(1);
	if ((dptr->idunion.mca_p1ml.file = 
			 fopen(dptr->idunion.mca_p1ml.fname, "w")) == NULL){
		perror(dptr->idunion.mca_p1ml.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************/
void
setup_mca_p1ml()
{
	idstruct	*dptr,*mptrL,*mptrl,*bptrL,*bptrl;
	int i,j,k,icount,nn;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != MCA_P1ML)
                        inerr("setup_mca_p1ml: expected id or keyword\n","");
        } else {
                printf("setup_mca_p1ml: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((mptrl = getstruct(tok)) == NULL) 
                inerr("setup_mca_p1ml: no such id\n","");
        if (mptrl->id_type != MCA_P1)
                inerr("setup_mca_p1ml: mca_p1l id expected\n","");
        gettok();
        if ((mptrL = getstruct(tok)) == NULL) 
                inerr("setup_mca_p1ml: no such id\n","");
        if (mptrL->id_type != MCA_P1)
                inerr("setup_mca_p1ml: mca_p1L id expected\n","");

        gettok();
        for(;;) {
         if(tok[0]==';') {
          break;
         } else if(!strcmp(tok,"RIGID")) {
          dptr->idunion.mca_p1ml.rigid=1;
          gettok();
         }
        }

        dptr->idunion.mca_p1ml.mcaL = mptrL;
        dptr->idunion.mca_p1ml.mcal = mptrl;

        printf("n1:              %d\n",mptrl->idunion.mca_p1.n1);
        printf("n123:            %d\n",mptrl->idunion.mca_p1.n123);
        printf("m1:              %d\n",mptrL->idunion.mca_p1.m1);
        printf("m123:            %d\n",mptrL->idunion.mca_p1.m123);

        dptr->idunion.mca_p1ml.dim_l_L123=
         mptrl->idunion.mca_p1.n1*mptrL->idunion.mca_p1.m123;
        dptr->idunion.mca_p1ml.dim_l123_L123=
         mptrl->idunion.mca_p1.n123*mptrL->idunion.mca_p1.m123;

        dptr->idunion.mca_p1ml.s_l_L123 = (_REAL *)
          get2(dptr->idunion.mca_p1ml.dim_l_L123,sizeof(_REAL));
        dptr->idunion.mca_p1ml.f_l_L123 = (_REAL *)
          get2(dptr->idunion.mca_p1ml.dim_l_L123,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p1ml.dim_l_L123;i++) {
         *dptr->idunion.mca_p1ml.s_l_L123++ =0.;
         *dptr->idunion.mca_p1ml.f_l_L123++ =0.;
        }
        dptr->idunion.mca_p1ml.s_l_L123-=
         dptr->idunion.mca_p1ml.dim_l_L123;
        dptr->idunion.mca_p1ml.f_l_L123-=
         dptr->idunion.mca_p1ml.dim_l_L123;

        dptr->idunion.mca_p1ml.s_l123_L123 = (_REAL *)
          get2(dptr->idunion.mca_p1ml.dim_l123_L123,sizeof(_REAL));
        dptr->idunion.mca_p1ml.f_l123_L123 = (_REAL *)
          get2(dptr->idunion.mca_p1ml.dim_l123_L123,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p1ml.dim_l123_L123;i++) {
         *dptr->idunion.mca_p1ml.s_l123_L123++ =0.;
         *dptr->idunion.mca_p1ml.f_l123_L123++ =0.;
        }
        dptr->idunion.mca_p1ml.s_l123_L123-=
         dptr->idunion.mca_p1ml.dim_l123_L123;
        dptr->idunion.mca_p1ml.f_l123_L123-=
         dptr->idunion.mca_p1ml.dim_l123_L123;

        printf("setup_mca_p1ml: boundaries of data.\n");

        bptrL=mptrL->idunion.mca_p1.bead;
        bptrl=mptrl->idunion.mca_p1.bead;

        dptr->idunion.mca_p1ml.Fb_l_L = rot_matrix
         (0,bptrl->idunion.bead.nbond-1,0,bptrL->idunion.bead.nbond-1);
        dptr->idunion.mca_p1ml.Fbm = rot_matrix
         (0,bptrl->idunion.bead.nbond-1,0,mptrL->idunion.mca_p1.m1-1);

        printf("dim_l_L123:      %d\n",dptr->idunion.mca_p1ml.dim_l_L123);
        printf("dim_l123_L123:   %d\n",dptr->idunion.mca_p1ml.dim_l123_L123);

}

/************************************************************************/
do_first_p1ml(idstruct *objptr,idstruct *serptr)
{
   _REAL *box=serptr->idunion.ser.set->boxc;
   idstruct *bptrL,*bptrl,*bptr,*mptrL,*mptrl;
   int i1,i2,i3,i1_0,i2_0,i3_0,i1_1,i2_1,i3_1,
       j1,j2,j3,j1_0,j2_0,j3_0,j1_1,j2_1,j3_1,
       i,j,k,l,ib,jb,ii,jj,kk,im,m1,m2,nbondL,nbondl,overlap,nmode,
       nlistL,nlistl;
   int *offset_l_l, *doffset_l;
   _REAL **r1, **r2;
   _REAL *rt1,*rt2,r,zi,zj,zr,cutsq,
         boxm1[3],x1,y1,z1,x2,y2,z2,dx,dy,dz,scp;
   bondstruct *lb,*lb_1,*lb_2,*lL,*ll,*lb1,*lb2;
   mod123struct *l123,*l123_1,*l123_2,*L123,*L123_1,*L123_2;
   _REAL xm,ym,zm,f[3][3],Ct[3],csum,D[3][3],b,c,weight;

   mptrL=objptr->idunion.mca_p1ml.mcaL;
   mptrl=objptr->idunion.mca_p1ml.mcal;
   bptrL=mptrL->idunion.mca_p1.bead;
   bptrl=mptrl->idunion.mca_p1.bead;
   bptr=bptrL;
   offset_l_l=bptrl->idunion.bead.l->offset_l;
   doffset_l=bptrl->idunion.bead.l->doffset;
   zr=bptrL->idunion.bead.zr;
   cutsq=bptrL->idunion.bead.cutsq;
   nbondL=bptrL->idunion.bead.nbond;
   nmode=mptrL->idunion.mca_p1.m1;
   nbondl=bptrl->idunion.bead.nbond;
   for(k=0;k<3;k++) {
    boxm1[k]=1./box[k];
   }
   b=1./l_sqrt(2.);

   ++objptr->idunion.mca_p1ml.cases;

   lL=bptrL->idunion.bead.l;
   ll=bptrl->idunion.bead.l;

   for (ib=0,lb1=ll;ib<nbondl;ib++,lb1++) {
    for (jb=0,lb2=lL;jb<nbondL;jb++,lb2++) {
#include "fbb.c"
     for(m1=0;m1<3;m1++) {
      for(m2=m1;m2<3;m2++) {
       f[m2][m1]=f[m1][m2];
      }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_p1ml.Fb_l_L[ib][jb].f[m1][m2]=f[m1][m2];
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
     for (jb=0;jb<nbondl;jb++) {
       for(m1=0;m1<3;m1++) {
        for(m2=0;m2<3;m2++) {
         f[m1][m2]+=
          objptr->idunion.mca_p1ml.Fb_l_L[i][jb].f[m1][m2] *
          mptrL->idunion.mca_p1.comb[jb+1][j+1];
        }
       }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_p1ml.Fbm[i][j].f[m1][m2]=f[m1][m2];
      }
     }
    }
   }

   lb=bptrl->idunion.bead.l;
   L123=mptrL->idunion.mca_p1.L123;
   for (i1=0,lb1=lb,i=0,kk=0;i1<nbondl;i1++,lb1++) {

   for (j1=0,j=0,L123_2=L123;j1<nmode;j1++) {
     for (j2=0;j2<nmode;j2++) {
      for (j3=j2;j3<nmode;j3++,j++,L123_2++,kk++) {

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
          objptr->idunion.mca_p1ml.Fbm[i1][j1].f[m1][l]  *
          L123_2->dphidm1[m2][l];
         f[m1][m2]+=
          objptr->idunion.mca_p1ml.Fbm[i1][j2].f[m1][l]  *
          L123_2->dphidm2[m2][l];
         f[m1][m2]+=
          objptr->idunion.mca_p1ml.Fbm[i1][j3].f[m1][l]  *
          L123_2->dphidm3[m2][l];
        }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mca_p1ml.f_l_L123 += scp/((_REAL)3);
     objptr->idunion.mca_p1ml.f_l_L123++;

     scp=0.;
     scp+=lb1->x*L123_2->phi[0];
     scp+=lb1->y*L123_2->phi[1];
     scp+=lb1->z*L123_2->phi[2];
     *objptr->idunion.mca_p1ml.s_l_L123 += scp/((_REAL)3);
     objptr->idunion.mca_p1ml.s_l_L123++;

     }
    }
   }

   }
   objptr->idunion.mca_p1ml.s_l_L123 -= objptr->idunion.mca_p1ml.dim_l_L123;
   objptr->idunion.mca_p1ml.f_l_L123 -= objptr->idunion.mca_p1ml.dim_l_L123;

   if(mptrl->idunion.mca_p1.sel_bonds)
   {

   l123=mptrl->idunion.mca_p1.l123;
   L123=mptrL->idunion.mca_p1.L123;
   i1_0=0;
   i1_1=nbondl;
   for (i1=i1_0,i=0,kk=0,l123_1=l123;i1<i1_1;i1++) {
    i2_0=i1-offset_l_l[i1];
    i2_1=i2_0+doffset_l[i1];
    for (i2=i2_0,ii=0;i2<i2_1;i2++,ii++) {
     i3_0=i2;
     i3_1=i3_0+doffset_l[i1]-ii;
     for (i3=i3_0;i3<i3_1;i3++,i++,l123_1++) {

   for (j1=0,j=0,L123_2=L123;j1<nmode;j1++) {
     for (j2=0;j2<nmode;j2++) {
      for (j3=j2;j3<nmode;j3++,j++,L123_2++,kk++) {

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
          objptr->idunion.mca_p1ml.Fbm[i1][j1].f[k][l]  *
          L123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm1[m1][k] *
          objptr->idunion.mca_p1ml.Fbm[i1][j2].f[k][l]  *
          L123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm1[m1][k] *
          objptr->idunion.mca_p1ml.Fbm[i1][j3].f[k][l]  *
          L123_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mca_p1ml.Fbm[i2][j1].f[k][l]  *
          L123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mca_p1ml.Fbm[i2][j2].f[k][l]  *
          L123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mca_p1ml.Fbm[i2][j3].f[k][l]  *
          L123_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mca_p1ml.Fbm[i3][j1].f[k][l]  *
          L123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mca_p1ml.Fbm[i3][j2].f[k][l]  *
          L123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mca_p1ml.Fbm[i3][j3].f[k][l]  *
          L123_2->dphidm3[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mca_p1ml.f_l123_L123 += scp/((_REAL)3);
     objptr->idunion.mca_p1ml.f_l123_L123++;

     scp=0.;
     scp+=l123_1->phi[0]*L123_2->phi[0];
     scp+=l123_1->phi[1]*L123_2->phi[1];
     scp+=l123_1->phi[2]*L123_2->phi[2];
     *objptr->idunion.mca_p1ml.s_l123_L123 += scp/((_REAL)3);
     objptr->idunion.mca_p1ml.s_l123_L123++;

     }
    }
   }

     }
    }
   }

   objptr->idunion.mca_p1ml.s_l123_L123 -= 
    objptr->idunion.mca_p1ml.dim_l123_L123;
   objptr->idunion.mca_p1ml.f_l123_L123 -= 
    objptr->idunion.mca_p1ml.dim_l123_L123;

   }

}

/************************************************************************/

void
write_p1ml(idstruct *dptr)
{
	int	i,j,k,l,m,i2,j2;
	int	states = dptr->idunion.mca_p1ml.cases;
	FILE 	*file = dptr->idunion.mca_p1ml.file,*dump;
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

	fact = 1./ ((_REAL)states);
        mptrL=dptr->idunion.mca_p1ml.mcaL;
        mptrl=dptr->idunion.mca_p1ml.mcal;
        bptrL=mptrL->idunion.mca_p1.bead;
        bptrl=mptrl->idunion.mca_p1.bead;

        printf("write_p1ml: number of states: %d\n",states);

        for (i=1; i<=mptrl->idunion.mca_p1.n1; i++) {
         for (j=1; j<=mptrL->idunion.mca_p1.m123; j++) {
           *dptr->idunion.mca_p1ml.s_l_L123++ *= fact;
           *dptr->idunion.mca_p1ml.f_l_L123++ *= fact;
         }
        }
        dptr->idunion.mca_p1ml.s_l_L123 -= dptr->idunion.mca_p1ml.dim_l_L123;
        dptr->idunion.mca_p1ml.f_l_L123 -= dptr->idunion.mca_p1ml.dim_l_L123;

        if(mptrl->idunion.mca_p1.sel_bonds!=0)
        {

        for (i=1; i<=mptrl->idunion.mca_p1.n123; i++) {
         for (j=1; j<=mptrL->idunion.mca_p1.m123; j++) {
           *dptr->idunion.mca_p1ml.s_l123_L123++ *= fact;
           *dptr->idunion.mca_p1ml.f_l123_L123++ *= fact;
         }
        }
        dptr->idunion.mca_p1ml.s_l123_L123 -= 
         dptr->idunion.mca_p1ml.dim_l123_L123;
        dptr->idunion.mca_p1ml.f_l123_L123 -= 
         dptr->idunion.mca_p1ml.dim_l123_L123;

        }

        if(dptr->idunion.mca_p1ml.dump==1) {
         dump=fopen("fslm.dump","w");
         fwrite(dptr->idunion.mca_p1ml.f_l_L123,
          dptr->idunion.mca_p1ml.dim_l_L123*sizeof(_REAL),1,dump);
         fwrite(dptr->idunion.mca_p1ml.s_l_L123,
          dptr->idunion.mca_p1ml.dim_l_L123*sizeof(_REAL),1,dump);
         fwrite(dptr->idunion.mca_p1ml.f_l123_L123,
          dptr->idunion.mca_p1ml.dim_l123_L123*sizeof(_REAL),1,dump);
         fwrite(dptr->idunion.mca_p1ml.s_l123_L123,
          dptr->idunion.mca_p1ml.dim_l123_L123*sizeof(_REAL),1,dump);
         fflush(dump);
         fclose(dump);
        } else if(dptr->idunion.mca_p1ml.dump==-1) {
         dump=fopen("fsLl.dump","r");
         fread(dptr->idunion.mca_p1ml.f_l_L123,
          dptr->idunion.mca_p1ml.dim_l_L123*sizeof(_REAL),1,dump);
         fread(dptr->idunion.mca_p1ml.s_l_L123,
          dptr->idunion.mca_p1ml.dim_l_L123*sizeof(_REAL),1,dump);
         fread(dptr->idunion.mca_p1ml.f_l123_L123,
          dptr->idunion.mca_p1ml.dim_l123_L123*sizeof(_REAL),1,dump);
         fread(dptr->idunion.mca_p1ml.s_l123_L123,
          dptr->idunion.mca_p1ml.dim_l123_L123*sizeof(_REAL),1,dump);
         fclose(dump);
        }

}
