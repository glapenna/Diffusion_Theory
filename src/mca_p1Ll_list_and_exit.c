#include "head.h"

// extern char tok[];
// char tok[];
// void *get2(int,int);

extern idstruct *startid(),*getstruct();
rotstruct **rot_matrix();

/************************************************************************/

void
init_mca_p1Ll()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(MCA_P1LL, id);

	dptr->idunion.mca_p1Ll.cases = 0;
	dptr->idunion.mca_p1Ll.dump = 1;
	dptr->idunion.mca_p1Ll.type = 1;
	dptr->idunion.mca_p1Ll.rigid = 0;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("MCA_P1LL: expected file name","");
	strcpy(dptr->idunion.mca_p1Ll.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.mca_p1Ll.fname, 0, &exists, "MCA_P1LL"))
		exit(1);
	if ((dptr->idunion.mca_p1Ll.file = 
			 fopen(dptr->idunion.mca_p1Ll.fname, "w")) == NULL){
		perror(dptr->idunion.mca_p1Ll.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************/
void
setup_mca_p1Ll()
{
	idstruct	*dptr,*mptrL,*mptrl,*bptrL,*bptrl;
	int i,j,k,icount,nn;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != MCA_P1LL)
                        inerr("setup_mca_p1Ll: expected id or keyword\n","");
        } else {
                printf("setup_mca_p1Ll: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((mptrL = getstruct(tok)) == NULL) 
                inerr("setup_mca_p1Ll: no such id\n","");
        if (mptrL->id_type != MCARIG)
                inerr("setup_mca_p1Ll: mcarigL id expected\n","");
        gettok();
        if ((mptrl = getstruct(tok)) == NULL) 
                inerr("setup_mca_p1Ll: no such id\n","");
        if (mptrl->id_type != MCARIG)
                inerr("setup_mca_p1Ll: mcarigl id expected\n","");

        gettok();
        for(;;) {
         if(tok[0]==';') {
          break;
         } else if(!strcmp(tok,"RIGID")) {
          dptr->idunion.mca_p1Ll.rigid=1;
          gettok();
         }
        }

        dptr->idunion.mca_p1Ll.mcaL = mptrL;
        dptr->idunion.mca_p1Ll.mcal = mptrl;

        dptr->idunion.mca_p1Ll.dim_L_l=
         mptrL->idunion.mcarig.n1*mptrl->idunion.mcarig.n1;
        dptr->idunion.mca_p1Ll.dim_L_l123=
         mptrL->idunion.mcarig.n1*mptrl->idunion.mcarig.n2;

/*
        dptr->idunion.mca_p1Ll.s_L_l = (_REAL *)
          get2(dptr->idunion.mca_p1Ll.dim_L_l,sizeof(_REAL));
        dptr->idunion.mca_p1Ll.f_L_l = (_REAL *)
          get2(dptr->idunion.mca_p1Ll.dim_L_l,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p1Ll.dim_L_l;i++) {
         *dptr->idunion.mca_p1Ll.s_L_l++ =0.;
         *dptr->idunion.mca_p1Ll.f_L_l++ =0.;
        }
        dptr->idunion.mca_p1Ll.s_L_l-=dptr->idunion.mca_p1Ll.dim_L_l;
        dptr->idunion.mca_p1Ll.f_L_l-=dptr->idunion.mca_p1Ll.dim_L_l;
*/

        dptr->idunion.mca_p1Ll.s_L_l123 = (_REAL *)
          get2(dptr->idunion.mca_p1Ll.dim_L_l123,sizeof(_REAL));
        dptr->idunion.mca_p1Ll.f_L_l123 = (_REAL *)
          get2(dptr->idunion.mca_p1Ll.dim_L_l123,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p1Ll.dim_L_l123;i++) {
         *dptr->idunion.mca_p1Ll.s_L_l123++ =0.;
         *dptr->idunion.mca_p1Ll.f_L_l123++ =0.;
        }
        dptr->idunion.mca_p1Ll.s_L_l123-=
         dptr->idunion.mca_p1Ll.dim_L_l123;
        dptr->idunion.mca_p1Ll.f_L_l123-=
         dptr->idunion.mca_p1Ll.dim_L_l123;

        printf("setup_mca_p1Ll: boundaries of data.\n");

        bptrL=mptrL->idunion.mcarig.bead;
        bptrl=mptrl->idunion.mcarig.bead;

        dptr->idunion.mca_p1Ll.Fb_L_l = rot_matrix
         (0,bptrL->idunion.bead.nbond-1,0,bptrl->idunion.bead.nbond-1);
        nn=mptrL->idunion.mcarig.n1+mptrl->idunion.mcarig.n1+mptrl->idunion.mcarig.n2;
/*
        dptr->idunion.mca_p1Ll.eval = (_REAL *)
         get2(nn,sizeof(_REAL));
        dptr->idunion.mca_p1Ll.evec = matrix( 1,nn, 1,nn);
*/

        printf("dim_L_l:      %d\n",dptr->idunion.mca_p1Ll.dim_L_l);
        printf("dim_L_l123:   %d\n",dptr->idunion.mca_p1Ll.dim_L_l123);

}

/************************************************************************/
do_first_p1Ll(idstruct *objptr,idstruct *serptr)
{
   _REAL *box=serptr->idunion.ser.set->boxc;
   idstruct *bptrL,*bptrl,*bptr,*mptrL,*mptrl;
   int i1,j1_0,j2_0,j3_0,j1,j2,j3,j1_1,j2_1,j3_1,
       i,j,k,l,ib,jb,jj,kk,im,m1,m2,nbondL,nbondl,overlap;
   int *offset_l, *doffset;
   _REAL **r1, **r2;
   _REAL *rt1,*rt2,r,zi,zj,zr,cutsq,
         boxm1[3],x1,y1,z1,x2,y2,z2,dx,dy,dz,scp;
   bondstruct *lb,*lL,*ll,*lb1,*lb2;
   mod123struct *l123,*l123_2;
   _REAL xm,ym,zm,f[3][3],Ct[3],csum,D[3][3],b,c;

   mptrL=objptr->idunion.mca_p1Ll.mcaL;
   mptrl=objptr->idunion.mca_p1Ll.mcal;
   bptrL=mptrL->idunion.mcarig.bead;
   bptrl=mptrl->idunion.mcarig.bead;
   bptr=bptrL;
   offset_l=bptrl->idunion.bead.l->offset_l;
   doffset=bptrl->idunion.bead.l->doffset;
   zr=bptrL->idunion.bead.zr;
   cutsq=bptrL->idunion.bead.cutsq;
   nbondL=bptrL->idunion.bead.nbond;
   nbondl=bptrl->idunion.bead.nbond;
   for(k=0;k<3;k++) {
    boxm1[k]=1./box[k];
   }
   b=1./l_sqrt(2.);

   ++objptr->idunion.mca_p1Ll.cases;

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
       objptr->idunion.mca_p1Ll.Fb_L_l[ib][jb].f[m1][m2]=f[m1][m2];
      }
     }
/*
     scp= lb1->x*lb2->x+lb1->y*lb2->y+lb1->z*lb2->z;
     *objptr->idunion.mca_p1Ll.s_L_l++ += scp/3.;
     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mca_p1Ll.f_L_l++ += scp/3.;
*/
    }
   }
/*
   objptr->idunion.mca_p1Ll.s_L_l -= objptr->idunion.mca_p1Ll.dim_L_l;
   objptr->idunion.mca_p1Ll.f_L_l -= objptr->idunion.mca_p1Ll.dim_L_l;
*/

   if(mptrl->idunion.mcarig.sel_bonds==0) return;

   lb=bptrL->idunion.bead.l;
   l123=mptrl->idunion.mcarig.l123;
   for (i1=0,lb1=lb,i=0,kk=0;i1<nbondL;i1++,lb1++) {

   j1_0=0;
   j1_1=nbondl;
   for (j1=j1_0,j=0,kk=0,l123_2=l123;j1<j1_1;j1++) {
    j2_0=j1-offset_l[j1];
    j2_1=j2_0+doffset[j1];
    for (j2=j2_0,jj=0;j2<j2_1;j2++,jj++) {
     j3_0=j2;
     j3_1=j3_0+doffset[j1]-jj;
     for (j3=j3_0;j3<j3_1;j3++,j++,l123_2++) {

/*
 printf("i1,j1,j2,j3= %d %d %d %d\n",i1+1,j1+1,j2+1,j3+1);
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
          objptr->idunion.mca_p1Ll.Fb_L_l[i1][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          rdelta(m1,k)*
          objptr->idunion.mca_p1Ll.Fb_L_l[i1][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          rdelta(m1,k)*
          objptr->idunion.mca_p1Ll.Fb_L_l[i1][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mca_p1Ll.f_L_l123 += scp/((_REAL)3);
     objptr->idunion.mca_p1Ll.f_L_l123++;

     scp=0.;
     scp+=lb1->x*l123_2->phi[0];
     scp+=lb1->y*l123_2->phi[1];
     scp+=lb1->z*l123_2->phi[2];
     *objptr->idunion.mca_p1Ll.s_L_l123 += scp/((_REAL)3);
     objptr->idunion.mca_p1Ll.s_L_l123++;

     }
    }
   }

   }
   objptr->idunion.mca_p1Ll.s_L_l123 -= 
    objptr->idunion.mca_p1Ll.dim_L_l123;
   objptr->idunion.mca_p1Ll.f_L_l123 -= 
    objptr->idunion.mca_p1Ll.dim_L_l123;

}

/************************************************************************/

void
write_p1Ll(idstruct *dptr)
{
	int	i,j,k,l,m,i2,j2;
	int	states = dptr->idunion.mca_p1Ll.cases;
	FILE 	*file = dptr->idunion.mca_p1Ll.file,*dump;
        idstruct *mptrL,*mptrl,*bptrL,*bptrl;

	_REAL fact;

        int nL,n1,n2,nn,nf,rank;
        _REAL **F,**S,**Fm1,*w,**tmp;
        _REAL sum,sum2,tol;
        _REAL *u,*v;
#ifdef IMSL
        d_complex *vi,*wi;
#else
        _REAL *vi;
#endif

	fact = 1./ ((_REAL)states);
        mptrL=dptr->idunion.mca_p1Ll.mcaL;
        mptrl=dptr->idunion.mca_p1Ll.mcal;
        bptrL=mptrL->idunion.mcarig.bead;
        bptrl=mptrl->idunion.mcarig.bead;

        printf("write_p1Ll: number of states: %d\n",states);

/*
        for (i=1; i<=mptrL->idunion.mcarig.n1; i++) {
         for (j=1; j<=mptrl->idunion.mcarig.n1; j++) {
           *dptr->idunion.mca_p1Ll.s_L_l++ *= fact;
           *dptr->idunion.mca_p1Ll.f_L_l++ *= fact;
         }
        }
        dptr->idunion.mca_p1Ll.s_L_l -= dptr->idunion.mca_p1Ll.dim_L_l;
        dptr->idunion.mca_p1Ll.f_L_l -= dptr->idunion.mca_p1Ll.dim_L_l;
*/

        if(mptrl->idunion.mcarig.sel_bonds!=0)
        {

        for (i=1; i<=mptrL->idunion.mcarig.n1; i++) {
         for (j=1; j<=mptrl->idunion.mcarig.n2; j++) {
           *dptr->idunion.mca_p1Ll.s_L_l123++ *= fact;
           *dptr->idunion.mca_p1Ll.f_L_l123++ *= fact;
         }
        }
        dptr->idunion.mca_p1Ll.s_L_l123 -= dptr->idunion.mca_p1Ll.dim_L_l123;
        dptr->idunion.mca_p1Ll.f_L_l123 -= dptr->idunion.mca_p1Ll.dim_L_l123;

        }

        if(dptr->idunion.mca_p1Ll.dump==1) {
         dump=fopen("fsILl.dump","w");
/*
         fwrite(dptr->idunion.mca_p1Ll.f_L_l,dptr->idunion.mca_p1Ll.dim_L_l*sizeof(_REAL),
                1,dump);
         fwrite(dptr->idunion.mca_p1Ll.s_L_l,dptr->idunion.mca_p1Ll.dim_L_l*sizeof(_REAL),
                1,dump);
*/
         fwrite(dptr->idunion.mca_p1Ll.f_L_l123,dptr->idunion.mca_p1Ll.dim_L_l123*sizeof(_REAL),
                1,dump);
         fwrite(dptr->idunion.mca_p1Ll.s_L_l123,dptr->idunion.mca_p1Ll.dim_L_l123*sizeof(_REAL),
                1,dump);
         fflush(dump);
         fclose(dump);
        } else if(dptr->idunion.mca_p1Ll.dump==-1) {
         dump=fopen("fsILl.dump","r");
         fread(dptr->idunion.mca_p1Ll.f_L_l,dptr->idunion.mca_p1Ll.dim_L_l*sizeof(_REAL),
                1,dump);
         fread(dptr->idunion.mca_p1Ll.s_L_l,dptr->idunion.mca_p1Ll.dim_L_l*sizeof(_REAL),
                1,dump);
         fread(dptr->idunion.mca_p1Ll.f_L_l123,dptr->idunion.mca_p1Ll.dim_L_l123*sizeof(_REAL),
                1,dump);
         fread(dptr->idunion.mca_p1Ll.s_L_l123,dptr->idunion.mca_p1Ll.dim_L_l123*sizeof(_REAL),
                1,dump);
         fclose(dump);
        }
    return;


/* invert F */
        nL = mptrL->idunion.mcarig.n1;
        n1 = mptrl->idunion.mcarig.n1;
        n2 = mptrl->idunion.mcarig.n2;
        nn = nL+n1+n2;
        F = matrix(1,nn,1,nn);
        S = matrix(1,nn,1,nn);

        for(i=1;i<=nL;i++) {
         for(j=i;j<=nL;j++) {
          F[i][j]=(*mptrL->idunion.mcarig.f1++);
          F[j][i]=F[i][j];
          S[i][j]=*mptrL->idunion.mcarig.s1++;
          S[j][i]=S[i][j];
         }
        }
        mptrL->idunion.mcarig.s1 -= mptrL->idunion.mcarig.dim1;
        mptrL->idunion.mcarig.f1 -= mptrL->idunion.mcarig.dim1;

        for(i=1;i<=nL;i++) {
         for(j=1;j<=n1;j++) {
          j2=nL+j;
          F[i][j2]=*dptr->idunion.mca_p1Ll.f_L_l++;
          F[j2][i]=F[i][j2];
          S[i][j2]=*dptr->idunion.mca_p1Ll.s_L_l++;
          S[j2][i]=S[i][j2];
         }
        }
        dptr->idunion.mca_p1Ll.s_L_l -= dptr->idunion.mca_p1Ll.dim_L_l;
        dptr->idunion.mca_p1Ll.f_L_l -= dptr->idunion.mca_p1Ll.dim_L_l;

        if(mptrl->idunion.mcarig.sel_bonds!=0) {

        for(i=1;i<=nL;i++) {
         for(j=1;j<=n2;j++) {
          j2=nL+n1+j;
          F[i][j2]=*dptr->idunion.mca_p1Ll.f_L_l123++;
          F[j2][i]=F[i][j2];
          S[i][j2]=*dptr->idunion.mca_p1Ll.s_L_l123++;
          S[j2][i]=S[i][j2];
         }
        }
        dptr->idunion.mca_p1Ll.s_L_l123 -= dptr->idunion.mca_p1Ll.dim_L_l123;
        dptr->idunion.mca_p1Ll.f_L_l123 -= dptr->idunion.mca_p1Ll.dim_L_l123;

        }

        for(i=1;i<=n1;i++) {
         for(j=i;j<=n1;j++) {
          i2=nL+i;
          j2=nL+j;
          F[i2][j2]=*mptrl->idunion.mcarig.f1++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*mptrl->idunion.mcarig.s1++;
          S[j2][i2]=S[i2][j2];
         }
        }
        mptrl->idunion.mcarig.s1 -= mptrl->idunion.mcarig.dim1;
        mptrl->idunion.mcarig.f1 -= mptrl->idunion.mcarig.dim1;

        if(mptrl->idunion.mcarig.sel_bonds!=0) {

        for(i=1;i<=n1;i++) {
         for(j=1;j<=n2;j++) {
          i2=nL+i;
          j2=nL+n1+j;
          F[i2][j2]=*mptrl->idunion.mcarig.f12++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*mptrl->idunion.mcarig.s12++;
          S[j2][i2]=S[i2][j2];
         }
        }
        mptrl->idunion.mcarig.s12 -= mptrl->idunion.mcarig.dim12;
        mptrl->idunion.mcarig.f12 -= mptrl->idunion.mcarig.dim12;

        for(i=1;i<=n2;i++) {
         for(j=i;j<=n2;j++) {
          i2=nL+n1+i;
          j2=nL+n1+j;
          F[i2][j2]=*mptrl->idunion.mcarig.f2++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*mptrl->idunion.mcarig.s2++;
          S[j2][i2]=S[i2][j2];
         }
        }
        mptrl->idunion.mcarig.s2 -= mptrl->idunion.mcarig.dim2;
        mptrl->idunion.mcarig.f2 -= mptrl->idunion.mcarig.dim2;

        }

/*
        fprintf(file,"S,F:\n");
        for(i=1;i<=nn;i++)
         for(j=1;j<=nn;j++)
          fprintf(file,"%d %d %e %e\n",i,j,S[i][j],F[i][j]);
        fflush(file);
*/
        

        for(i=1;i<=nn;i++) {
         if(S[i][i]<ZERO) {
          printf("write_p1Ll:\n");
          printf("S element is wrong: S(%d,%d) = %e\n",i,i,S[i][i]);
         }
        }

        Fm1 = matrix(1,nn,1,nn);
#ifdef IMSL
        tol = TOLSVD;
        u = vector(0,nn*nn-1);
        if(dptr->idunion.mca_p1Ll.rigid==0) {
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
        if(dptr->idunion.mca_p1Ll.rigid==0) {
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
        if(dptr->idunion.mca_p1Ll.rigid==0) {
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
        if(dptr->idunion.mca_p1Ll.rigid==0) {
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

        if(dptr->idunion.mca_p1Ll.rigid==0) {

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
        if(dptr->idunion.mca_p1Ll.rigid==0) {
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
        if(dptr->idunion.mca_p1Ll.rigid==0) {
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
          dptr->idunion.mca_p1Ll.evec[i][j] = vi[k].re;
#else
          dptr->idunion.mca_p1Ll.evec[i][j] = vi[k];
#endif
         }
        }
        eigsrt(w,dptr->idunion.mca_p1Ll.evec,nn);
        for(i=1;i<=nn;i++) {
         /* dptr->idunion.mca_p1Ll.eval[i-1] = 1./w[i]; */
         dptr->idunion.mca_p1Ll.eval[i-1] = w[i];
        }
/*
        free_vector(w,1,nn);
        free(vi);
*/

/*
        fprintf(file,"evec:\n");
        for(i=1;i<=bptr->idunion.bead.nbond;i++) {
         for(j=1;j<=bptr->idunion.bead.nbond;j++) {
          fprintf(file,"%d %d %e\n",i,j,dptr->idunion.mca_p1Ll.evec[i][j]);
         }
        }
        fflush(file);
*/

        if(dptr->idunion.mca_p1Ll.type==1) {
         nf=nn;
        } else {
         nf=dptr->idunion.mca_p1Ll.sel_modes;
        }
        if(dptr->idunion.mca_p1Ll.rigid==1) {
         if(dptr->idunion.mca_p1Ll.type==1) {
          nf=3;
         } else {
          nf=dptr->idunion.mca_p1Ll.sel_modes;
         }
        }

        dptr->idunion.mca_p1Ll.nf = nf;
/*
        dptr->idunion.mca_p1Ll.dim1=
         dptr->idunion.mca_p1Ll.nf*(dptr->idunion.mca_p1Ll.nf+1)/2;
*/

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mca_p1Ll.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mca_p1Ll.evec[k][i]*tmp[k][j];
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
          dptr->idunion.mca_p1Ll.evec[i][j]/=sum;
         }
        }
        if(k!=nf) {
         printf("write_p1Ll: missed zero in Ct*S*C.\n");
        }
/*
#ifdef IMSL

#endif

#ifdef EISPACK
        adjust_ev(nn,Fm1,tmp,dptr->idunion.mca_p1Ll.evec);
#endif
*/

/* transform F in Ct*F*C and S in Ct*S*C */
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=F[i][k]*dptr->idunion.mca_p1Ll.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mca_p1Ll.evec[k][i]*tmp[k][j];
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
           sum+=S[i][k]*dptr->idunion.mca_p1Ll.evec[k][j];
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
           sum+=dptr->idunion.mca_p1Ll.evec[k][i]*tmp[k][j];
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
