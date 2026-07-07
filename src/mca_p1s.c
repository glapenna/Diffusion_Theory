#include "head.h"

// void *get2(int,int);


/************************************************************************/
do_first_p1s(idstruct *objptr,idstruct *serptr)
{
   _REAL *box=serptr->idunion.ser.set->boxc;
   idstruct *bptr,*mptr;
   int i1,i2,i3,j1_0,j2_0,j3_0,j1,j2,j3,kp1,kp2,kp3,
       i,j,k,l,ib,jb,kk,im,m1,m2,nbondm,nmode,nmode1,nmode2,overlap;
   _REAL **r1, **r2;
   _REAL *rt1,*rt2,r,zi,zj,zr,cutsq,
         boxm1[3],x1,y1,z1,x2,y2,z2,dx,dy,dz,scp;
   bondstruct *lb,*lb1,*lb2;
   modstruct *p,*p1,*p2,*p3;
   mod123struct *l123,*l123_1,*l123_2;
   _REAL xm,ym,zm,f[3][3],Ct[3],csum,D[3][3],rd1;

   mptr=objptr->idunion.mca_p1.mca;
   bptr=mptr->idunion.mcarig.bead;
   zr=bptr->idunion.bead.zr;
   cutsq=bptr->idunion.bead.cutsq;
   nmode=objptr->idunion.mca_p1.n0;
   nmode1=objptr->idunion.mca_p1.m1;
   nmode2=objptr->idunion.mca_p1.m2;
   nbondm=bptr->idunion.bead.nbond;
   for(k=0;k<3;k++) {
    boxm1[k]=1./box[k];
   }

   ++objptr->idunion.mca_p1.cases;

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

   p=objptr->idunion.mca_p1.m;
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

   if(objptr->idunion.mca_p1.type>1) {
   l123=objptr->idunion.mca_p1.l123;
   p=objptr->idunion.mca_p1.m;
   for (i1=0;i1<nmode1;i1++) {
    for (i2=0;i2<nmode1;i2++) {
     for (i3=i2;i3<nmode1;i3++) {
      ip1=sorted123[i1][1];
      ip2=sorted123[i2][2];
      ip3=sorted123[i3][3];
      scp=p[ip2].m[0]*p[ip3].m[0]+p[ip2].m[1]*p[ip3].m[1]+p[ip2].m[2]*p[ip3].m[2];
      for(k=0;k<3;k++) {
       l123->phi[k]=p[ip1].m[k]*scp;
       for(kk=0;kk<3;kk++) {
        l123->dphidm1[k][kk]=rdelta(k,kk)*scp;
        l123->dphidm2[k][kk]=p[ip1].m[k]*p[ip3].m[kk];
        l123->dphidm3[k][kk]=p[ip1].m[k]*p[ip2].m[kk];
       }
      }
      l123++;
     }
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
       objptr->idunion.mca_p1.Fb[ib][jb].f[m1][m2]=f[m1][m2];
       objptr->idunion.mca_p1.Fb[jb][ib].f[m2][m1]=f[m1][m2];
      }
     }
    }
   }

   for (i=0;i<nmode;i++) {
    ip1=sorted[i];
    for (j=0;j<nmode;j++) {
     jp1=sorted[j];
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
          mptr->idunion.mcarig.evec[ib+1][ip1+1] *
          objptr->idunion.mca_p1.Fb[ib][jb].f[m1][m2] *
          mptr->idunion.mcarig.evec[jb+1][jp1+1];
        }
       }
      }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_p1.Fm[i][j].f[m1][m2]=f[m1][m2];
      }
     }
    }
   }

   p=objptr->idunion.mca_p1.m;
   for (i1=0,kk=0;i1<nmode;i1++) {
    ip1=sorted[i1];
    for (j1=i1;j1<nmode;j1++) {
     jp1=sorted[j1];

     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       f[m1][m2]=objptr->idunion.mca_p1.Fm[i1][j1].f[m1][m2];
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mca_p1.f1 += scp/((_REAL)3);
     objptr->idunion.mca_p1.f1++;

     scp=0.;
     scp+=p[ip1].m[0]*p[jp1].m[0];
     scp+=p[ip1].m[1]*p[jp1].m[1];
     scp+=p[ip1].m[2]*p[jp1].m[2];
     *objptr->idunion.mca_p1.s1 += scp/((_REAL)3);
     objptr->idunion.mca_p1.s1++;

     kk++;
    }
   }
   objptr->idunion.mca_p1.s1 -= objptr->idunion.mca_p1.dim1;
   objptr->idunion.mca_p1.f1 -= objptr->idunion.mca_p1.dim1;

   if(objptr->idunion.mca_p1.type==1) return;

   p=objptr->idunion.mca_p1.m;
   l123=objptr->idunion.mca_p1.l123;
sono arrivato qui
   for (i1=0,i=0,kk=0,p1=p;i1<nmode;i1++,p1++) {

   for (j1=0,j=0,l123_2=l123;j1<nmode1;j1++) {
     for (j2=0;j2<nmode1;j2++) {
      for (j3=j2;j3<nmode1;j3++,j++,l123_2++,kk++) {

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
          objptr->idunion.mca_p1.Fm[i1][j1].f[m1][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          objptr->idunion.mca_p1.Fm[i1][j2].f[m1][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          objptr->idunion.mca_p1.Fm[i1][j3].f[m1][l]  *
          l123_2->dphidm3[m2][l];
        }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mca_p1.f12 += scp/((_REAL)3);
     objptr->idunion.mca_p1.f12++;

     scp=0.;
     scp+=p1->m[0]*l123_2->phi[0];
     scp+=p1->m[1]*l123_2->phi[1];
     scp+=p1->m[2]*l123_2->phi[2];
     *objptr->idunion.mca_p1.s12 += scp/((_REAL)3);
     objptr->idunion.mca_p1.s12++;

     }
    }
   }

   }
   objptr->idunion.mca_p1.s12 -= objptr->idunion.mca_p1.dim12;
   objptr->idunion.mca_p1.f12 -= objptr->idunion.mca_p1.dim12;

   l123=objptr->idunion.mca_p1.l123;
   for (i1=0,i=0,kk=0,l123_1=l123;i1<nmode1;i1++) {
    for (i2=0;i2<nmode1;i2++) {
     for (i3=i2;i3<nmode1;i3++,i++,l123_1++) {
     j1_0=i1;
     j2_0=i2;
     j3_0=i3;
   for (j1=j1_0,j=i,l123_2=l123_1;j1<nmode1;j1++,j3_0=0,j2_0=0) {
    for (j2=j2_0;j2<nmode1;j2++,j3_0=j2) {
     for (j3=j3_0;j3<nmode1;j3++,j++,l123_2++,kk++) {

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
          objptr->idunion.mca_p1.Fm[i1][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm1[m1][k] *
          objptr->idunion.mca_p1.Fm[i1][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm1[m1][k] *
          objptr->idunion.mca_p1.Fm[i1][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mca_p1.Fm[i2][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mca_p1.Fm[i2][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm2[m1][k] *
          objptr->idunion.mca_p1.Fm[i2][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mca_p1.Fm[i3][j1].f[k][l]  *
          l123_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mca_p1.Fm[i3][j2].f[k][l]  *
          l123_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123_1->dphidm3[m1][k] *
          objptr->idunion.mca_p1.Fm[i3][j3].f[k][l]  *
          l123_2->dphidm3[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=f[1][1];
     scp+=f[2][2];
     *objptr->idunion.mca_p1.f2 += scp/((_REAL)3);
     objptr->idunion.mca_p1.f2++;

     scp=0.;
     scp+=l123_1->phi[0]*l123_2->phi[0];
     scp+=l123_1->phi[1]*l123_2->phi[1];
     scp+=l123_1->phi[2]*l123_2->phi[2];
     *objptr->idunion.mca_p1.s2 += scp/((_REAL)3);
     objptr->idunion.mca_p1.s2++;

     }
    }
   }
     }
    }
   }
   objptr->idunion.mca_p1.s2 -= objptr->idunion.mca_p1.dim2;
   objptr->idunion.mca_p1.f2 -= objptr->idunion.mca_p1.dim2;

}

/************************************************************************/
void
write_first_p1(idstruct *dptr)
{
	int	i,j,i2,j2,k,l,m;
	int	states = dptr->idunion.mca_p1.cases;
	FILE 	*file = dptr->idunion.mca_p1.file,*dump;
        idstruct *mptr,*bptr;

	_REAL fact,tol,factns,sum,sum2;
        int rank,nbondm;

        int n1,n123,nn,nf;
        _REAL **S;
        _REAL **F,**Fm1,**tmp,csum,ctmp;
        _REAL *ui,*v,*w;
#ifdef IMSL
        d_complex *wi,*vi;
#else
        _REAL *vi;
#endif

        n1=dptr->idunion.mca_p1.n1;
        n123=dptr->idunion.mca_p1.n123;
        mptr=dptr->idunion.mca_p1.mca;
        bptr=mptr->idunion.mcarig.bead;
	fact = 1./ ((_REAL)states);
        factns = bptr->idunion.bead.factns;

        printf("write_first_p1: number of states: %d\n",states);
        for (i=1; i<=n1; i++) {
         for (j=i; j<=n1; j++) {
           *dptr->idunion.mca_p1.s1++ *= fact;
           *dptr->idunion.mca_p1.f1++ *= fact;
         }
        }
        dptr->idunion.mca_p1.s1 -= dptr->idunion.mca_p1.dim1;
        dptr->idunion.mca_p1.f1 -= dptr->idunion.mca_p1.dim1;

        if(dptr->idunion.mca_p1.type>1) {
        for (i=1; i<=n1; i++) {
         for (j=1; j<=n123; j++) {
           *dptr->idunion.mca_p1.s12++ *= fact;
           *dptr->idunion.mca_p1.f12++ *= fact;
         }
        }
        dptr->idunion.mca_p1.s12 -= dptr->idunion.mca_p1.dim12;
        dptr->idunion.mca_p1.f12 -= dptr->idunion.mca_p1.dim12;
        for (i=1; i<=n123; i++) {
         for (j=i; j<=n123; j++) {
           *dptr->idunion.mca_p1.s2++ *= fact;
           *dptr->idunion.mca_p1.f2++ *= fact;
         }
        }
        dptr->idunion.mca_p1.s2 -= dptr->idunion.mca_p1.dim2;
        dptr->idunion.mca_p1.f2 -= dptr->idunion.mca_p1.dim2;
        }

/* invert F */
        if(dptr->idunion.mca_p1.type==1) {
         nn = n1;
        } else if(dptr->idunion.mca_p1.type==2) {
         nn = n1+n123;
        } else {
         printf("write_first_p1: type not found.\n");
         exit(1);
        }

        F = matrix(1,nn,1,nn);
        S = matrix(1,nn,1,nn);

        for(i=1;i<=n1;i++) {
         for(j=i;j<=n1;j++) {
          F[i][j]=*dptr->idunion.mca_p1.f1++;
          F[j][i]=F[i][j];
          S[i][j]=*dptr->idunion.mca_p1.s1++;
          S[j][i]=S[i][j];
         }
        }
        dptr->idunion.mca_p1.s1 -= dptr->idunion.mca_p1.dim1;
        dptr->idunion.mca_p1.f1 -= dptr->idunion.mca_p1.dim1;

        if(dptr->idunion.mca_p1.type>1) {
        for(i=1;i<=n1;i++) {
         for(j=1;j<=n123;j++) {
          j2=n1+j;
          F[i][j2]=*dptr->idunion.mca_p1.f12++;
          F[j2][i]=F[i][j2];
          S[i][j2]=*dptr->idunion.mca_p1.s12++;
          S[j2][i]=S[i][j2];
         }
        }
        dptr->idunion.mca_p1.s12 -= dptr->idunion.mca_p1.dim12;
        dptr->idunion.mca_p1.f12 -= dptr->idunion.mca_p1.dim12;

        for(i=1;i<=n123;i++) {
         for(j=i;j<=n123;j++) {
          i2=n1+i;
          j2=n1+j;
          F[i2][j2]=*dptr->idunion.mca_p1.f2++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*dptr->idunion.mca_p1.s2++;
          S[j2][i2]=S[i2][j2];
         }
        }
        dptr->idunion.mca_p1.s2 -= dptr->idunion.mca_p1.dim2;
        dptr->idunion.mca_p1.f2 -= dptr->idunion.mca_p1.dim2;
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
          printf("write_first:\n");
          printf("S element is wrong: S(%d,%d) = %e\n",i,i,S[i][i]);
         }
        }

        fflush(stdout);
        Fm1 = matrix(1,nn,1,nn);
#ifdef IMSL
        tol=TOLSVD;
        ui = (_REAL *)get2(nn*nn,sizeof(_REAL));
        if(dptr->idunion.mca_p1.rigid==0) {
         for(i=1,k=0;i<=nn;i++)
          for(j=1;j<=nn;j++,k++)
           ui[k]=S[i][j];
        } else {
         for(i=1,k=0;i<=nn;i++)
          for(j=1;j<=nn;j++,k++)
           ui[k]=F[i][j];
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
        if(dptr->idunion.mca_p1.rigid==0) {
         rank=invert(S,nn,nn,Fm1);
        } else {
         rank=invert(F,nn,nn,Fm1);
        }
#endif
        fprintf(file,"write_first_p1: rank of F= %d\n",rank); fflush(file);

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

        if(dptr->idunion.mca_p1.rigid==0) {

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          csum=0.;
          for(k=1;k<=nn;k++) {
           csum+=Fm1[i][k]*F[k][j];
          }
          tmp[i][j]=csum;
         }
        }

        } else {

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
/*
        for(i=1;i<=nn;i++) {
         tol=1./S[i][i];
         csum=tol*Fm1[i][i];
         fprintf(file,"tau(%d)= %e\n",i,factns*csum);
        }
        fflush(file);
*/
        }

        w=vector(1,nn);
#ifdef IMSL
        for(i=1,l=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,l++) {
          ui[l]=tmp[i][j];
         }
        }
        wi = imsl_d_eig_gen(nn,ui,IMSL_VECTORS,&vi,0);
        if(dptr->idunion.mca_p1.rigid==0) {
         for(i=0;i<nn;i++) {
          w[i+1]=1./wi[i].re;
         }
        } else {
         for(i=0;i<nn;i++) {
          w[i+1]=wi[i].re;
         }
        }
        free(ui);
        free(wi);
#else
        vi=(_REAL *)get2(nn*nn,sizeof(_REAL));
#ifdef EISPACK
        eigvec_ev(nn,tmp,vi,w);
        if(dptr->idunion.mca_p1.rigid==0) {
         for(i=0;i<nn;i++) {
          w[i+1]=1./w[i+1];
         }
        }
#else
        eigvec_nr(tmp,nn,w,vi);
#endif
#endif

        fprintf(file,"Fm1*S eigenvalues:\n");
        nf=0;
        for(i=1;i<=nn;i++) {
         fprintf(file,"w(%d)= %e\n",i,w[i]);
         if(Nabs(w[i])>=ZERO) {
          nf++;
         }
        }
        fflush(file);

        if(nf!=dptr->idunion.mca_p1.n1) {
         printf("write_first_p1: WARNING: nf=%d is different from n1=%d.\n",
                 nf,n1);
        }
        fflush(file);

        dptr->idunion.mca_p1.nf = nf;
        nf=nn;
        dptr->idunion.mca_p1.eval= (_REAL *)get2(nf,sizeof(_REAL));
        dptr->idunion.mca_p1.evec= matrix(1,nn,1,nf);

        nf=0;
        for(i=1;i<=nn;i++) {
         dptr->idunion.mca_p1.eval[i-1]=1./w[i];
        }
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
#ifdef IMSL
          dptr->idunion.mca_p1.evec[i][j]=vi[k].re;
#else
          dptr->idunion.mca_p1.evec[i][j]=vi[k];
#endif
         }
        }
        eigsrt(w,dptr->idunion.mca_p1.evec,nn);
        for(i=1;i<=nn;i++) {
         dptr->idunion.mca_p1.eval[i-1]=1./w[i];
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
          dptr->idunion.mca_p1.eval[nf-1]=1./wi[i-1].re;
          vt=&vi[k];
          for(j=1;j<=nn;j++) {
           dptr->idunion.mca_p1.evec[j][nf]=(*vt).re;
           vt+=nn;
          }
         }
        }
*/

/*
        fprintf(file,"evec:\n");
        for(j=1;j<=3;j++) {
         for(i=1;i<=nn;i++) {
          fprintf(file,"%d %d %e\n",i,j,dptr->idunion.mca_p1.evec[i][j]);
         }
        }
*/

        nf=nn;
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mca_p1.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mca_p1.evec[k][i]*tmp[k][j];
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
           sum+=dptr->idunion.mca_p1.evec[i][k]*sum2;
          }
          v[m]=sum;
         }
        }
        for(i=1,k=0;i<=nf;i++) {
         for(j=1;j<=nf;j++,k++) {
          dptr->idunion.mca_p1.evec[i][j]=v[k];
         }
        }
*/

/* normalize eigenvectors */
        for(j=1,k=0;j<=nn;j++) {
         if(Nabs(Fm1[j][j]) > ZERO) {
          k++;
          sum=Fm1[j][j];
         } else {
          sum=1.;
         }
         sum=l_sqrt(sum);
         for(i=1;i<=nn;i++) {
          dptr->idunion.mca_p1.evec[i][j]/=sum;
         }
        }
        if(k!=nf) {
         printf("write_first: missed zero in Ct*S*C.\n");
        }
/*
        sum2=0.;
        for(j=1,k=0;j<=nn;j++) {
         if(Nabs(Fm1[j][j]) > ZERO) {
          k++;
          sum=Fm1[j][j];
         } else {
          sum=1.;
         }
         sum2+=sum;
        }
        sum2=l_sqrt(sum2);
        if(k!=nf) {
         printf("write_first: missed zero in Ct*S*C.\n");
        }
        for(j=1,k=0;j<=nn;j++) {
         for(i=1;i<=nn;i++) {
          dptr->idunion.mca_p1.evec[i][j]/=sum2;
         }
        }
*/

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mca_p1.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        sum2=0.;
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mca_p1.evec[k][i]*tmp[k][j];
          }
          S[i][j]=sum;
         }
         sum2+=S[i][i];
        }
        printf("sum_i(psi_i*psi_i)= %e\n",sum2);
/*
        fprintf(file,"Ct*S1*C:\n");
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nf;j++) {
          fprintf(file,"%d %d %e\n",i,j,S[i][j]);
         }
        }
        fflush(file);
*/

/*
        fprintf(file,"evec2:\n");
        for(j=1;j<=nf;j++) {
         for(i=1;i<=nn;i++) {
          fprintf(file,"%d %d %e\n",i,j,dptr->idunion.mca_p1.evec[i][j]);
         }
        }
*/

        free_matrix(F,1,nn,1,nn);
        free_matrix(S,1,nn,1,nn);
        free_matrix(tmp,1,nn,1,nn);

/*
        if(dptr->idunion.mca_p1.p1_to_p2) {
         dptr->idunion.mca_p1.npsi2=nn*(nn+1)/2;
         dptr->idunion.mca_p1.psi2=
          (rot2struct *)get2(nn*(nn+1)/2,sizeof(rot2struct));
        } 
*/
}
