#include "head.h"

// void *get2(int,int);

/************************************************************************/

void
solve_mateq(_REAL *s,_REAL *f,int n, int dim, int states,
            _REAL statesw, _REAL **evec, _REAL *eval,
            int rigid, _REAL factns, FILE *file, int *nf)
{
	int	i,j,k,l,m,nn,rank,nft;
	_REAL fact,tol,
              **F,**S,**Fm1,*w,**tmp,
              sum,sum1,sum2;
        _REAL *u,*v;
#ifdef IMSL
        d_complex *vi,*wi;
#else
        _REAL *vi;
#endif

	/* fact = 1./ ((_REAL)states); */
	fact = 1./ statesw;

        printf("solve_mateq: number of states: %d\n",states);
        printf("solve_mateq: matrix order:     %d\n",n);
        fflush(stdout);

        nn = n;
        F = matrix(1,nn,1,nn);
        S = matrix(1,nn,1,nn);

        for (i=1; i<=n; i++) {
         for (j=i; j<=n; j++) {
           *s *= fact;
           *f *= fact;
           S[i][j]=*s;
           S[j][i]=S[i][j];
           F[i][j]=*f;
           F[j][i]=F[i][j];
           s++;
           f++;
         }
        }
        s -= dim;
        f -= dim;

/* invert F */

/*
*/
        fprintf(file,"S,F:\n");
        for(i=1;i<=nn;i++)
         for(j=1;j<=nn;j++)
          fprintf(file,"%d %d %e %e\n",i,j,S[i][j],F[i][j]);
        fflush(file);

        for(i=1;i<=nn;i++) {
         if(S[i][i]<ZERO) {
          printf("solve_mateq:\n");
          printf("S element is wrong: S(%d,%d) = %e\n",i,i,S[i][i]);
         }
        }

/*
        eval = (_REAL *) get2(nn,sizeof(_REAL));
        evec = matrix(1,nn,1,nn);
*/

        Fm1 = matrix(1,nn,1,nn);
        tmp = matrix(1,nn,1,nn);
/*
#ifdef IMSL
        tol = TOLSVD;
        u = vector(0,nn*nn-1);
        if(rigid==0) {
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
        if(rigid==0) {
         rank=invert(F,nn,nn,Fm1);
        } else {
         printf("not implemented...\n");
        }
#endif
        printf("solve_mat: rank of F: %d\n",rank);

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
          printf("solve_mateq:\n");
          printf("SF^-1S diagonal element is wrong: F(%d,%d) = %e\n",
           i,i,Fm1[i][i]);
         }
         if(Nabs(S[i][i])>ZERO) {
          Fm1[i][i] /= S[i][i];
         }
        }

        for(i=1;i<=nn;i++) {
         printf("solve_mateq: tau(%d)= %e ",i,1.*Fm1[i][i]);
         printf("%e\n",factns*Fm1[i][i]);
        }
        fflush(stdout);
*/

#ifdef IMSL
        tol = TOLSVD;
        u = vector(0,nn*nn-1);
        if(rigid==0) {
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
        if(rigid==0) {
         rank=invert(S,nn,nn,Fm1);
        } else {
         rank=invert(F,nn,nn,Fm1);
        }
#endif
        printf("solve_mateq: rank of S: %d\n",rank);

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

        if(rigid==0) {

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

/* diagonalize A = S-1*F or A = F^-1*S if rigid */
        
        w=vector(1,nn);
#ifdef IMSL
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
          u[k]=tmp[i][j];
         }
        }
        wi = imsl_d_eig_gen(nn,u,IMSL_VECTORS,&vi,0);
/*
        if(rigid==0) {
         for(i=0;i<nn;i++) {
          w[i+1] = 1./wi[i].re;
         }
        } else {
         for(i=0;i<nn;i++) {
          w[i+1] = wi[i].re;
         }
        }
*/
        if(rigid==0) {
         for(i=0;i<nn;i++) {
          w[i+1] = wi[i].re;
         }
        } else {
         for(i=0;i<nn;i++) {
          w[i+1] = 1./wi[i].re;
         }
        }
        free(wi);
        free_vector(u,0,nn*nn-1);
#else
        vi = (_REAL *)get2(nn*nn,sizeof(_REAL));
#ifdef EISPACK
        eigvec_ev(nn,tmp,vi,w);
/*
        if(rigid==0) {
         for(i=0;i<nn;i++) {
          w[i+1] = 1./w[i+1];
         }
        }
*/
#else
        eigvec_nr(tmp,nn,w,vi);
#endif
#endif
        /* printf("w1= %e\n",w[1]); fflush(stdout); exit(0); */

        for(i=1,sum=0.,fact=0.,nft=0;i<=nn;i++) {
         sum+=w[i];
         if(Nabs(w[i])>ZERO)
         {
          nft++;
          /* printf("solve_mateq: lambda2(%d)= %e\n",i,1./w[i]); */
          printf("solve_mateq: lambda2(%d)= %e\n",i,w[i]);
          fact+= w[i];
         }
        }
        printf("\n"); fflush(stdout);

        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
#ifdef IMSL
          evec[i][j] = vi[k].re;
#else
          evec[i][j] = vi[k];
#endif
         }
        }
        eigsrt(w,evec,nn);
        for(i=1;i<=nn;i++) {
         /* this is for write_acff */
         /* eval[i-1] = w[i]; */
         /* this is for write_acff2 */
         /* eval[i-1] = 1./w[i]; */
         eval[i-1] = w[i];
        }
/*
        free_vector(w,1,nn);
        free(vi);
*/

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=evec[k][i]*tmp[k][j];
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
          evec[i][j]/=sum;
         }
        }
        if(k!=nft) {
         printf("solve_mateq: missed zero in Ct*S*C.\n");
        }

        printf("solve_mateq: nf= %d\n",nft); fflush(stdout);
        /* *nf = nft; */
        *nf = nn;

	{ FILE *ftmp;
        ftmp=fopen("evec_m1","w");	 
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          fprintf(ftmp,"%d %d %e\n",i,j,evec[i][j]);
         }
        }
	fclose(ftmp);
        }

        return;

/* transform F in Ct*F*C and S in Ct*S*C */
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=F[i][k]*evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=evec[k][i]*tmp[k][j];
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
           sum+=S[i][k]*evec[k][j];
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
           sum+=evec[k][i]*tmp[k][j];
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
