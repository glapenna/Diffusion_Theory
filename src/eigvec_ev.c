#include "head.h"

// void *get2(int,int);

#ifdef EISPACK
#define MAXEIGVEC 2000
_REAL u_ev[MAXEIGVEC][MAXEIGVEC],v_ev[MAXEIGVEC][MAXEIGVEC],wi_ev[MAXEIGVEC];

void eigvec_ev(int n,_REAL **u,_REAL *v,_REAL *w)
{
 int i,j,k,l,max;
/* _REAL **a,sum; */

 if(n>MAXEIGVEC) {
  printf("eigvec: too large matrices for static version of eigvec.\n");
  printf("MAXEIGVEC= %d\n",MAXEIGVEC);
  exit(1);
 }
 max=MAXEIGVEC;
 
/* a=matrix(1,n,1,n); */
 for(i=1;i<=n;i++) {
  for(j=1;j<=n;j++) {
   u_ev[j][i]=u[i][j];
/*   a[i][j]=u[i][j]; */
  }
 }

 eigvec_(&max,&n,&u_ev[1][1],&w[1],&wi_ev[1],&v_ev[1][1]);
 for(i=1;i<=n;i++) {
  for(j=1;j<=n;j++) {
   u[i][j]=v_ev[j][i];
  }
 }
 eigsrt(w,u,n);
 for(i=1,k=0;i<=n;i++) {
  for(j=1;j<=n;j++,k++) {
   v[k]=u[i][j];
  }
 }
/* test if Ct*A*C=D */
/*
 for(j=1;j<=n;j++) {
  for(i=1;i<=n;i++) {
   for(k=1,sum=0.;k<=n;k++) {
    sum+=a[i][k]*u[k][j];
   }
   printf("%d %d %e %e\n",i,j,sum,w[j]*u[i][j]);
  }
 }
 free_matrix(a,1,n,1,n);
*/

}

void adjust_ev(int n, _REAL **Fm1, _REAL **tmp, _REAL **evec)
{
  _REAL *w,*vi,sum,sum2;
  int i,j,k,l,m;

  vi=(_REAL *)get2(n*n,sizeof(_REAL));
  w=vector(1,n);
  eigvec_ev(n,Fm1,vi,w);

  for(i=1,k=0;i<=n;i++) {
   if(Nabs(w[i]) > ZERO) {
    w[i]=1./l_sqrt(w[i]);
    for(j=1;j<=n;j++,k++) {
     tmp[i][j]=vi[k];
    }
   } else {
    w[i]=1.;
    for(j=1;j<=n;j++,k++) {
     tmp[i][j]=0.;
    }
    tmp[i][i]=1.;
   }
  }
  for(i=1,m=0;i<=n;i++) {
   for(j=1;j<=n;j++,m++) {
    for(k=1,sum=0.;k<=n;k++) {
     for(l=1,sum2=0.;l<=n;l++) {
      sum2+=tmp[k][l]*w[l]*tmp[j][l];
     }
     sum+=evec[i][k]*sum2;
    }
    vi[m]=sum;
   }
  }
  for(i=1,k=0;i<=n;i++) {
   for(j=1;j<=n;j++,k++) {
    evec[i][j]=vi[k];
   }
  }

  free(vi);
  free_vector(w,1,n);
}

#endif
