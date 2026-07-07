#include "head.h"
#define MAXDEG 2

void ortho(_REAL *eval,_REAL **evec, _REAL **metr, int n, int nf)
{
 int i,j,k,*deg,**ideg;
 _REAL v1,v2;

 eval--;
 deg=ivector(1,n);
 ideg=imatrix(1,n,1,MAXDEG);

 for(i=1;i<=n;i++) {
  deg[i]=1;
  for(j=1;j<=MAXDEG;j++) {
   ideg[i][j]=i;
  }
 }
 for(i=1;i<=n-1;i++) {
  for(j=i+1;j<=n;j++) {
   if( Nabs(eval[i]-eval[j]) < ZERO && deg[j]!=0) {
    deg[i]++; deg[j]=0; 
    if(deg[i]>MAXDEG) {
     printf("ortho: ideg dimension exceeded.\n");
     exit(1);
    }
    ideg[i][deg[i]]=j;
    ideg[j][1]=-1;
   }
  }
 }

/*
 for(i=1;i<=n;i++) {
  printf("deg(%d)= ",i);
  for(j=1;j<=deg[i];j++) {
   printf("%d ",ideg[i][j]);
  }
  printf("\n");
 }
*/

/*
 for(i=1;i<=n;i++) {
  if(deg[i]==2) {
   k=ideg[i][deg[i]];
   for(j=1;j<=n;j++) {
    v1=evec[j][i];
    v2=evec[j][k];
    v2-=v1*metr[i][k];
    evec[j][k]=v2;
   }
  }
 }
*/
 for(i=1,nf=0;i<=n;i++) {
  if(deg[i]==0) {
   eval[i]=1e+09;
   nf++;
  }
 }
 nf=n-nf;

}

