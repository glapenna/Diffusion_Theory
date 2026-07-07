/************************************************************************/
void
write_acff_stress(idstruct *dptr)
{
 extern parmstruct *prm;
 _REAL **proj=dptr->idunion.acff.proj;
 _REAL *eval;
 _REAL **evec;
 _REAL tau,fact,sum,norm;
 _REAL factns;
 int i,j,k,l,nf,mm;
 FILE *file=dptr->idunion.acff.file;

/* fact=1./((_REAL)dptr->idunion.acff.cases); */
 fact=1./dptr->idunion.acff.casesw;
 mm=5;

 if(dptr->idunion.acff.model->id_type==MCA_P2) {
  nf = 1;
  factns = modptr->idunion.mca_p2.bead->idunion.bead.factns;
  eval=modptr->idunion.mca_p2.eval;
  evec=modptr->idunion.mca_p2.evec;
  mm=5;
 } else {
  printf("write_acff_stress: model not implemented yet.\n");
  fflush(stdout);
  exit(1);
 }
 norm=1./((_REAL)mm);

 printf("write_acff_stress: number of basis functions, nbase= %d\n",dptr->idunion.acff.nbase);
 printf("write_acff_stress: nuber of functions, nf= %d\n",nf);
 printf("write_acff_stress: number of 2nd-rank components, mm= %d\n",mm);
 fflush(stdout);

 for(j=1;j<=dptr->idunion.acff.nbase;j++) {
  proj[1][j]*=fact;
  fprintf(file,"%e\n",proj[1][j]);
 }
 fclose(file);

/*
 printf("projections:\n");
 for(j=1;j<=dptr->idunion.acff.nbase;j++) {
  printf("proj(%d)= %e\n",j,proj[1][j]);
 }
 printf("eigenvectors:\n");
 for(j=1;j<=dptr->idunion.acff.nbase;j++) {
  printf("eval(%d)= %e\n",j,eval[j-1]);
  for(i=1;i<=modptr->idunion.mca_p1l.n1+modptr->idunion.mca_p1l.n123m;i++) {
   printf("evec(%d,%d)= %e\n",i,j,evec[i][j]);
  }
 }
*/

 sum=0.;
 for(j=1;j<=dptr->idunion.acff.nbase;j++) {
  for(k=1;k<=dptr->idunion.acff.nbase;k++) {
   sum += proj[1][j]*proj[1][k]*evec[j][i]*evec[k][i];
  }
 }
 sum/=norm;
 printf("ACF0= %e\n",sum);
 fflush(stdout);
}
