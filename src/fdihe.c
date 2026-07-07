/* this file is included in mca_dihe.c */
{ rdiadstruct D[1];
  int k1,k2;
  _REAL rm1,rm2,rmod,rsq,Rm1,zim1,zij,zsq,aa,bb;
  
   f=0.;

   for(k1=0;k1<4;k1++)
   {
    for(k2=0;k2<4;k2++)
    {
     if(lb1->id[k1] == lb2->id[k2]) {
      zim1 = 1./(lb1->zeta[k1]);
      f += zim1*lb1->fx[k1]*lb2->fx[k2];
      f += zim1*lb1->fy[k1]*lb2->fy[k2];
      f += zim1*lb1->fz[k1]*lb2->fz[k2];
     } else {
      zi=lb1->zeta[k1];
      zj=lb2->zeta[k2];
      dx=lb1->x[k1] - lb2->x[k2];
      dy=lb1->y[k1] - lb2->y[k2];
      dz=lb1->z[k1] - lb2->z[k2];
#ifdef PBC
      dx -= box[0]*rint(dx*boxm1[0]);
      dy -= box[1]*rint(dy*boxm1[1]);
      dz -= box[2]*rint(dz*boxm1[2]);
#endif
#include "hydrod.c"
      f += D->xx*lb1->fx[k1]*lb2->fx[k2];
      f += D->xy*lb1->fx[k1]*lb2->fy[k2];
      f += D->xz*lb1->fx[k1]*lb2->fz[k2];
      f += D->xy*lb1->fy[k1]*lb2->fx[k2];
      f += D->yy*lb1->fy[k1]*lb2->fy[k2];
      f += D->yz*lb1->fy[k1]*lb2->fz[k2];
      f += D->xz*lb1->fz[k1]*lb2->fx[k2];
      f += D->yz*lb1->fz[k1]*lb2->fy[k2];
      f += D->zz*lb1->fz[k1]*lb2->fz[k2];
     }
    }
   }

}
