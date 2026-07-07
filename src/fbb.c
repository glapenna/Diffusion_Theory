/* this file is included in mcarig_m1.c */
{ rdiadstruct D[1];
  _REAL rm1,rm2,rmod,rsq,Rm1,zim1,zij,zsq,aa,bb;
  
     for(m1=0;m1<3;m1++)
      for(m2=m1;m2<3;m2++)
       f[m1][m2]=0.;

     dx=lb[ib].x1 - lb[jb].x1;
     dy=lb[ib].y1 - lb[jb].y1;
     dz=lb[ib].z1 - lb[jb].z1;
#ifdef PBC
     dx -= box[0]*rint(dx*boxm1[0]);
     dy -= box[1]*rint(dy*boxm1[1]);
     dz -= box[2]*rint(dz*boxm1[2]);
#endif
     if(lb[ib].id1 == lb[jb].id1) {
      zim1 = 1./(lb[ib].zeta1);
      f[0][0] += zim1;
      f[1][1] += zim1;
      f[2][2] += zim1;
     } else {
      zi=lb[ib].zeta1;
      zj=lb[jb].zeta1;
#include "hydrod.c"
      f[0][0] += D->xx;
      f[0][1] += D->xy;
      f[0][2] += D->xz;
      f[1][1] += D->yy;
      f[1][2] += D->yz;
      f[2][2] += D->zz;
     }
     dx=lb[ib].x1 - lb[jb].x2;
     dy=lb[ib].y1 - lb[jb].y2;
     dz=lb[ib].z1 - lb[jb].z2;
#ifdef PBC
     dx -= box[0]*rint(dx*boxm1[0]);
     dy -= box[1]*rint(dy*boxm1[1]);
     dz -= box[2]*rint(dz*boxm1[2]);
#endif
     if(lb[ib].id1 == lb[jb].id2) {
      zim1 = 1./(lb[ib].zeta1);
      f[0][0] -= zim1;
      f[1][1] -= zim1;
      f[2][2] -= zim1;
     } else {
      zi=lb[ib].zeta1;
      zj=lb[jb].zeta2;
#include "hydrod.c"
      f[0][0] -= D->xx;
      f[0][1] -= D->xy;
      f[0][2] -= D->xz;
      f[1][1] -= D->yy;
      f[1][2] -= D->yz;
      f[2][2] -= D->zz;
     }
     dx=lb[ib].x2 - lb[jb].x1;
     dy=lb[ib].y2 - lb[jb].y1;
     dz=lb[ib].z2 - lb[jb].z1;
#ifdef PBC
     dx -= box[0]*rint(dx*boxm1[0]);
     dy -= box[1]*rint(dy*boxm1[1]);
     dz -= box[2]*rint(dz*boxm1[2]);
#endif
     if(lb[ib].id2 == lb[jb].id1) {
      zim1 = 1./(lb[ib].zeta2);
      f[0][0] -= zim1;
      f[1][1] -= zim1;
      f[2][2] -= zim1;
     } else {
      zi=lb[ib].zeta2;
      zj=lb[jb].zeta1;
#include "hydrod.c"
      f[0][0] -= D->xx;
      f[0][1] -= D->xy;
      f[0][2] -= D->xz;
      f[1][1] -= D->yy;
      f[1][2] -= D->yz;
      f[2][2] -= D->zz;
     }
     dx=lb[ib].x2 - lb[jb].x2;
     dy=lb[ib].y2 - lb[jb].y2;
     dz=lb[ib].z2 - lb[jb].z2;
#ifdef PBC
     dx -= box[0]*rint(dx*boxm1[0]);
     dy -= box[1]*rint(dy*boxm1[1]);
     dz -= box[2]*rint(dz*boxm1[2]);
#endif
     if(lb[ib].id2 == lb[jb].id2) {
      zim1 = 1./(lb[ib].zeta2);
      f[0][0] += zim1;
      f[1][1] += zim1;
      f[2][2] += zim1;
     } else {
      zi=lb[ib].zeta2;
      zj=lb[jb].zeta2;
#include "hydrod.c"
      f[0][0] += D->xx;
      f[0][1] += D->xy;
      f[0][2] += D->xz;
      f[1][1] += D->yy;
      f[1][2] += D->yz;
      f[2][2] += D->zz;
     }

}
