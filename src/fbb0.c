/* this file is included in mcarig_m1.c */
{ _REAL zim1;
  
     for(m1=0;m1<3;m1++)
      for(m2=0;m2<3;m2++)
       f[m1][m2]=0.;

     if(lb1->id1 == lb2->id1) {
      zim1 = 1.;
      f[0][0] += zim1;
      f[1][1] += zim1;
      f[2][2] += zim1;
     }
     if(lb1->id1 == lb2->id2) {
      zim1 = 1.;
      f[0][0] -= zim1;
      f[1][1] -= zim1;
      f[2][2] -= zim1;
     }
     if(lb1->id2 == lb2->id1) {
      zim1 = 1.;
      f[0][0] -= zim1;
      f[1][1] -= zim1;
      f[2][2] -= zim1;
     }
     if(lb1->id2 == lb2->id2) {
      zim1 = 1.;
      f[0][0] += zim1;
      f[1][1] += zim1;
      f[2][2] += zim1;
     }
}
