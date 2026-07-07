/* this file is included in fbb.c and in rotdiff.c */
/* computes the hydrodynamic interaction */
{ _REAL zit,zjt,r4,r6,f1,f2,screen;

      D->xx = 0.;
      D->xy = 0.;
      D->xz = 0.;
      D->yy = 0.;
      D->yz = 0.;
      D->zz = 0.;

      zit=zi;
      zjt=zj;
      overlap=0;
      rsq = dx*dx+dy*dy+dz*dz;
//      if(rsq<=cutsq)
      if(rsq>cutsq)
      {

      rmod = l_sqrt(rsq);
      /* weight=rmod; */
      zsq = (zit+zjt)*(zit+zjt);
      rm2 = 1./rsq;
      rm1 = 1./rmod;
      dx *= rm1; dy *= rm1; dz *= rm1;

#ifdef HYDRO_ROTNEPRAGER
       if(rsq>=zsq) {
        Rm1 = rm1*zr;
        zij = 0.5*(zit+zjt);
/*        zij = (zit*zit+zjt*zjt)*rm2; */
        zij = 2.*zij*zij*rm2;
        D->xx = dx*dx;
        D->xy = dx*dy;
        D->xz = dx*dz;
        D->yy = dy*dy;
        D->yz = dy*dz;
        D->zz = dz*dz;
        D->xx = Rm1*(1. + D->xx + zij*(1./3. -D->xx));
        D->xy = Rm1*(     D->xy + zij*(      -D->xy));
        D->xz = Rm1*(     D->xz + zij*(      -D->xz));
        D->yy = Rm1*(1. + D->yy + zij*(1./3. -D->yy));
        D->yz = Rm1*(     D->yz + zij*(      -D->yz));
        D->zz = Rm1*(1. + D->zz + zij*(1./3. -D->zz));
       } else {
        overlap=1;
        Rm1 = 4.*zr/3.;
        zij = 0.5*(zit+zjt);
        aa = 1./zij;
        bb = 3.*rmod/zij/32.;
        D->xx = Rm1*aa*( 1.-3.*bb + bb*dx*dx);
        D->xy = Rm1*aa*(            bb*dx*dy);
        D->xz = Rm1*aa*(            bb*dx*dz);
        D->yy = Rm1*aa*( 1.-3.*bb + bb*dy*dy);
        D->yz = Rm1*aa*(            bb*dy*dz);
        D->zz = Rm1*aa*( 1.-3.*bb + bb*dz*dz);
/*
        D->xx = 0.;
        D->xy = 0.;
        D->xz = 0.;
        D->yy = 0.;
        D->yz = 0.;
        D->zz = 0.;
*/
       }
#endif
#ifdef HYDRO_OSEEN
       if(bptr->idunion.bead.screen!=0.) 
        screen = exp(-rmod/bptr->idunion.bead.screen);
       else screen=1.;
       Rm1 = screen*rm1*zr;
       D->xx = Rm1*(1. + dx*dx);
       D->xy = Rm1*(     dx*dy);
       D->xz = Rm1*(     dx*dz);
       D->yy = Rm1*(1. + dy*dy);
       D->yz = Rm1*(     dy*dz);
       D->zz = Rm1*(1. + dz*dz);
#endif
#ifdef HYDRO_OTTINGER
       zij = 0.5*(zit+zjt);
       aa = 2.*zij/l_sqrt(3.);
       aa *= aa;
       bb = rsq + aa;
       bb = 1./bb;
       bb = bb*bb*bb;
       Rm1 = rm1*bb*zr;
       r4 = rsq*rsq;
       r6 = r4*rsq;
       f1 = r6 + (7./2.)*aa*r4 + (9./2.)*aa*aa*rsq;
       f2 = r6 + (3./2.)*aa*r4 - (3./2.)*aa*aa*rsq;
       D->xx = dx*dx;
       D->xy = dx*dy;
       D->xz = dx*dz;
       D->yy = dy*dy;
       D->yz = dy*dz;
       D->zz = dz*dz;
       D->xx = Rm1*(f1 + f2*D->xx);
       D->xy = Rm1*(     f2*D->xy);
       D->xz = Rm1*(     f2*D->xz);
       D->yy = Rm1*(f1 + f2*D->yy);
       D->yz = Rm1*(     f2*D->yz);
       D->zz = Rm1*(f1 + f2*D->zz);
#endif
      } /* endif cutsq */
}
