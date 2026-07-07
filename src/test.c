{
 _REAL d,dx0,dy0,dz0,dx,dy,dz,
       xt[8],yt[8],zt[8],uxt[6],uyt[6],uzt[6],dmin,dmax,dper,Rf,Rc,
       dt,RcRf,RcRf2,dper2;
 int imin,k,min,max;
/*  possibile specificare un delta solo per
     l'ultimo passo  */
/* if(sigma<=0.001) return(1); */
 
 dx0=x-xa; dy0=y-ya; dz0=z-za;
 d=dx0*dx0+dy0*dy0+dz0*dz0;
 igood=0;
 if (d<=sum2) {
  if(d<diff2) {
   igood=-1;
  } else {
   *rt++ = x+cellh; *rt++ = y+cellh; *rt++ = z-cellh;
   *rt++ = x-cellh; *rt++ = y+cellh; *rt++ = z-cellh;
   *rt++ = x-cellh; *rt++ = y-cellh; *rt++ = z-cellh;
   *rt++ = x+cellh; *rt++ = y-cellh; *rt++ = z-cellh;
   *rt++ = x+cellh; *rt++ = y+cellh; *rt++ = z+cellh;
   *rt++ = x-cellh; *rt++ = y+cellh; *rt++ = z+cellh;
   *rt++ = x-cellh; *rt++ = y-cellh; *rt++ = z+cellh;
   *rt++ = x+cellh; *rt++ = y-cellh; *rt++ = z+cellh;
   rt-=24;
   for(k=0,imin=0;k<8;k++) {
    dx=*rt++ -xa; dy=*rt++ -ya; dz=*rt++ -za;
    dt=dx*dx+dy*dy+dz*dz;
    if(dt<sigma2) imin++;
   }
   rt-=24;
/*
   xt[0]=x+cellh; yt[0]=y+cellh; zt[0]=z-cellh;
   xt[1]=x-cellh; yt[1]=y+cellh; zt[1]=z-cellh;
   xt[2]=x-cellh; yt[2]=y-cellh; zt[2]=z-cellh;
   xt[3]=x+cellh; yt[3]=y-cellh; zt[3]=z-cellh;
   xt[4]=x+cellh; yt[4]=y+cellh; zt[4]=z+cellh;
   xt[5]=x-cellh; yt[5]=y+cellh; zt[5]=z+cellh;
   xt[6]=x-cellh; yt[6]=y-cellh; zt[6]=z+cellh;
   xt[7]=x+cellh; yt[7]=y-cellh; zt[7]=z+cellh;
   for(k=0,imin=0;k<8;k++) {
    dx=xt[k]-xa; dy=yt[k]-ya; dz=zt[k]-za;
    dt=dx*dx+dy*dy+dz*dz;
    if(dt<sigma2) imin++;
   }
*/
   if(imin==0) {
    xt[0]=x-cellh; yt[0]=y; zt[0]=z;
    xt[1]=x+cellh; yt[1]=y; zt[1]=z;
    xt[2]=x; yt[2]=y-cellh; zt[2]=z;
    xt[3]=x; yt[3]=y+cellh; zt[3]=z;
    xt[4]=x; yt[4]=y; zt[4]=z-cellh;
    xt[5]=x; yt[5]=y; zt[5]=z+cellh;
    uxt[0]=-1.; uyt[0]= 0.; uzt[0]= 0.;
    uxt[1]=+1.; uyt[1]= 0.; uzt[1]= 0.;
    uxt[2]= 0.; uyt[2]=-1.; uzt[2]= 0.;
    uxt[3]= 0.; uyt[3]=+1.; uzt[3]= 0.;
    uxt[4]= 0.; uyt[4]= 0.; uzt[4]=-1.;
    uxt[5]= 0.; uyt[5]= 0.; uzt[5]= 1.;
    for(k=0,dmin=10000.;k<6;k++) {
     dx=xt[k]-xa; dy=yt[k]-ya; dz=zt[k]-za;
     dt = Nabs(dx*uxt[k]+dy*uyt[k]+dz*uzt[k]);
     if(dt<dmin) { 
      min=k;
      dmin=dt;
     }
    }
    dper=dmin-cellh;
    dper2=dper*dper;
    Rc = l_sqrt(sigma2 - dper2);
    Rf = cellh;
    RcRf = Rc+Rf;
    RcRf2 = RcRf*RcRf;
    dper = dx0*uxt[min]+dy0*uyt[min]+dz0*uzt[min];
    dper2=dper*dper;
    dt = d - dper2;
    if(dt<=RcRf2) {
     igood=0;
    } else {
     igood=1;
    }
   } else if(imin==8) {
    igood=-1;
   }
  }
 } else {
  igood=1;
 }
}
