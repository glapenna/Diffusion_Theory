#include "head.h"

// extern char tok[];
// char tok[];
void *get2(int,int);
idstruct *startid(),*getstruct();

void setup_asa()
{
 int i;
 char id[TOKLEN];
 idstruct *dptr,*grptr,*grptrc;

 dptr=startid(ASA,id);

 dptr->idunion.asa.f = 1.;
 dptr->idunion.asa.close = 0;
 dptr->idunion.asa.probe = 0.;
 dptr->idunion.asa.res = 0.05;
 dptr->idunion.asa.nmax = 10;
 dptr->idunion.asa.groupc = NULL;
 dptr->idunion.asa.output = 0;

 gettok();
 if ((grptr = getstruct(tok)) == NULL)
         inerr("setup_asa: no such id: ", tok);
 if (grptr->id_type != GROUP)
         inerr("setup_asa: wrong type of id: ", tok);
 dptr->idunion.asa.group=grptr;

 for (;;) {
     gettok();
      if (!strcmp(tok, ";")) {
       break;
      } else if (!strcmp(tok, "FACT")) {
      gettok();
      if (!tokREAL(&dptr->idunion.asa.f))
       inerr("setup_asa: expected f", "");
      } else if (!strcmp(tok, "PROBE")) {
      gettok();
      if (!tokREAL(&dptr->idunion.asa.probe))
       inerr("setup_asa: expected probe", "");
      } else if (!strcmp(tok, "RES")) {
      gettok();
      if (!tokREAL(&dptr->idunion.asa.res))
       inerr("setup_asa: expected res", "");
      } else if (!strcmp(tok, "NMAX")) {
      gettok();
      if (!tokint(&dptr->idunion.asa.nmax))
       inerr("setup_asa: expected nmax", "");
      } else if (!strcmp(tok, "CLOSE")) {
       dptr->idunion.asa.close=1;
       gettok();
       if ((grptrc = getstruct(tok)) == NULL)
         inerr("setup_asa: no such id: ", tok);
       if (grptrc->id_type != GROUP)
         inerr("setup_asa: wrong type of id: ", tok);
       dptr->idunion.asa.groupc=grptrc;
      } else if (!strcmp(tok, "OUTPUT")) {
       dptr->idunion.asa.output=1;
       gettok();
       if ((strstr(tok, ".pdb")) == NULL)
         inerr("setup_asa: wrong filename: must has suffix .pdb:", tok);
       dptr->idunion.asa.fpdb=fopen(tok,"w");
      } else break;
 }
}


typedef struct child {
 _REAL x,y,z;
 int   lab[8];
} childstruct;

int test(_REAL x,_REAL y,_REAL z,_REAL sigmac,
         _REAL xa,_REAL ya,_REAL za,_REAL sigma,_REAL cell);

do_asa(idstruct *dptr)
{
 int i,j,k,l,m,j3,kk,nx,ny,nz;
 _REAL safe,cell,sigmac,sigma,x,y,z,xa,ya,za,xc,yc,zc;
 int ical,cntrl,igood,n_good,nn_good,pass;
 childstruct *good,*new_good;
 int *cubo;
/* definisco i formati output dell PDB di controllo
   della superficie                                 */
 char pdbform[]="ATOM  %5d  O   WAT %5d    %8.3f%8.3f%8.3f%8.3f\nTER\n";
 char pdberr[]="FAILED! Resolution %5.3f\n";
 idstruct *gptr,*gptrc;
 gptr=dptr->idunion.asa.group;
 gptrc=dptr->idunion.asa.groupc;

 safe=10.;
 cell=2.;
 
 gptr->idunion.group.xmin=1000000.;
 gptr->idunion.group.xmax=-1000000.;
 gptr->idunion.group.ymin=1000000.;
 gptr->idunion.group.ymax=-1000000.;
 gptr->idunion.group.zmin=1000000.;
 gptr->idunion.group.zmax=-1000000.;
 for(l=0;l<gptr->idunion.group.natoms;l++) {
  j3=gptr->idunion.group.AtomXlist[l];
  x=gptr->idunion.group.set->crd[j3++];
  y=gptr->idunion.group.set->crd[j3++];
  z=gptr->idunion.group.set->crd[j3];
  if(x<gptr->idunion.group.xmin) gptr->idunion.group.xmin=x;
  if(x>gptr->idunion.group.xmax) gptr->idunion.group.xmax=x;
  if(y<gptr->idunion.group.ymin) gptr->idunion.group.ymin=y;
  if(y>gptr->idunion.group.ymax) gptr->idunion.group.ymax=y;
  if(z<gptr->idunion.group.zmin) gptr->idunion.group.zmin=z;
  if(z>gptr->idunion.group.zmax) gptr->idunion.group.zmax=z;
 }
 gptr->idunion.group.xmin-=safe; gptr->idunion.group.xmax+=safe;
 gptr->idunion.group.ymin-=safe; gptr->idunion.group.ymax+=safe;
 gptr->idunion.group.zmin-=safe; gptr->idunion.group.zmax+=safe;

 nx = (int)((gptr->idunion.group.xmax-gptr->idunion.group.xmin)/cell);
 ny = (int)((gptr->idunion.group.ymax-gptr->idunion.group.ymin)/cell);
 nz = (int)((gptr->idunion.group.zmax-gptr->idunion.group.zmin)/cell);
 gptr->idunion.group.xmax = gptr->idunion.group.xmin+((_REAL)nx)*cell;
 gptr->idunion.group.ymax = gptr->idunion.group.ymin+((_REAL)ny)*cell;
 gptr->idunion.group.zmax = gptr->idunion.group.zmin+((_REAL)nz)*cell;

 printf("group %s extents:\n",gptr->id_name);
 printf("xmin= %f --- xmax= %f in %d steps\n",
  gptr->idunion.group.xmin,
  gptr->idunion.group.xmax,nx);
 printf("ymin= %f --- ymax= %f in %d steps\n",
  gptr->idunion.group.ymin,
  gptr->idunion.group.ymax,ny);
 printf("zmin= %f --- zmax= %f in %d steps\n",
  gptr->idunion.group.zmin,
  gptr->idunion.group.zmax,nz);
 
 cubo = (int *)get2(nx*ny*nz,sizeof(int));

 ical=0;
 sigmac = 0.5*cell*l_sqrt(3.);
 for(i=0;i<nx*ny*nz;i++,cubo++) {
  *cubo=gptr->idunion.group.natoms;
 }
 cubo-=nx*ny*nz;
 for(i=0,x=gptr->idunion.group.xmin+0.5*cell;i<nx;i++,x+=cell) {
  for(j=0,y=gptr->idunion.group.ymin+0.5*cell;j<ny;j++,y+=cell) {
   for(k=0,z=gptr->idunion.group.zmin+0.5*cell;k<nz;k++,z+=cell,cubo++) {
    for(l=0;l<gptr->idunion.group.natoms;l++) {
     j3=gptr->idunion.group.AtomXlist[l];
     sigma = gptr->idunion.group.vdw[l]+dptr->idunion.asa.probe;
     xa=gptr->idunion.group.set->crd[j3++];
     ya=gptr->idunion.group.set->crd[j3++];
     za=gptr->idunion.group.set->crd[j3];
     igood=test(x,y,z,sigmac,xa,ya,za,sigma,cell);
     if(igood>0) {
      *cubo -= 1;
     } else if(igood<0) {
      *cubo = -1;
      break;
     }
    }
    if(*cubo>0) ical++;
   }
  }
 }
 cubo-=nx*ny*nz;
 n_good = ical;
 good = (childstruct *)get2(n_good,sizeof(childstruct));

 ical=0;
 for(i=0,x=gptr->idunion.group.xmin+0.5*cell;i<nx;i++,x+=cell) {
  for(j=0,y=gptr->idunion.group.ymin+0.5*cell;j<ny;j++,y+=cell) {
   for(k=0,z=gptr->idunion.group.zmin+0.5*cell;k<nz;k++,z+=cell,cubo++) {
    if(*cubo>0) {
     good->x = x;
     good->y = y;
     good->z = z;
     good++;
     ical++;
    }
   }
  }
 }
 good-=n_good;
 cubo-=nx*ny*nz;
 if(ical!=n_good) {
  printf("asa: ical different from n_good.\n");
  exit(1);
 }
 free(cubo);

 pass=0;
 while(0.5*cell>=dptr->idunion.asa.res && pass<dptr->idunion.asa.nmax ) {
  
  pass++;
  cell=cell*0.5;
  
  sigmac=cell*0.5*l_sqrt(3.);
  ical=0;
  for(m=0;m<n_good;m++,good++) {
   cntrl=0;
   xc=good->x; yc=good->y; zc=good->z;
   for(i=0;i<8;i++) {
    good->lab[i]=gptr->idunion.group.natoms;
   }
   for(i=0,x=xc-0.5*cell,kk=0;i<2;i++,x+=cell) {
    for(j=0,y=yc-0.5*cell;j<2;j++,y+=cell) {
     for(k=0,z=zc-0.5*cell;k<2;k++,z+=cell,kk++) {
      for(l=0;l<gptr->idunion.group.natoms;l++) {
       j3=gptr->idunion.group.AtomXlist[l];
       sigma = gptr->idunion.group.vdw[l]+dptr->idunion.asa.probe;
       xa=gptr->idunion.group.set->crd[j3++];
       ya=gptr->idunion.group.set->crd[j3++];
       za=gptr->idunion.group.set->crd[j3];
       igood=test(x,y,z,sigmac,xa,ya,za,sigma,cell);
       if(igood>0) {
        good->lab[kk] -=1;
       } else if(igood<0) {
        good->lab[kk] = -1;
        break;
       }
      }
     }
    }
   }
   for(i=0,kk=0;i<2;i++) {
    for(j=0;j<2;j++) {
     for(k=0;k<2;k++,kk++) {
      if(good->lab[kk]>0) {
       ical++;
       cntrl++;
      }
     }
    }
   }
/*
   if(cntrl==0) {
    printf("asa: cntrl zero in pass %d.\n",pass);
    exit(1);
   }
*/
  }
  good-=n_good;
  nn_good=ical;
/*  printf("allocating new_good of nn_good = %d\n",nn_good); fflush(stdout); */
  new_good = (childstruct *)get2(nn_good,sizeof(childstruct));
  for(m=0;m<n_good;m++,good++) {
   xc=good->x; yc=good->y; zc=good->z;
   for(i=0,x=xc-0.5*cell,kk=0;i<2;i++,x+=cell) {
    for(j=0,y=yc-0.5*cell;j<2;j++,y+=cell) {
     for(k=0,z=zc-0.5*cell;k<2;k++,z+=cell,kk++) {
      if(good->lab[kk]>0) {
       new_good->x = x;
       new_good->y = y;
       new_good->z = z;
       new_good++;
      }
     }
    }
   }
  }
  good-=n_good;
  new_good-=nn_good;
  free(good);
  n_good=nn_good;
  good = (childstruct *)get2(n_good,sizeof(childstruct));
  for(m=0;m<n_good;m++,good++,new_good++) {
   good->x = new_good->x;
   good->y = new_good->y;
   good->z = new_good->z;
  }
  good-=n_good;
  new_good-=n_good;
  free(new_good);

 }

 printf("final pass= %d\n",pass);
 printf("number of points= %d\n",n_good);
 printf("surface(Ang**2) = %f\n",dptr->idunion.asa.f*((_REAL)n_good));
 printf("r eq group %s + r probe (Ang) = %f\n",
  gptr->id_name,
  l_sqrt(dptr->idunion.asa.f*((_REAL)n_good)/(4.*_PI)));
 fflush(stdout);

 cell=cell*0.5;
 sigmac=cell*0.5*l_sqrt(3.);
 
 if(dptr->idunion.asa.close) {
 pass++;
 ical=0;
 for(m=0;m<n_good;m++,good++) {
   cntrl=0;
   x=good->x; y=good->y; z=good->z;
   good->lab[0]=gptrc->idunion.group.natoms;
   for(l=0;l<gptrc->idunion.group.natoms;l++) {
    j3=gptrc->idunion.group.AtomXlist[l];
    sigma = gptrc->idunion.group.vdw[l]+dptr->idunion.asa.probe;
    xa=gptrc->idunion.group.set->crd[j3++];
    ya=gptrc->idunion.group.set->crd[j3++];
    za=gptrc->idunion.group.set->crd[j3];
    igood=test(x,y,z,sigmac,xa,ya,za,sigma,cell);
    if(igood>0) {
     good->lab[0] -=1;
    } else if(igood<0) {
     good->lab[0] = -1;
     break;
    }
   }
   if(good->lab[0]>=0) {
    ical++;
    cntrl++;
   }
/*
   if(cntrl==0) {
    printf("asa: cntrl zero in pass %d.\n",pass);
    exit(1);
   }
*/
 }
 good-=n_good;
 nn_good=ical;
 new_good = (childstruct *)get2(nn_good,sizeof(childstruct));
 for(m=0;m<n_good;m++,good++) {
  x=good->x; y=good->y; z=good->z;
  if(good->lab[0]>=0) {
   new_good->x = x;
   new_good->y = y;
   new_good->z = z;
   new_good++;
  }
 }
 good-=n_good;
 new_good-=nn_good;
 free(good);
 n_good=nn_good;
 good = (childstruct *)get2(n_good,sizeof(childstruct));
 for(m=0;m<n_good;m++,good++,new_good++) {
  good->x = new_good->x;
  good->y = new_good->y;
  good->z = new_good->z;
 }
 good-=n_good;
 new_good-=n_good;

 }

 printf("final pass= %d\n",pass);
 printf("ASA: number of points= %d\n",n_good);
 printf("ASA: surface(Ang**2) = %f\n",dptr->idunion.asa.f*((_REAL)n_good));
 printf("ASA: r eq group %s + r probe (Ang) = %f\n",
  gptr->id_name,
  l_sqrt(dptr->idunion.asa.f*((_REAL)n_good)/(4.*_PI)));
 fflush(stdout);

/* scrive file PDB di controllo della superficie */
 if (dptr->idunion.asa.output) {
   if (cell<0.04) {
     fprintf(dptr->idunion.asa.fpdb,pdberr,cell);
   } else 
   {
     for(m=0;m<n_good;m++,good++) {
       xc=good->x; yc=good->y; zc=good->z;
       fprintf(dptr->idunion.asa.fpdb,pdbform,m,m,xc,yc,zc,sigmac);
     }
     good-=n_good;
   }
   rewind(dptr->idunion.asa.fpdb);
   fclose(dptr->idunion.asa.fpdb);
 }

 free(good);
 free(new_good);
}

int test(_REAL x,_REAL y,_REAL z,_REAL sigmac,
         _REAL xa,_REAL ya,_REAL za,_REAL sigma,_REAL cell)
{
 _REAL d,sum,diff,dx0,dy0,dz0,dx,dy,dz,delta,
       xt[8],yt[8],zt[8],uxt[6],uyt[6],uzt[6],dmin,dmax,dper,Rf,Rc,
       dt;
 int igood,imin,k,min,max;
/*  possibile specificare un delta solo per
     l'ultimo passo  */
 delta=0.00*cell;
 
 dx0=x-xa; dy0=y-ya; dz0=z-za;
 d=l_sqrt(dx0*dx0+dy0*dy0+dz0*dz0);
 sum=sigma+sigmac;
 diff=Nabs(sigma-sigmac);
 igood=0;
 if (d<=sum) {
  if(d<diff) {
   igood=-1;
   return(igood);
  } else {
   xt[0]=x+cell*0.5; yt[0]=y+cell*0.5; zt[0]=z-cell*0.5;
   xt[1]=x-cell*0.5; yt[1]=y+cell*0.5; zt[1]=z-cell*0.5;
   xt[2]=x-cell*0.5; yt[2]=y-cell*0.5; zt[2]=z-cell*0.5;
   xt[3]=x+cell*0.5; yt[3]=y-cell*0.5; zt[3]=z-cell*0.5;
   xt[4]=x+cell*0.5; yt[4]=y+cell*0.5; zt[4]=z+cell*0.5;
   xt[5]=x-cell*0.5; yt[5]=y+cell*0.5; zt[5]=z+cell*0.5;
   xt[6]=x-cell*0.5; yt[6]=y-cell*0.5; zt[6]=z+cell*0.5;
   xt[7]=x+cell*0.5; yt[7]=y-cell*0.5; zt[7]=z+cell*0.5;
   for(k=0,imin=0;k<8;k++) {
    dx=xt[k]-xa; dy=yt[k]-ya; dz=zt[k]-za;
    dt=l_sqrt(dx*dx+dy*dy+dz*dz);
    if(dt<sigma) imin++;
   }
   if(imin==0) {
    xt[0]=x-cell*0.5; yt[0]=y; zt[0]=z;
    xt[1]=x+cell*0.5; yt[1]=y; zt[1]=z;
    xt[2]=x; yt[2]=y-0.5*cell; zt[2]=z;
    xt[3]=x; yt[3]=y+0.5*cell; zt[3]=z;
    xt[4]=x; yt[4]=y; zt[4]=z-0.5*cell;
    xt[5]=x; yt[5]=y; zt[5]=z+0.5*cell;
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
    dper=dmin-0.5*cell;
    Rc = l_sqrt(sigma*sigma - dper*dper);
    Rf = 0.5*cell;
    dper = dx0*uxt[min]+dy0*uyt[min]+dz0*uzt[min];
    dt = l_sqrt(d*d - dper*dper);
    if(dt<=Rc+Rf) {
     igood=0;
    } else {
     igood=1;
    }
    return(igood);
   } else if(imin==8) {
    igood=-1;
    return(igood);
   }
   igood=0;
   return(igood);
  }
 } else {
  if(d>sum+delta) {
   igood=1;
   return(igood);
  }
 }
 return(igood);
}
