#include "head.h"

// extern char             tok[];
// char             tok[];
extern idstruct         *startid(),*getstruct();
extern parmstruct *prm;

int *ivector();
short *ivector2();
void Maps();

/* this function is used in Map to give the cell index */
int ICELL( int ix,  int iy,  int iz , int Mx, int My, int Mz)
{
  int index;
  index= 1 + (ix-1+Mx)%Mx + ((iy-1+My)%My)*Mx  + ((iz-1+Mz)%Mz)*Mx*My;
  return (index);
}

/*** Map and List setup *****************************************/
void init_map()
{
 idstruct	*dptr,*serptr;
 char		id[TOKLEN];

 dptr = startid(MAP, id);

 dptr->idunion.map.Nx = 0;
 dptr->idunion.map.Ny = 0;
 dptr->idunion.map.Nz = 0;
 dptr->idunion.map.Mx = 3;
 dptr->idunion.map.My = 3;
 dptr->idunion.map.Mz = 3;
 dptr->idunion.map.NN = 342;
 dptr->idunion.map.ncell = 0;

 dptr->idunion.map.List = ivector(1,prm->Nsit);

 gettok();
 if ((serptr = getstruct(tok)) != NULL) {
         if (serptr->id_type != STREAM)
                 inerr("setup_map: expected stream id\n","");
 } else {
         printf("setup_map: id %s not found in idlist\n", tok);
 }
 dptr->idunion.map.set = serptr->idunion.ser.set;

 gettok();
 if (tok[0] != ';')
 {
  printf("error reading setup_map options.\n"); 
  exit(1);
 }
}

void setup_map(idstruct *dptr)
{
 int ix,iy,iz,imap,k;
 _REAL d,dx,dy,dz;
 void extents(_REAL *r, _REAL *box, parmstruct *prm, idstruct *ptr);

 extents(dptr->idunion.map.set->crd,dptr->idunion.map.set->boxc,
         prm,dptr);

  dptr->idunion.map.rcut = 10.;
  dptr->idunion.map.cell = 2.;
  dptr->idunion.map.cellx = dptr->idunion.map.cell;
  dptr->idunion.map.celly = dptr->idunion.map.cell;
  dptr->idunion.map.cellz = dptr->idunion.map.cell;
  dptr->idunion.map.cellxm1 = 1./dptr->idunion.map.cellx;
  dptr->idunion.map.cellym1 = 1./dptr->idunion.map.celly;
  dptr->idunion.map.cellzm1 = 1./dptr->idunion.map.cellz;
  dptr->idunion.map.Nx=
   (int)(dptr->idunion.map.box[0]/ dptr->idunion.map.cell);
  dptr->idunion.map.Ny=
   (int)(dptr->idunion.map.box[1]/ dptr->idunion.map.cell);
  dptr->idunion.map.Nz=
   (int)(dptr->idunion.map.box[2]/ dptr->idunion.map.cell);
  dptr->idunion.map.ncell = dptr->idunion.map.Nx*
                            dptr->idunion.map.Ny*
                            dptr->idunion.map.Nz;

  printf("Nx= %d, Ny= %d, Nz= %d, ncell= %d\n",
   dptr->idunion.map.Nx, 
   dptr->idunion.map.Ny, 
   dptr->idunion.map.Nz,dptr->idunion.map.ncell);
  printf("Mx= %d, My= %d, Mz= %d\n",
   dptr->idunion.map.Mx, dptr->idunion.map.My, dptr->idunion.map.Mz);
  fflush(stdout);

  imap=0;
  k=0;

  for (ix=-dptr->idunion.map.Mx;ix<=dptr->idunion.map.Mx;ix++) {
   for (iy=-dptr->idunion.map.My;iy<=dptr->idunion.map.My;iy++) {
    for (iz=-dptr->idunion.map.Mz;iz<=dptr->idunion.map.Mz;iz++,k++) {
     if(!(ix==0 && iy==0 && iz==0)) 
     {
     dx=dptr->idunion.map.cellx*(_REAL)ix;
     dy=dptr->idunion.map.celly*(_REAL)iy;
     dz=dptr->idunion.map.cellz*(_REAL)iz;
     d=l_sqrt(dx*dx+dy*dy+dz*dz);
     d-=dptr->idunion.map.cell*l_sqrt(3.);
     if(d<=dptr->idunion.map.rcut)
     {
      imap++;
     }
     }
    }
   }
  }

  dptr->idunion.map.NN = imap;
  printf("NN= %d\n", dptr->idunion.map.NN);
  dptr->idunion.map.mapsize= dptr->idunion.map.NN * dptr->idunion.map.ncell;
  dptr->idunion.map.offset = ivector(0,3*dptr->idunion.map.NN-1);
  dptr->idunion.map.Map = ivector(1,dptr->idunion.map.mapsize);
  dptr->idunion.map.Head = ivector(1,dptr->idunion.map.ncell);

  imap=0;
  k=0;

  for (ix=-dptr->idunion.map.Mx;ix<=dptr->idunion.map.Mx;ix++) {
   for (iy=-dptr->idunion.map.My;iy<=dptr->idunion.map.My;iy++) {
    for (iz=-dptr->idunion.map.Mz;iz<=dptr->idunion.map.Mz;iz++,k++) {
     if(!(ix==0 && iy==0 && iz==0)) 
     {
     dx=dptr->idunion.map.cellx*(_REAL)ix;
     dy=dptr->idunion.map.celly*(_REAL)iy;
     dz=dptr->idunion.map.cellz*(_REAL)iz;
     d=l_sqrt(dx*dx+dy*dy+dz*dz);
     d-=dptr->idunion.map.cell*l_sqrt(3.);
     if(d<=dptr->idunion.map.rcut)
     {
      imap++;
      *dptr->idunion.map.offset++ = ix;
      *dptr->idunion.map.offset++ = iy;
      *dptr->idunion.map.offset++ = iz;
     }
     }
    }
   }
  }
  dptr->idunion.map.offset -= 3*dptr->idunion.map.NN;

/* create neighbors list: this needs to be done once only in the program if constant V */
  Maps( dptr,
       dptr->idunion.map.Nx,dptr->idunion.map.Ny,dptr->idunion.map.Nz,
       dptr->idunion.map.Map,dptr->idunion.map.NN,
       dptr->idunion.map.Mx,dptr->idunion.map.My,dptr->idunion.map.Mz);

 printf("esco da Maps.\n"); fflush(stdout);
}


/*****************************************************************
 MAPS: sets up an indexed spatial decomposition map  
 for use in direct_sum force calcs
*****************************************************************/

void Maps( idstruct *dptr, int Nx, int Ny, int Nz,int *Map, int NN, 
           int Mx, int My, int Mz ) 
{
  int ICELL(int ix, int iy, int iz , int Nx, int Ny, int Nz);
  int ix,iy,iz,imap,k,kk;
  int x,y,z;
 
  for (z=1;z<=Nz;z++) {
   for (y=1;y<=Ny;y++) {
    for (x=1;x<=Nx;x++)   {
     k=1;
     imap = (ICELL(x,y,z,Nx,Ny,Nz) - 1 ) * NN;
     /* first set */
     for (kk=1;kk<=NN;kk++,k++) {
      ix=*dptr->idunion.map.offset++;
      iy=*dptr->idunion.map.offset++;
      iz=*dptr->idunion.map.offset++;
      Map[ imap + k] =  ICELL( ix+x, iy+y, iz+z, Nx,Ny,Nz);
     }
     dptr->idunion.map.offset-=3*NN;
    }
   }
  }
}


/****************************************************************/
/* extents: puts the molecule in a box in 0->Lx... */
void extents(_REAL *r, _REAL *box, parmstruct *prm, idstruct *dptr)
{
 _REAL xmax,ymax,zmax,xmin,ymin,zmin,xC,yC,zC,x,y,z,boxx,boxy,boxz,
       boxxh,boxyh,boxzh;
 int i,quoz,resto,ibox;

 xmax=-10000.; xmin=10000.;
 ymax=-10000.; ymin=10000.;
 zmax=-10000.; zmin=10000.;
 xC=0.; yC=0.; zC=0.;

 for(i=0;i<prm->Natom;i++)
 {
   x=*r++; y=*r++; z=*r++;
   xC+=x; yC+=y; zC+=z;
   if(x>xmax) xmax=x;
   if(x<xmin) xmin=x;
   if(y>ymax) ymax=y;
   if(y<ymin) ymin=y;
   if(z>zmax) zmax=z;
   if(z<zmin) zmin=z;
 }
 r-=3*prm->Natom;
 xC/=prm->Natom;
 yC/=prm->Natom;
 zC/=prm->Natom;
 boxx=xmax-xmin+12.;
 boxy=ymax-ymin+12.;
 boxz=zmax-zmin+12.;
 boxxh=boxx*0.5; boxyh=boxy*0.5; boxzh=boxz*0.5;

 for(i=0;i<prm->Natom;i++)
 {
   x=*r++; y=*r++; z=*r++; r-=3;
   x-=xC; y-=yC; z-=zC;
   x+=boxxh; y+=boxyh; z+=boxzh;
   *r++ = x; *r++ = y; *r++ = z;
 }
 r-=3*prm->Natom;
 ibox=(int) (boxx);
 quoz=ibox/2;
 resto=ibox%2;
 if(resto>0) quoz+=1;
 ibox=quoz*2;
 boxx=(_REAL)ibox;
 ibox=(int) (boxy);
 quoz=ibox/2;
 resto=ibox%2;
 if(resto>0) quoz+=1;
 ibox=quoz*2;
 boxy=(_REAL)ibox;
 ibox=(int) (boxz);
 quoz=ibox/2;
 resto=ibox%2;
 if(resto>0) quoz+=1;
 ibox=quoz*2;
 boxz=(_REAL)ibox;

 dptr->idunion.map.box[0]=boxx;
 dptr->idunion.map.box[1]=boxy;
 dptr->idunion.map.box[2]=boxz;
 printf("new box: %f %f %f\n",dptr->idunion.map.box[0],
                             dptr->idunion.map.box[1],
                             dptr->idunion.map.box[2]);
}

/****************************************************************/
/* Link: this sets up the head of the chain and linked list */
void Link(idstruct *dptr, idstruct *serptr)
{

  _REAL *r=serptr->idunion.ser.set->crd;
  parmstruct *prm=serptr->idunion.ser.prm;

  int icell, i, Mx=dptr->idunion.map.Nx,
                My=dptr->idunion.map.Ny,
                Mz=dptr->idunion.map.Nz;

  /* first zero head of chain array */
  for (icell=1;icell<=dptr->idunion.map.ncell;icell++) {
    dptr->idunion.map.Head[icell]= -1;
  }
  for (i=0;i<prm->Nsit;i++) {
    dptr->idunion.map.List[i]= -1;
  }
  
  /* sort all atoms */
  for (i=0; i< prm->Nsit;i++){
    icell= 1 + (int)( (r[0]*dptr->idunion.map.cellxm1) ) + 
               (int)( (r[1]*dptr->idunion.map.cellym1) ) * Mx +  
               (int)( (r[2]*dptr->idunion.map.cellzm1) ) * Mx * My;
/*
  printf("Link:rx... %e %e %e\n",r[0],r[1],r[2]);
  printf("Link:ix... %d\n",icell);
*/
    r+=3;
    dptr->idunion.map.List[i] = dptr->idunion.map.Head[icell];
    dptr->idunion.map.Head[icell]= i;
/*
    printf("Link:cell,i,l,h: %d %d %d %d\n",
     icell, i,dptr->idunion.map.List[i], dptr->idunion.map.Head[icell]);
*/
  }
  r-=3*(prm->Nsit);
  
/*
  {
   int ii,kk,jp,jv,j;
   for(i=1,kk=0,jp=0,jv=0;i<=dptr->idunion.map.ncell;i++) {
    ii=dptr->idunion.map.Head[i];
    if(ii>-1) {
     kk++;
     jp++;
     j = dptr->idunion.map.List[ii];
     while(j>-1) {
      kk++;
      j = dptr->idunion.map.List[j];
     }
    } else jv++;
   }
   printf("kk= %d jp= %d jv= %d\n",kk,jp,jv); fflush(stdout);
  }
*/

} /* end link */
/*****************************************************************************/
