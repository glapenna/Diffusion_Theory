#include "head.h"

// extern char tok[];
// char tok[];
// void *get2(int,int);
extern idstruct         *startid();
extern parmstruct *prm;

void setup_asaf()
{
 idstruct *dptr, *mapptr;
 char id[TOKLEN];
 int i;

 dptr = startid(ASAF, id);

 dptr->idunion.asaf.fsurf = 0.;
 dptr->idunion.asaf.probe = 0.;
 dptr->idunion.asaf.res = 0.05;
 dptr->idunion.asaf.nmax = 10;
 dptr->idunion.asaf.cell = 2.;

 gettok();
 tokREAL(&dptr->idunion.asaf.probe,tok);
 gettok();
 tokREAL(&dptr->idunion.asaf.res,tok);
 gettok();
 tokint(&dptr->idunion.asaf.nmax,tok);
 gettok();
 tokREAL(&dptr->idunion.asaf.fsurf,tok);
 gettok();
 tokREAL(&dptr->idunion.asaf.cell,tok);

 for (;;) {
  gettok();
  if (tok[0]==';') {
   break;
  }
 }

 dptr->idunion.asaf.sigma=(_REAL *)get2(prm->Natom,sizeof(_REAL));
 dptr->idunion.asaf.vdw=(_REAL *)get2(prm->Natom,sizeof(_REAL));
 dptr->idunion.asaf.vdw2=(_REAL *)get2(prm->Natom,sizeof(_REAL));
 dptr->idunion.asaf.sum2=(_REAL *)get2(prm->Natom,sizeof(_REAL));
 dptr->idunion.asaf.diff2=(_REAL *)get2(prm->Natom,sizeof(_REAL));
 for(i=0;i<prm->Natom;i++)
 {
  printf ("i= %d\n",i); fflush(stdout);
  dptr->idunion.asaf.vdw[i]=2.0;
  dptr->idunion.asaf.sigma[i]=0.;
  dptr->idunion.asaf.vdw[i]+=dptr->idunion.asaf.probe;
  dptr->idunion.asaf.vdw2[i]=dptr->idunion.asaf.vdw[i]*dptr->idunion.asaf.vdw[i];
  dptr->idunion.asaf.sigma[i]*=dptr->idunion.asaf.fsurf;
 }

}


typedef struct child {
 _REAL x,y,z;
 int nlist;
 int *list;
 struct child  *childs;
} childstruct;

surf_mol(idstruct *dptr, idstruct *serptr, idstruct *mapptr)
{
 void print_cube();
 _REAL *ri=serptr->idunion.ser.set->crd;
 _REAL cell=dptr->idunion.asaf.cell,cellh,sqrt_three;
 int i,j,k,l,m,kk,n,nx,ny,nz,atom,maxlist;
 _REAL safe,sigmac,sigma,x,y,z,xa,ya,za,xc,yc,zc,
       xmin,xmax,ymin,ymax,zmin,zmax,f,deltaG,*rt;
 int ical,cntrl,igood,n_good,nn_good,pass,
     icell,jcell,iat,nabor;
 int *cubo_list;
 childstruct *good,*new_good;
 int *cubo;
 _REAL sum2,diff2,sigma2;
 _REAL time0,time1,time;

 sqrt_three=l_sqrt(3.);
 xmin=0.;
 ymin=0.;
 zmin=0.;
 nx=mapptr->idunion.map.Nx;
 ny=mapptr->idunion.map.Ny;
 nz=mapptr->idunion.map.Nz;
 maxlist=50;
 printf("creo spazio per cubo...\n"); fflush(stdout);
 cubo = (int *)get2(nx*ny*nz,sizeof(int));
 cubo_list = (int *)get2(maxlist*nx*ny*nz,sizeof(int));
 printf("ho creato spazio per cubo.\n"); fflush(stdout);

 ical=0;
 cellh=0.5*cell;
 sigmac = cellh*sqrt_three;
 for(i=0;i<prm->Natom;i++) {
  sigma=dptr->idunion.asaf.vdw[i];
  sum2=sigma+sigmac;
  sum2*=sum2;
  dptr->idunion.asaf.sum2[i]=sum2;
  diff2=sigma-sigmac;
  diff2*=diff2;
  dptr->idunion.asaf.diff2[i]=diff2;
 }

 for(i=0;i<nx*ny*nz;i++,cubo++) {
  for(j=0;j<maxlist;j++)
  {
   *cubo_list++ = -1;
  }
  *cubo=0;
 }
 cubo-=nx*ny*nz;
 cubo_list-=maxlist*nx*ny*nz;

 for(i=1;i<=nx*ny*nz;i++,cubo++) {
  iat=mapptr->idunion.map.Head[i];
  while(iat>-1)
  {
   (*cubo)++;
   iat=mapptr->idunion.map.List[iat];
  }
  for (nabor=1; nabor<=mapptr->idunion.map.NN; nabor ++)
  {
   jcell= mapptr->idunion.map.Map[ (mapptr->idunion.map.NN)*(i-1) + nabor];
   iat = mapptr->idunion.map.Head[jcell];
   while(iat>-1)
   {
    (*cubo)++;
    iat=mapptr->idunion.map.List[iat];
   }
  }
 }
 cubo-=nx*ny*nz;

/*
 icell=1;
 for(k=0,z=zmin+cellh;k<nz;k++,z+=cell)
 {
  for(j=0,y=ymin+cellh;j<ny;j++,y+=cell)
  {
   for(i=0,x=xmin+cellh;i<nx;i++,x+=cell,cubo++,icell++)
   {
    if(*cubo>0)
    {
     iat=*cubo;
     print_cube(&icell,&x,&y,&z,&iat);
    }
   }
  }
 }
 cubo-=nx*ny*nz;
*/

 icell=1;
 for(k=0,z=zmin+cellh;k<nz;k++,z+=cell)
 {
  for(j=0,y=ymin+cellh;j<ny;j++,y+=cell)
  {
   for(i=0,x=xmin+cellh;i<nx;i++,x+=cell,cubo++,cubo_list+=maxlist,icell++)
   {
    kk=0;
    iat = mapptr->idunion.map.Head[icell];
    while(iat>-1)
    {
/*
 if(icell==0) {
  printf("icell,iat= %d %d\n",icell,iat); fflush(stdout);
 }
*/
     xa=ri[3*iat]; ya=ri[3*iat+1]; za=ri[3*iat+2];
     sigma=dptr->idunion.asaf.vdw[iat];
     sigma2=dptr->idunion.asaf.vdw2[iat];
     sum2=dptr->idunion.asaf.sum2[iat];
     diff2=dptr->idunion.asaf.diff2[iat];
#include "test.c"
     /* igood=test(x,y,z,sigmac,xa,ya,za,sigma,cell,cellh); */
     /* igood=test2(x,y,z,sigmac,xa,ya,za,sigma,cell); */
     if(igood>0) 
     {
      *cubo -= 1;
     }
     else if(igood==0) 
     {
      cubo_list[kk] = iat;
      kk++;
      if(kk>maxlist)
      {
       printf("maxlist exceeded.\n");
       exit(1);
      }
     }
     else if(igood<0) 
     {
      *cubo = -1;
      break;
     }
     iat=mapptr->idunion.map.List[iat];
    }
    for (nabor=1; nabor<=mapptr->idunion.map.NN; nabor ++)
    {
     jcell= mapptr->idunion.map.Map[ (mapptr->idunion.map.NN)*(icell-1) + nabor];
     iat = mapptr->idunion.map.Head[jcell];
     while(iat>-1)
     {
/*
 if(icell==0) {
  printf("nabor,jcell,iat= %d %d %d\n",nabor,jcell,iat); fflush(stdout);
 }
*/
      xa=ri[3*iat]; ya=ri[3*iat+1]; za=ri[3*iat+2];
      sigma=dptr->idunion.asaf.vdw[iat];
      sigma2=dptr->idunion.asaf.vdw2[iat];
      sum2=dptr->idunion.asaf.sum2[iat];
      diff2=dptr->idunion.asaf.diff2[iat];
#include "test.c"
      /* igood=test(x,y,z,sigmac,xa,ya,za,sigma,cell,cellh); */
      /* igood=test2(x,y,z,sigmac,xa,ya,za,sigma,cell); */
      if(igood>0) 
      {
       *cubo -= 1;
      }
      else if(igood==0) 
      {
       cubo_list[kk] = iat;
       kk++;
       if(kk>maxlist)
       {
        printf("maxlist exceeded.\n");
        exit(1);
       }
      } 
      else if(igood<0) 
      {
       *cubo = -1;
       break;
      }
      iat=mapptr->idunion.map.List[iat];
     }
    }
    if(*cubo>0) ical++;
   }
  }
 }
 cubo-=nx*ny*nz;
 cubo_list-=maxlist*nx*ny*nz;

 printf("ho finito di creare cubo.\n"); fflush(stdout);
/*
 icell=1;
 for(k=0,z=zmin+cellh;k<nz;k++,z+=cell)
 {
  for(j=0,y=ymin+cellh;j<ny;j++,y+=cell)
  {
   for(i=0,x=xmin+cellh;i<nx;i++,x+=cell,cubo++,icell++)
   {
    if(*cubo>0)
    {
     iat=*cubo;
     print_cube(&icell,&x,&y,&z,&iat);
    }
   }
  }
 }
 cubo-=nx*ny*nz;
*/

 n_good = ical;
 good = (childstruct *)get2(n_good,sizeof(childstruct));

 ical=0;
 kk=0;
 for(k=0,z=zmin+cellh;k<nz;k++,z+=cell) 
 {
  for(j=0,y=ymin+cellh;j<ny;j++,y+=cell) 
  {
   for(i=0,x=xmin+cellh;i<nx;i++,x+=cell,
        cubo++,cubo_list+=maxlist,kk++) 
   {
    if(*cubo>0) {
     good->x = x;
     good->y = y;
     good->z = z;
     good->nlist = *cubo;
     good->list = (int *)get2(good->nlist,sizeof(int));
 /* printf("cubo(%d,%d,%d)= %d\n",i,j,k,*cubo); fflush(stdout); */
    /* print_cube(&kk,&x,&y,&z); */
     for(l=0;l<good->nlist;l++)
     {
      good->list[l] = cubo_list[l];
/*
      printf("kk,ical,x,y,z,nlist,l= %d %d %f %f %f %d %d\n",
       kk,ical,x,y,z,good->nlist,good->list[l]);
*/
     }
     good++;
     ical++;
    }
   }
  }
 }
 good-=n_good;
 cubo-=nx*ny*nz;
 cubo_list-=maxlist*nx*ny*nz;
 free(cubo);
 free(cubo_list);
 if(ical!=n_good) {
  printf("asa: ical different from n_good.\n");
  exit(1);
 }

 printf("entro nel ciclo di pass.\n"); fflush(stdout);
 time0=clock();
 time1=clock();
 time=time1-time0;
 printf("time= %e\n",time); fflush(stdout);
 time0=time1;
 pass=0;
/*
 while(cellh>=dptr->idunion.asaf.res && pass<dptr->idunion.asaf.nmax ) {
*/
 while( pass<dptr->idunion.asaf.nmax ) {
  
  pass++;
  cell=cell*0.5;
  cellh=cell*0.5;
  sigmac=cellh*sqrt_three;
  for(i=0;i<prm->Natom;i++) {
   sigma=dptr->idunion.asaf.vdw[i];
   sum2=sigma+sigmac;
   sum2*=sum2;
   dptr->idunion.asaf.sum2[i]=sum2;
   diff2=sigma-sigmac;
   diff2*=diff2;
   dptr->idunion.asaf.diff2[i]=diff2;
  }

  ical=0;
  for(m=0;m<n_good;m++,good++) {
   cntrl=0;
   xc=good->x; yc=good->y; zc=good->z;
   good->childs = (childstruct *)get2(8,sizeof(childstruct));
   for(i=0;i<8;i++) {
    good->childs[i].nlist = 0;
   }

   for(i=0;i<8;i++) {
    good->childs[i].nlist = good->nlist;
    good->childs[i].list = (int *)get2(good->nlist,sizeof(int));
    for(j=0;j<good->nlist;j++) {
     good->childs[i].list[j] = good->list[j];
    }
   }

   for(k=0,z=zc-cellh,kk=0;k<2;k++,z+=cell) 
   {
    for(j=0,y=yc-cellh;j<2;j++,y+=cell) 
    {
     for(i=0,x=xc-cellh;i<2;i++,x+=cell,kk++) 
     {
      n=0;
      for(l=0;l<good->nlist;l++) 
      {
       atom=good->list[l];
       xa=ri[3*atom]; ya=ri[3*atom+1]; za=ri[3*atom+2];
       sigma=dptr->idunion.asaf.vdw[atom];
       sigma2=dptr->idunion.asaf.vdw2[atom];
       sum2=dptr->idunion.asaf.sum2[atom];
       diff2=dptr->idunion.asaf.diff2[atom];
#include "test.c"
       /* igood=test(x,y,z,sigmac,xa,ya,za,sigma,cell,cellh); */
       /* igood=test2(x,y,z,sigmac,xa,ya,za,sigma,cell); */
       if(igood>0) 
       {
        good->childs[kk].nlist -=1;
       } 
       else if(igood==0) 
       {
        good->childs[kk].list[n]=atom;
        n++;
       }
       else if(igood<0) 
       {
        good->childs[kk].nlist = -1;
        break;
       }
      }
     }
    }
   }
   for(k=0,kk=0;k<2;k++) {
    for(j=0;j<2;j++) {
     for(i=0;i<2;i++,kk++) {
      if(good->childs[kk].nlist>0) {
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
   for(k=0,z=zc-cellh,kk=0;k<2;k++,z+=cell) {
    for(j=0,y=yc-cellh;j<2;j++,y+=cell) {
     for(i=0,x=xc-cellh;i<2;i++,x+=cell,kk++) {
      if(good->childs[kk].nlist>0) {
       new_good->x = x;
       new_good->y = y;
       new_good->z = z;
       new_good->nlist = good->childs[kk].nlist;
       new_good->list = (int *)get2(good->childs[kk].nlist,sizeof(int));
       for(l=0;l<good->childs[kk].nlist;l++) {
        new_good->list[l]=good->childs[kk].list[l];
       }
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
   good->nlist = new_good->nlist;
   good->list = (int *)get2(new_good->nlist,sizeof(int));
   for(l=0;l<new_good->nlist;l++) {
    good->list[l]=new_good->list[l];
   }
  }
  good-=n_good;
  new_good-=n_good;
  free(new_good);

 time1=clock();
 time=time1-time0;
 printf("time= %e\n",time); fflush(stdout);
 time0=time1;

 } /* endwhile */

 printf("Final pass= %d\n",pass);
 printf("Final cell= %e\n",cell);
 printf("Number of points= %d\n",n_good);
 printf("Molecular surface (Ang.^2)= %f\n",dptr->idunion.asaf.fsurf*n_good);

 deltaG=0.;
 for(m=0;m<n_good;m++,good++) {
/* scommentare per scrivere pdb della superficie */
/*  print_cube(&m,good);  */
  for(l=0;l<good->nlist;l++) {
   j=good->list[l];
   deltaG+=dptr->idunion.asaf.sigma[j]/(good->nlist);
  }
 }
 good-=n_good;

 printf("Ricetta 3: DeltaG (kJ/mol)= %f\n",deltaG);
 fflush(stdout);

 free(good);

}

void print_cube(int *i, childstruct *good)
{
 int k;
 char pdbform2[]="ATOM  %5d  XE  XE  %5d    %8.3f%8.3f%8.3f%8.3f\nTER\n";

 k=*i;
 printf(pdbform2,k+1,(*i)+1,good->x,
                            good->y,
                            good->z);

}
