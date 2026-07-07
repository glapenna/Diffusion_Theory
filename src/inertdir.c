#include "head.h"

FILE *outdir=NULL;

void
inertdir(idstruct *serptr)
{
/*************************************************************************
 Evaluate inertia, shape tensors, director and order.
 Molecules must be non-truncated by PBC!!! Apply sputinbox first!
*************************************************************************/

	_REAL   *r = serptr->idunion.ser.set->crd;
	_REAL   *box = serptr->idunion.ser.set->boxc;
	int     nmol = serptr->idunion.ser.prm->Nmol;
	_REAL   facmol,facat;
	int     natom = serptr->idunion.ser.prm->Natom;
	int     nsit = serptr->idunion.ser.prm->Nsit;
	_REAL   *mass = serptr->idunion.ser.prm->Masses;

	int i,j,ia,im,ia3;
        _REAL x,y,z,m1,xcm,ycm,zcm,masstot,ax,ay,az,xi,yi,zi,proj,dirl[3];
        _REAL **In,**evec,*eval,**Qzz,**Qyy,**Qxx,**S;
        _REAL **matrix(),*vector();
        void eigsrt(),jacobi();
        void euler(_REAL uz[],_REAL ut[],_REAL ux[],
                   _REAL uy[],_REAL ang[]);
        void proquat(_REAL q1[],_REAL q2[],_REAL q3[]);
	_REAL ux[3],uy[3],uz[3],ut[3],utmp[4],q[4],qt[4],
	      qtmp[4],ang[3],beh,alefp,alefm,
	      Asph,Ksq_num,Ksq_den,Ksq,trace,average,S_num,S_den,S_asph,
	      Delta_num,Delta_den,Delta,x1,y1,z1;

        printf("iset= %d\n",serptr->idunion.ser.iset);
        if(serptr->idunion.ser.iset==0)
	{
         In = matrix(1,3,1,3);
         S = matrix(1,3,1,3);
         Qzz = matrix(1,3,1,3);
         Qyy = matrix(1,3,1,3);
         Qxx = matrix(1,3,1,3);
         evec = matrix(1,3,1,3);
         eval = vector(1,3);
	}
        facmol = 1./(_REAL)(nmol);
        facat = 1./(_REAL)(natom);
        dirl[0]=0.;
        dirl[1]=0.;
        dirl[2]=0.;

        for(i=1;i<=3;i++) for(j=i;j<=3;j++) Qzz[i][j] = 0.;
        for(i=1;i<=3;i++) for(j=i;j<=3;j++) Qyy[i][j] = 0.;
        for(i=1;i<=3;i++) for(j=i;j<=3;j++) Qxx[i][j] = 0.;

        for(im=0; im<nmol; im++) 
	{
/* computes the center of mass */
         xcm=0.; ycm=0.; zcm=0.; masstot=0.;
        /* for(ia=0; ia<21; ia++,r+=3) {} */
         for(ia=0,ia3=0; ia<natom; ia++,ia3+=3) 
	 {
          x1 = r[ia3];
          y1 = r[ia3+1];
          z1 = r[ia3+2];
	  m1=mass[ia];
// if ignoring masses
	  m1=1.;
          masstot += m1;
          xcm += m1*x1;
          ycm += m1*y1;
          zcm += m1*z1;
/*      printf("mass= %f",mass[ia]);*/
         }
         xcm/=masstot;
         ycm/=masstot;
         zcm/=masstot;

/* computes the inertia and shape tensors */
         for(i=1;i<=3;i++) 
	 {
	  for(j=i;j<=3;j++) 
	  {
	   In[i][j] = 0.;
	   S[i][j] = 0.;
	  }
	 }
         for(ia=0,ia3=0; ia<natom; ia++,ia3+=3) 
	 {
          x1 = r[ia3];
          y1 = r[ia3+1];
          z1 = r[ia3+2];
          x = x1-xcm;
          y = y1-ycm;
          z = z1-zcm;
	  m1=mass[ia];
// if ignoring mass
	  m1=1.;
          In[1][1] += m1*(y*y+z*z);
          In[1][2] -= m1*x*y;
          In[1][3] -= m1*x*z;
          In[2][2] += m1*(x*x+z*z);
          In[2][3] -= m1*y*z;
          In[3][3] += m1*(x*x+y*y);
          S[1][1] += x*x;
          S[1][2] += x*y;
          S[1][3] += x*z;
          S[2][2] += y*y;
          S[2][3] += y*z;
          S[3][3] += z*z;
         }
         In[2][1]=In[1][2];
         In[3][1]=In[1][3];
         In[3][2]=In[2][3];
         S[2][1]=S[1][2];
         S[3][1]=S[1][3];
         S[3][2]=S[2][3];
         for(i=1;i<=3;i++) 
	 {
	  for(j=i;j<=3;j++) 
	  {
	   S[i][j] *= facat;
	  }
	 }
/* diagonalise the Shape tensor */
         jacobi(S,3,eval,evec);
/* diagonalise the Inertia tensor */
         eigsrt(eval,evec,3);
	 Asph=eval[3]-0.5*(eval[1]+eval[2]);
	 Ksq_num=(eval[1]*eval[2]+eval[2]*eval[3]+eval[3]*eval[1]);
	 Ksq_den=(eval[3]+eval[2]+eval[1])*(eval[3]+eval[2]+eval[1]);
	 Ksq=1.-3.*Ksq_num/Ksq_den;
      printf ("Shape1= %e %e %e\n",eval[1],eval[2],eval[3]);
      printf ("Shape2= %e %e\n",Asph,Ksq);
/* The Shape tensor of Dima04 is the same of S */
	 trace=(eval[1]+eval[2]+eval[3]);
	 average=trace/3.;
	 Delta_num = 1.5*((eval[1]-average)*(eval[1]-average)+
	                  (eval[2]-average)*(eval[2]-average)+
	                  (eval[3]-average)*(eval[3]-average));
	 Delta_den = trace*trace;
	 Delta = Delta_num/Delta_den;
	 S_num = 27.*(eval[1]-average)*(eval[2]-average)*(eval[3]-average);
	 S_den = Delta_den*trace;
	 S_asph = S_num/S_den;
      printf ("Dima04= %e %e\n",Delta,S_asph);
/* allign the z axis of Inertia with the atomic order in the mol */
         jacobi(In,3,eval,evec);
         eigsrt(eval,evec,3);
      printf ("Inert= %e %e %e\n",eval[1],eval[2],eval[3]);
         xi=r[3*natom-3] - xcm;
         yi=r[3*natom-2] - ycm;
         zi=r[3*natom-1] - zcm;
         proj = evec[1][3]*xi+evec[2][3]*yi+evec[3][3]*zi;
         if (proj < 0.) 
	 {
          evec[1][3] = -evec[1][3];
          evec[2][3] = -evec[2][3];
          evec[3][3] = -evec[3][3];
          evec[1][1] = -evec[1][1];
          evec[2][1] = -evec[2][1];
          evec[3][1] = -evec[3][1];
         }

         ax=evec[1][3];
         ay=evec[2][3];
         az=evec[3][3];
         dirl[0]+=ax;
         dirl[1]+=ay;
         dirl[2]+=az;

         uz[0]=evec[1][3];
         uz[1]=evec[2][3];
         uz[2]=evec[3][3];
         ut[0]=evec[1][1];
         ut[1]=evec[2][1];
         ut[2]=evec[3][1];
         euler(uz,ut,ux,uy,ang);
         beh = 0.5*ang[1];
         alefp = 0.5*(ang[0]+ang[2]);
         alefm = 0.5*(ang[0]-ang[2]);
         q[0] =  cos(beh)*cos(alefp);
         q[1] =  sin(beh)*cos(alefm);
         q[2] =  sin(beh)*sin(alefm);
         q[3] =  cos(beh)*sin(alefp);
         qt[0]= q[0];
         qt[1]=-q[1];
         qt[2]=-q[2];
         qt[3]=-q[3];

         ut[0]=0.;
         for(ia=0,ia3=0;ia<natom;ia++,ia3+=3)
         {
          ut[1] = r[ia3] -xcm;
          ut[2] = r[ia3+1] -ycm;
          ut[3] = r[ia3+2] -zcm;
          proquat(q,ut,qtmp);
          proquat(qtmp,qt,utmp);
          r[ia3]=utmp[1];
          r[ia3+1]=utmp[2];
          r[ia3+2]=utmp[3];
         }
         ang[0]=3.0*_PI/4.0;
         ang[1]=acos(sqrt(3.0)/3.0);
         ang[2]=0.0;
         beh = 0.5*ang[1];
         alefp = 0.5*(ang[0]+ang[2]);
         alefm = 0.5*(ang[0]-ang[2]);
         q[0] =  cos(beh)*cos(alefp);
         q[1] =  sin(beh)*cos(alefm);
         q[2] =  sin(beh)*sin(alefm);
         q[3] =  cos(beh)*sin(alefp);
         qt[0]= q[0];
         qt[1]=-q[1];
         qt[2]=-q[2];
         qt[3]=-q[3];
         ut[0]=0.;
         for(ia=0,ia3=0;ia<natom;ia++,ia3+=3)
         {
          ut[1] = r[ia3];
          ut[2] = r[ia3+1];
          ut[3] = r[ia3+2];
          proquat(qt,ut,qtmp);
          proquat(qtmp,q,utmp);
          r[ia3]=utmp[1];
          r[ia3+1]=utmp[2];
          r[ia3+2]=utmp[3];
         }
/*
      printf ("a3:\n");
      printf("%e %e %e\n",ax,ay,az);
*/
         Qzz[1][1]+=ax*ax;
         Qzz[1][2]+=ax*ay;
         Qzz[1][3]+=ax*az;
         Qzz[2][2]+=ay*ay;
         Qzz[2][3]+=ay*az;
         Qzz[3][3]+=az*az;

         ax=evec[1][2];
         ay=evec[2][2];
         az=evec[3][2];
         Qyy[1][1]+=ax*ax;
         Qyy[1][2]+=ax*ay;
         Qyy[1][3]+=ax*az;
         Qyy[2][2]+=ay*ay;
         Qyy[2][3]+=ay*az;
         Qyy[3][3]+=az*az;
         ax=evec[1][3];
         ay=evec[2][3];
         az=evec[3][3];
         Qxx[1][1]+=ax*ax;
         Qxx[1][2]+=ax*ay;
         Qxx[1][3]+=ax*az;
         Qxx[2][2]+=ay*ay;
         Qxx[2][3]+=ay*az;
         Qxx[3][3]+=az*az;

         r+=3*natom;
        }
        r-=3*nsit;

        dirl[0] *= facmol;
        dirl[1] *= facmol;
        dirl[2] *= facmol;
      printf ("P1= %+12.6e %+12.6e %+12.6e\n",dirl[0],dirl[1],dirl[2]);
        ax=1./l_sqrt(dirl[0]*dirl[0]+dirl[1]*dirl[1]+dirl[2]*dirl[2]);
        dirl[0] *= ax;
        dirl[1] *= ax;
        dirl[2] *= ax;
        Qzz[1][1] = 0.5*(3.*facmol*Qzz[1][1]-1.);
        Qzz[1][2] = 0.5*(3.*facmol*Qzz[1][2]);
        Qzz[1][3] = 0.5*(3.*facmol*Qzz[1][3]);
        Qzz[2][2] = 0.5*(3.*facmol*Qzz[2][2]-1.);
        Qzz[2][3] = 0.5*(3.*facmol*Qzz[2][3]);
        Qzz[3][3] = 0.5*(3.*facmol*Qzz[3][3]-1.);
        Qyy[1][1] = 0.5*(3.*facmol*Qyy[1][1]-1.);
        Qyy[1][2] = 0.5*(3.*facmol*Qyy[1][2]);
        Qyy[1][3] = 0.5*(3.*facmol*Qyy[1][3]);
        Qyy[2][2] = 0.5*(3.*facmol*Qyy[2][2]-1.);
        Qyy[2][3] = 0.5*(3.*facmol*Qyy[2][3]);
        Qyy[3][3] = 0.5*(3.*facmol*Qyy[3][3]-1.);
        Qxx[1][1] = 0.5*(3.*facmol*Qxx[1][1]-1.);
        Qxx[1][2] = 0.5*(3.*facmol*Qxx[1][2]);
        Qxx[1][3] = 0.5*(3.*facmol*Qxx[1][3]);
        Qxx[2][2] = 0.5*(3.*facmol*Qxx[2][2]-1.);
        Qxx[2][3] = 0.5*(3.*facmol*Qxx[2][3]);
        Qxx[3][3] = 0.5*(3.*facmol*Qxx[3][3]-1.);

        for(i=1;i<=3;i++) for(j=i+1;j<=3;j++) Qzz[j][i] = Qzz[i][j];
        for(i=1;i<=3;i++) for(j=i+1;j<=3;j++) Qyy[j][i] = Qyy[i][j];
        for(i=1;i<=3;i++) for(j=i+1;j<=3;j++) Qxx[j][i] = Qxx[i][j];


/* diagonalise the Saupe matrix */
      jacobi(Qzz,3,eval,evec);
      eigsrt(eval,evec,3);
/* store the director */
      for(i=1;i<=3;i++) 
       serptr->idunion.ser.set->dir[i-1] = evec[i][1];
      for(i=1;i<=3;i++) 
       serptr->idunion.ser.set->dir[i-1] = dirl[i-1];

      for(i=1;i<=3;i++) for(j=1;j<=3;j++) 
       serptr->idunion.ser.set->rotdir[i-1][j-1] = evec[i][j];
      serptr->idunion.ser.set->P2 = eval[1];

      printf ("ordering matrix eigenvalues:\n");
      printf ("Sxx= %+12.6e  Syy= %+12.6e  Szz= %+12.6e\n",
               eval[1], eval[2], eval[3]);
      printf ("director:\n");
      printf ("nx= %+12.6e  ny= %+12.6e  nz= %+12.6e\n",
       serptr->idunion.ser.set->dir[0],
       serptr->idunion.ser.set->dir[1],
       serptr->idunion.ser.set->dir[2]);

/*
      {
       float xt[3];
       extern FILE *outdir;
       if(outdir==NULL) {
        outdir=fopen("dir.out","w");
        if(outdir==NULL) {
         printf("error opening dir.out\n");
         exit(1);
        }
       }
       xt[0] = (float) serptr->idunion.ser.set->dir[0];
       xt[1] = (float) serptr->idunion.ser.set->dir[1];
       xt[2] = (float) serptr->idunion.ser.set->dir[2];
       printf("writing director on to %d\n",outdir);
       fwrite(xt,3*sizeof(float),1,outdir);
       fflush(outdir);
      }
*/
/*
      printf ("P2: %e\n",serptr->idunion.ser.set->P2);
*/

}


void
orientmol(_REAL *crd,int natom)
{
   void euler(_REAL uz[],_REAL ut[],_REAL ux[],
        _REAL uy[],_REAL ang[]);
   void proquat(_REAL q1[],_REAL q2[],_REAL q3[]);

   int j1,j2,j3,i,j,k;
   _REAL scp,ux[3],uy[3],uz[3],ut[3],xc,yc,zc,
              ang[3],beh,alefp,alefm,
              q[4],qt[4],qtmp[4],utmp[4],
              ux0[]={0.,1.,0.,0.},
              uy0[]={0.,0.,1.,0.},
              uz0[]={0.,0.,0.,1.};

   j1=27;
   j2=28;
   j3=29;
   j1--; j2--; j3--;
   j1*=3; j2*=3; j3*=3;

   xc=crd[j2];
   yc=crd[j2+1];
   zc=crd[j2+2];
   uz[0]=crd[j3]-crd[j2];
   uz[1]=crd[j3+1]-crd[j2+1];
   uz[2]=crd[j3+2]-crd[j2+2];
   ut[0]=crd[j1]-crd[j2];
   ut[1]=crd[j1+1]-crd[j2+1];
   ut[2]=crd[j1+2]-crd[j2+2];
   scp=1./l_sqrt(uz[0]*uz[0]+uz[1]*uz[1]+uz[2]*uz[2]);
   for(i=0;i<3;i++) uz[i]*=scp;
   scp=1./l_sqrt(ut[0]*ut[0]+ut[1]*ut[1]+ut[2]*ut[2]);
   for(i=0;i<3;i++) ut[i]*=scp;
   euler(uz,ut,ux,uy,ang);

   beh = 0.5*ang[1];
   alefp = 0.5*(ang[0]+ang[2]);
   alefm = 0.5*(ang[0]-ang[2]);
   q[0] =  cos(beh)*cos(alefp);
   q[1] =  sin(beh)*cos(alefm);
   q[2] =  sin(beh)*sin(alefm);
   q[3] =  cos(beh)*sin(alefp);
   qt[0]= q[0];
   qt[1]=-q[1];
   qt[2]=-q[2];
   qt[3]=-q[3];

   ut[0]=0.;
   for(i=0;i<natom;i++)
   {
    ut[1]=crd[0]-xc;
    ut[2]=crd[1]-yc;
    ut[3]=crd[2]-zc;
    proquat(q,ut,qtmp);
    proquat(qtmp,qt,utmp);
    crd[0]=utmp[1];
    crd[1]=utmp[2];
    crd[2]=utmp[3];
    crd+=3;
  }
  crd-=3*natom;

}
