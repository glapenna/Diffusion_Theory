#include "head.h"

/**** extern ****/
// extern char             tok[];
// char             tok[];
// void *get2(int,int);
extern idstruct         *startid(),*getstruct();

/************************************************************************/
void init_rig()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(RIG, id);

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("RIG: expected file name","");
	strcpy(dptr->idunion.rig.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.rig.fname, 0, &exists, "RIG"))
		exit(1);
	if ((dptr->idunion.rig.file = 
				  fopen(dptr->idunion.rig.fname, "w")) == NULL){
		perror(dptr->idunion.rig.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

void setup_rig()
{
	idstruct	*dptr,*bptr;
	int i,j,k;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != RIG)
                        inerr("setup_rig: expected id or keyword\n","");
        } else {
                printf("setup_rig: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((bptr = getstruct(tok)) == NULL) 
                inerr("setup_rig: no such id\n","");
        if (bptr->id_type != BEAD)
                inerr("setup_rig: bead id expected\n","");
        gettok();

        dptr->idunion.rig.bead = bptr;
        dptr->idunion.rig.type = 0;
}

void rotdiff(idstruct *dptr, idstruct *serptr)
{
   extern parmstruct *prm;
   _REAL *box=serptr->idunion.ser.set->boxc;
   _REAL *boxm1=serptr->idunion.ser.set->boxm1;
   idstruct *bptr=dptr->idunion.rig.bead;
   int i,j,k,l,i3,j3,nn,N,rank,max,overlap;
   _REAL **r1, **r2;
   _REAL *rt1,*rt2,*z,zi,zj,zij,zim1,rm1,rsq,rm2,rmod,aa,bb,zsq,Rm1,zav,zr,
         cutsq;
   _REAL x1,y1,z1,x2,y2,z2,dx,dy,dz,f[3][3];
   rdiadstruct D[1];
   _REAL **Q,**u,**uu,*ut,*v,*w,**S,**SS,tol,sum,zsimm;
   _REAL **Kt,**KOr,**KOc,**KRr,*ROR,**KRc,**St,**Stt,**Tt;
   FILE *file=dptr->idunion.rig.file;
   _REAL Rm1av,Rm1avz,zi3,zi2,RKR,RKRZ,weight;

   _REAL C1[4][4],D1[4][4],D2[6][6],C2[6][6],lambda1[4],lambda2[6];
   int n1,n2;

   _REAL Dxx,Dyy,Dzz,Dp,Dm,Dtr,Delta;
   void euler();
   _REAL dirx[3],diry[3],dirz[3],vx[3],vy[3],ang[3],betah,
         alphap,alpham;

  printf("sono in rotdiff.\n"); fflush(stdout);

   zav=bptr->idunion.bead.zav;
   z=bptr->idunion.bead.z;
   zr=bptr->idunion.bead.zr;
   cutsq=bptr->idunion.bead.cutsq;
   r1=bptr->idunion.bead.r;
   r2=bptr->idunion.bead.r;

   printf("rotdiff: building Q matrix...\n"); fflush(stdout);

   nn=3*bptr->idunion.bead.nbead;
   Q = matrix(1,nn,1,nn);
   S = matrix(1,nn,1,nn);
   for(i=1;i<=nn;i++) {
    for(j=1;j<=nn;j++) {
     Q[i][j] = 0.;
     S[i][j] = 0.;
    }
   }

   Rm1av=0.;
   Rm1avz=0.;
   zi3=0.;
   zi2=0.;
   for (i=1,i3=3;i<=bptr->idunion.bead.nbead;i++,i3+=3,r1++) {
    rt1=*r1;
    zi = z[i-1];
    zi3+=zi*zi*zi;
    zi2+=zi*zi;
    zim1 = 1./zi;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
    for (j=i,j3=i3,r2=r1;j<=bptr->idunion.bead.nbead;j++,j3+=3,r2++) {

    zj=z[j-1];
    rt2=*r2;
    x2=rt2[0]; y2=rt2[1]; z2=rt2[2];

    dx=x1-x2;
    dy=y1-y2;
    dz=z1-z2;
#ifdef PBC
     dx -= box[0]*rint(dx*boxm1[0]);
     dy -= box[1]*rint(dy*boxm1[1]);
     dz -= box[2]*rint(dz*boxm1[2]);
#endif
     if(i == j) {
      f[0][0] = zim1; f[1][1] = zim1; f[2][2] = zim1;
      f[0][1] = 0.; f[0][2] = 0.; f[1][2] = 0.;
     } else {
#include "hydrod.c"
      f[0][0] = D->xx;
      f[0][1] = D->xy;
      f[0][2] = D->xz;
      f[1][1] = D->yy;
      f[1][2] = D->yz;
      f[2][2] = D->zz;
      if(overlap) printf("spheres %d and %d overlap\n",i,j);
      Rm1av+=2.*zi*zj*rm1;
      Rm1avz+=2.*zi*zi*zj*zj*rm1;
     }

     zsimm=l_sqrt(zi*zj);
     Q[i3-2][j3-2] = zsimm*f[0][0];
     Q[i3-2][j3-1] = zsimm*f[0][1];
     Q[i3-2][j3  ] = zsimm*f[0][2];
     Q[i3-1][j3-1] = zsimm*f[1][1];
     Q[i3-1][j3  ] = zsimm*f[1][2];
     Q[i3  ][j3  ] = zsimm*f[2][2];
     Q[i3-1][j3-2] = zsimm*f[0][1];
     Q[i3  ][j3-2] = zsimm*f[0][2];
     Q[i3  ][j3-1] = zsimm*f[1][2];
     Q[j3-2][i3-2] = Q[i3-2][j3-2];
     Q[j3-2][i3-1] = Q[i3-2][j3-1];
     Q[j3-2][i3  ] = Q[i3-2][j3  ];
     Q[j3-1][i3-2] = Q[i3-1][j3-2];
     Q[j3-1][i3-1] = Q[i3-1][j3-1];
     Q[j3-1][i3  ] = Q[i3-1][j3  ];
     Q[j3  ][i3-2] = Q[i3  ][j3-2];
     Q[j3  ][i3-1] = Q[i3  ][j3-1];
     Q[j3  ][i3  ] = Q[i3  ][j3  ];
    } /* enddo j */
   } /* enddo i */

   RKR = ((_REAL)bptr->idunion.bead.nbead)*zav/
          (1. + Rm1av/((_REAL)bptr->idunion.bead.nbead)/zav);
   RKRZ = zi2*zi2/ (zi3 + Rm1avz);
   fprintf(file,"RKR (Ang.)= %e, RKRZ = %e\n",prm->sig0r*RKR,prm->sig0r*RKRZ);

   r1=bptr->idunion.bead.r;
   r2=bptr->idunion.bead.r;

   printf("rotdiff: inverting Q matrix...\n"); fflush(stdout);

/*
*/
#ifdef IMSL
 tol=TOLSVD;
 ut = (_REAL *)get2(nn*nn,sizeof(_REAL));
 for(i=1,k=0;i<=nn;i++) {
  for(j=1;j<=nn;j++,k++) {
   ut[k]=Q[i][j];
  }
 }
 w = imsl_d_lin_svd_gen(nn,nn,ut,
       IMSL_RANK, tol, &rank,
       IMSL_INVERSE,&v,0);
 free(w);
 for(i=1,k=0;i<=nn;i++) {
  for(j=1;j<=nn;j++,k++) {
   S[i][j]=v[k];
  }
 }
 free(v);
 free(ut);
#else
 rank=invert(Q,nn,nn,S);
#endif

 printf("rotdiff: rank of u: %d\n",rank); fflush(stdout);

 free_matrix(Q,1,nn,1,nn);

   printf("rotdiff: building K matrices...\n"); fflush(stdout);

 N=3;
 Kt = matrix(1,N,1,N);
 KOc = matrix(1,N,1,N);
 KRc = matrix(1,N,1,N);
 KOr = matrix(1,N,1,N);
 KRr = matrix(1,N,1,N);
 St = matrix(1,N,1,N);
 Tt = matrix(1,N,1,N);
 ROR = vector(1,N);
 for(i=1;i<=3;i++) {
  for(j=1;j<=3;j++) {
   Kt[i][j] = 0.;
   KOc[i][j] = 0.;
   KRc[i][j] = 0.;
   KOr[i][j] = 0.;
   KRr[i][j] = 0.;
  }
 }
 for(i=1,i3=3;i<=bptr->idunion.bead.nbead;i++,i3+=3,r1++) {
  rt1=*r1;
  x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
  for(j=1,j3=3;j<=bptr->idunion.bead.nbead;j++,j3+=3,r2++) {
   zim1=l_sqrt(z[i-1]*z[j-1]);
   rt2=*r2;
   x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
/*
  printf("r1: %e %e %e\n",x1,y1,z1);
  printf("r2: %e %e %e\n",x2,y2,z2);
*/
   St[1][1] = S[i3-2][j3-2];
   St[1][2] = S[i3-2][j3-1];
   St[1][3] = S[i3-2][j3  ];
   St[2][1] = S[i3-1][j3-2];
   St[2][2] = S[i3-1][j3-1];
   St[2][3] = S[i3-1][j3  ];
   St[3][1] = S[i3  ][j3-2];
   St[3][2] = S[i3  ][j3-1];
   St[3][3] = S[i3  ][j3  ];
   for(k=1;k<=3;k++)
    for(l=1;l<=3;l++)
     Kt[k][l] += zim1*St[k][l];
   Tt[1][1] = y1*St[3][1] - z1*St[2][1];
   Tt[1][2] = y1*St[3][2] - z1*St[2][2];
   Tt[1][3] = y1*St[3][3] - z1*St[2][3];
   Tt[2][1] = z1*St[1][1] - x1*St[3][1];
   Tt[2][2] = z1*St[1][2] - x1*St[3][2];
   Tt[2][3] = z1*St[1][3] - x1*St[3][3];
   Tt[3][1] = x1*St[2][1] - y1*St[1][1];
   Tt[3][2] = x1*St[2][2] - y1*St[1][2];
   Tt[3][3] = x1*St[2][3] - y1*St[1][3];
   for(k=1;k<=3;k++)
    for(l=1;l<=3;l++)
     KOc[k][l] += zim1*Tt[k][l];
   St[1][1] = Tt[1][2]*z2 -Tt[1][3]*y2;
   St[1][2] = Tt[1][3]*x2 -Tt[1][1]*z2;
   St[1][3] = Tt[1][1]*y2 -Tt[1][2]*x2;
   St[2][1] = Tt[2][2]*z2 -Tt[2][3]*y2;
   St[2][2] = Tt[2][3]*x2 -Tt[2][1]*z2;
   St[2][3] = Tt[2][1]*y2 -Tt[2][2]*x2;
   St[3][1] = Tt[3][2]*z2 -Tt[3][3]*y2;
   St[3][2] = Tt[3][3]*x2 -Tt[3][1]*z2;
   St[3][3] = Tt[3][1]*y2 -Tt[3][2]*x2;
   for(k=1;k<=3;k++)
    for(l=1;l<=3;l++)
     KOr[k][l] -= zim1*St[k][l];
  }
  r2-=bptr->idunion.bead.nbead;
 }
 r1-=bptr->idunion.bead.nbead;

 free_matrix(S,1,nn,1,nn);

/*
 printf("KOc:\n");
 for(i=1;i<=3;i++) {
  for(j=1;j<=3;j++) printf("%f ",KOc[i][j]);
  printf("\n");
 }
*/

 u = matrix(1,N,1,N);

 u[1][1] =  Kt[2][2]+Kt[3][3];
 u[1][2] = -Kt[1][2];
 u[1][3] = -Kt[1][3];
 u[2][1] = -Kt[2][1];
 u[2][2] =  Kt[1][1]+Kt[3][3];
 u[2][3] = -Kt[2][3];
 u[3][1] = -Kt[3][1];
 u[3][2] = -Kt[3][2];
 u[3][3] =  Kt[1][1]+Kt[2][2];

#ifdef IMSL
 tol=TOLSVD;
 ut = (_REAL *)get2(N*N,sizeof(_REAL));
 for(i=1,k=0;i<=N;i++)
  for(j=1;j<=N;j++,k++)
   ut[k]=u[i][j];
 w = imsl_d_lin_svd_gen(N,N,ut,
       IMSL_RANK, tol, &rank,
       IMSL_INVERSE,&v,0);
 free(w);
 for(i=1,k=0;i<=N;i++) {
  for(j=1;j<=N;j++,k++) {
   St[i][j]=v[k];
  }
 }
 free(v);
 free(ut);
#else
 rank=invert(u,N,N,St);
#endif

 printf("rotdiff: rank of u: %d\n",rank);

 for(i=1;i<=3;i++) {
  ROR[i]=0.;
 }
 r1=bptr->idunion.bead.r;
 for(i=0;i<nn/3;i++,r1++) {
  rt1=*r1;
  x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
  ROR[1]+=x1; ROR[2]+=y1; ROR[3]+=z1;
 }
 for(i=1;i<=3;i++) {
  ROR[i]/=(_REAL)nn;
 }
 r1=bptr->idunion.bead.r;

 w = vector(1,N);
 w[1]=KOc[3][2]-KOc[2][3];
 w[2]=KOc[1][3]-KOc[3][1];
 w[3]=KOc[2][1]-KOc[1][2];
/*
*/
 for(i=1;i<=3;i++) {
  for(j=1,sum=0.;j<=3;j++) {
   sum+=St[i][j]*w[j];
  }
  ROR[i]=sum;
 }
 for(i=1;i<=3;i++) {
  for(j=1;j<=3;j++) {
   KRr[i][j]=KOr[i][j];
  }
 }

 x1=ROR[1]; y1=ROR[2]; z1=ROR[3];
 fprintf(file,"ROR = %e %e %e\n",x1,y1,z1);
 Tt[1][1] = y1*Kt[3][1] - z1*Kt[2][1];
 Tt[1][2] = y1*Kt[3][2] - z1*Kt[2][2];
 Tt[1][3] = y1*Kt[3][3] - z1*Kt[2][3];
 Tt[2][1] = z1*Kt[1][1] - x1*Kt[3][1];
 Tt[2][2] = z1*Kt[1][2] - x1*Kt[3][2];
 Tt[2][3] = z1*Kt[1][3] - x1*Kt[3][3];
 Tt[3][1] = x1*Kt[2][1] - y1*Kt[1][1];
 Tt[3][2] = x1*Kt[2][2] - y1*Kt[1][2];
 Tt[3][3] = x1*Kt[2][3] - y1*Kt[1][3];
/*
 for(i=1;i<=3;i++)
  for(j=1;j<=3;j++)
   KRc[i][j] = KOc[i][j]-Tt[i][j];
*/
 St[1][1] = Tt[1][2]*z1 - Tt[1][3]*y1;
 St[1][2] = Tt[1][3]*x1 - Tt[1][1]*z1;
 St[1][3] = Tt[1][1]*y1 - Tt[1][2]*x1;
 St[2][1] = Tt[2][2]*z1 - Tt[2][3]*y1;
 St[2][2] = Tt[2][3]*x1 - Tt[2][1]*z1;
 St[2][3] = Tt[2][1]*y1 - Tt[2][2]*x1;
 St[3][1] = Tt[3][2]*z1 - Tt[3][3]*y1;
 St[3][2] = Tt[3][3]*x1 - Tt[3][1]*z1;
 St[3][3] = Tt[3][1]*y1 - Tt[3][2]*x1;
 for(i=1;i<=3;i++)
  for(j=1;j<=3;j++)
   KRr[i][j] -= St[i][j];
 Tt[1][1] = KOc[1][2]*z1 - KOc[1][3]*y1;
 Tt[1][2] = KOc[1][3]*x1 - KOc[1][1]*z1;
 Tt[1][3] = KOc[1][1]*y1 - KOc[1][2]*x1;
 Tt[2][1] = KOc[2][2]*z1 - KOc[2][3]*y1;
 Tt[2][2] = KOc[2][3]*x1 - KOc[2][1]*z1;
 Tt[2][3] = KOc[2][1]*y1 - KOc[2][2]*x1;
 Tt[3][1] = KOc[3][2]*z1 - KOc[3][3]*y1;
 Tt[3][2] = KOc[3][3]*x1 - KOc[3][1]*z1;
 Tt[3][3] = KOc[3][1]*y1 - KOc[3][2]*x1;
 for(i=1;i<=3;i++)
  for(j=1;j<=3;j++)
   KRr[i][j] += Tt[i][j];
/*
 Tt[1][1] = y1*KOc[3][1] - z1*KOc[2][1];
 Tt[1][2] = z1*KOc[1][1] - x1*KOc[3][1];
 Tt[1][3] = x1*KOc[2][1] - y1*KOc[1][1];
 Tt[2][1] = y1*KOc[3][2] - z1*KOc[2][2];
 Tt[2][2] = z1*KOc[1][2] - x1*KOc[3][2];
 Tt[2][3] = x1*KOc[2][2] - y1*KOc[1][2];
 Tt[3][1] = y1*KOc[3][3] - z1*KOc[2][3];
 Tt[3][2] = z1*KOc[1][3] - x1*KOc[3][3];
 Tt[3][3] = x1*KOc[2][3] - y1*KOc[1][3];
 for(i=1;i<=3;i++)
  for(j=1;j<=3;j++)
   KRr[i][j] += Tt[i][j];
*/
 Tt[1][1] = y1*KOc[1][3] - z1*KOc[1][2];
 Tt[1][2] = y1*KOc[2][3] - z1*KOc[2][2];
 Tt[1][3] = y1*KOc[3][3] - z1*KOc[3][2];
 Tt[2][1] = z1*KOc[1][1] - x1*KOc[1][3];
 Tt[2][2] = z1*KOc[2][1] - x1*KOc[2][3];
 Tt[2][3] = z1*KOc[3][1] - x1*KOc[3][3];
 Tt[3][1] = x1*KOc[1][2] - y1*KOc[1][1];
 Tt[3][2] = x1*KOc[2][2] - y1*KOc[2][1];
 Tt[3][3] = x1*KOc[3][2] - y1*KOc[3][1];
 for(i=1;i<=3;i++)
  for(j=1;j<=3;j++)
   KRr[i][j] -= Tt[i][j];

/*
 printf("KRr:\n");
 for(i=1;i<=3;i++) {
  for(j=1;j<=3;j++) printf("%f ",KRr[i][j]);
  printf("\n");
 }
 printf("KRc:\n");
 for(i=1;i<=3;i++) {
  for(j=1;j<=3;j++) printf("%f ",KRc[i][j]);
  printf("\n");
 }
*/

 free_matrix(St,1,N,1,N);
 free_matrix(u,1,N,1,N);
 free_vector(w,1,N);
 N=6;
 Stt = matrix(1,N,1,N);
 uu = matrix(1,N,1,N);

 for(i=1;i<=3;i++) {
  for(j=1;j<=3;j++) {
   uu[i][j]=KRr[i][j];
   uu[i][3+j]=KRc[i][j];
   uu[3+i][j]=KRc[j][i];
   uu[3+i][3+j]=Kt[i][j];
  }
 }

#ifdef IMSL
 tol=TOLSVD;
 ut = (_REAL *)get2(N*N,sizeof(_REAL));
 for(i=1,k=0;i<=N;i++)
  for(j=1;j<=N;j++,k++)
   ut[k]=uu[i][j];
 w = imsl_d_lin_svd_gen(N,N,ut,
       IMSL_RANK, tol, &rank,
       IMSL_INVERSE,&v,0);
 free(w);
 for(i=1,k=0;i<=N;i++) {
  for(j=1;j<=N;j++,k++) {
   Stt[i][j]=v[k];
  }
 }
 free(v);
 free(ut);
#else
 rank=invert(uu,N,N,Stt);
#endif

 printf("rotdiff: rank of u: %d\n",rank);
 fprintf(file,"D:\n");
 for(i=1;i<=6;i++) {
  for(j=1;j<=6;j++) fprintf(file,"%f ",Stt[i][j]);
  fprintf(file,"\n");
 }

 w=vector(1,N);
 SS=matrix(1,N,1,N);
 jacobi(Stt,N,w,SS);
/*
*/
 for(i=1;i<=3;i++) {
  if(Nabs(w[i])<=ZERO) w[i]=1.e+08;
 }
 N=3;
 eigsrt(w,SS,N);

 fprintf(file,"D eigenvalues:\n");
 for(i=1;i<=6;i++) {
  fprintf(file,"%e ",w[i]);
 }
 fprintf(file,"\n"); fflush(file);

 fprintf(file,"DR= %e\n",bptr->idunion.bead.sigma*(w[1]+w[2]+w[3])/3.);
 fprintf(file,"DT= %e\n",bptr->idunion.bead.sigma*prm->sig0r*prm->sig0r*((w[4]+w[5]+w[6])/3.)
              *1.e-20);
 fprintf(file,"Req (Ang.)= %e\n",prm->sig0r/((w[4]+w[5]+w[6])/3.));

 for(j=1;j<=3;j++) {
  for(i=1,sum=0.;i<=3;i++) {
   sum+=SS[i][j]*SS[i][j];
  }
  sum=l_sqrt(sum);
  for(i=1;i<=3;i++) {
   serptr->idunion.ser.set->rotdir[i-1][j-1]=SS[i][j]/sum;
  }
 }
 N=6;
 free_matrix(SS,1,N,1,N);
 free_matrix(Stt,1,N,1,N);
 N=3;
 free_matrix(Kt,1,N,1,N);
 free_matrix(KOc,1,N,1,N);
 free_matrix(KRc,1,N,1,N);
 free_matrix(KOr,1,N,1,N);
 free_matrix(KRr,1,N,1,N);
 free_matrix(Tt,1,N,1,N);
 free_vector(ROR,1,N);

 N=6;
 free_vector(w,1,N);
 fprintf(file,"directors:\n");
 for(i=0;i<3;i++) {
  fprintf(file,"%f ",serptr->idunion.ser.set->rotdir[i][0]);
  serptr->idunion.ser.set->dir[i]=
               serptr->idunion.ser.set->rotdir[i][0];
  dirz[i]= serptr->idunion.ser.set->rotdir[i][0];
 }
 fprintf(file,"\n");
 for(i=0;i<3;i++) {
  fprintf(file,"%f ",serptr->idunion.ser.set->rotdir[i][1]);
  diry[i]= -serptr->idunion.ser.set->rotdir[i][1];
 }
 fprintf(file,"\n");
 for(i=0;i<3;i++) {
  fprintf(file,"%f ",serptr->idunion.ser.set->rotdir[i][2]);
  dirx[i]= serptr->idunion.ser.set->rotdir[i][2];
 }
 fprintf(file,"\n"); fflush(file);

 Dxx=w[3];
 Dyy=w[2];
 Dzz=w[1];
 euler(dirz,diry,vx,vy,ang);
 fprintf(file,
  "alpha= %e, beta= %e, gamma= %e\n",TODEG*ang[0],TODEG*ang[1],TODEG*ang[2]);
 alphap=0.5*(ang[0]+ang[2]);
 alpham=0.5*(ang[0]-ang[2]);
 betah=0.5*ang[1];
 dptr->idunion.rig.qdir[0]=cos(betah)*cos(alphap);
 dptr->idunion.rig.qdir[1]=sin(betah)*cos(alpham);
 dptr->idunion.rig.qdir[2]=sin(betah)*sin(alpham);
 dptr->idunion.rig.qdir[3]=cos(betah)*sin(alphap);

 fprintf(file,"\n rank 1 G rates (Favro):\n");
 fprintf(file,"lambda_1= %e\n",Dxx+Dyy);
 fprintf(file,"lambda_2= %e\n",Dxx+Dzz);
 fprintf(file,"lambda_3= %e\n",Dyy+Dzz);
 fflush(file);

 Dtr = (Dxx+Dyy+Dzz)/3.;
 Delta = l_sqrt(Dxx*Dxx+Dyy*Dyy+Dzz*Dzz-Dxx*Dyy-Dxx*Dzz-Dyy*Dzz);
 fprintf(file,"\n rank 2 G rates (Favro):\n");
 fprintf(file,"lambda_1= %e\n",3.*(Dtr+Dzz));
 fprintf(file,"lambda_2= %e\n",3.*(Dtr+Dxx));
 fprintf(file,"lambda_3= %e\n",3.*(Dtr+Dyy));
 fprintf(file,"lambda_4= %e\n",6.*Dtr+2.*Delta);
 fprintf(file,"lambda_5= %e\n",6.*Dtr-2.*Delta);
 fflush(file);

 Dp = (Dxx+Dyy)/2.;
 Dm = (Dxx-Dyy)/2.;

 n1=3;
 D1[1][1] = 2*Dp+Dzz-Dp;
 D1[1][2] = Dm;
 D1[1][3] = 0.;
 D1[2][1] = D1[1][2];
 D1[2][2] = 2*Dp+Dzz-Dp;
 D1[2][3] = 0.;
 D1[3][1] = D1[1][3];
 D1[3][2] = D1[2][3];
 D1[3][3] = 2*Dp;

/*
*/
 fprintf(file,"\n Diffusion operator on rank 1 basis:\n");
 for(i=1;i<=n1;i++) {
  for(j=1;j<=n1;j++) {
   fprintf(file,"%e ",D1[i][j]);
  }
  fprintf(file,"\n");
 }
 fprintf(file,"\n"); fflush(file);

/*
*/
 { 
   lambda1[1]=Dyy+Dzz;
   lambda1[2]=Dxx+Dzz;
   lambda1[3]=Dxx+Dyy;

   C1[1][1]=-1./l_sqrt(2.);
   C1[2][1]= 1./l_sqrt(2.);
   C1[3][1]= 0.;

   C1[1][2]= 1./l_sqrt(2.);
   C1[2][2]= 1./l_sqrt(2.);
   C1[3][2]= 0.;

   C1[1][3]= 0.;
   C1[2][3]= 0.;
   C1[3][3]= 1.;
 }

  dptr->idunion.rig.eigenvalues1[1]=lambda1[1];
  dptr->idunion.rig.eigenvalues1[3]=lambda1[2];
  dptr->idunion.rig.eigenvalues1[2]=lambda1[3];
  dptr->idunion.rig.eigenstates1[1][1]=C1[1][1];
  dptr->idunion.rig.eigenstates1[3][1]=C1[2][1];
  dptr->idunion.rig.eigenstates1[2][1]=C1[3][1];
  dptr->idunion.rig.eigenstates1[1][2]=C1[1][3];
  dptr->idunion.rig.eigenstates1[3][2]=C1[2][3];
  dptr->idunion.rig.eigenstates1[2][2]=C1[3][3];
  dptr->idunion.rig.eigenstates1[1][3]=C1[1][2];
  dptr->idunion.rig.eigenstates1[3][3]=C1[2][2];
  dptr->idunion.rig.eigenstates1[2][3]=C1[3][2];

 n2=5;
 D2[1][1] = 6*Dp+4.*(Dzz-Dp);
 D2[1][2] = 0.;
 D2[1][3] = l_sqrt(6.)*Dm;
 D2[1][4] = 0.;
 D2[1][5] = 0.;
 D2[2][1] = D2[1][2];
 D2[2][2] = 6*Dp+4.*(Dzz-Dp);
 D2[2][3] = l_sqrt(6.)*Dm;
 D2[2][4] = 0.;
 D2[2][5] = 0.;
 D2[3][1] = D2[1][3];
 D2[3][2] = D2[2][3];
 D2[3][3] = 6*Dp;
 D2[3][4] = 0.;
 D2[3][5] = 0.;
 D2[4][1] = D2[1][4];
 D2[4][2] = D2[2][4];
 D2[4][3] = D2[3][4];
 D2[4][4] = 6.*Dp+Dzz-Dp;
 D2[4][5] = 3.*Dm;
 D2[5][1] = D2[1][5];
 D2[5][2] = D2[2][5];
 D2[5][3] = D2[3][5];
 D2[5][4] = D2[4][5];
 D2[5][5] = 6*Dp+Dzz-Dp;

/*
*/
 fprintf(file,"\n Diffusion operator on rank 2 basis:\n");
 for(i=1;i<=n2;i++) {
  for(j=1;j<=n2;j++) {
   fprintf(file,"%e ",D2[i][j]);
  }
  fprintf(file,"\n");
 }
 fprintf(file,"\n"); fflush(file);

/*
*/
 { _REAL a,b,N,delta,thetaI,sum;
   _REAL tmp1[6][6],tmp2[6][6];
   delta = l_sqrt((Dxx-Dyy)*(Dxx-Dyy)+(Dzz-Dxx)*(Dzz-Dyy));
   a = l_sqrt(3.)*(Dxx-Dyy);
   b = 2.*Dzz - Dxx-Dyy +2.*delta;
   N = 2.*l_sqrt(delta*b);

   lambda2[1]=3.*(Dtr+Dzz);
   lambda2[2]=6.*Dtr+2.*delta;
   lambda2[3]=6.*Dtr-2.*delta;
   lambda2[4]=3.*(Dtr+Dyy);
   lambda2[5]=3.*(Dtr+Dxx);

   C2[1][1]=-1./l_sqrt(2.);
   C2[2][1]= 1./l_sqrt(2.);
   C2[3][1]= 0.;
   C2[4][1]= 0.;
   C2[5][1]= 0.;

   C2[1][2]= b/(l_sqrt(2.)*N);
   C2[2][2]= b/(l_sqrt(2.)*N);
   C2[3][2]= a/N;
   C2[4][2]= 0.;
   C2[5][2]= 0.;

   C2[1][3]= -a/(l_sqrt(2.)*N);
   C2[2][3]= -a/(l_sqrt(2.)*N);
   C2[3][3]=  b/N;
   C2[4][3]= 0.;
   C2[5][3]= 0.;

   C2[1][4]= 0.;
   C2[2][4]= 0.;
   C2[3][4]= 0.;
   C2[4][4]=-1./l_sqrt(2.);
   C2[5][4]= 1./l_sqrt(2.);

   C2[1][5]= 0.;
   C2[2][5]= 0.;
   C2[3][5]= 0.;
   C2[4][5]= 1./l_sqrt(2.);
   C2[5][5]= 1./l_sqrt(2.);

   sum=0.;
   sum+= (a*a/(N*N))/lambda2[2];
   sum+= (b*b/(N*N))/lambda2[3];
   fprintf(file,"Pecora's F00 tau= %e\n",sum);
   sum=0.;
   sum+= (l_sqrt(2.)*a*b/(N*N))/lambda2[2];
   sum-= (l_sqrt(2.)*a*b/(N*N))/lambda2[3];
   fprintf(file,"Pecora's F20 tau= %e\n",sum);
 }

  dptr->idunion.rig.eigenvalues2[1]=lambda2[1];
  dptr->idunion.rig.eigenvalues2[5]=lambda2[2];
  dptr->idunion.rig.eigenvalues2[3]=lambda2[3];
  dptr->idunion.rig.eigenvalues2[2]=lambda2[4];
  dptr->idunion.rig.eigenvalues2[4]=lambda2[5];
  dptr->idunion.rig.eigenstates2[1][1]=C2[1][1];
  dptr->idunion.rig.eigenstates2[5][1]=C2[2][1];
  dptr->idunion.rig.eigenstates2[3][1]=C2[3][1];
  dptr->idunion.rig.eigenstates2[2][1]=C2[4][1];
  dptr->idunion.rig.eigenstates2[4][1]=C2[5][1];
  dptr->idunion.rig.eigenstates2[1][2]=C2[1][4];
  dptr->idunion.rig.eigenstates2[5][2]=C2[2][4];
  dptr->idunion.rig.eigenstates2[3][2]=C2[3][4];
  dptr->idunion.rig.eigenstates2[2][2]=C2[4][4];
  dptr->idunion.rig.eigenstates2[4][2]=C2[5][4];
  dptr->idunion.rig.eigenstates2[1][3]=C2[1][3];
  dptr->idunion.rig.eigenstates2[5][3]=C2[2][3];
  dptr->idunion.rig.eigenstates2[3][3]=C2[3][3];
  dptr->idunion.rig.eigenstates2[2][3]=C2[4][3];
  dptr->idunion.rig.eigenstates2[4][3]=C2[5][3];
  dptr->idunion.rig.eigenstates2[1][4]=C2[1][5];
  dptr->idunion.rig.eigenstates2[5][4]=C2[2][5];
  dptr->idunion.rig.eigenstates2[3][4]=C2[3][5];
  dptr->idunion.rig.eigenstates2[2][4]=C2[4][5];
  dptr->idunion.rig.eigenstates2[4][4]=C2[5][5];
  dptr->idunion.rig.eigenstates2[1][5]=C2[1][2];
  dptr->idunion.rig.eigenstates2[5][5]=C2[2][2];
  dptr->idunion.rig.eigenstates2[3][5]=C2[3][2];
  dptr->idunion.rig.eigenstates2[2][5]=C2[4][2];
  dptr->idunion.rig.eigenstates2[4][5]=C2[5][2];

}

void do_acff_rig(idstruct *acfptr, idstruct *serptr)
{
   extern parmstruct *prm;
   idstruct *modptr=acfptr->idunion.acff.model;
   _REAL *modsq=acfptr->idunion.acff.modsq;
   _REAL **relax_Sz=acfptr->idunion.acff.relax_Sz;
   _REAL **relax_Sxy=acfptr->idunion.acff.relax_Sxy;
   _REAL **relax_Iz_Sz=acfptr->idunion.acff.relax_Iz_Sz;
   idstruct *bptr;
   _REAL *box=serptr->idunion.ser.set->boxc;
   _REAL *qd=modptr->idunion.rig.qdir;

   _REAL factns,A[6][6],omega,taub1[4],taub2[6],*taur_p1,*taur_p2,
         c,d,J_0,J_S,J_I,J_ImS,J_IpS;
   int i,j,k,n1,n2;
   _REAL **r1, **r2;
   _REAL *rt1,*rt2;
   _REAL x1,y1,z1,x2,y2,z2,r[3],rsq,ql[4],qld[4],qdt[4],tau,tau_p1,tau_p2;
   void r2q(), proquat(), qeuler();
   fcomplex d0,d1,d0cc,d1cc,Dt,
     dd1,dd2,dd3,dd4,dd5,dd6,dd7,dd8,dd9,dd10,
     Pak,D[6][6];
   _REAL sqrt_six,sqrt_two,msqrt_two;
   FILE *file=acfptr->idunion.acff.file;
   void genquatch1();
   _REAL rtmp[3];
 
   bptr=modptr->idunion.rig.bead;
   factns = bptr->idunion.bead.factns;
   sqrt_six=l_sqrt(6.);
   sqrt_two=l_sqrt(2.);
   msqrt_two = -sqrt_two;
   n1=3;
   n2=5;
/*
  printf("lambda:\n");
  for(j=1;j<=n1;j++) {
   printf("%e ",modptr->idunion.rig.eigenvalues1[j]);
  }
  printf("\n"); fflush(stdout);
  printf("cc:\n");
  for(j=1;j<=n1;j++) {
   for(k=1;k<=n1;k++) {
    printf("%e ",modptr->idunion.rig.eigenstates1[j][k]);
   }
   printf("\n");
  }
  printf("\n");
*/

   for(k=1;k<=n1;k++) {
    taub1[k] = factns/modptr->idunion.rig.eigenvalues1[k];
   }
   for(k=1;k<=n2;k++) {
    taub2[k] = factns/modptr->idunion.rig.eigenvalues2[k];
   }
   qdt[0]=qd[0]; for(j=1;j<4;j++) qdt[j]= -qd[j];

   taur_p1=vector(1,acfptr->idunion.acff.nbond);
   taur_p2=vector(1,acfptr->idunion.acff.nbond);
   r1=acfptr->idunion.acff.r1;
   r2=acfptr->idunion.acff.r2;
   for (i=0;i<acfptr->idunion.acff.nbond;i++,r1++,r2++) {


    rt1=*r1; rt2=*r2;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
/* add an extra atom with some construction... */
/*
    genquatch1(&rt1[-3],&rt1[0],&rt1[3],rtmp);
    x2=rtmp[0]; y2=rtmp[1]; z2=rtmp[2];
*/
    x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
    r[0]=x2-x1; r[1]=y2-y1; r[2]=z2-z1;
    rsq=r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    rsq=1./rsq;
    modsq[i+1]+=rsq*rsq*rsq;
    /* modsq[i+1]=1./pow(1.02,6.); */
    r2q(r,ql);
    proquat(ql,qdt,qld);

/*
    {
     _REAL qlt[4],rt[4],rtt[4],qldt[4],qtmp[4];
     rt[0]=0.; rt[1]=0.; rt[2]=0.; rt[3]=1.;
     qlt[0]=ql[0]; for(j=1;j<4;j++) qlt[j]=-ql[j];
     proquat(qlt,rt,qtmp);
     proquat(qtmp,ql,rtt);
     printf("legame L->E = %e %e %e\n",rtt[1],rtt[2],rtt[3]);

     qldt[0]=qld[0]; for(j=1;j<4;j++) qldt[j]=-qld[j];
     proquat(qldt,rt,qtmp);
     proquat(qtmp,qld,rtt);
     proquat(qdt,rtt,qtmp);
     proquat(qtmp,qd,rtt);
     printf("legame M->E L->M = %e %e %e\n",rtt[1],rtt[2],rtt[3]);
    }
*/
    { _REAL cb,cbh,ca,sa,cg,sg,sb,a0,a1,a2,p2;

     qeuler(qld,&cbh,&ca,&sa,&cg,&sg);
     cb=2.*cbh*cbh-1.;
     sb=l_sqrt(1.-cb*cb);
     printf("orient(%d)= %f %f %f %f %f\n",i+1,cb,ca,sa,cg,sg);

/*
     a0 = cb*cb;
     a1 = sb*sb;
     tau_p1=0.;
     tau_p1+=a0/modptr->idunion.rig.eigenvalues1[2];
     tau_p1+=a1/modptr->idunion.rig.eigenvalues1[3];
     fprintf(file,"%e ",tau_p1);

     p2 = 0.5*(3.*cb*cb-1.);
     a0 = p2*p2;
     a1 = 3.*cb*sb*cb*sb;
     a2 = 3.*sb*sb*sb*sb/4.;
     tau_p2=0.;
     tau_p2+=a0/modptr->idunion.rig.eigenvalues2[3];
     tau_p2+=a1/modptr->idunion.rig.eigenvalues2[4];
     tau_p2+=a2/modptr->idunion.rig.eigenvalues2[5];
     fprintf(file,"%e ",tau_p2);
*/
    }

/*  d0 = Complex(qdl[0],-qdl[3]); d1 = Complex(qdl[2],-qdl[1]); */
    d0 = Complex(qld[0],-qld[3]); d1 = Complex(qld[2],-qld[1]);
    d0cc = Conjg(d0); d1cc = Conjg(d1);
    dd1 = Cmul(d0,d0);
    dd2 = Cmul(d0,d1);
    dd3 = Cmul(d0,d1cc);
    dd4 = Cmul(d0,d0cc);
    dd5 = Cmul(d1,d1);
    dd6 = Cmul(d1,d1cc);
    dd7 = Cmul(d1,d0cc);
    dd8 = Cmul(d1cc,d1cc);
    dd9 = Cmul(d1cc,d0cc);
    dd10 = Cmul(d0cc,d0cc);

/*
    D[1][1] = dd1;
*/
    D[1][2] = RCmul(msqrt_two,dd3);
/*
    D[1][3] = dd8;
    D[2][1] = RCmul(sqrt_two,dd2);
*/
    D[2][2] = Csub(dd4,dd6);
/*
    D[2][3] = RCmul(msqrt_two,dd9);
    D[3][1] = dd5;
*/
    D[3][2] = RCmul(sqrt_two,dd7);
/*
    D[3][3] = dd10;
*/

    tau_p1=0.;
    tau=0.;
    for(k=1;k<=n1;k++) {
     Pak = Complex(0.,0.);
     for(j=1;j<=n1;j++) {
      Dt = RCmul(modptr->idunion.rig.eigenstates1[j][k],D[j][2]);
      Pak = Cadd(Pak,Dt);
     }
     A[1][k] = Pak.r*Pak.r+Pak.i*Pak.i;
     tau += A[1][k]/modptr->idunion.rig.eigenvalues1[k];
    }
    /* fprintf(file,"%e ",factns*tau); */
    tau_p1+=tau;

    tau=0.;
    for(k=1;k<=n1;k++) {
     Pak = Complex(0.,0.);
     for(j=1;j<=n1;j++) {
      Dt = RCmul(modptr->idunion.rig.eigenstates1[j][k],D[j][2]);
      Pak = Cadd(Pak,Dt);
     }
     A[3][k] = Pak.r*Pak.r+Pak.i*Pak.i;
     tau += A[3][k]/modptr->idunion.rig.eigenvalues1[k];
    }
    /* fprintf(file,"%e ",factns*tau); */
    tau_p1+=tau;

    tau=0.;
    for(k=1;k<=n1;k++) {
     Pak = Complex(0.,0.);
     for(j=1;j<=n1;j++) {
      Dt = RCmul(modptr->idunion.rig.eigenstates1[j][k],D[j][2]);
      Pak = Cadd(Pak,Dt);
     }
     A[2][k] = Pak.r*Pak.r+Pak.i*Pak.i;
     tau += A[2][k]/modptr->idunion.rig.eigenvalues1[k];
    }
    /* fprintf(file,"%e ",factns*tau); */
    tau_p1+=tau;
    tau_p1 /= (_REAL)(2*1+1);
    taur_p1[i+1] = factns*tau_p1;
/*
    fprintf(file,"%e ",1.*tau_p1);
    fprintf(file,"%e ",factns*tau_p1);
*/

/*
    printf("\n");
    for(j=1;j<=n1;j++) {
     for(k=1;k<=n1;k++) {
      printf("A: %e ",A[j][k]);
     }
     printf("\n");
    }
*/
/*
    D[1][1] = Complex(0.,0.);
    Dt = Cmul(dd1,dd1);
    D[1][1] = Cadd(D[1][1],Dt);

    D[1][2] = Complex(0.,0.);
    Dt = Cmul(dd3,dd1);
    Dt = RCmul(-2.,Dt);
    D[1][2] = Cadd(D[1][2],Dt);
*/

    D[1][3] = Complex(0.,0.);
    Dt = Cmul(dd8,dd1);
    Dt = RCmul(sqrt_six,Dt);
    D[1][3] = Cadd(D[1][3],Dt);

/*
    D[1][4] = Complex(0.,0.);
    Dt = Cmul(dd8,dd3);
    Dt = RCmul(-2.0,Dt);
    D[1][4] = Cadd(D[1][4],Dt);

    D[1][5] = Complex(0.,0.);
    Dt = Cmul(dd8,dd8);
    D[1][5] = Cadd(D[1][5],Dt);

    D[2][1] = Complex(0.,0.);
    Dt = Cmul(dd2,dd1);
    Dt = RCmul(2.,Dt);
    D[2][1] = Cadd(D[2][1],Dt);

    D[2][2] = Complex(0.,0.);
    Dt = Cmul(dd6,dd1);
    Dt = RCmul(-3.,Dt);
    D[2][2] = Cadd(D[2][2],Dt);
    Dt = Cmul(dd4,dd1);
    D[2][2] = Cadd(D[2][2],Dt);
*/

    D[2][3] = Complex(0.,0.);
    Dt = Cmul(dd8,dd2);
    Dt = RCmul(sqrt_six,Dt);
    D[2][3] = Cadd(D[2][3],Dt);
    Dt = Cmul(dd9,dd1);
    Dt = RCmul(-sqrt_six,Dt);
    D[2][3] = Cadd(D[2][3],Dt);

/*
    D[2][4] = Complex(0.,0.);
    Dt = Cmul(dd8,dd6);
    Dt = RCmul(-1.,Dt);
    D[2][4] = Cadd(D[2][4],Dt);
    Dt = Cmul(dd9,dd3);
    Dt = RCmul(3.,Dt);
    D[2][4] = Cadd(D[2][4],Dt);

    D[2][5] = Complex(0.,0.);
    Dt = Cmul(dd9,dd8);
    Dt = RCmul(-2.,Dt);
    D[2][5] = Cadd(D[2][5],Dt);

    D[3][1] = Complex(0.,0.);
    Dt = Cmul(dd5,dd1);
    Dt = RCmul(sqrt_six,Dt);
    D[3][1] = Cadd(D[3][1],Dt);

    D[3][2] = Complex(0.,0.);
    Dt = Cmul(dd6,dd2);
    Dt = RCmul(-sqrt_six,Dt);
    D[3][2] = Cadd(D[3][2],Dt);
    Dt = Cmul(dd7,dd1);
    Dt = RCmul(sqrt_six,Dt);
    D[3][2] = Cadd(D[3][2],Dt);
*/

    D[3][3] = Complex(0.,0.);
    Dt = Cmul(dd8,dd5);
    D[3][3] = Cadd(D[3][3],Dt);
    Dt = Cmul(dd9,dd2);
    Dt = RCmul(-4.,Dt);
    D[3][3] = Cadd(D[3][3],Dt);
    Dt = Cmul(dd10,dd1);
    D[3][3] = Cadd(D[3][3],Dt);

/*
    D[3][4] = Complex(0.,0.);
    Dt = Cmul(dd9,dd6);
    Dt = RCmul(sqrt_six,Dt);
    D[3][4] = Cadd(D[3][4],Dt);
    Dt = Cmul(dd10,dd3);
    Dt = RCmul(-sqrt_six,Dt);
    D[3][4] = Cadd(D[3][4],Dt);

    D[3][5] = Complex(0.,0.);
    Dt = Cmul(dd10,dd8);
    Dt = RCmul(sqrt_six,Dt);
    D[3][5] = Cadd(D[3][5],Dt);

    D[4][1] = Complex(0.,0.);
    Dt = Cmul(dd5,dd2);
    Dt = RCmul(2.,Dt);
    D[4][1] = Cadd(D[4][1],Dt);

    D[4][2] = Complex(0.,0.);
    Dt = Cmul(dd6,dd5);
    Dt = RCmul(-1.,Dt);
    D[4][2] = Cadd(D[4][2],Dt);
    Dt = Cmul(dd7,dd2);
    Dt = RCmul(3.,Dt);
    D[4][2] = Cadd(D[4][2],Dt);
*/

    D[4][3] = Complex(0.,0.);
    Dt = Cmul(dd9,dd5);
    Dt = RCmul(-sqrt_six,Dt);
    D[4][3] = Cadd(D[4][3],Dt);
    Dt = Cmul(dd10,dd2);
    Dt = RCmul(sqrt_six,Dt);
    D[4][3] = Cadd(D[4][3],Dt);

/*
    D[4][4] = Complex(0.,0.);
    Dt = Cmul(dd10,dd6);
    Dt = RCmul(-3.,Dt);
    D[4][4] = Cadd(D[4][4],Dt);
    Dt = Cmul(dd10,dd4);
    D[4][4] = Cadd(D[4][4],Dt);

    D[4][5] = Complex(0.,0.);
    Dt = Cmul(dd10,dd9);
    Dt = RCmul(-2.,Dt);
    D[4][5] = Cadd(D[4][5],Dt);

    D[5][1] = Complex(0.,0.);
    Dt = Cmul(dd5,dd5);
    D[5][1] = Cadd(D[5][1],Dt);

    D[5][2] = Complex(0.,0.);
    Dt = Cmul(dd7,dd5);
    Dt = RCmul(2.,Dt);
    D[5][2] = Cadd(D[5][2],Dt);
*/

    D[5][3] = Complex(0.,0.);
    Dt = Cmul(dd10,dd5);
    Dt = RCmul(sqrt_six,Dt);
    D[5][3] = Cadd(D[5][3],Dt);

/*
    D[5][4] = Complex(0.,0.);
    Dt = Cmul(dd10,dd7);
    Dt = RCmul(2.,Dt);
    D[5][4] = Cadd(D[5][4],Dt);

    D[5][5] = Complex(0.,0.);
    Dt = Cmul(dd10,dd10);
    D[5][5] = Cadd(D[5][5],Dt);
*/

    tau_p2=0.;
    tau=0.;
    for(k=1;k<=n2;k++) {
     Pak = Complex(0.,0.);
     for(j=1;j<=n2;j++) {
      Dt = RCmul(modptr->idunion.rig.eigenstates2[j][k],D[j][3]);
      Pak = Cadd(Pak,Dt);
     }
     A[1][k] = Pak.r*Pak.r+Pak.i*Pak.i;
     tau += A[1][k]/modptr->idunion.rig.eigenvalues2[k];
    }
    /* fprintf(file,"%e ",factns*tau); */
    tau_p2+=tau;

    tau=0.;
    for(k=1;k<=n2;k++) {
     Pak = Complex(0.,0.);
     for(j=1;j<=n2;j++) {
      Dt = RCmul(modptr->idunion.rig.eigenstates2[j][k],D[j][3]);
      Pak = Cadd(Pak,Dt);
     }
     A[5][k] = Pak.r*Pak.r+Pak.i*Pak.i;
     tau += A[5][k]/modptr->idunion.rig.eigenvalues2[k];
    }
    /* fprintf(file,"%e ",factns*tau); */
    tau_p2+=tau;

    tau=0.;
    for(k=1;k<=n2;k++) {
     Pak = Complex(0.,0.);
     for(j=1;j<=n2;j++) {
      Dt = RCmul(modptr->idunion.rig.eigenstates2[j][k],D[j][3]);
      Pak = Cadd(Pak,Dt);
     }
     A[2][k] = Pak.r*Pak.r+Pak.i*Pak.i;
     tau += A[2][k]/modptr->idunion.rig.eigenvalues2[k];
    }
    /* fprintf(file,"%e ",factns*tau); */
    tau_p2+=tau;

    tau=0.;
    for(k=1;k<=n2;k++) {
     Pak = Complex(0.,0.);
     for(j=1;j<=n2;j++) {
      Dt = RCmul(modptr->idunion.rig.eigenstates2[j][k],D[j][3]);
      Pak = Cadd(Pak,Dt);
     }
     A[4][k] = Pak.r*Pak.r+Pak.i*Pak.i;
     tau += A[4][k]/modptr->idunion.rig.eigenvalues2[k];
    }
    /* fprintf(file,"%e ",factns*tau); */
    tau_p2+=tau;

    tau=0.;
    for(k=1;k<=n2;k++) {
     Pak = Complex(0.,0.);
     for(j=1;j<=n2;j++) {
      Dt = RCmul(modptr->idunion.rig.eigenstates2[j][k],D[j][3]);
      Pak = Cadd(Pak,Dt);
     }
     A[3][k] = Pak.r*Pak.r+Pak.i*Pak.i;
     tau += A[3][k]/modptr->idunion.rig.eigenvalues2[k];
    }
    /* fprintf(file,"%e ",factns*tau); */
    tau_p2+=tau;
    tau_p2 /= (_REAL)(2*2+1);
    taur_p2[i+1]=factns*tau_p2;

/*
    printf("\n");
    for(j=1;j<=5;j++) {
     for(k=1;k<=5;k++) {
      printf("A: %e ",A[j][k]);
     }
     printf("\n");
    }
*/
/*
   { _REAL taum,taue,Ssq;
     FILE *filet;
     filet=fopen("j.tmp","w");
     taum=3.2;
     taue=0.340;
     Ssq = 0.77;
     for(k=0;k<acfptr->idunion.acff.nomega;k++) {
      omega=acfptr->idunion.acff.omegaH[i][k]-acfptr->idunion.acff.omegaN[i][k];
      omega *= omega;
      J0 = (2./5.)*(Ssq*taum/(1.+omega*taum*taum) + 
                    (1.-Ssq)*taue );
      if(i==0) {
       fprintf(filet,"%e %e\n",l_sqrt(omega),J0);
      }
      omega=acfptr->idunion.acff.omegaN[i][k];
      omega *= omega;
      J1 = (2./5.)*(Ssq*taum/(1.+omega*taum*taum) + 
                    (1.-Ssq)*taue );
      omega=acfptr->idunion.acff.omegaH[i][k]+acfptr->idunion.acff.omegaN[i][k];
      omega *= omega;
      J2 = (2./5.)*(Ssq*taum/(1.+omega*taum*taum) + 
                    (1.-Ssq)*taue );
      relax = (1./4.)*(J0+3.*J1+6.*J2);
      relax *= acfptr->idunion.acff.gNgH[i];
      relax *= acfptr->idunion.acff.modsq[i+1];

      omega=acfptr->idunion.acff.omegaN[i][k];
      omega *= omega;
      relax += acfptr->idunion.acff.csa[i]*omega*omega*J1;
      relax = 1000./relax;
      fprintf(file,"%e ",relax);
     }
    }
*/
    d = acfptr->idunion.acff.gNgH[i]/20.;
    d *= modsq[i+1];

    if(acfptr->idunion.acff.type==1) {
     for(k=0;k<acfptr->idunion.acff.nomega;k++) {
      omega=acfptr->idunion.acff.omegaN[i][k];
      omega *= omega;
      c = acfptr->idunion.acff.csa[i]*omega;
      omega=0.;
      omega *= omega;
      J_0=0.;
      for(j=1;j<=5;j++) {
       J_0 += 2.*A[3][j]*taub2[j]/(1.+taub2[j]*taub2[j]*omega);
      }
      omega=acfptr->idunion.acff.omegaH[i][k]-acfptr->idunion.acff.omegaN[i][k];
      omega *= omega;
      J_ImS=0.;
      for(j=1;j<=5;j++) {
       J_ImS += 2.*A[3][j]*taub2[j]/(1.+taub2[j]*taub2[j]*omega);
      }
      omega=acfptr->idunion.acff.omegaN[i][k];
      omega *= omega;
      J_S=0.;
      for(j=1;j<=5;j++) {
       J_S += 2.*A[4][j]*taub2[j]/(1.+taub2[j]*taub2[j]*omega);
      }
      omega=acfptr->idunion.acff.omegaH[i][k];
      omega *= omega;
      J_I=0.;
      for(j=1;j<=5;j++) {
       J_I += 2.*A[4][j]*taub2[j]/(1.+taub2[j]*taub2[j]*omega);
      }
      omega=acfptr->idunion.acff.omegaH[i][k]+acfptr->idunion.acff.omegaN[i][k];
      omega *= omega;
      J_IpS=0.;
      for(j=1;j<=5;j++) {
       J_IpS += 2.*A[5][j]*taub2[j]/(1.+taub2[j]*taub2[j]*omega);
      }

      relax_Sz[i+1][k+1] = d*J_ImS +(3.*d+c)*J_S + 6.*d*J_IpS;
      relax_Sxy[i+1][k+1] = (6.*d+2*c)*J_0/3. +d*J_ImS/2. +(3.*d+c)*J_S/2.
                       +3.*d*J_I +3.*d*J_IpS;
      relax_Iz_Sz[i+1][k+1] = -d*J_ImS + 6.*d*J_IpS;
     }
    }

   } /* enddo bonds */
   r1-=acfptr->idunion.acff.nbond;
   r2-=acfptr->idunion.acff.nbond;

   for(i=1;i<=acfptr->idunion.acff.nbond;i++) {
    fprintf(file,"%d %e %e ",i,taur_p1[i],taur_p2[i]);
    for(k=1;k<=acfptr->idunion.acff.nomega;k++) {
     fprintf(file,"%e ",relax_Sz[i][k]);
    }
    fprintf(file,"\n");
   }
   for(i=1;i<=acfptr->idunion.acff.nbond;i++) {
    fprintf(file,"%d %e %e ",i,taur_p1[i],taur_p2[i]);
    for(k=1;k<=acfptr->idunion.acff.nomega;k++) {
     fprintf(file,"%e ",relax_Sxy[i][k]);
    }
    fprintf(file,"\n");
   }
   for(i=1;i<=acfptr->idunion.acff.nbond;i++) {
    fprintf(file,"%d %e %e ",i,taur_p1[i],taur_p2[i]);
    for(k=1;k<=acfptr->idunion.acff.nomega;k++) {
     fprintf(file,"%e ",relax_Iz_Sz[i][k]);
    }
    fprintf(file,"\n");
   }

   free_vector(taur_p1,1,acfptr->idunion.acff.nbond);
   free_vector(taur_p2,1,acfptr->idunion.acff.nbond);
   fflush(file);

}
