#include "head.h"
#include <time.h>

#include <openacc.h>

// extern char tok[];
// char tok[];
// void *get2(int,int);

extern idstruct *startid(),*getstruct();
rotstruct **rot_matrix();
rot1struct **rot1_matrix();

_REAL **Fb_xx,**Fb_xy,**Fb_xz,**Fb_yy,**Fb_yz,**Fb_zz,*C1,*C2;

/************************************************************************/

void
init_mca_p2p1()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(MCA_P2P1, id);

	dptr->idunion.mca_p2p1.cases = 0;
	dptr->idunion.mca_p2p1.casesw = 0.;
	dptr->idunion.mca_p2p1.dump = 1;
	dptr->idunion.mca_p2p1.type = 2;
	dptr->idunion.mca_p2p1.rigid = 0;
	dptr->idunion.mca_p2p1.m1 = 3;
	dptr->idunion.mca_p2p1.m2 = 0;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("MCA_P2P1: expected file name","");
	strcpy(dptr->idunion.mca_p2p1.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.mca_p2p1.fname, 0, &exists, "MCA_P2P1"))
		exit(1);
	if ((dptr->idunion.mca_p2p1.file = 
			 fopen(dptr->idunion.mca_p2p1.fname, "w")) == NULL){
		perror(dptr->idunion.mca_p2p1.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
#ifdef _OPENACC
{
   int mygpu, myrealgpu, num_devices;
   acc_device_t my_device_type;

   my_device_type = acc_device_nvidia;
   mygpu=0;
//   if(argc == 1) mygpu = 0; else mygpu = atoi(argv[1]);
   acc_set_device_type(my_device_type) ;
   num_devices = acc_get_num_devices(my_device_type) ;
   fprintf(stderr,"Number of devices available: %d \n",num_devices);
   acc_set_device_num(mygpu,my_device_type);
   fprintf(stderr,"Trying to use GPU: %d \n",mygpu);
   myrealgpu = acc_get_device_num(my_device_type);
   fprintf(stderr,"Actually I am using GPU: %d \n",myrealgpu);
   if(mygpu != myrealgpu) {
     fprintf(stderr,"I cannot use the requested GPU: %d\n",mygpu);
     exit(1);
   }
}
#endif
        return;
}

/************************************************************************/
void
setup_mca_p2p1()
{
	idstruct	*dptr,*bptr,*mptr;
	int i,j,k,icount;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != MCA_P2P1)
                        inerr("setup_mca_p2p1: expected id or keyword\n","");
        } else {
                printf("setup_mca_p2p1: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((mptr = getstruct(tok)) == NULL) 
                inerr("setup_mca_p2p1: no such id\n","");
        if (mptr->id_type != MCARIG)
                inerr("setup_mca_p2p1: mcarig id expected\n","");
        gettok();
        if ((bptr = getstruct(tok)) == NULL) 
                inerr("setup_mca_p2p1: no such id\n","");
        if (bptr->id_type != BEAD)
                inerr("setup_mca_p2p1: bead id expected\n","");

        dptr->idunion.mca_p2p1.mca = mptr;
        dptr->idunion.mca_p2p1.bead = bptr;
        dptr->idunion.mca_p2p1.m1=mptr->idunion.mcarig.nf;
        dptr->idunion.mca_p2p1.m2=0;

        gettok();
        for(;;) {
         if(tok[0]==';') {
          break;
         } else if(!strcmp(tok,"READFSII")) {
          dptr->idunion.mca_p2p1.dump=-1;
          gettok();
         } else if(!strcmp(tok,"RM2-I")) {
          dptr->idunion.mca_p2p1.type=1;
          gettok();
         } else if(!strcmp(tok,"RM2-II")) {
          dptr->idunion.mca_p2p1.type=2;
          gettok();
         } else if(!strcmp(tok,"RM2-III")) {
          dptr->idunion.mca_p2p1.type=3;
          gettok();
         } else if(!strcmp(tok,"RM2-IIm")) {
          gettok();
          tokint(&dptr->idunion.mca_p2p1.m1,tok);
          gettok();
          tokint(&dptr->idunion.mca_p2p1.m2,tok);
	  if(dptr->idunion.mca_p2p1.m2==0)
	  {
           dptr->idunion.mca_p2p1.type=1;
          }
	  else if(dptr->idunion.mca_p2p1.m2>0)
	  {
           dptr->idunion.mca_p2p1.type=2;
	  }
          gettok();
         } else if(!strcmp(tok,"RIGID")) {
          dptr->idunion.mca_p2p1.rigid=1;
          gettok();
         } else {
          printf("setup_mca_p2p1: option not found.\n");
          exit(1);
         }
        }

 printf("setup_p2p1: type= %d\n",dptr->idunion.mca_p2p1.type); fflush(stdout);

        dptr->idunion.mca_p2p1.n0=dptr->idunion.mca_p2p1.m1
                                 +dptr->idunion.mca_p2p1.m2;
        dptr->idunion.mca_p2p1.n12=dptr->idunion.mca_p2p1.n0*
         (dptr->idunion.mca_p2p1.n0+1)/2;
        dptr->idunion.mca_p2p1.n1234=dptr->idunion.mca_p2p1.m1*
         (dptr->idunion.mca_p2p1.m1+1)/2
                                 *dptr->idunion.mca_p2p1.m1*
         (dptr->idunion.mca_p2p1.m1+1)/2;
        dptr->idunion.mca_p2p1.n123456=dptr->idunion.mca_p2p1.m1*
         (dptr->idunion.mca_p2p1.m1+1)/2
                                 *dptr->idunion.mca_p2p1.m1*
         (dptr->idunion.mca_p2p1.m1+1)/2
                                 *dptr->idunion.mca_p2p1.m1*
         (dptr->idunion.mca_p2p1.m1+1)/2;

        dptr->idunion.mca_p2p1.dim_I_I=
         dptr->idunion.mca_p2p1.n12*(dptr->idunion.mca_p2p1.n12+1)/2;
        dptr->idunion.mca_p2p1.dim_I_II=
         dptr->idunion.mca_p2p1.n12*dptr->idunion.mca_p2p1.n1234;
        dptr->idunion.mca_p2p1.dim_I_III=
         dptr->idunion.mca_p2p1.n12*dptr->idunion.mca_p2p1.n123456;
        dptr->idunion.mca_p2p1.dim_II_II=
         dptr->idunion.mca_p2p1.n1234*(dptr->idunion.mca_p2p1.n1234+1)/2;
        dptr->idunion.mca_p2p1.dim_II_III=
         dptr->idunion.mca_p2p1.n1234*dptr->idunion.mca_p2p1.n123456;
        dptr->idunion.mca_p2p1.dim_III_III=
         dptr->idunion.mca_p2p1.n123456*(dptr->idunion.mca_p2p1.n123456+1)/2;

        dptr->idunion.mca_p2p1.s_I_I = (_REAL *)
          get2(dptr->idunion.mca_p2p1.dim_I_I,sizeof(_REAL));
        dptr->idunion.mca_p2p1.f_I_I = (_REAL *)
          get2(dptr->idunion.mca_p2p1.dim_I_I,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p2p1.dim_I_I;i++) {
         *dptr->idunion.mca_p2p1.s_I_I++ =0.;
         *dptr->idunion.mca_p2p1.f_I_I++ =0.;
        }
        dptr->idunion.mca_p2p1.s_I_I-=dptr->idunion.mca_p2p1.dim_I_I;
        dptr->idunion.mca_p2p1.f_I_I-=dptr->idunion.mca_p2p1.dim_I_I;

        if(dptr->idunion.mca_p2p1.type>1) {

        dptr->idunion.mca_p2p1.s_I_II = (_REAL *)
          get2(dptr->idunion.mca_p2p1.dim_I_II,sizeof(_REAL));
        dptr->idunion.mca_p2p1.f_I_II = (_REAL *)
          get2(dptr->idunion.mca_p2p1.dim_I_II,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p2p1.dim_I_II;i++) {
         *dptr->idunion.mca_p2p1.s_I_II++ =0.;
         *dptr->idunion.mca_p2p1.f_I_II++ =0.;
        }
        dptr->idunion.mca_p2p1.s_I_II-=dptr->idunion.mca_p2p1.dim_I_II;
        dptr->idunion.mca_p2p1.f_I_II-=dptr->idunion.mca_p2p1.dim_I_II;

        dptr->idunion.mca_p2p1.s_II_II = (_REAL *)
          get2(dptr->idunion.mca_p2p1.dim_II_II,sizeof(_REAL));
        dptr->idunion.mca_p2p1.f_II_II = (_REAL *)
          get2(dptr->idunion.mca_p2p1.dim_II_II,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p2p1.dim_II_II;i++) {
         *dptr->idunion.mca_p2p1.s_II_II++ =0.;
         *dptr->idunion.mca_p2p1.f_II_II++ =0.;
        }
        dptr->idunion.mca_p2p1.s_II_II-=dptr->idunion.mca_p2p1.dim_II_II;
        dptr->idunion.mca_p2p1.f_II_II-=dptr->idunion.mca_p2p1.dim_II_II;

        }

        if(dptr->idunion.mca_p2p1.type>2) {

        dptr->idunion.mca_p2p1.s_I_III = (_REAL *)
          get2(dptr->idunion.mca_p2p1.dim_I_III,sizeof(_REAL));
        dptr->idunion.mca_p2p1.f_I_III = (_REAL *)
          get2(dptr->idunion.mca_p2p1.dim_I_III,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p2p1.dim_I_III;i++) {
         *dptr->idunion.mca_p2p1.s_I_III++ =0.;
         *dptr->idunion.mca_p2p1.f_I_III++ =0.;
        }
        dptr->idunion.mca_p2p1.s_I_III-=dptr->idunion.mca_p2p1.dim_I_III;
        dptr->idunion.mca_p2p1.f_I_III-=dptr->idunion.mca_p2p1.dim_I_III;

        dptr->idunion.mca_p2p1.s_II_III = (_REAL *)
          get2(dptr->idunion.mca_p2p1.dim_II_III,sizeof(_REAL));
        dptr->idunion.mca_p2p1.f_II_III = (_REAL *)
          get2(dptr->idunion.mca_p2p1.dim_II_III,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p2p1.dim_II_III;i++) {
         *dptr->idunion.mca_p2p1.s_II_III++ =0.;
         *dptr->idunion.mca_p2p1.f_II_III++ =0.;
        }
        dptr->idunion.mca_p2p1.s_II_III-=dptr->idunion.mca_p2p1.dim_II_III;
        dptr->idunion.mca_p2p1.f_II_III-=dptr->idunion.mca_p2p1.dim_II_III;

        dptr->idunion.mca_p2p1.s_III_III = (_REAL *)
          get2(dptr->idunion.mca_p2p1.dim_III_III,sizeof(_REAL));
        dptr->idunion.mca_p2p1.f_III_III = (_REAL *)
          get2(dptr->idunion.mca_p2p1.dim_III_III,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_p2p1.dim_III_III;i++) {
         *dptr->idunion.mca_p2p1.s_III_III++ =0.;
         *dptr->idunion.mca_p2p1.f_III_III++ =0.;
        }
        dptr->idunion.mca_p2p1.s_III_III-=dptr->idunion.mca_p2p1.dim_III_III;
        dptr->idunion.mca_p2p1.f_III_III-=dptr->idunion.mca_p2p1.dim_III_III;

        }

        dptr->idunion.mca_p2p1.m = (modstruct *)
          get2(dptr->idunion.mca_p2p1.n0,sizeof(modstruct));
        dptr->idunion.mca_p2p1.l12 = (mod12struct *)
          get2(dptr->idunion.mca_p2p1.n12,sizeof(mod12struct));
        if(dptr->idunion.mca_p2p1.type>1) {
        dptr->idunion.mca_p2p1.l1234 = (mod1234struct *)
          get2(dptr->idunion.mca_p2p1.n1234,sizeof(mod1234struct));
        }
        if(dptr->idunion.mca_p2p1.type>2) {
        dptr->idunion.mca_p2p1.l123456 = (mod123456struct *)
          get2(dptr->idunion.mca_p2p1.n123456,sizeof(mod123456struct));
        }

        Fb_xx = matrix(0,bptr->idunion.bead.nbond-1,0,bptr->idunion.bead.nbond-1);
        Fb_xy = matrix(0,bptr->idunion.bead.nbond-1,0,bptr->idunion.bead.nbond-1);
        Fb_xz = matrix(0,bptr->idunion.bead.nbond-1,0,bptr->idunion.bead.nbond-1);
        Fb_yy = matrix(0,bptr->idunion.bead.nbond-1,0,bptr->idunion.bead.nbond-1);
        Fb_yz = matrix(0,bptr->idunion.bead.nbond-1,0,bptr->idunion.bead.nbond-1);
        Fb_zz = matrix(0,bptr->idunion.bead.nbond-1,0,bptr->idunion.bead.nbond-1);
        C1 = (_REAL*)get2(bptr->idunion.bead.nbond,sizeof(_REAL));
        C2 = (_REAL*)get2(bptr->idunion.bead.nbond,sizeof(_REAL));

        dptr->idunion.mca_p2p1.Fm = rot_matrix
         (0,dptr->idunion.mca_p2p1.n0-1,0,dptr->idunion.mca_p2p1.n0-1);

        printf("setup_mca_p2p1:\n");
        printf("m1:        %d\n",dptr->idunion.mca_p2p1.m1);
        printf("m2:        %d\n",dptr->idunion.mca_p2p1.m2);
        printf("n0:        %d\n",dptr->idunion.mca_p2p1.n0);
        printf("n12:       %d\n",dptr->idunion.mca_p2p1.n12);
        printf("n1234:     %d\n",dptr->idunion.mca_p2p1.n1234);
        printf("n123456:   %d\n",dptr->idunion.mca_p2p1.n123456);
        printf("dim_I_I:      %d\n",dptr->idunion.mca_p2p1.dim_I_I);
        printf("dim_I_II:     %d\n",dptr->idunion.mca_p2p1.dim_I_II);
        printf("dim_I_III:    %d\n",dptr->idunion.mca_p2p1.dim_I_III);
        printf("dim_II_II:    %ld\n",dptr->idunion.mca_p2p1.dim_II_II);
        printf("dim_II_III:   %ld\n",dptr->idunion.mca_p2p1.dim_II_III);
        printf("dim_III_III:  %ld\n",dptr->idunion.mca_p2p1.dim_III_III);

	return;
}

/************************************************************************/
do_first_p2p1(idstruct *objptr,idstruct *serptr)
{
   _REAL *box=serptr->idunion.ser.set->boxc;
   idstruct *bptr,*mptr;
   int i1,i2,i3,i4,j1_0,j2_0,j3_0,j4_0,j1,j2,j3,j4,
       i5,i6,j5_0,j6_0,j5,j6,
       i,j,k,l,m,ib,jb,jb0,kk,im,m1,m2,nbondm,nmode,nmode1,nmode2,overlap;
   _REAL **r1, **r2;
   _REAL *rt1,*rt2,r,zi,zj,zr,cutsq,
         boxm1[3],x1,y1,z1,x2,y2,z2,dx,dy,dz,scp,scp34,scp56;
   bondstruct *lb;
   modstruct *p,*p1,*p2,*p3,*p4,*p5,*p6;
   mod12struct *l12,*l12_1,*l12_2;
   mod1234struct *l1234,*l1234_1,*l1234_2;
   mod123456struct *l123456,*l123456_1,*l123456_2;
   _REAL xm,ym,zm,f[5][5],Ct[3],csum,D[3][3],b,c,weight,ft;

   mptr=objptr->idunion.mca_p2p1.mca;
   bptr=mptr->idunion.mcarig.bead;
   zr=bptr->idunion.bead.zr;
   cutsq=bptr->idunion.bead.cutsq;
   nmode=objptr->idunion.mca_p2p1.n0;
   nmode1=objptr->idunion.mca_p2p1.m1;
   nmode2=objptr->idunion.mca_p2p1.m2;
   nbondm=bptr->idunion.bead.nbond;
   for(k=0;k<3;k++) {
    boxm1[k]=1./box[k];
   }
   if(serptr->idunion.ser.iweight)
   {
    weight=serptr->idunion.ser.w->idunion.serw.w;
   } else {
    weight=1.;
   }
   b=l_sqrt(3./8.);
   c=0.5;
/*
   b=0.5;
   c=1./l_sqrt(6.);
*/

   ++objptr->idunion.mca_p2p1.cases;
   objptr->idunion.mca_p2p1.casesw+=weight;

   lb=bptr->idunion.bead.l;
   r1=bptr->idunion.bead.r1;
   r2=bptr->idunion.bead.r2;
   for (i=0;i<nbondm;i++,r1++,r2++) {
    rt1=*r1; rt2=*r2;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
    x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
    dx=x2-x1; dy=y2-y1; dz=z2-z1;
    lb[i].id1 = (int)(*r1);
    lb[i].id2 = (int)(*r2);
    lb[i].x1 = x1; lb[i].y1 = y1; lb[i].z1 = z1;
    lb[i].x2 = x2; lb[i].y2 = y2; lb[i].z2 = z2;
   }

   p=objptr->idunion.mca_p2p1.m;
   for (j=0;j<nmode;j++) {
    xm=0.; ym=0.; zm=0.;
    r1=bptr->idunion.bead.r1;
    r2=bptr->idunion.bead.r2;
    for (i=0;i<nbondm;i++,r1++,r2++) {
     rt1=*r1; rt2=*r2;
     x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
     x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
     dx=x2-x1; dy=y2-y1; dz=z2-z1;
#ifdef PBC
     dx -= box[0]*rint(dx*boxm1[0]);
     dy -= box[1]*rint(dy*boxm1[1]);
     dz -= box[2]*rint(dz*boxm1[2]);
#endif
     xm += mptr->idunion.mcarig.evec[i+1][j+1]*dx;
     ym += mptr->idunion.mcarig.evec[i+1][j+1]*dy;
     zm += mptr->idunion.mcarig.evec[i+1][j+1]*dz;
    }
    p->mod=l_sqrt(xm*xm+ym*ym+zm*zm);
    p->m[0]=xm;
    p->m[1]=ym;
    p->m[2]=zm;
    p++;
   }

   p=objptr->idunion.mca_p2p1.m;
   l12=objptr->idunion.mca_p2p1.l12;
   for (i1=0,p1=p;i1<nmode;i1++,p1++) {
    for (i2=i1,p2=p1;i2<nmode;i2++,p2++) {
     scp = p1->m[0]*p2->m[0]+p1->m[1]*p2->m[1]+p1->m[2]*p2->m[2];
     l12->phi[0]=c*(3.*p1->m[2]*p2->m[2] - scp);
     l12->dphidm1[0][0]=-c*p2->m[0];
     l12->dphidm1[0][1]=-c*p2->m[1];
     l12->dphidm1[0][2]= 2.*c*p2->m[2];
     l12->dphidm2[0][0]=-c*p1->m[0];
     l12->dphidm2[0][1]=-c*p1->m[1];
     l12->dphidm2[0][2]= 2.*c*p1->m[2];
     l12->phi[1]=b*(p1->m[0]*p2->m[0] - p1->m[1]*p2->m[1]);
     l12->dphidm1[1][0]= b*p2->m[0];
     l12->dphidm1[1][1]=-b*p2->m[1];
     l12->dphidm1[1][2]= 0.;
     l12->dphidm2[1][0]= b*p1->m[0];
     l12->dphidm2[1][1]=-b*p1->m[1];
     l12->dphidm2[1][2]= 0.;
     l12->phi[2]=b*(p1->m[0]*p2->m[1] + p1->m[1]*p2->m[0]);
     l12->dphidm1[2][0]= b*p2->m[1];
     l12->dphidm1[2][1]= b*p2->m[0];
     l12->dphidm1[2][2]= 0.;
     l12->dphidm2[2][0]= b*p1->m[1];
     l12->dphidm2[2][1]= b*p1->m[0];
     l12->dphidm2[2][2]= 0.;
     l12->phi[3]=b*(p1->m[0]*p2->m[2] + p1->m[2]*p2->m[0]);
     l12->dphidm1[3][0]= b*p2->m[2];
     l12->dphidm1[3][1]= 0.;
     l12->dphidm1[3][2]= b*p2->m[0];
     l12->dphidm2[3][0]= b*p1->m[2];
     l12->dphidm2[3][1]= 0.;
     l12->dphidm2[3][2]= b*p1->m[0];
     l12->phi[4]=b*(p1->m[1]*p2->m[2] + p1->m[2]*p2->m[1]);
     l12->dphidm1[4][0]= 0.;
     l12->dphidm1[4][1]= b*p2->m[2];
     l12->dphidm1[4][2]= b*p2->m[1];
     l12->dphidm2[4][0]= 0.;
     l12->dphidm2[4][1]= b*p1->m[2];
     l12->dphidm2[4][2]= b*p1->m[1];
     l12++;
    }
   }

   if(objptr->idunion.mca_p2p1.type>1) {

   l1234=objptr->idunion.mca_p2p1.l1234;
   l12=objptr->idunion.mca_p2p1.l12;
   p=objptr->idunion.mca_p2p1.m;
   for (i1=0;i1<nmode1;i1++) {
    for (i2=i1;i2<nmode1;i2++) {
     m=(2*nmode-i1+1)*i1/2+(i2-i1);
     for (i3=0;i3<nmode1;i3++) {
      for (i4=i3;i4<nmode1;i4++) {
       scp=p[i3].m[0]*p[i4].m[0]+p[i3].m[1]*p[i4].m[1]+p[i3].m[2]*p[i4].m[2];
       for(k=0;k<5;k++) {
        l1234->phi[k]=l12[m].phi[k]*scp;
        for(kk=0;kk<3;kk++) {
         l1234->dphidm1[k][kk]=l12[m].dphidm1[k][kk]*scp;
         l1234->dphidm2[k][kk]=l12[m].dphidm2[k][kk]*scp;
         l1234->dphidm3[k][kk]=l12[m].phi[k]*p[i4].m[kk];
         l1234->dphidm4[k][kk]=l12[m].phi[k]*p[i3].m[kk];
        }
       }
       l1234++;
      }
     }
    }
   }

   }

   if(objptr->idunion.mca_p2p1.type>2) {

   l123456=objptr->idunion.mca_p2p1.l123456;
   l12=objptr->idunion.mca_p2p1.l12;
   p=objptr->idunion.mca_p2p1.m;
   for (i1=0,p1=p;i1<nmode1;i1++,p1++) {
    for (i2=i1,p2=p1;i2<nmode1;i2++,p2++,l12++) {
     for (i3=0,p3=p;i3<nmode1;i3++,p3++) {
      for (i4=i3,p4=p3;i4<nmode1;i4++,p4++) {
       scp34=p3->m[0]*p4->m[0]+p3->m[1]*p4->m[1]+p3->m[2]*p4->m[2];
       for (i5=0,p5=p;i5<nmode1;i5++,p5++) {
        for (i6=i5,p6=p5;i6<nmode1;i6++,p6++) {
         scp56=p5->m[0]*p6->m[0]+p5->m[1]*p6->m[1]+p5->m[2]*p6->m[2];

       for(k=0;k<5;k++) {
        l123456->phi[k]=l12->phi[k]*scp34*scp56;
        for(kk=0;kk<3;kk++) {
         l123456->dphidm1[k][kk]=l12->dphidm1[k][kk]*scp34*scp56;
         l123456->dphidm2[k][kk]=l12->dphidm2[k][kk]*scp34*scp56;
         l123456->dphidm3[k][kk]=l12->phi[k]*p4->m[kk]*scp56;
         l123456->dphidm4[k][kk]=l12->phi[k]*p3->m[kk]*scp56;
         l123456->dphidm5[k][kk]=l12->phi[k]*scp34*p6->m[kk];
         l123456->dphidm6[k][kk]=l12->phi[k]*scp34*p5->m[kk];
        }
       }

         l123456++;
        }
       }
      }
     }
    }
   }

   }
       
// to exit without computing F and S matrices 
//   return;
/* BEGIN area for GPU */
// #pragma acc kernels
{
// #pragma acc parallel loop collapse(2)
   for (ib=0;ib<nbondm;ib++) {
    for (jb=0;jb<nbondm;jb++) {
#include "fbb.c"
     Fb_xx[ib][jb] = f[0][0];
     Fb_xy[ib][jb] = f[0][1];
     Fb_xz[ib][jb] = f[0][2];
     Fb_yy[ib][jb] = f[1][1];
     Fb_yz[ib][jb] = f[1][2];
     Fb_zz[ib][jb] = f[2][2];
    }
   }
}
/* END area for GPU */
   printf("Fb0_xx= %e\n",Fb_xx[0][0]); fflush(stdout);

   for (i=0;i<nmode;i++) {
    for (j=0;j<nmode;j++) {
     for (ib=0;ib<nbondm;ib++) {
      C1[ib]=mptr->idunion.mcarig.evec[ib+1][i+1];
      C2[ib]=mptr->idunion.mcarig.evec[ib+1][j+1];
     }
     ft=0.;
#pragma acc parallel loop collapse(2) reduction(+:ft)
     for (ib=0;ib<nbondm;ib++) {
      for (jb=0;jb<nbondm;jb++) {
          ft += C1[ib] * Fb_xx[ib][jb] * C2[jb];
      }
     }
     f[0][0]=ft;
     ft=0.;
#pragma acc parallel loop collapse(2) reduction(+:ft)
     for (ib=0;ib<nbondm;ib++) {
      for (jb=0;jb<nbondm;jb++) {
          ft += C1[ib] * Fb_xy[ib][jb] * C2[jb];
      }
     }
     f[0][1]=ft;
     ft=0.;
#pragma acc parallel loop collapse(2) reduction(+:ft)
     for (ib=0;ib<nbondm;ib++) {
      for (jb=0;jb<nbondm;jb++) {
          ft += C1[ib] * Fb_xz[ib][jb] * C2[jb];
      }
     }
     f[0][2]=ft;
     ft=0.;
#pragma acc parallel loop collapse(2) reduction(+:ft)
     for (ib=0;ib<nbondm;ib++) {
      for (jb=0;jb<nbondm;jb++) {
          ft += C1[ib] * Fb_yy[ib][jb] * C2[jb];
      }
     }
     f[1][1]=ft;
     ft=0.;
#pragma acc parallel loop collapse(2) reduction(+:ft)
     for (ib=0;ib<nbondm;ib++) {
      for (jb=0;jb<nbondm;jb++) {
          ft += C1[ib] * Fb_yz[ib][jb] * C2[jb];
      }
     }
     f[1][2]=ft;
     ft=0.;
#pragma acc parallel loop collapse(2) reduction(+:ft)
     for (ib=0;ib<nbondm;ib++) {
      for (jb=0;jb<nbondm;jb++) {
          ft += C1[ib] * Fb_zz[ib][jb] * C2[jb];
      }
     }
     f[2][2]=ft;

     for(m1=0;m1<3;m1++) {
      for(m2=m1;m2<3;m2++) {
       objptr->idunion.mca_p2p1.Fm[i][j].f[m1][m2]=f[m1][m2];
       objptr->idunion.mca_p2p1.Fm[i][j].f[m2][m1]=f[m1][m2];
      }
     }
    }
   }
   printf("Fm0= %e\n",objptr->idunion.mca_p2p1.Fm[0][0].f[0][0]); fflush(stdout);

   l12=objptr->idunion.mca_p2p1.l12;
   for (i1=0,i=0,kk=0,l12_1=l12;i1<nmode;i1++) {
    for (i2=i1;i2<nmode;i2++,i++,l12_1++) {
     j1_0=i1;
     j2_0=i2;

   for (j1=j1_0,j=i,l12_2=l12_1;j1<nmode;j1++,j2_0=j1) {
    for (j2=j2_0;j2<nmode;j2++,j++,l12_2++) {

     for(m1=0;m1<5;m1++) {
      for(m2=0;m2<5;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<5;m1++) {
      {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j1].f[k][l]  *
          l12_2->dphidm1[m2][l];
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j2].f[k][l]  *
          l12_2->dphidm2[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j1].f[k][l]  *
          l12_2->dphidm1[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j2].f[k][l]  *
          l12_2->dphidm2[m2][l];
        }
       }
  
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     scp+=2.*f[3][3];
     scp+=2.*f[4][4];
     *objptr->idunion.mca_p2p1.f_I_I += weight*scp/((_REAL)5);
     objptr->idunion.mca_p2p1.f_I_I++;

     scp=0.;
     scp+=l12_1->phi[0]*l12_2->phi[0];
     scp+=2.*l12_1->phi[1]*l12_2->phi[1];
     scp+=2.*l12_1->phi[2]*l12_2->phi[2];
     scp+=2.*l12_1->phi[3]*l12_2->phi[3];
     scp+=2.*l12_1->phi[4]*l12_2->phi[4];
     *objptr->idunion.mca_p2p1.s_I_I += weight*scp/((_REAL)5);
     objptr->idunion.mca_p2p1.s_I_I++;

     kk++;
    }
   }
    }
   }
   objptr->idunion.mca_p2p1.s_I_I -= objptr->idunion.mca_p2p1.dim_I_I;
   objptr->idunion.mca_p2p1.f_I_I -= objptr->idunion.mca_p2p1.dim_I_I;

   if(objptr->idunion.mca_p2p1.type == 1) return;

   l12=objptr->idunion.mca_p2p1.l12;
   l1234=objptr->idunion.mca_p2p1.l1234;
   for (i1=0,i=0,kk=0,l12_1=l12;i1<nmode;i1++) {
    for (i2=i1;i2<nmode;i2++,i++) {

   for (j1=0,j=0,l1234_2=l1234;j1<nmode1;j1++) {
    for (j2=j1;j2<nmode1;j2++) {
     for (j3=0;j3<nmode1;j3++) {
      for (j4=j3;j4<nmode1;j4++,j++,l1234_2++,kk++) {

     for(m1=0;m1<5;m1++) {
      for(m2=0;m2<5;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<5;m1++) {
       {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l12[i].dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l12[i].dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l12[i].dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l12[i].dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
         f[m1][m2]+=
          l12[i].dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l12[i].dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l12[i].dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l12[i].dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     scp+=2.*f[3][3];
     scp+=2.*f[4][4];
     *objptr->idunion.mca_p2p1.f_I_II += weight*scp/((_REAL)5);
     objptr->idunion.mca_p2p1.f_I_II++;

     scp=0.;
     scp+=l12[i].phi[0]*l1234_2->phi[0];
     scp+=2.*l12[i].phi[1]*l1234_2->phi[1];
     scp+=2.*l12[i].phi[2]*l1234_2->phi[2];
     scp+=2.*l12[i].phi[3]*l1234_2->phi[3];
     scp+=2.*l12[i].phi[4]*l1234_2->phi[4];
     *objptr->idunion.mca_p2p1.s_I_II += weight*scp/((_REAL)5);
     objptr->idunion.mca_p2p1.s_I_II++;

      }
     }
    }
   }
    }
   }
   objptr->idunion.mca_p2p1.s_I_II -= objptr->idunion.mca_p2p1.dim_I_II;
   objptr->idunion.mca_p2p1.f_I_II -= objptr->idunion.mca_p2p1.dim_I_II;

   l1234=objptr->idunion.mca_p2p1.l1234;
   for (i1=0,i=0,kk=0,l1234_1=l1234;i1<nmode1;i1++) {
    for (i2=i1;i2<nmode1;i2++) {
     for (i3=0;i3<nmode1;i3++) {
      for (i4=i3;i4<nmode1;i4++,i++,l1234_1++) {
     j1_0=i1;
     j2_0=i2;
     j3_0=i3;
     j4_0=i4;
   for (j1=j1_0,j=i,l1234_2=l1234_1;j1<nmode1;j1++,j2_0=j1) {
    for (j2=j2_0;j2<nmode1;j2++,j4_0=0,j3_0=0) {
     for (j3=j3_0;j3<nmode1;j3++,j4_0=j3) {
      for (j4=j4_0;j4<nmode1;j4++,j++,l1234_2++) {

     for(m1=0;m1<5;m1++) {
      for(m2=0;m2<5;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<5;m1++) {
       {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j1].f[k][l]  *
          l1234_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j2].f[k][l]  *
          l1234_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j3].f[k][l]  *
          l1234_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j4].f[k][l]  *
          l1234_2->dphidm4[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     scp+=2.*f[3][3];
     scp+=2.*f[4][4];
     *objptr->idunion.mca_p2p1.f_II_II += weight*scp/((_REAL)5);
     objptr->idunion.mca_p2p1.f_II_II++;

     scp=0.;
     scp+=l1234_1->phi[0]*l1234_2->phi[0];
     scp+=2.*l1234_1->phi[1]*l1234_2->phi[1];
     scp+=2.*l1234_1->phi[2]*l1234_2->phi[2];
     scp+=2.*l1234_1->phi[3]*l1234_2->phi[3];
     scp+=2.*l1234_1->phi[4]*l1234_2->phi[4];
     *objptr->idunion.mca_p2p1.s_II_II += weight*scp/((_REAL)5);
     objptr->idunion.mca_p2p1.s_II_II++;

       kk++;
      }
     }
    }
   }
      }
     }
    }
   }
   objptr->idunion.mca_p2p1.s_II_II -= objptr->idunion.mca_p2p1.dim_II_II;
   objptr->idunion.mca_p2p1.f_II_II -= objptr->idunion.mca_p2p1.dim_II_II;

   if(objptr->idunion.mca_p2p1.type == 2) return;

   l12=objptr->idunion.mca_p2p1.l12;
   l123456=objptr->idunion.mca_p2p1.l123456;
   for (i1=0,i=0,kk=0,l12_1=l12;i1<nmode;i1++) {
    for (i2=i1;i2<nmode;i2++,i++,l12_1++) {

   for (j1=0,j=0,l123456_2=l123456;j1<nmode1;j1++) {
    for (j2=j1;j2<nmode1;j2++) {
     for (j3=0;j3<nmode1;j3++) {
      for (j4=j3;j4<nmode1;j4++) {
       for (j5=0;j5<nmode1;j5++) {
        for (j6=j5;j6<nmode1;j6++,j++,l123456_2++,kk++) {

     for(m1=0;m1<5;m1++) {
      for(m2=0;m2<5;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<5;m1++) {
       {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j1].f[k][l]  *
          l123456_2->dphidm1[m2][l];
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j2].f[k][l]  *
          l123456_2->dphidm2[m2][l];
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j3].f[k][l]  *
          l123456_2->dphidm3[m2][l];
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j4].f[k][l]  *
          l123456_2->dphidm4[m2][l];
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j5].f[k][l]  *
          l123456_2->dphidm5[m2][l];
         f[m1][m2]+=
          l12_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j6].f[k][l]  *
          l123456_2->dphidm6[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j1].f[k][l]  *
          l123456_2->dphidm1[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j2].f[k][l]  *
          l123456_2->dphidm2[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j3].f[k][l]  *
          l123456_2->dphidm3[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j4].f[k][l]  *
          l123456_2->dphidm4[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j5].f[k][l]  *
          l123456_2->dphidm5[m2][l];
         f[m1][m2]+=
          l12_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j6].f[k][l]  *
          l123456_2->dphidm6[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     scp+=2.*f[3][3];
     scp+=2.*f[4][4];
     *objptr->idunion.mca_p2p1.f_I_III += weight*scp/((_REAL)5);
     objptr->idunion.mca_p2p1.f_I_III++;

     scp=0.;
     scp+=l12_1->phi[0]*l123456_2->phi[0];
     scp+=2.*l12_1->phi[1]*l123456_2->phi[1];
     scp+=2.*l12_1->phi[2]*l123456_2->phi[2];
     scp+=2.*l12_1->phi[3]*l123456_2->phi[3];
     scp+=2.*l12_1->phi[4]*l123456_2->phi[4];
     *objptr->idunion.mca_p2p1.s_I_III += weight*scp/((_REAL)5);
     objptr->idunion.mca_p2p1.s_I_III++;

        }
       }
      }
     }
    }
   }
    }
   }
   objptr->idunion.mca_p2p1.s_I_III -= objptr->idunion.mca_p2p1.dim_I_III;
   objptr->idunion.mca_p2p1.f_I_III -= objptr->idunion.mca_p2p1.dim_I_III;

   l1234=objptr->idunion.mca_p2p1.l1234;
   l123456=objptr->idunion.mca_p2p1.l123456;
   for (i1=0,i=0,kk=0,l1234_1=l1234;i1<nmode1;i1++) {
    for (i2=i1;i2<nmode1;i2++) {
     for (i3=0;i3<nmode1;i3++) {
      for (i4=i3;i4<nmode1;i4++,i++,l1234_1++) {

   for (j1=0,j=0,l123456_2=l123456;j1<nmode1;j1++) {
    for (j2=j1;j2<nmode1;j2++) {
     for (j3=0;j3<nmode1;j3++) {
      for (j4=j3;j4<nmode1;j4++) {
       for (j5=0;j5<nmode1;j5++) {
        for (j6=j5;j6<nmode1;j6++,j++,l123456_2++,kk++) {

     for(m1=0;m1<5;m1++) {
      for(m2=0;m2<5;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<5;m1++) {
       {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j1].f[k][l]  *
          l123456_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j2].f[k][l]  *
          l123456_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j3].f[k][l]  *
          l123456_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j4].f[k][l]  *
          l123456_2->dphidm4[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j5].f[k][l]  *
          l123456_2->dphidm5[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j6].f[k][l]  *
          l123456_2->dphidm6[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j1].f[k][l]  *
          l123456_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j2].f[k][l]  *
          l123456_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j3].f[k][l]  *
          l123456_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j4].f[k][l]  *
          l123456_2->dphidm4[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j5].f[k][l]  *
          l123456_2->dphidm5[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j6].f[k][l]  *
          l123456_2->dphidm6[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j1].f[k][l]  *
          l123456_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j2].f[k][l]  *
          l123456_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j3].f[k][l]  *
          l123456_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j4].f[k][l]  *
          l123456_2->dphidm4[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j5].f[k][l]  *
          l123456_2->dphidm5[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j6].f[k][l]  *
          l123456_2->dphidm6[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j1].f[k][l]  *
          l123456_2->dphidm1[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j2].f[k][l]  *
          l123456_2->dphidm2[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j3].f[k][l]  *
          l123456_2->dphidm3[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j4].f[k][l]  *
          l123456_2->dphidm4[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j5].f[k][l]  *
          l123456_2->dphidm5[m2][l];
         f[m1][m2]+=
          l1234_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j6].f[k][l]  *
          l123456_2->dphidm6[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     scp+=2.*f[3][3];
     scp+=2.*f[4][4];
     *objptr->idunion.mca_p2p1.f_II_III += weight*scp/((_REAL)5);
     objptr->idunion.mca_p2p1.f_II_III++;

     scp=0.;
     scp+=l1234_1->phi[0]*l123456_2->phi[0];
     scp+=2.*l1234_1->phi[1]*l123456_2->phi[1];
     scp+=2.*l1234_1->phi[2]*l123456_2->phi[2];
     scp+=2.*l1234_1->phi[3]*l123456_2->phi[3];
     scp+=2.*l1234_1->phi[4]*l123456_2->phi[4];
     *objptr->idunion.mca_p2p1.s_II_III += weight*scp/((_REAL)5);
     objptr->idunion.mca_p2p1.s_II_III++;

        }
       }
      }
     }
    }
   }
      }
     }
    }
   }
   objptr->idunion.mca_p2p1.s_II_III -= objptr->idunion.mca_p2p1.dim_II_III;
   objptr->idunion.mca_p2p1.f_II_III -= objptr->idunion.mca_p2p1.dim_II_III;

   l123456=objptr->idunion.mca_p2p1.l123456;
   for (i1=0,i=0,kk=0,l123456_1=l123456;i1<nmode1;i1++) {
    for (i2=i1;i2<nmode1;i2++) {
     for (i3=0;i3<nmode1;i3++) {
      for (i4=i3;i4<nmode1;i4++) {
       for (i5=0;i5<nmode1;i5++) {
        for (i6=i5;i6<nmode1;i6++,i++,l123456_1++) {
     j1_0=i1;
     j2_0=i2;
     j3_0=i3;
     j4_0=i4;
     j5_0=i5;
     j6_0=i6;
   for (j1=j1_0,j=i,l123456_2=l123456_1;j1<nmode1;j1++,j2_0=j1) {
    for (j2=j2_0;j2<nmode1;j2++,j4_0=0,j3_0=0) {
     for (j3=j3_0;j3<nmode1;j3++,j4_0=j3) {
      for (j4=j4_0;j4<nmode1;j4++,j5_0=0,j6_0=0) {
       for (j5=j5_0;j5<nmode1;j5++,j6_0=j5) {
        for (j6=j6_0;j6<nmode1;j6++,j++,l123456_2++) {

     for(m1=0;m1<5;m1++) {
      for(m2=0;m2<5;m2++) {
       f[m1][m2]=0.;
      }
     }

     for(m1=0;m1<5;m1++) {
       {
       m2=m1;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         f[m1][m2]+=
          l123456_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j1].f[k][l]  *
          l123456_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j2].f[k][l]  *
          l123456_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j3].f[k][l]  *
          l123456_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j4].f[k][l]  *
          l123456_2->dphidm4[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j5].f[k][l]  *
          l123456_2->dphidm5[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm1[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i1][j6].f[k][l]  *
          l123456_2->dphidm6[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j1].f[k][l]  *
          l123456_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j2].f[k][l]  *
          l123456_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j3].f[k][l]  *
          l123456_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j4].f[k][l]  *
          l123456_2->dphidm4[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j5].f[k][l]  *
          l123456_2->dphidm5[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm2[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i2][j6].f[k][l]  *
          l123456_2->dphidm6[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j1].f[k][l]  *
          l123456_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j2].f[k][l]  *
          l123456_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j3].f[k][l]  *
          l123456_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j4].f[k][l]  *
          l123456_2->dphidm4[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j5].f[k][l]  *
          l123456_2->dphidm5[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm3[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i3][j6].f[k][l]  *
          l123456_2->dphidm6[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j1].f[k][l]  *
          l123456_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j2].f[k][l]  *
          l123456_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j3].f[k][l]  *
          l123456_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j4].f[k][l]  *
          l123456_2->dphidm4[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j5].f[k][l]  *
          l123456_2->dphidm5[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm4[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i4][j6].f[k][l]  *
          l123456_2->dphidm6[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm5[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i5][j1].f[k][l]  *
          l123456_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm5[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i5][j2].f[k][l]  *
          l123456_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm5[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i5][j3].f[k][l]  *
          l123456_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm5[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i5][j4].f[k][l]  *
          l123456_2->dphidm4[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm5[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i5][j5].f[k][l]  *
          l123456_2->dphidm5[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm5[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i5][j6].f[k][l]  *
          l123456_2->dphidm6[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm6[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i6][j1].f[k][l]  *
          l123456_2->dphidm1[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm6[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i6][j2].f[k][l]  *
          l123456_2->dphidm2[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm6[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i6][j3].f[k][l]  *
          l123456_2->dphidm3[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm6[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i6][j4].f[k][l]  *
          l123456_2->dphidm4[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm6[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i6][j5].f[k][l]  *
          l123456_2->dphidm5[m2][l];
         f[m1][m2]+=
          l123456_1->dphidm6[m1][k] *
          objptr->idunion.mca_p2p1.Fm[i6][j6].f[k][l]  *
          l123456_2->dphidm6[m2][l];
        }
       }
      }
     }

     scp=0.;
     scp+=f[0][0];
     scp+=2.*f[1][1];
     scp+=2.*f[2][2];
     scp+=2.*f[3][3];
     scp+=2.*f[4][4];
     *objptr->idunion.mca_p2p1.f_III_III += weight*scp/((_REAL)5);
     objptr->idunion.mca_p2p1.f_III_III++;

     scp=0.;
     scp+=l123456_1->phi[0]*l123456_2->phi[0];
     scp+=2.*l123456_1->phi[1]*l123456_2->phi[1];
     scp+=2.*l123456_1->phi[2]*l123456_2->phi[2];
     scp+=2.*l123456_1->phi[3]*l123456_2->phi[3];
     scp+=2.*l123456_1->phi[4]*l123456_2->phi[4];
     *objptr->idunion.mca_p2p1.s_III_III += weight*scp/((_REAL)5);
     objptr->idunion.mca_p2p1.s_III_III++;

       kk++;
        }
       }
      }
     }
    }
   }
        }
       }
      }
     }
    }
   }
   objptr->idunion.mca_p2p1.s_III_III -= objptr->idunion.mca_p2p1.dim_III_III;
   objptr->idunion.mca_p2p1.f_III_III -= objptr->idunion.mca_p2p1.dim_III_III;

   return;
}

/************************************************************************/
void
write_first_p2p1(idstruct *dptr)
{
	int	i,j,i2,j2,k,l,m,bytes;
	int	states = dptr->idunion.mca_p2p1.cases;
	FILE 	*file = dptr->idunion.mca_p2p1.file,*dump;
        idstruct *mptr,*bptr;

	_REAL fact;
        int n12,n1234,n123456,nn,nf;

        n12=dptr->idunion.mca_p2p1.n12;
        n1234=dptr->idunion.mca_p2p1.n1234;
        n123456=dptr->idunion.mca_p2p1.n123456;
        mptr=dptr->idunion.mca_p2p1.mca;
        bptr=mptr->idunion.mcarig.bead;
	/* fact = 1./ ((_REAL)states); */
	fact = 1./ dptr->idunion.mca_p2p1.casesw;

        printf("write_first_p2p1: number of states: %d\n",states);

        for (i=1; i<=n12; i++) {
         for (j=i; j<=n12; j++) {
           *dptr->idunion.mca_p2p1.s_I_I++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1.s_I_I -= dptr->idunion.mca_p2p1.dim_I_I;
        for (i=1; i<=n12; i++) {
         for (j=i; j<=n12; j++) {
           *dptr->idunion.mca_p2p1.f_I_I++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1.f_I_I -= dptr->idunion.mca_p2p1.dim_I_I;

        if(dptr->idunion.mca_p2p1.type>1) {

        for (i=1; i<=n12; i++) {
         for (j=1; j<=n1234; j++) {
           *dptr->idunion.mca_p2p1.s_I_II++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1.s_I_II -= dptr->idunion.mca_p2p1.dim_I_II;
        for (i=1; i<=n12; i++) {
         for (j=1; j<=n1234; j++) {
           *dptr->idunion.mca_p2p1.f_I_II++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1.f_I_II -= dptr->idunion.mca_p2p1.dim_I_II;
        for (i=1; i<=n1234; i++) {
         for (j=i; j<=n1234; j++) {
           *dptr->idunion.mca_p2p1.s_II_II++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1.s_II_II -= dptr->idunion.mca_p2p1.dim_II_II;
        for (i=1; i<=n1234; i++) {
         for (j=i; j<=n1234; j++) {
           *dptr->idunion.mca_p2p1.f_II_II++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1.f_II_II -= dptr->idunion.mca_p2p1.dim_II_II;

        }

        if(dptr->idunion.mca_p2p1.type>2) {

        for (i=1; i<=n12; i++) {
         for (j=1; j<=n123456; j++) {
           *dptr->idunion.mca_p2p1.s_I_III++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1.s_I_III -= dptr->idunion.mca_p2p1.dim_I_III;
        for (i=1; i<=n12; i++) {
         for (j=1; j<=n123456; j++) {
           *dptr->idunion.mca_p2p1.f_I_III++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1.f_I_III -= dptr->idunion.mca_p2p1.dim_I_III;
        for (i=1; i<=n1234; i++) {
         for (j=1; j<=n123456; j++) {
           *dptr->idunion.mca_p2p1.s_II_III++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1.s_II_III -= dptr->idunion.mca_p2p1.dim_II_III;
        for (i=1; i<=n1234; i++) {
         for (j=1; j<=n123456; j++) {
           *dptr->idunion.mca_p2p1.f_II_III++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1.f_II_III -= dptr->idunion.mca_p2p1.dim_II_III;
        for (i=1; i<=n123456; i++) {
         for (j=i; j<=n123456; j++) {
           *dptr->idunion.mca_p2p1.s_III_III++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1.s_III_III -= dptr->idunion.mca_p2p1.dim_III_III;
        for (i=1; i<=n123456; i++) {
         for (j=i; j<=n123456; j++) {
           *dptr->idunion.mca_p2p1.f_III_III++ *= fact;
         }
        }
        dptr->idunion.mca_p2p1.f_III_III -= dptr->idunion.mca_p2p1.dim_III_III;

        }

        if(dptr->idunion.mca_p2p1.dump==1) {
         for (i=1; i<=n12; i++) 
         {
          for (j=i; j<=n12; j++) 
          {
           printf("F(%d,%d)= %e\n",i,j,*dptr->idunion.mca_p2p1.f_I_I++);
          }
         }
         dptr->idunion.mca_p2p1.f_I_I -= dptr->idunion.mca_p2p1.dim_I_I;

         dump=fopen("fs.dump","w");
         bytes=sizeof(_REAL)*dptr->idunion.mca_p2p1.dim_I_I;
         fwrite(&bytes,sizeof(int),1,dump);
         fwrite(dptr->idunion.mca_p2p1.f_I_I,dptr->idunion.mca_p2p1.dim_I_I*sizeof(_REAL),
                1,dump);
         fwrite(&bytes,sizeof(int),1,dump);
         fwrite(&bytes,sizeof(int),1,dump);
         fwrite(dptr->idunion.mca_p2p1.s_I_I,dptr->idunion.mca_p2p1.dim_I_I*sizeof(_REAL),
                1,dump);
         fwrite(&bytes,sizeof(int),1,dump);
         if(dptr->idunion.mca_p2p1.type>1)
         {
          bytes=sizeof(_REAL)*dptr->idunion.mca_p2p1.dim_I_II;
          fwrite(&bytes,sizeof(int),1,dump);
          fwrite(dptr->idunion.mca_p2p1.f_I_II,dptr->idunion.mca_p2p1.dim_I_II*sizeof(_REAL),
                 1,dump);
          fwrite(&bytes,sizeof(int),1,dump);
          fwrite(&bytes,sizeof(int),1,dump);
          fwrite(dptr->idunion.mca_p2p1.s_I_II,dptr->idunion.mca_p2p1.dim_I_II*sizeof(_REAL),
                 1,dump);
          fwrite(&bytes,sizeof(int),1,dump);
          bytes=sizeof(_REAL)*dptr->idunion.mca_p2p1.dim_II_II;
          fwrite(&bytes,sizeof(int),1,dump);
          fwrite(dptr->idunion.mca_p2p1.f_II_II,dptr->idunion.mca_p2p1.dim_II_II*sizeof(_REAL),
                 1,dump);
          fwrite(&bytes,sizeof(int),1,dump);
          fwrite(&bytes,sizeof(int),1,dump);
          fwrite(dptr->idunion.mca_p2p1.s_II_II,dptr->idunion.mca_p2p1.dim_II_II*sizeof(_REAL),
                 1,dump);
          fwrite(&bytes,sizeof(int),1,dump);
         }
         if(dptr->idunion.mca_p2p1.type>2)
         {
          printf("write_first_p2p1: 3rd order dump not implemented yet.\n");
	  exit(1);
	 }
         fflush(dump);
         fclose(dump);
        }

        printf("S/F matrix elements written as binary fortran.\n");
        /*
        */
        printf("Now code stops before relaxation modes.\n");
        return;
}
