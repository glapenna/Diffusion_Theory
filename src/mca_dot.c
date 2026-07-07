#include "head.h"
#include <time.h>

// extern char tok[];
// char tok[];
// void *get2(int,int);

extern idstruct *startid(),*getstruct();
rotstruct **rot_matrix();
rot1struct **rot1_matrix();

/************************************************************************/

void
init_mca_dot()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(MCA_DOT, id);

	dptr->idunion.mca_dot.cases = 0;
	dptr->idunion.mca_dot.casesw = 0.;
	dptr->idunion.mca_dot.dump = 0;
	dptr->idunion.mca_dot.type = 1;
	dptr->idunion.mca_dot.list = 0;
	dptr->idunion.mca_dot.n1 = 3;
	dptr->idunion.mca_dot.n12 = 0;
	dptr->idunion.mca_dot.m1 = 0;
	dptr->idunion.mca_dot.m12 = 0;
	dptr->idunion.mca_dot.sel_bonds = 0;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("MCA_DOT: expected file name","");
	strcpy(dptr->idunion.mca_dot.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.mca_dot.fname, 0, &exists, "MCA_DOT"))
		exit(1);
	if ((dptr->idunion.mca_dot.file = 
			 fopen(dptr->idunion.mca_dot.fname, "w")) == NULL){
		perror(dptr->idunion.mca_dot.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************/
void
setup_mca_dot()
{
	idstruct	*dptr,*bptr;
	int i,j,k,icount,nt,offset_r;
        char fname_comb[TOKLEN];
        FILE *file_comb,*ftmp;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != MCA_DOT)
                        inerr("setup_mca_dot: expected id or keyword\n","");
        } else {
                printf("setup_mca_dot: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((bptr = getstruct(tok)) == NULL) 
                inerr("setup_mca_dot: no such id\n","");
        if (bptr->id_type != BEAD)
                inerr("setup_mca_dot: bead id expected\n","");

        dptr->idunion.mca_dot.m1=0;

        gettok();
        for(;;) {
         if(tok[0]==';') {
          break;
         } else if(!strcmp(tok,"NBONDS")) {
          gettok();
          tokint(&dptr->idunion.mca_dot.sel_bonds,tok);
          gettok();
         } else if(!strcmp(tok,"LIST")) {
          dptr->idunion.mca_dot.list=1;
          gettok();
         } else if(!strcmp(tok,"RIGID")) {
          dptr->idunion.mca_dot.rigid=1;
          gettok();
         } else if(!strcmp(tok,"NMODES")) {
          gettok();
          tokint(&dptr->idunion.mca_dot.m1,tok);
          dptr->idunion.mca_dot.type=2;
          gettok();
         } else if(!strcmp(tok,"FILE_COMB")) {
          gettok();
	  strcpy(fname_comb, tok);
	  if ((file_comb = fopen(fname_comb, "r")) == NULL){
           perror(fname_comb);
           exit(1);
	  }
          gettok();
         } else {
          printf("setup_mca_dot: option not found.\n");
          exit(1);
         }
        }

        dptr->idunion.mca_dot.bead = bptr;
        dptr->idunion.mca_dot.n1=bptr->idunion.bead.nbond;

        nt=0;
        if(dptr->idunion.mca_dot.sel_bonds!=0)
        {

        if(dptr->idunion.mca_dot.list==0)
        {

         bptr->idunion.bead.l->offset_l = (int *)
           get2(dptr->idunion.mca_dot.n1,sizeof(int));
         bptr->idunion.bead.l->doffset = (int *)
           get2(dptr->idunion.mca_dot.n1,sizeof(int));

         nt=0;
         for(i=0;i<bptr->idunion.bead.nbond;i++) {
          bptr->idunion.bead.l->offset_l[i]= 0;
                             offset_r   = dptr->idunion.mca_dot.sel_bonds -1;
          bptr->idunion.bead.l->doffset[i]=  offset_r
                                         +bptr->idunion.bead.l->offset_l[i] +1;
 
          if(i+offset_r > bptr->idunion.bead.nbond-1)
          {
                             offset_r   -=i+1
                                         +dptr->idunion.mca_dot.sel_bonds-1
                                         -bptr->idunion.bead.nbond;
          }
          bptr->idunion.bead.l->doffset[i]=  offset_r
                                        +bptr->idunion.bead.l->offset_l[i] +1;
/*
         printf("offset(%d)= %d %d\n",
          i,bptr->idunion.bead.l->offset_l[i],bptr->idunion.bead.l->doffset[i]);
         fflush(stdout);
*/

          nt+=bptr->idunion.bead.l->doffset[i]*
              (bptr->idunion.bead.l->doffset[i]+1)/2;


         }
        } else {

          ftmp=fopen("j1j2.selected","r");
          nt=0;
          while( !feof(ftmp) ) 
          {
           fscanf(ftmp,"%*s %*d %*d %*d %*d ");
           nt++;
          }
          printf("setup_mca_dot: reading %d entries from file j1j2j3.selected\n",nt);
          fflush(stdout);
          dptr->idunion.mca_dot.nlist=nt;
          dptr->idunion.mca_dot.j1=(int *)get2(nt,sizeof(int));
          dptr->idunion.mca_dot.j2=(int *)get2(nt,sizeof(int));
          rewind(ftmp);
          nt=0;
          while( !feof(ftmp) ) 
          {
           fscanf(ftmp,"%*s %*d %d %d ",
              &dptr->idunion.mca_dot.j1[nt],
              &dptr->idunion.mca_dot.j2[nt]);
           dptr->idunion.mca_dot.j1[nt]--;
           dptr->idunion.mca_dot.j2[nt]--;
           nt++;
          }
          if(nt!=dptr->idunion.mca_dot.nlist)
          {
           printf("setup_mca_dot: error reading file j1j2.selected.\n");
           exit(1);
          }
          fclose(ftmp);

        }
        }
        dptr->idunion.mca_dot.n12=nt;
        if(dptr->idunion.mca_dot.type>1)
        {
         dptr->idunion.mca_dot.m12= dptr->idunion.mca_dot.m1*
          (dptr->idunion.mca_dot.m1+1)/2;
        }

        dptr->idunion.mca_dot.dim_1_l12=
         dptr->idunion.mca_dot.n12;
        dptr->idunion.mca_dot.dim_1_L12=
         dptr->idunion.mca_dot.m12;
        dptr->idunion.mca_dot.dim_l12_l12=
         dptr->idunion.mca_dot.n12*(dptr->idunion.mca_dot.n12+1)/2;
        dptr->idunion.mca_dot.dim_l12_L12=
         dptr->idunion.mca_dot.n12*dptr->idunion.mca_dot.m12;
        dptr->idunion.mca_dot.dim_L12_L12=
         dptr->idunion.mca_dot.m12*(dptr->idunion.mca_dot.m12+1)/2;

        if(dptr->idunion.mca_dot.sel_bonds!=0) {

        dptr->idunion.mca_dot.s_1_l12 = (_REAL *)
          get2(dptr->idunion.mca_dot.dim_1_l12,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_dot.dim_1_l12;i++) {
         *dptr->idunion.mca_dot.s_1_l12++ =0.;
        }
        dptr->idunion.mca_dot.s_1_l12-=dptr->idunion.mca_dot.dim_1_l12;

        dptr->idunion.mca_dot.s_l12_l12 = (_REAL *)
          get2(dptr->idunion.mca_dot.dim_l12_l12,sizeof(_REAL));
        dptr->idunion.mca_dot.f_l12_l12 = (_REAL *)
          get2(dptr->idunion.mca_dot.dim_l12_l12,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_dot.dim_l12_l12;i++) {
         *dptr->idunion.mca_dot.s_l12_l12++ =0.;
         *dptr->idunion.mca_dot.f_l12_l12++ =0.;
        }
        dptr->idunion.mca_dot.s_l12_l12-=dptr->idunion.mca_dot.dim_l12_l12;
        dptr->idunion.mca_dot.f_l12_l12-=dptr->idunion.mca_dot.dim_l12_l12;
	}

        if(dptr->idunion.mca_dot.type>1) {

        dptr->idunion.mca_dot.s_1_L12 = (_REAL *)
          get2(dptr->idunion.mca_dot.dim_1_L12,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_dot.dim_1_L12;i++) {
         *dptr->idunion.mca_dot.s_1_L12++ =0.;
        }
        dptr->idunion.mca_dot.s_1_L12-=dptr->idunion.mca_dot.dim_1_L12;

        if(dptr->idunion.mca_dot.sel_bonds!=0) {

        dptr->idunion.mca_dot.s_l12_L12 = (_REAL *)
          get2(dptr->idunion.mca_dot.dim_l12_L12,sizeof(_REAL));
        dptr->idunion.mca_dot.f_l12_L12 = (_REAL *)
          get2(dptr->idunion.mca_dot.dim_l12_L12,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_dot.dim_l12_L12;i++) {
         *dptr->idunion.mca_dot.s_l12_L12++ =0.;
         *dptr->idunion.mca_dot.f_l12_L12++ =0.;
        }
        dptr->idunion.mca_dot.s_l12_L12-=dptr->idunion.mca_dot.dim_l12_L12;
        dptr->idunion.mca_dot.f_l12_L12-=dptr->idunion.mca_dot.dim_l12_L12;

        }

        dptr->idunion.mca_dot.s_L12_L12 = (_REAL *)
          get2(dptr->idunion.mca_dot.dim_L12_L12,sizeof(_REAL));
        dptr->idunion.mca_dot.f_L12_L12 = (_REAL *)
          get2(dptr->idunion.mca_dot.dim_L12_L12,sizeof(_REAL));
        for(i=0;i<dptr->idunion.mca_dot.dim_L12_L12;i++) {
         *dptr->idunion.mca_dot.s_L12_L12++ =0.;
         *dptr->idunion.mca_dot.f_L12_L12++ =0.;
        }
        dptr->idunion.mca_dot.s_L12_L12-=dptr->idunion.mca_dot.dim_L12_L12;
        dptr->idunion.mca_dot.f_L12_L12-=dptr->idunion.mca_dot.dim_L12_L12;

        }


        dptr->idunion.mca_dot.Fb = rot_matrix
         (0,bptr->idunion.bead.nbond-1,0,bptr->idunion.bead.nbond-1);

        if(dptr->idunion.mca_dot.sel_bonds) 
	{
         dptr->idunion.mca_dot.l12 = (dotstruct *)
          get2(dptr->idunion.mca_dot.n12,sizeof(dotstruct));
        }
        if(dptr->idunion.mca_dot.type>1) 
	{
        dptr->idunion.mca_dot.m = (modstruct *)
          get2(dptr->idunion.mca_dot.m1,sizeof(modstruct));
        dptr->idunion.mca_dot.L12 = (dotstruct *)
          get2(dptr->idunion.mca_dot.m12,sizeof(dotstruct));

        dptr->idunion.mca_dot.comb = matrix
         (1,bptr->idunion.bead.nbond,1,dptr->idunion.mca_dot.m1);
        for(j=1;j<=dptr->idunion.mca_dot.m1;j++)
        {
         for(i=1;i<=bptr->idunion.bead.nbond;i++)
         {
          fscanf(file_comb,"%*d %*d %le ",&dptr->idunion.mca_dot.comb[i][j]);  
         }
        }
        fclose(file_comb); 
        dptr->idunion.mca_dot.Fbm = rot_matrix
         (0,bptr->idunion.bead.nbond-1,0,dptr->idunion.mca_dot.m1-1);
        dptr->idunion.mca_dot.Fm = rot_matrix
         (0,dptr->idunion.mca_dot.m1-1,0,dptr->idunion.mca_dot.m1-1);
        }

        printf("setup_mca_dot: n1:    %d\n",dptr->idunion.mca_dot.n1);
        printf("setup_mca_dot: n12:   %d\n",dptr->idunion.mca_dot.n12);
        printf("setup_mca_dot: m1:    %d\n",dptr->idunion.mca_dot.m1);
        printf("setup_mca_dot: m12:   %d\n",dptr->idunion.mca_dot.m12);
        printf("setup_mca_dot: dim_1_l12:      %d\n",
         dptr->idunion.mca_dot.dim_1_l12);
        printf("setup_mca_dot: dim_1_L12:      %d\n",
         dptr->idunion.mca_dot.dim_1_L12);
        printf("setup_mca_dot: dim_l12_l12:    %d\n",
         dptr->idunion.mca_dot.dim_l12_l12);
        printf("setup_mca_dot: dim_l12_L12:    %d\n",
         dptr->idunion.mca_dot.dim_l12_L12);
        printf("setup_mca_dot: dim_L12_L12:    %d\n",
         dptr->idunion.mca_dot.dim_L12_L12);
}

/************************************************************************/
do_first_dot(idstruct *objptr,idstruct *serptr)
{
   _REAL *box=serptr->idunion.ser.set->boxc;
   idstruct *bptr;
   int i1,i2,i1_0,i2_0,j1_0,j2_0,j1,j2,
       i1_1,i2_1,j1_1,j2_1,
       i,j,k,l,ib,jb,kk,im,m1,m2,nbond,nmode,overlap,
	       *offset_l,*doffset;
   _REAL **r1, **r2;
   _REAL *rt1,*rt2,r,zi,zj,zr,cutsq,
         boxm1[3],x1,y1,z1,x2,y2,z2,dx,dy,dz,scp;
   bondstruct *lb,*lb1,*lb2,*lb3;
   modstruct *p,*p1,*p2,*p3;
   dotstruct *l12,*l12_1,*l12_2,*L12,*L12_1,*L12_2;
   _REAL xm,ym,zm,f[3][3],weight;

   bptr=objptr->idunion.mca_dot.bead;
   zr=bptr->idunion.bead.zr;
   cutsq=bptr->idunion.bead.cutsq;
   nbond=objptr->idunion.mca_dot.n1;
   nmode=objptr->idunion.mca_dot.m1;
   offset_l=bptr->idunion.bead.l->offset_l;
   doffset=bptr->idunion.bead.l->doffset;

   for(k=0;k<3;k++) {
    boxm1[k]=1./box[k];
   }

   if(serptr->idunion.ser.iweight)
   {            
    weight=serptr->idunion.ser.w->idunion.serw.w;
   } else {
    weight=1.;
   }

   ++objptr->idunion.mca_dot.cases;
   objptr->idunion.mca_dot.casesw+=weight;

   lb=bptr->idunion.bead.l;
   r1=bptr->idunion.bead.r1;
   r2=bptr->idunion.bead.r2;
   for (i=0;i<nbond;i++,r1++,r2++) {
    rt1=*r1; rt2=*r2;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
    x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
    dx=x2-x1; dy=y2-y1; dz=z2-z1;
    lb->id1 = (int)(*r1);
    lb->id2 = (int)(*r2);
    lb->x1 = x1; lb->y1 = y1; lb->z1 = z1;
    lb->x2 = x2; lb->y2 = y2; lb->z2 = z2;
#ifdef PBC
    dx -= box[0]*rint(dx*boxm1[0]);
    dy -= box[1]*rint(dy*boxm1[1]);
    dz -= box[2]*rint(dz*boxm1[2]);
#endif 
    lb->x = dx;
    lb->y = dy;
    lb->z = dz;
    lb++;
   }

   if(objptr->idunion.mca_dot.sel_bonds)
   {

   l12=objptr->idunion.mca_dot.l12;
   lb=bptr->idunion.bead.l;
   for (i1=0,lb1=lb;i1<nbond;i1++,lb1++) {
    for (i2=i1,lb2=lb1;i2<i1+doffset[i1];i2++,lb2++) {
       scp=lb1->x*lb2->x+lb1->y*lb2->y+lb1->z*lb2->z;
       l12->phi=scp;
       l12->dphidm1[0]=lb2->x;
       l12->dphidm1[1]=lb2->y;
       l12->dphidm1[2]=lb2->z;
       l12->dphidm2[0]=lb1->x;
       l12->dphidm2[1]=lb1->x;
       l12->dphidm2[2]=lb1->x;

/*
  printf("i= %d,%d\n",i1+1,i2+1);
  printf("l12= %e\n",l12->phi);
  fflush(stdout);
*/
      l12++;
    }
   }

   }

   if(objptr->idunion.mca_dot.type>1)
   {

   p=objptr->idunion.mca_dot.m;
   for (j=0;j<nmode;j++) {
    xm=0.; ym=0.; zm=0.;
    r1=bptr->idunion.bead.r1;
    r2=bptr->idunion.bead.r2;
    for (i=0;i<nbond;i++,r1++,r2++) {
     rt1=*r1; rt2=*r2;
     x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
     x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
     dx=x2-x1; dy=y2-y1; dz=z2-z1;
     xm += objptr->idunion.mca_dot.comb[i+1][j+1]*dx;
     ym += objptr->idunion.mca_dot.comb[i+1][j+1]*dy;
     zm += objptr->idunion.mca_dot.comb[i+1][j+1]*dz;
    }
    p->mod=l_sqrt(xm*xm+ym*ym+zm*zm);
    p->m[0]=xm;
    p->m[1]=ym;
    p->m[2]=zm;
    p++;
   }

   L12=objptr->idunion.mca_dot.L12;
   p=objptr->idunion.mca_dot.m;
   for (i1=0,p1=p;i1<nmode;i1++,p1++) {
    for (i2=i1,p2=p1;i2<nmode;i2++,p2++) {
      scp=p1->m[0]*p2->m[0]+p1->m[1]*p2->m[1]+p1->m[2]*p2->m[2];
      L12->phi=scp;
      L12->dphidm1[0]=p2->m[0];
      L12->dphidm1[1]=p2->m[1];
      L12->dphidm1[2]=p2->m[2];
      L12->dphidm2[0]=p1->m[0];
      L12->dphidm2[1]=p1->m[1];
      L12->dphidm2[2]=p1->m[2];
      L12++;
    }
   }

   }

   lb=bptr->idunion.bead.l;
   for (ib=0,kk=0,lb1=lb;ib<nbond;ib++,lb1++) {
    for (jb=ib,lb2=lb1;jb<nbond;jb++,lb2++,kk++) {
#include "fbb.c"

     objptr->idunion.mca_dot.Fb[ib][jb].f[0][0]=f[0][0];
     objptr->idunion.mca_dot.Fb[ib][jb].f[0][1]=f[0][1];
     objptr->idunion.mca_dot.Fb[ib][jb].f[0][2]=f[0][2];
     objptr->idunion.mca_dot.Fb[ib][jb].f[1][0]=f[0][1];
     objptr->idunion.mca_dot.Fb[ib][jb].f[1][1]=f[1][1];
     objptr->idunion.mca_dot.Fb[ib][jb].f[1][2]=f[1][2];
     objptr->idunion.mca_dot.Fb[ib][jb].f[2][0]=f[0][2];
     objptr->idunion.mca_dot.Fb[ib][jb].f[2][1]=f[1][2];
     objptr->idunion.mca_dot.Fb[ib][jb].f[2][2]=f[2][2];

     objptr->idunion.mca_dot.Fb[jb][ib].f[0][0]=f[0][0];
     objptr->idunion.mca_dot.Fb[jb][ib].f[0][1]=f[0][1];
     objptr->idunion.mca_dot.Fb[jb][ib].f[0][2]=f[0][2];
     objptr->idunion.mca_dot.Fb[jb][ib].f[1][0]=f[0][1];
     objptr->idunion.mca_dot.Fb[jb][ib].f[1][1]=f[1][1];
     objptr->idunion.mca_dot.Fb[jb][ib].f[1][2]=f[1][2];
     objptr->idunion.mca_dot.Fb[jb][ib].f[2][0]=f[0][2];
     objptr->idunion.mca_dot.Fb[jb][ib].f[2][1]=f[1][2];
     objptr->idunion.mca_dot.Fb[jb][ib].f[2][2]=f[2][2];

    }
   }

   if(objptr->idunion.mca_dot.sel_bonds)
   {

   l12=objptr->idunion.mca_dot.l12;
   i1_0=0;
   i1_1=nbond;
   for (i1=i1_0,i=0,kk=0,l12_1=l12;i1<i1_1;i1++) {
    i2_0=i1;
    i2_1=i2_0+doffset[i1];
    for (i2=i2_0;i2<i2_1;i2++,l12_1++,i++) {

     scp=l12_1->phi;
     *objptr->idunion.mca_dot.s_1_l12 += weight*scp;
     objptr->idunion.mca_dot.s_1_l12++;

    }
   }
   objptr->idunion.mca_dot.s_1_l12 -= objptr->idunion.mca_dot.dim_1_l12;

   l12=objptr->idunion.mca_dot.l12;
   i1_0=0;
   i1_1=nbond;
   for (i1=i1_0,i=0,kk=0,l12_1=l12;i1<i1_1;i1++) {
    i2_0=i1;
    i2_1=i2_0+doffset[i1];
    for (i2=i2_0;i2<i2_1;i2++,l12_1++,i++) {

     j1_0=i1;
     j2_0=i2;
     j1_1=i1_1;
     j2_1=i2_1;

   for (j1=j1_0,j=i,l12_2=l12_1;j1<j1_1;j1++,j2_0=j1,j2_1=j2_0+doffset[j1])
   {
    for (j2=j2_0;j2<j2_1;j2++,l12_2++,j++,kk++)
    {

       scp = 0.;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         scp +=
          l12_1->dphidm1[k] *
          objptr->idunion.mca_dot.Fb[i1][j1].f[k][l]  *
          l12_2->dphidm1[l];
         scp +=
          l12_1->dphidm1[k] *
          objptr->idunion.mca_dot.Fb[i1][j1].f[k][l]  *
          l12_2->dphidm2[l];
         scp +=
          l12_1->dphidm2[k] *
          objptr->idunion.mca_dot.Fb[i1][j1].f[k][l]  *
          l12_2->dphidm1[l];
         scp +=
          l12_1->dphidm2[k] *
          objptr->idunion.mca_dot.Fb[i1][j2].f[k][l]  *
          l12_2->dphidm2[l];

        }
       }

     *objptr->idunion.mca_dot.f_l12_l12 += weight*scp;
     objptr->idunion.mca_dot.f_l12_l12++;

     scp=l12_1->phi*l12_2->phi;
     *objptr->idunion.mca_dot.s_l12_l12 += weight*scp;
     objptr->idunion.mca_dot.s_l12_l12++;

    }
   }
    }
   }
   objptr->idunion.mca_dot.s_l12_l12 -= objptr->idunion.mca_dot.dim_l12_l12;
   objptr->idunion.mca_dot.f_l12_l12 -= objptr->idunion.mca_dot.dim_l12_l12;

// printf("kk(l123_l123)= %d\n",kk);

   }

   if(objptr->idunion.mca_dot.type==1) return;

   for (i=0;i<nmode;i++) {
    for (j=0;j<nmode;j++) {
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       f[m1][m2]=0.;
      }
     }
     for (ib=0;ib<nbond;ib++) {
      for (jb=0;jb<nbond;jb++) {
       for(m1=0;m1<3;m1++) {
        for(m2=0;m2<3;m2++) {
         f[m1][m2]+=
          objptr->idunion.mca_dot.comb[ib+1][i+1] *
          objptr->idunion.mca_dot.Fb[ib][jb].f[m1][m2] *
          objptr->idunion.mca_dot.comb[jb+1][j+1];
        }
       }
      }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_dot.Fm[i][j].f[m1][m2]=f[m1][m2];
      }
     }
    }
   }
   for (i=0;i<nbond;i++) {
    for (j=0;j<nmode;j++) {
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       f[m1][m2]=0.;
      }
     }
     for (jb=0;jb<nbond;jb++) {
       for(m1=0;m1<3;m1++) {
        for(m2=0;m2<3;m2++) {
         f[m1][m2]+=
          objptr->idunion.mca_dot.Fb[i][jb].f[m1][m2] *
          objptr->idunion.mca_dot.comb[jb+1][j+1];
        }
       }
     }
     for(m1=0;m1<3;m1++) {
      for(m2=0;m2<3;m2++) {
       objptr->idunion.mca_dot.Fbm[i][j].f[m1][m2]=f[m1][m2];
      }
     }
    }
   }

   if(objptr->idunion.mca_dot.sel_bonds)
   {

   l12=objptr->idunion.mca_dot.l12;
   L12=objptr->idunion.mca_dot.L12;
   i1_0=0;
   i1_1=nbond;
   for (i1=i1_0,i=0,kk=0,l12_1=l12;i1<i1_1;i1++) {
    i2_0=i1;
    i2_1=i2_0+doffset[i1];
    for (i2=i2_0;i2<i2_1;i2++,l12_1++,i++) {

   for (j1=0,j=0,L12_2=L12;j1<nmode;j1++) {
     for (j2=j1;j2<nmode;j2++,L12_2++,j++,kk++) {

       scp = 0.;      
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         scp +=
          l12_1->dphidm1[k] *
          objptr->idunion.mca_dot.Fbm[i1][j1].f[k][l]  *
          L12_2->dphidm1[l];
         scp +=
          l12_1->dphidm1[k] *
          objptr->idunion.mca_dot.Fbm[i1][j1].f[k][l]  *
          L12_2->dphidm2[l];
         scp +=
          l12_1->dphidm2[k] *
          objptr->idunion.mca_dot.Fbm[i1][j1].f[k][l]  *
          L12_2->dphidm1[l];
         scp +=
          l12_1->dphidm2[k] *
          objptr->idunion.mca_dot.Fbm[i1][j1].f[k][l]  *
          L12_2->dphidm2[l];

        }
       }

     *objptr->idunion.mca_dot.f_l12_L12 += weight*scp;
     objptr->idunion.mca_dot.f_l12_L12++;

     scp=l12_1->phi*L12_2->phi;
     *objptr->idunion.mca_dot.s_l12_L12 += weight*scp;
     objptr->idunion.mca_dot.s_l12_L12++;

    }
   }

    }
   }

   objptr->idunion.mca_dot.s_l12_L12 -= objptr->idunion.mca_dot.dim_l12_L12;
   objptr->idunion.mca_dot.f_l12_L12 -= objptr->idunion.mca_dot.dim_l12_L12;

// printf("kk(l123_L123)= %d\n",kk);

   }

   L12=objptr->idunion.mca_dot.L12;
   for (i1=0,i=0,kk=0,L12_1=L12;i1<nmode;i1++) {
    for (i2=i1;i2<nmode;i2++,L12_1++,i++) {

     scp=L12_1->phi;
     *objptr->idunion.mca_dot.s_1_L12 += weight*scp;
     objptr->idunion.mca_dot.s_1_L12++;

    }
   }
   objptr->idunion.mca_dot.s_1_L12 -= objptr->idunion.mca_dot.dim_1_L12;

   L12=objptr->idunion.mca_dot.L12;
   for (i1=0,i=0,kk=0,L12_1=L12;i1<nmode;i1++) {
    for (i2=i1;i2<nmode;i2++,L12_1++,i++) {
     j1_0=i1;
     j2_0=i2;
   for (j1=j1_0,j=i,L12_2=L12_1;j1<nmode;j1++,j2_0=j1) {
    for (j2=j2_0;j2<nmode;j2++,L12_2++,j++,kk++) {

       scp = 0.;
       for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
         scp +=
          L12_1->dphidm1[k] *
          objptr->idunion.mca_dot.Fm[i1][j1].f[k][l]  *
          L12_2->dphidm1[l];
         scp +=
          L12_1->dphidm1[k] *
          objptr->idunion.mca_dot.Fm[i1][j1].f[k][l]  *
          L12_2->dphidm2[l];
         scp +=
          L12_1->dphidm2[k] *
          objptr->idunion.mca_dot.Fm[i1][j1].f[k][l]  *
          L12_2->dphidm1[l];
         scp +=
          L12_1->dphidm2[k] *
          objptr->idunion.mca_dot.Fm[i1][j1].f[k][l]  *
          L12_2->dphidm2[l];

        }
       }

     *objptr->idunion.mca_dot.f_L12_L12 += weight*scp;
     objptr->idunion.mca_dot.f_L12_L12++;

     scp=L12_1->phi*L12_2->phi;
     *objptr->idunion.mca_dot.s_L12_L12 += weight*scp;
     objptr->idunion.mca_dot.s_L12_L12++;

    }
   }
    }
   }
   objptr->idunion.mca_dot.s_L12_L12 -= objptr->idunion.mca_dot.dim_L12_L12;
   objptr->idunion.mca_dot.f_L12_L12 -= objptr->idunion.mca_dot.dim_L12_L12;

 printf("kk(L123_L123)= %d\n",kk);

}

/************************************************************************/
void
write_first_dot(idstruct *dptr)
{
	int	i,j,i2,j2,k,l,m;
	int	states = dptr->idunion.mca_dot.cases;
	_REAL   statesw = dptr->idunion.mca_dot.casesw;
	FILE 	*file = dptr->idunion.mca_dot.file,*dump;
        idstruct *bptr;

	_REAL fact,tol,factns,sum,sum2;
        int rank,nbondm;

        int n1,n12,m1,m12,nn,nf;
        _REAL **S;
        _REAL **F,**Fm1,**tmp,csum,ctmp;
        _REAL *ui,*v,*w;
#ifdef IMSL
        d_complex *wi,*vi;
#else
        _REAL *vi;
#endif

        n1=dptr->idunion.mca_dot.n1;
        n12=dptr->idunion.mca_dot.n12;
        m1=dptr->idunion.mca_dot.m1;
        m12=dptr->idunion.mca_dot.m12;
        bptr=dptr->idunion.mca_dot.bead;
	fact = 1./ statesw;
        factns = bptr->idunion.bead.factns;

        printf("write_first_dot: number of states: %d\n",states);

        if(dptr->idunion.mca_dot.sel_bonds)
        {

        for (j=1; j<=n12; j++) {
          *dptr->idunion.mca_dot.s_1_l12++ *= fact;
        }
        dptr->idunion.mca_dot.s_1_l12 -= dptr->idunion.mca_dot.dim_1_l12;

        for (i=1; i<=n12; i++) {
         for (j=i; j<=n12; j++) {
           *dptr->idunion.mca_dot.s_l12_l12++ *= fact;
           *dptr->idunion.mca_dot.f_l12_l12++ *= fact;
         }
        }
        dptr->idunion.mca_dot.s_l12_l12 -= dptr->idunion.mca_dot.dim_l12_l12;
        dptr->idunion.mca_dot.f_l12_l12 -= dptr->idunion.mca_dot.dim_l12_l12;

        }

        if(dptr->idunion.mca_dot.type>1) {

        if(dptr->idunion.mca_dot.sel_bonds) {

        for (i=1; i<=n12; i++) {
         for (j=1; j<=m12; j++) {
           *dptr->idunion.mca_dot.s_l12_L12++ *= fact;
           *dptr->idunion.mca_dot.f_l12_L12++ *= fact;
         }
        }
        dptr->idunion.mca_dot.s_l12_L12 -= dptr->idunion.mca_dot.dim_l12_L12;
        dptr->idunion.mca_dot.f_l12_L12 -= dptr->idunion.mca_dot.dim_l12_L12;

        }

        for (j=1; j<=m12; j++) {
          *dptr->idunion.mca_dot.s_1_L12++ *= fact;
        }
        dptr->idunion.mca_dot.s_1_L12 -= dptr->idunion.mca_dot.dim_1_L12;

        for (i=1; i<=m12; i++) {
         for (j=i; j<=m12; j++) {
           *dptr->idunion.mca_dot.s_L12_L12++ *= fact;
           *dptr->idunion.mca_dot.f_L12_L12++ *= fact;
         }
        }
        dptr->idunion.mca_dot.s_L12_L12 -= dptr->idunion.mca_dot.dim_L12_L12;
        dptr->idunion.mca_dot.f_L12_L12 -= dptr->idunion.mca_dot.dim_L12_L12;

        }

/* invert F */
        nn=1;
        if(dptr->idunion.mca_dot.sel_bonds)
         nn += n12;
        if(dptr->idunion.mca_dot.type==2)
         nn += m12;
        printf("write_first_dot: nn= %d\n",nn); fflush(stdout);

        F = matrix(1,nn,1,nn);
        S = matrix(1,nn,1,nn);

	S[1][1]=1.;
	F[1][1]=0.;

        if(dptr->idunion.mca_dot.sel_bonds)
        {

        for (j=1; j<=n12; j++) {
           i2=1;
           j2=1+j;
           S[i2][j2]=*dptr->idunion.mca_dot.s_1_l12;
           S[j2][i2]=S[i2][j2];
           F[i2][j2]=0.;
           F[j2][i2]=0.;
           dptr->idunion.mca_dot.s_1_l12++;
        }
        dptr->idunion.mca_dot.s_1_l12 -= dptr->idunion.mca_dot.dim_1_l12;

        for (i=1; i<=n12; i++) {
         for (j=i; j<=n12; j++) {
           i2=1+i;
           j2=1+j;
           S[i2][j2]=*dptr->idunion.mca_dot.s_l12_l12;
           S[j2][i2]=S[i2][j2];
           F[i2][j2]=*dptr->idunion.mca_dot.f_l12_l12;
           F[j2][i2]=F[i2][j2];
           dptr->idunion.mca_dot.s_l12_l12++;
           dptr->idunion.mca_dot.f_l12_l12++;
         }
        }
        dptr->idunion.mca_dot.s_l12_l12 -= dptr->idunion.mca_dot.dim_l12_l12;
        dptr->idunion.mca_dot.f_l12_l12 -= dptr->idunion.mca_dot.dim_l12_l12;

        }

        if(dptr->idunion.mca_dot.type>1) {

        if(dptr->idunion.mca_dot.sel_bonds) {

        for(i=1;i<=n12;i++) {
         for(j=1;j<=m12;j++) {
          i2=1+i;
          j2=1+n12+j;
          F[i2][j2]=*dptr->idunion.mca_dot.f_l12_L12++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*dptr->idunion.mca_dot.s_l12_L12++;
          S[j2][i2]=S[i2][j2];
         }
        }
        dptr->idunion.mca_dot.s_l12_L12 -= dptr->idunion.mca_dot.dim_l12_L12;
        dptr->idunion.mca_dot.f_l12_L12 -= dptr->idunion.mca_dot.dim_l12_L12;

        }

        for(j=1;j<=m12;j++) {
          i2=1;
          j2=1+n12+j;
          F[i2][j2]=0.;
          F[j2][i2]=0.;
          S[i2][j2]=*dptr->idunion.mca_dot.s_1_L12++;
          S[j2][i2]=S[i2][j2];
        }
        dptr->idunion.mca_dot.s_1_L12 -= dptr->idunion.mca_dot.dim_1_L12;

        for(i=1;i<=m12;i++) {
         for(j=i;j<=m12;j++) {
          i2=1+n12+i;
          j2=1+n12+j;
          F[i2][j2]=*dptr->idunion.mca_dot.f_L12_L12++;
          F[j2][i2]=F[i2][j2];
          S[i2][j2]=*dptr->idunion.mca_dot.s_L12_L12++;
          S[j2][i2]=S[i2][j2];
         }
        }
        dptr->idunion.mca_dot.s_L12_L12 -= dptr->idunion.mca_dot.dim_L12_L12;
        dptr->idunion.mca_dot.f_L12_L12 -= dptr->idunion.mca_dot.dim_L12_L12;

        }

/*
*/
        fprintf(file,"S,F:\n");
        for(i=1;i<=nn;i++)
         for(j=1;j<=nn;j++)
          fprintf(file,"%d %d %e %e\n",i,j,S[i][j],F[i][j]);
        fflush(file);

        for(i=1;i<=nn;i++) {
         if(S[i][i]<ZERO) {
          printf("write_first:\n");
          printf("S element is wrong: S(%d,%d) = %e\n",i,i,S[i][i]);
         }
        }

        fflush(stdout);
        Fm1 = matrix(1,nn,1,nn);
        tmp = matrix(1,nn,1,nn);

#ifdef IMSL
        tol=TOLSVD;
        ui = (_REAL *)get2(nn*nn,sizeof(_REAL));
        for(i=1,k=0;i<=nn;i++)
         for(j=1;j<=nn;j++,k++)
          ui[k]=F[i][j];
        w = imsl_d_lin_svd_gen(nn,nn,ui,
               IMSL_RANK,tol,&rank,
               IMSL_INVERSE,&v,0);
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
          Fm1[i][j]=v[k];
         }
        }
        free(v);
        free(w);
#else
        rank=invert(F,nn,nn,Fm1);
#endif
        fprintf(file,"write_first_p1: rank of F= %d\n",rank); fflush(file);
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          csum=0.;
          for(k=1;k<=nn;k++) {
           csum+=Fm1[i][k]*S[k][j];
          }
          tmp[i][j]=csum;
         }
        }
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          csum=0.;
          for(k=1;k<=nn;k++) {
           csum+=S[i][k]*tmp[k][j];
          }
          Fm1[i][j]=csum;
         }
        }
/*
*/
        for(i=1;i<=nn;i++) {
         tol=1./S[i][i];
         csum=tol*Fm1[i][i];
         fprintf(file,"tau(%d)= %e\n",i,factns*csum);
        }
        fflush(file);

#ifdef IMSL
        tol=TOLSVD;
        ui = (_REAL *)get2(nn*nn,sizeof(_REAL));
        if(dptr->idunion.mca_dot.rigid==0) {
         for(i=1,k=0;i<=nn;i++)
          for(j=1;j<=nn;j++,k++)
           ui[k]=S[i][j];
        } else {
         for(i=1,k=0;i<=nn;i++)
          for(j=1;j<=nn;j++,k++)
           ui[k]=F[i][j];
        }
        w = imsl_d_lin_svd_gen(nn,nn,ui,
               IMSL_RANK,tol,&rank,
               IMSL_INVERSE,&v,0);
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
          Fm1[i][j]=v[k];
         }
        }
        free(v);
        free(w);
#else
        if(dptr->idunion.mca_dot.rigid==0) {
         rank=invert(S,nn,nn,Fm1);
        } else {
         rank=invert(F,nn,nn,Fm1);
        }
#endif
        fprintf(file,"write_first_p1: rank of S= %d\n",rank); fflush(file);

/*
        fprintf(file,"Fm1*F:\n");
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          csum=0.;
          for(k=1;k<=nn;k++) {
           csum+=Fm1[i][k]*F[k][j];
          }
          fprintf(file,"%d %d %e\n",i,j,csum);
         }
        }
        fflush(file);
*/

        if(dptr->idunion.mca_dot.rigid==0) {

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          csum=0.;
          for(k=1;k<=nn;k++) {
           csum+=Fm1[i][k]*F[k][j];
          }
          tmp[i][j]=csum;
         }
        }

        } else {

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nn;j++) {
          csum=0.;
          for(k=1;k<=nn;k++) {
           csum+=Fm1[i][k]*S[k][j];
          }
          tmp[i][j]=csum;
         }
        }
        }

        w=vector(1,nn);
#ifdef IMSL
        for(i=1,l=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,l++) {
          ui[l]=tmp[i][j];
         }
        }
        wi = imsl_d_eig_gen(nn,ui,IMSL_VECTORS,&vi,0);
        if(dptr->idunion.mca_dot.rigid==0) {
         for(i=0;i<nn;i++) {
          w[i+1]=1./wi[i].re;
         }
        } else {
         for(i=0;i<nn;i++) {
          w[i+1]=wi[i].re;
         }
        }
        free(ui);
        free(wi);
#else
        vi=(_REAL *)get2(nn*nn,sizeof(_REAL));
#ifdef EISPACK
        eigvec_ev(nn,tmp,vi,w);
        if(dptr->idunion.mca_dot.rigid==0) {
         for(i=0;i<nn;i++) {
          w[i+1]=1./w[i+1];
         }
        }
#else
        eigvec_nr(tmp,nn,w,vi);
#endif
#endif

        fprintf(file,"Fm1*S eigenvalues:\n");
        nf=0;
        for(i=1;i<=nn;i++) {
         fprintf(file,"w(%d)= %e\n",i,w[i]);
         if(Nabs(w[i])>=ZERO) {
          nf++;
         }
        }
        fflush(file);

        if(nf!=dptr->idunion.mca_dot.n1) {
         printf("write_first_p1: WARNING: nf=%d is different from n1=%d.\n",
                 nf,n1);
        }
        fflush(file);

        dptr->idunion.mca_dot.nf = nf;
        nf=nn;
        dptr->idunion.mca_dot.eval= (_REAL *)get2(nf,sizeof(_REAL));
        dptr->idunion.mca_dot.evec= matrix(1,nn,1,nf);

        nf=0;
        for(i=1;i<=nn;i++) {
         dptr->idunion.mca_dot.eval[i-1]=1./w[i];
        }
        for(i=1,k=0;i<=nn;i++) {
         for(j=1;j<=nn;j++,k++) {
#ifdef IMSL
          dptr->idunion.mca_dot.evec[i][j]=vi[k].re;
#else
          dptr->idunion.mca_dot.evec[i][j]=vi[k];
#endif
         }
        }
        eigsrt(w,dptr->idunion.mca_dot.evec,nn);
        for(i=1;i<=nn;i++) {
         dptr->idunion.mca_dot.eval[i-1]=1./w[i];
        }
/*
        free_vector(w,1,nn);
        free(vi);
*/
/*
        for(i=1,k=0;i<=nn;i++,k++) {
         if(Nabs(wi[i-1].re)>=ZERO)
         {
          nf++;
          dptr->idunion.mca_dot.eval[nf-1]=1./wi[i-1].re;
          vt=&vi[k];
          for(j=1;j<=nn;j++) {
           dptr->idunion.mca_dot.evec[j][nf]=(*vt).re;
           vt+=nn;
          }
         }
        }
*/

/*
        fprintf(file,"evec:\n");
        for(j=1;j<=3;j++) {
         for(i=1;i<=nn;i++) {
          fprintf(file,"%d %d %e\n",i,j,dptr->idunion.mca_dot.evec[i][j]);
         }
        }
*/

        nf=nn;
        for(i=1;i<=nn;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mca_dot.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nn;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mca_dot.evec[k][i]*tmp[k][j];
          }
          Fm1[i][j]=sum;
         }
        }
/*
        fprintf(file,"Ct*S*C:\n");
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nf;j++) {
          fprintf(file,"%d %d %e\n",i,j,Fm1[i][j]);
         }
        }
        fflush(file);
*/

/*
        for(i=1,k=0;i<=nf;i++) {
         for(j=1;j<=nf;j++,k++) {
          ui[k]=Fm1[i][j];
         }
        }
        w = imsl_d_eig_sym(nf,ui,IMSL_VECTORS,&v,0);
        for(i=1,k=0;i<=nf;i++) {
         if(Nabs(w[i-1]) > ZERO)
         {
          w[i-1]=1./l_sqrt(w[i-1]);
          for(j=1;j<=nf;j++,k++) {
           tmp[i][j]=v[k];
          }
         } else {
          w[i-1]=1.;
          for(j=1;j<=nf;j++,k++) {
           tmp[i][j]=0.;
          }
          tmp[i][i]=1.;
         }
        }
        for(i=1,m=0;i<=nf;i++) {
         for(j=1;j<=nf;j++,m++) {
          for(k=1,sum=0.;k<=nf;k++) {
           for(l=1,sum2=0.;l<=nf;l++) {
            sum2+=tmp[k][l]*w[l-1]*tmp[j][l];
           }
           sum+=dptr->idunion.mca_dot.evec[i][k]*sum2;
          }
          v[m]=sum;
         }
        }
        for(i=1,k=0;i<=nf;i++) {
         for(j=1;j<=nf;j++,k++) {
          dptr->idunion.mca_dot.evec[i][j]=v[k];
         }
        }
*/

/* normalize eigenvectors */
        for(j=1,k=0;j<=nn;j++) {
         if(Nabs(Fm1[j][j]) > ZERO) {
          k++;
          sum=Fm1[j][j];
         } else {
          sum=1.;
         }
         sum=l_sqrt(sum);
         for(i=1;i<=nn;i++) {
          dptr->idunion.mca_dot.evec[i][j]/=sum;
         }
        }
        if(k!=nf) {
         printf("write_first: missed zero in Ct*S*C.\n");
        }
/*
        sum2=0.;
        for(j=1,k=0;j<=nn;j++) {
         if(Nabs(Fm1[j][j]) > ZERO) {
          k++;
          sum=Fm1[j][j];
         } else {
          sum=1.;
         }
         sum2+=sum;
        }
        sum2=l_sqrt(sum2);
        if(k!=nf) {
         printf("write_first: missed zero in Ct*S*C.\n");
        }
        for(j=1,k=0;j<=nn;j++) {
         for(i=1;i<=nn;i++) {
          dptr->idunion.mca_dot.evec[i][j]/=sum2;
         }
        }
*/

        for(i=1;i<=nn;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=S[i][k]*dptr->idunion.mca_dot.evec[k][j];
          }
          tmp[i][j]=sum;
         }
        }
        sum2=0.;
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nf;j++) {
          for(k=1,sum=0.;k<=nn;k++) {
           sum+=dptr->idunion.mca_dot.evec[k][i]*tmp[k][j];
          }
          S[i][j]=sum;
         }
         sum2+=S[i][i];
        }
        printf("sum_i(psi_i*psi_i)= %e\n",sum2);
/*
        fprintf(file,"Ct*S1*C:\n");
        for(i=1;i<=nf;i++) {
         for(j=1;j<=nf;j++) {
          fprintf(file,"%d %d %e\n",i,j,S[i][j]);
         }
        }
        fflush(file);
*/

/*
        fprintf(file,"evec2:\n");
        for(j=1;j<=nf;j++) {
         for(i=1;i<=nn;i++) {
          fprintf(file,"%d %d %e\n",i,j,dptr->idunion.mca_dot.evec[i][j]);
         }
        }
*/

        free_matrix(F,1,nn,1,nn);
        free_matrix(S,1,nn,1,nn);
        free_matrix(tmp,1,nn,1,nn);

/*
        if(dptr->idunion.mca_dot.p1_to_p2) {
         dptr->idunion.mca_dot.npsi2=nn*(nn+1)/2;
         dptr->idunion.mca_dot.psi2=
          (rot2struct *)get2(nn*(nn+1)/2,sizeof(rot2struct));
        } 
*/
}
