#include "head.h"

// extern char		tok[];
// char		tok[];
// void *get2(int,int);

/**** extern ****/
extern idstruct		*startid(),*getstruct();
extern _REAL **matrix(),*vector();

/************************************************************************/

void
init_ccff()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(CCFF, id);

	dptr->idunion.ccff.type = 1;
	dptr->idunion.ccff.cases = 0;
	dptr->idunion.ccff.dump = 0;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("CCFF: expected file name","");
	strcpy(dptr->idunion.ccff.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.ccff.fname, 0, &exists, "CCFF"))
		exit(1);
	if ((dptr->idunion.ccff.file = 
				  fopen(dptr->idunion.ccff.fname, "w")) == NULL){
		perror(dptr->idunion.ccff.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************/

void
setup_ccff()
{
	idstruct	*dptr,*serptr,*bptr,*modptr;
        extern parmstruct *prm;
	int i,j,k,icount,ifound,mm;
        int resnum1,resnum2,resnum3,resnum4,resnumt;
        char atnam1[2*ATOMTOKLEN], atnam2[2*ATOMTOKLEN],
             atnam3[2*ATOMTOKLEN], atnam4[2*ATOMTOKLEN],*atnamt;
        _REAL **r1t, **r2t, **r3t, **r4t;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != CCFF)
                        inerr("setupccff: expected id or keyword\n","");
        } else {
                printf("setupccff: id %s not found in idlist\n", tok);
        }

        gettok();
        if ((bptr = getstruct(tok)) != NULL) {
                if (bptr->id_type != BEAD)
                        inerr("setupccff: expected id or keyword\n","");
        } else {
                printf("setupccff: id %s not found in idlist\n", tok);
        }

        gettok();
        if ((serptr = getstruct(tok)) == NULL) 
                inerr("setupccff: no such id\n","");
        if (serptr->id_type != STREAM)
                inerr("setupccff: stream id expected\n","");

        gettok();
        if ((modptr = getstruct(tok)) == NULL) 
                inerr("setupccff: no such id\n","");
        if (   modptr->id_type != MCARIG 
            && modptr->id_type != RIG
            && modptr->id_type != MCA_P2P1
            && modptr->id_type != MCA_P2P1MM
            && modptr->id_type != MCA_P1
            && modptr->id_type != MCA_P2P1L
            && modptr->id_type != MCA_P1L
            && modptr->id_type != MCA_RL )
                inerr("setupccff: some model id expected\n","");

        dptr->idunion.ccff.model = modptr;
        dptr->idunion.ccff.prm = serptr->idunion.ser.prm;

        r1t = (_REAL **) get2(dptr->idunion.ccff.prm->Nsit,sizeof(long));
        r2t = (_REAL **) get2(dptr->idunion.ccff.prm->Nsit,sizeof(long));
        r3t = (_REAL **) get2(dptr->idunion.ccff.prm->Nsit,sizeof(long));
        r4t = (_REAL **) get2(dptr->idunion.ccff.prm->Nsit,sizeof(long));

        icount=0;
        ifound=0;
        for(;;) {
         gettok();
         if (tok[0] == ';') {
          break;

         } else if(!strcmp(tok,"READFSI")) {
          dptr->idunion.ccff.dump=-1;
         } else {
          icount++;

          if(icount>1) {
           r1t++;
           r2t++;
           r3t++;
           r4t++;
          }
          sprintf(atnam1,"%s    ",tok);
          gettok();
          tokint(&resnum1,tok);
          /* printf("atnam1: %s %d\n",atnam1,resnum1); */
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           atnamt = &serptr->idunion.ser.prm->AtomNames[ATOMTOKLEN*i];
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           if (!strncmp(atnam1,atnamt,ATOMTOKLEN)&&resnum1==resnumt) {
            *r1t = serptr->idunion.ser.set->crd;
            /* printf("atom found:\n"); */
            ifound++;
            /* printf("atom1: %s %d\n",atnam1,resnum1); */
            /* printf("atomt: %s %d\n",atnamt,resnumt); */
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;

          gettok();
          sprintf(atnam2,"%s    ",tok);
          gettok();
          tokint(&resnum2,tok);
          /* printf("atnam2: %s %d\n",atnam2,resnum2); */
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           atnamt = &serptr->idunion.ser.prm->AtomNames[ATOMTOKLEN*i];
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           if (!strncmp(atnam2,atnamt,ATOMTOKLEN)&&resnum2==resnumt) {
            *r2t = serptr->idunion.ser.set->crd;
            /* printf("atom found:\n"); */
            ifound++;
            /* printf("atom2: %s %d\n",atnam2,resnum2); */
            /* printf("atomt: %s %d\n",atnamt,resnumt); */
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;

          gettok();
          sprintf(atnam3,"%s    ",tok);
          gettok();
          tokint(&resnum3,tok);
          /* printf("atnam3: %s %d\n",atnam3,resnum3); */
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           atnamt = &serptr->idunion.ser.prm->AtomNames[ATOMTOKLEN*i];
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           if (!strncmp(atnam3,atnamt,ATOMTOKLEN)&&resnum3==resnumt) {
            *r3t = serptr->idunion.ser.set->crd;
            /* printf("atom found:\n"); */
            ifound++;
            /* printf("atom3: %s %d\n",atnam3,resnum3); */
            /* printf("atomt: %s %d\n",atnamt,resnumt); */
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;

          gettok();
          sprintf(atnam4,"%s    ",tok);
          gettok();
          tokint(&resnum4,tok);
          /* printf("atnam4: %s %d\n",atnam4,resnum4); */
          for(i=0;i<serptr->idunion.ser.prm->Nsit;i++) {
           atnamt = &serptr->idunion.ser.prm->AtomNames[ATOMTOKLEN*i];
           resnumt=serptr->idunion.ser.prm->ResNum[i];
           if (!strncmp(atnam4,atnamt,ATOMTOKLEN)&&resnum4==resnumt) {
            *r4t = serptr->idunion.ser.set->crd;
            /* printf("atom found:\n"); */
            ifound++;
            /* printf("atom4: %s %d\n",atnam4,resnum4); */
            /* printf("atomt: %s %d\n",atnamt,resnumt); */
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*serptr->idunion.ser.prm->Nsit;

         }

        }
        r1t-=icount-1;
        r2t-=icount-1;
        r3t-=icount-1;
        r4t-=icount-1;

        if(icount!=ifound/4) {
         printf("setup_ccff: atom not found.\n");
         exit(1);
        }
        dptr->idunion.ccff.nbond = icount;

        dptr->idunion.ccff.r1 = (_REAL **)
            get2(icount,sizeof(long));
        dptr->idunion.ccff.r2 = (_REAL **)
            get2(icount,sizeof(long));
        dptr->idunion.ccff.r3 = (_REAL **)
            get2(icount,sizeof(long));
        dptr->idunion.ccff.r4 = (_REAL **)
            get2(icount,sizeof(long));
        for(i=1;i<=icount;i++) {
         *dptr->idunion.ccff.r1++ = *r1t++;
         *dptr->idunion.ccff.r2++ = *r2t++;
         *dptr->idunion.ccff.r3++ = *r3t++;
         *dptr->idunion.ccff.r4++ = *r4t++;
        }
        r1t-=icount;
        r2t-=icount;
        r3t-=icount;
        r4t-=icount;
        free(r1t);
        free(r2t);
        free(r3t);
        free(r4t);
        dptr->idunion.ccff.r1-=icount;
        dptr->idunion.ccff.r2-=icount;
        dptr->idunion.ccff.r3-=icount;
        dptr->idunion.ccff.r4-=icount;

   if(modptr->id_type==MCARIG) {

    dptr->idunion.ccff.nbase = modptr->idunion.mcarig.bead->idunion.bead.nbond;
    dptr->idunion.ccff.rtmp =
     (_REAL *)get2(3*dptr->idunion.ccff.nbase,sizeof(_REAL));
    dptr->idunion.ccff.b1tmp =
     (_REAL *)get2(3*dptr->idunion.ccff.nbond,sizeof(_REAL));
    dptr->idunion.ccff.b2tmp =
     (_REAL *)get2(3*dptr->idunion.ccff.nbond,sizeof(_REAL));
    mm=5;

   } else if (modptr->id_type==MCA_P2P1 || 
              modptr->id_type==MCA_P2P1MM || 
              modptr->id_type==MCA_P2P1L) {

    dptr->idunion.ccff.nbase = 5;
    dptr->idunion.ccff.b1tmp =
     (_REAL *)get2(5*dptr->idunion.ccff.nbond,sizeof(_REAL));
    dptr->idunion.ccff.b2tmp =
     (_REAL *)get2(5*dptr->idunion.ccff.nbond,sizeof(_REAL));
    mm=5;

   } else if (modptr->id_type==MCA_P1 ||
              modptr->id_type==MCA_P1L) {

    dptr->idunion.ccff.nbase = 3;
    dptr->idunion.ccff.b1tmp =
     (_REAL *)get2(3*dptr->idunion.ccff.nbond,sizeof(_REAL));
    dptr->idunion.ccff.b2tmp =
     (_REAL *)get2(3*dptr->idunion.ccff.nbond,sizeof(_REAL));
    mm=3;

   } else if (modptr->id_type==MCA_RL) {

    dptr->idunion.ccff.nbase = 5;
    dptr->idunion.ccff.b1tmp =
     (_REAL *)get2(5*dptr->idunion.ccff.nbond,sizeof(_REAL));
    dptr->idunion.ccff.b2tmp =
     (_REAL *)get2(5*dptr->idunion.ccff.nbond,sizeof(_REAL));
    mm=5;

   } else if (modptr->id_type==RIG) {

    dptr->idunion.ccff.nbase = 5;
    dptr->idunion.ccff.b1tmp =
     (_REAL *)get2(5*dptr->idunion.ccff.nbond,sizeof(_REAL));
    dptr->idunion.ccff.b2tmp =
     (_REAL *)get2(5*dptr->idunion.ccff.nbond,sizeof(_REAL));
    mm=1;

   }

   dptr->idunion.ccff.tmp =
     vector(1,dptr->idunion.ccff.nbond);
   dptr->idunion.ccff.modsq =
     vector(1,dptr->idunion.ccff.nbond);
   for(i=1;i<=dptr->idunion.ccff.nbond;i++) {
    dptr->idunion.ccff.modsq[i]=0.;
   }
   dptr->idunion.ccff.proj1 =
     matrix(1,dptr->idunion.ccff.nbond,1,dptr->idunion.ccff.nbase);
   for(i=1;i<=dptr->idunion.ccff.nbond;i++)
    for(j=1;j<=dptr->idunion.ccff.nbase;j++)
     dptr->idunion.ccff.proj1[i][j]=0.;
   dptr->idunion.ccff.proj2 =
     matrix(1,dptr->idunion.ccff.nbond,1,dptr->idunion.ccff.nbase);
   for(i=1;i<=dptr->idunion.ccff.nbond;i++)
    for(j=1;j<=dptr->idunion.ccff.nbase;j++)
     dptr->idunion.ccff.proj2[i][j]=0.;

   if(dptr->idunion.ccff.type>0) {
 
    dptr->idunion.ccff.domega = 2.*_PI*0.005;
    dptr->idunion.ccff.nt = 1000;
    dptr->idunion.ccff.dtime = _PI/
        ((_REAL)dptr->idunion.ccff.nt*dptr->idunion.ccff.domega);
    dptr->idunion.ccff.dtime /= bptr->idunion.bead.factns;
     /* dptr->idunion.ccff.p2t = vector(1,dptr->idunion.ccff.nt); */
    dptr->idunion.ccff.domega *= 2.;
   }
}

/************************************************************************/
void
adjust_ccff(idstruct *dptr)
{
   idstruct	*modptr;
   int i,j;

/*
   free_matrix(dptr->idunion.ccff.proj1,
     1,dptr->idunion.ccff.nbond,1,dptr->idunion.ccff.nbase);
   free_matrix(dptr->idunion.ccff.proj2,
     1,dptr->idunion.ccff.nbond,1,dptr->idunion.ccff.nbase);
*/
   modptr=dptr->idunion.ccff.model;

   if(modptr->id_type==MCA_P2P1) {
    if(modptr->idunion.mca_p2p1.type==1) {
     dptr->idunion.ccff.nbase = modptr->idunion.mca_p2p1.n12;
    } else if(modptr->idunion.mca_p2p1.type==2) {
     dptr->idunion.ccff.nbase = modptr->idunion.mca_p2p1.n12
                               +modptr->idunion.mca_p2p1.n1234;
    } else if(modptr->idunion.mca_p2p1.type==3) {
     dptr->idunion.ccff.nbase = modptr->idunion.mca_p2p1.n12
                               +modptr->idunion.mca_p2p1.n1234
                               +modptr->idunion.mca_p2p1.n123456;
    }
   } else if(modptr->id_type==MCA_P2P1MM) {
    dptr->idunion.ccff.nbase = modptr->idunion.mca_p2p1mm.n12
                              +modptr->idunion.mca_p2p1mm.n1234;
   } else if(modptr->id_type==MCA_P1) {
    if(modptr->idunion.mca_p1.type==1) {
     dptr->idunion.ccff.nbase = modptr->idunion.mca_p1.n1;
    } else if(modptr->idunion.mca_p1.type==2) {
     dptr->idunion.ccff.nbase = modptr->idunion.mca_p1.n1
                               +modptr->idunion.mca_p1.n123;
    }
   } else if(modptr->id_type==MCA_P1L) {
    dptr->idunion.ccff.nbase = modptr->idunion.mca_p1l.n0m
                              +modptr->idunion.mca_p1l.n0l;
   } else if(modptr->id_type==MCA_P2P1L) {
    if(modptr->idunion.mca_p2p1l.type==1) {
     dptr->idunion.ccff.nbase = modptr->idunion.mca_p2p1l.n12;
    } else if(modptr->idunion.mca_p2p1l.type==2) {
    dptr->idunion.ccff.nbase = modptr->idunion.mca_p2p1l.n12
                              +modptr->idunion.mca_p2p1l.n1234;
    }
   } else if(modptr->id_type==MCA_RL) {
    dptr->idunion.ccff.nbase = modptr->idunion.mca_rl.mcar->idunion.mca_p2p1.n12
                              +modptr->idunion.mca_rl.mcal->idunion.mca_p2p1.n12
                              +modptr->idunion.mca_rl.mcar->idunion.mca_p2p1.n1234
                              +modptr->idunion.mca_rl.mcal->idunion.mca_p2p1.n1234;
   }
   dptr->idunion.ccff.proj1 =
     matrix(1,dptr->idunion.ccff.nbond,1,dptr->idunion.ccff.nbase);
   for(i=1;i<=dptr->idunion.ccff.nbond;i++)
    for(j=1;j<=dptr->idunion.ccff.nbase;j++)
      dptr->idunion.ccff.proj1[i][j]=0.;
   dptr->idunion.ccff.proj2 =
     matrix(1,dptr->idunion.ccff.nbond,1,dptr->idunion.ccff.nbase);
   for(i=1;i<=dptr->idunion.ccff.nbond;i++)
    for(j=1;j<=dptr->idunion.ccff.nbase;j++)
      dptr->idunion.ccff.proj2[i][j]=0.;

   for(i=1;i<=dptr->idunion.ccff.nbond;i++)
    dptr->idunion.ccff.modsq[i]=0.;

   printf("new nbase= %d\n",dptr->idunion.ccff.nbase);

}

/************************************************************************/
void
do_ccff1(idstruct *dptr, idstruct *serptr)
{
   idstruct *modptr=dptr->idunion.ccff.model;
   _REAL *box=serptr->idunion.ser.set->boxc;
   _REAL *rt=dptr->idunion.ccff.rtmp;
   _REAL *b1t=dptr->idunion.ccff.b1tmp;
   _REAL *b2t=dptr->idunion.ccff.b2tmp;
   int i,j,j2,k,mm;
   _REAL **r1, **r2;
   _REAL *rt1,*rt2;
   _REAL x1,y1,z1,x2,y2,z2,dx,dy,dz,scp,lsq;
   void genquatch1();
   _REAL rtmp[3];

   mm=3;
   ++dptr->idunion.ccff.cases;

   r1=dptr->idunion.ccff.r1;
   r2=dptr->idunion.ccff.r2;
   for (i=0;i<dptr->idunion.ccff.nbond;i++,r1++,r2++) {
    rt1=*r1; rt2=*r2;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
/* add an extra atom with some construction... */
/*
    genquatch1(&rt2[-3],&rt2[0],&rt2[3],rtmp);
    x2=rtmp[0]; y2=rtmp[1]; z2=rtmp[2];
*/
    x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
    dx=x1-x2;
    dy=y1-y2;
    dz=z1-z2;
    scp=l_sqrt(dx*dx+dy*dy+dz*dz);
    *b1t++ = dx/scp;
    *b1t++ = dy/scp;
    *b1t++ = dz/scp;
   }
   b1t-=3*dptr->idunion.ccff.nbond;
   r1-=dptr->idunion.ccff.nbond;
   r2-=dptr->idunion.ccff.nbond;

   r1=dptr->idunion.ccff.r3;
   r2=dptr->idunion.ccff.r4;
   for (i=0;i<dptr->idunion.ccff.nbond;i++,r1++,r2++) {
    rt1=*r1; rt2=*r2;
    x1=rt1[0]; y1=rt1[1]; z1=rt1[2];
/* add an extra atom with some construction... */
/*
    genquatch1(&rt2[-3],&rt2[0],&rt2[3],rtmp);
    x2=rtmp[0]; y2=rtmp[1]; z2=rtmp[2];
*/
    x2=rt2[0]; y2=rt2[1]; z2=rt2[2];
    dx=x1-x2;
    dy=y1-y2;
    dz=z1-z2;
    scp=l_sqrt(dx*dx+dy*dy+dz*dz);
    *b2t++ = dx/scp;
    *b2t++ = dy/scp;
    *b2t++ = dz/scp;
   }
   b2t-=3*dptr->idunion.ccff.nbond;
   r1-=dptr->idunion.ccff.nbond;
   r2-=dptr->idunion.ccff.nbond;

   for (i=1;i<=dptr->idunion.ccff.nbond;i++,b1t+=3,b2t+=3) {
    dptr->idunion.ccff.modsq[i]+=b1t[0]*b2t[0]+b1t[1]*b2t[1]+b1t[2]*b2t[2];
    for (j=1;j<=modptr->idunion.mca_p1.n1;j++) {
     scp=0.;
     scp+= modptr->idunion.mca_p1.m->m[0]*b1t[0];
     scp+= modptr->idunion.mca_p1.m->m[1]*b1t[1];
     scp+= modptr->idunion.mca_p1.m->m[2]*b1t[2];
     dptr->idunion.ccff.proj1[i][j] += scp/3.;
     scp=0.;
     scp+= modptr->idunion.mca_p1.m->m[0]*b2t[0];
     scp+= modptr->idunion.mca_p1.m->m[1]*b2t[1];
     scp+= modptr->idunion.mca_p1.m->m[2]*b2t[2];
     dptr->idunion.ccff.proj2[i][j] += scp/3.;
     modptr->idunion.mca_p1.m++;
    }
    modptr->idunion.mca_p1.m-=modptr->idunion.mca_p1.n1;
    if(modptr->idunion.mca_p1.type>1) {
    for (j=1;j<=modptr->idunion.mca_p1.n123;j++) {
     j2=modptr->idunion.mca_p1.n1+j;
     scp=0.;
     scp+= modptr->idunion.mca_p1.l123->phi[0]*b1t[0];
     scp+= modptr->idunion.mca_p1.l123->phi[1]*b1t[1];
     scp+= modptr->idunion.mca_p1.l123->phi[2]*b1t[2];
     dptr->idunion.ccff.proj1[i][j2] += scp/3.;
     scp=0.;
     scp+= modptr->idunion.mca_p1.l123->phi[0]*b2t[0];
     scp+= modptr->idunion.mca_p1.l123->phi[1]*b2t[1];
     scp+= modptr->idunion.mca_p1.l123->phi[2]*b2t[2];
     dptr->idunion.ccff.proj2[i][j2] += scp/3.;
     modptr->idunion.mca_p1.l123++;
    }
    modptr->idunion.mca_p1.l123-=modptr->idunion.mca_p1.n123;
    }
   }
   b1t-=3*dptr->idunion.ccff.nbond;
   b2t-=3*dptr->idunion.ccff.nbond;
}

/************************************************************************/
void
write_ccff2(idstruct *dptr)
{
 extern parmstruct *prm;
 idstruct *modptr=dptr->idunion.ccff.model;
 _REAL *modsq=dptr->idunion.ccff.modsq;
 _REAL **proj1=dptr->idunion.ccff.proj1;
 _REAL **proj2=dptr->idunion.ccff.proj2;
 _REAL dtime=dptr->idunion.ccff.dtime;
 _REAL *eval;
 _REAL **evec;
 _REAL fact,sum,**cc,norm,time,*tau;
 _REAL factns;
 int i,j,k,l,nf,mm;
 FILE *file=dptr->idunion.ccff.file;

 fact=1./((_REAL)dptr->idunion.ccff.cases);

 mm=5;
 norm=1./((_REAL)mm);
 if(dptr->idunion.ccff.model->id_type==MCA_P2P1) {
  nf = modptr->idunion.mca_p2p1.nf;
  factns = modptr->idunion.mca_p2p1.bead->idunion.bead.factns;
  eval=modptr->idunion.mca_p2p1.eval;
  evec=modptr->idunion.mca_p2p1.evec;
  /* dptr->idunion.ccff.nbase=modptr->idunion.mca_p2p1.n12; */
 } else if(dptr->idunion.ccff.model->id_type==MCA_P2P1MM) {
  nf = modptr->idunion.mca_p2p1mm.nf;
  factns = modptr->idunion.mca_p2p1mm.bead->idunion.bead.factns;
  eval=modptr->idunion.mca_p2p1mm.eval;
  evec=modptr->idunion.mca_p2p1mm.evec;
 } else if(dptr->idunion.ccff.model->id_type==MCA_P2P1L) {
  nf = modptr->idunion.mca_p2p1l.nf;
  factns = modptr->idunion.mca_p2p1l.beadm->idunion.bead.factns;
  eval=modptr->idunion.mca_p2p1l.eval;
  evec=modptr->idunion.mca_p2p1l.evec;
 } else if(dptr->idunion.ccff.model->id_type==MCA_RL) {
  nf = modptr->idunion.mca_rl.nf;
  factns = modptr->idunion.mca_rl.mcar->idunion.mca_p2p1.bead->idunion.bead.factns;
  eval=modptr->idunion.mca_rl.eval;
  evec=modptr->idunion.mca_rl.evec;
 } else if(dptr->idunion.ccff.model->id_type==MCA_P1) {
  nf = modptr->idunion.mca_p1.nf;
  factns = modptr->idunion.mca_p1.bead->idunion.bead.factns;
  eval=modptr->idunion.mca_p1.eval;
  evec=modptr->idunion.mca_p1.evec;
  mm=3;
 } else if(dptr->idunion.ccff.model->id_type==MCA_P1L) {
  nf = modptr->idunion.mca_p1l.nf;
  factns = modptr->idunion.mca_p1l.beadm->idunion.bead.factns;
  eval=modptr->idunion.mca_p1l.eval;
  evec=modptr->idunion.mca_p1l.evec;
  mm=3;
 }
 norm=1./((_REAL)mm);

 for(i=1;i<=dptr->idunion.ccff.nbond;i++) {
  modsq[i] *= fact;
  for(j=1;j<=dptr->idunion.ccff.nbase;j++) {
   proj1[i][j]*=fact;
   proj2[i][j]*=fact;
  }
 }

 printf("write_ccff2: nbond= %d\n",dptr->idunion.ccff.nbond);
 printf("write_ccff2: nbase= %d\n",dptr->idunion.ccff.nbase);
 printf("write_ccff2: nf=    %d\n",nf);
 printf("write_ccff2: dtime= %f\n",dptr->idunion.ccff.dtime);
 fflush(stdout);

/*
 printf("projections:\n");
 for(i=1;i<=dptr->idunion.ccff.nbond;i++) {
  printf("modsq(%d)= %e\n",i,modsq[i]);
  for(j=1;j<=dptr->idunion.ccff.nbase;j++) {
   printf("proj1(%d,%d)= %e\n",i,j,proj1[i][j]);
  }
 }
 printf("eigenvectors:\n");
 for(j=1;j<=dptr->idunion.ccff.nbase;j++) {
  printf("eval(%d)= %e\n",j,eval[j-1]);
  for(i=1;i<=modptr->idunion.mca_p1.n1+modptr->idunion.mca_p1.n123;i++) {
   printf("evec(%d,%d)= %e\n",i,j,evec[i][j]);
  }
 }
*/
 cc = matrix(1,dptr->idunion.ccff.nbond,1,dptr->idunion.ccff.nbase);

 for(l=1;l<=dptr->idunion.ccff.nbond;l++) {
  sum=0.;
  for(i=1;i<=nf;i++) {
   cc[l][i] = 0.; 
   for(j=1;j<=dptr->idunion.ccff.nbase;j++) {
    for(k=1;k<=dptr->idunion.ccff.nbase;k++) {
     cc[l][i] += proj1[l][j]*proj2[l][k]*evec[j][i]*evec[k][i]/modsq[l]/norm;
    }
   }
   sum+=cc[l][i];
  }
  sum=1./sum;
  for(i=1;i<=nf;i++) {
   cc[l][i]*=sum;
  }
 }
 if(dptr->idunion.ccff.type==1) {
  tau=vector(1,dptr->idunion.ccff.nbond);
  for(i=1;i<=dptr->idunion.ccff.nbond;i++) {
   for(j=1,sum=0.;j<=nf;j++) {
    sum += cc[i][j]/eval[j-1];
   }
   tau[i]=sum*factns;
  } /* enddo i */
  for(i=1;i<=dptr->idunion.ccff.nbond;i++) {
   printf("tauccf(%d)= %e\n",i,tau[i]);
  }
  for(k=1,time=0.;k<=dptr->idunion.ccff.nt;k++,time+=dtime) {
   fprintf(file,"%e ",time*factns);
   for(i=1;i<=dptr->idunion.ccff.nbond;i++) {
    for(j=1,sum=0.;j<=nf;j++) {
     sum += cc[i][j]*exp(-eval[j-1]*time);
    }
    fprintf(file,"%e ",sum);
   } /* enddo i */
   fprintf(file,"\n");
  } /* enddo k */


 } /* endif type */

}
