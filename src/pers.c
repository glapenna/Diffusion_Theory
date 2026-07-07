#include "head.h"

// extern char		tok[];
// char		tok[];
// void *get2(int,int);

/**** extern ****/
extern idstruct		*startid(),*getstruct();

/************************************************************************/

void
initpers()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(PERS, id);

	dptr->idunion.pers.center = 1;
	dptr->idunion.pers.distmax = 10;
	dptr->idunion.pers.cases = 0;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("PERS: expected file name","");
	strcpy(dptr->idunion.pers.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.pers.fname, 0, &exists, "PERS"))
		exit(1);
	if ((dptr->idunion.pers.file = 
				  fopen(dptr->idunion.pers.fname, "w")) == NULL){
		perror(dptr->idunion.pers.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************/

void
setuppers()
{
	idstruct	*dptr,*serptr;
        extern parmstruct *prm;
	int i,j,k;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != PERS)
                        inerr("setuppers: expected id or keyword\n","");
        } else {
                printf("setuppers: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((serptr = getstruct(tok)) == NULL) 
                inerr("setuppers: no such id\n","");
        if (serptr->id_type != STREAM)
                inerr("setuppers: stream id expected\n","");

        for(;;) {
	        gettok();
        	if (tok[0] == ';') {
                 break;
                } else if (!strcmp(tok, "CENTER")) {
                        gettok();
                        if (!tokint(&dptr->idunion.pers.center))
                         inerr("PERS: expected integer center of distr.", "");
                } else if (!strcmp(tok, "DIST")) {
                        gettok();
                        if (!tokint(&dptr->idunion.pers.distmax))
                         inerr("PERS: expected integer distance in distr.", "");
                } else if (!strcmp(tok, "OFFSET")) {
                        gettok();
                        if (!tokint(&dptr->idunion.pers.offset))
                         inerr("PERS: expected integer offset in center.", "");
		} else inerr("setuppers: option not found", "");
	}

       dptr->idunion.pers.prm = prm;
/*
       dptr->idunion.pers.distr = (_REAL *)
        get2(dptr->idunion.pers.distmax,sizeof(_REAL));
       for(i=0;i<dptr->idunion.pers.distmax;i++)
        dptr->idunion.pers.distr[i]=0.;
*/
       dptr->idunion.pers.nn=prm->Nres;
       dptr->idunion.pers.distr = (_REAL *)
        get2(dptr->idunion.pers.nn,sizeof(_REAL));
       dptr->idunion.pers.count = (int *)
        get2(dptr->idunion.pers.nn,sizeof(int));
       for(i=0;i<dptr->idunion.pers.nn;i++)
        dptr->idunion.pers.distr[i]=0.;
        dptr->idunion.pers.count[i]=0;
}


/************************************************************************/

void
do_pers(idstruct *objptr,idstruct *serptr)
{
   parmstruct   *prm = objptr->idunion.pers.prm;
   _REAL   *r = serptr->idunion.ser.set->crd;
   _REAL   *box = serptr->idunion.ser.set->boxc;
   int center=objptr->idunion.pers.center;
   int i,im,k,atom1,atom2,res;
   _REAL *r1,*r2,dr[3],u0[3],u[3],mod,cos_theta;
   void images(_REAL *dist, idstruct *serptr);

   for(im=0;im<prm->Nmol;im++)
   {

   objptr->idunion.pers.cases++;

   atom1=center;
   atom2=center+1;
/*
   if( !((atom1%5==4) || (atom2%5==4)) )
   {
*/

   r1 = r+3*(atom1-1);
   r2 = r+3*(atom2-1);
   dr[0]=r1[0]-r2[0];
   dr[1]=r1[1]-r2[1];
   dr[2]=r1[2]-r2[2];
   images(dr,serptr);
   for(k=0,mod=0.;k<3;k++) mod+=dr[k]*dr[k];
   mod=1./l_sqrt(mod);
   for(k=0;k<3;k++) u0[k]=dr[k]*mod;

   for(i=0;i<objptr->idunion.pers.distmax;i++) 
   {
    atom1+=1;
    atom2+=1;
/*
    if( !((atom1%5==4) || (atom2%5==4)) )
    {
*/

    r1 = r+3*(atom1-1);
    r2 = r+3*(atom2-1);
    dr[0]=r1[0]-r2[0];
    dr[1]=r1[1]-r2[1];
    dr[2]=r1[2]-r2[2];
    images(dr,serptr);
    for(k=0,mod=0.;k<3;k++) mod+=dr[k]*dr[k];
    mod=1./l_sqrt(mod);
    for(k=0;k<3;k++) u[k]=dr[k]*mod;
    cos_theta = u[0]*u0[0]+u[1]*u0[1]+u[2]*u0[2];
    objptr->idunion.pers.distr[i]+=cos_theta;

/*
    }
*/
   }

/*
   }
*/

    r+=3*prm->Natom;
   }
   r-=3*(prm->Natom*prm->Nmol);

}

/************************************************************************/
do_pers_flory(idstruct *objptr,idstruct *serptr)
{
   parmstruct   *prm = objptr->idunion.pers.prm;
   _REAL   *r = serptr->idunion.ser.set->crd;
   _REAL   *box = serptr->idunion.ser.set->boxc;
   int center=objptr->idunion.pers.center;
   int i,im,k,atom1,atom2,nat_per_res;
   _REAL *r1,*r2,dr[3],u0[3],h[3],mod,cos_theta;
   void images(_REAL *dist, idstruct *serptr);

   for(im=0;im<prm->Nmol;im++)
   {

   nat_per_res = prm->Natom/prm->Nres;
   objptr->idunion.pers.cases++;

   r2 = r+3*(prm->Natom-1);
   r1 = r;
   h[0]=r2[0]-r1[0];
   h[1]=r2[1]-r1[1];
   h[2]=r2[2]-r1[2];
   /* images(dr,serptr); */

   atom1=center;
   atom2=center+nat_per_res;

   for(i=0;i<prm->Nres-1;i++)
   {

   r1 = r+3*(atom1-1);
   r2 = r+3*(atom2-1);
   dr[0]=r2[0]-r1[0];
   dr[1]=r2[1]-r1[1];
   dr[2]=r2[2]-r1[2];
   /* images(dr,serptr); */
   for(k=0,mod=0.;k<3;k++) mod+=dr[k]*dr[k];
   mod=1./l_sqrt(mod);
   for(k=0;k<3;k++) u0[k]=dr[k]*mod;

   cos_theta = u0[0]*h[0]+u0[1]*h[1]+u0[2]*h[2];
   objptr->idunion.pers.distr[i]+=cos_theta;

   r+=3*nat_per_res;

   }

   }
   r-=3*(prm->Nres-1)*(nat_per_res);

}

/************************************************************************/

void
do_pers_vb(idstruct *objptr,idstruct *serptr)
{
   parmstruct   *prm = objptr->idunion.pers.prm;
   _REAL   *r = serptr->idunion.ser.set->crd;
   _REAL   *box = serptr->idunion.ser.set->boxc;
   int center=objptr->idunion.pers.center;
   int i,im,k,atom1,atom2,res;
   _REAL *r1,*r2,dr[3],u0[3],u[3],mod,cos_theta;
   void images(_REAL *dist, idstruct *serptr);

   for(im=0;im<prm->Nmol;im++)
   {

   objptr->idunion.pers.cases++;

   atom1=prm->Ipres[center-1];
   atom2=prm->Ipres[center];
   atom1+=objptr->idunion.pers.offset;
   atom2+=objptr->idunion.pers.offset;
   r1 = r+3*(atom1-1);
   r2 = r+3*(atom2-1);
   dr[0]=r1[0]-r2[0];
   dr[1]=r1[1]-r2[1];
   dr[2]=r1[2]-r2[2];
   images(dr,serptr);
   for(k=0,mod=0.;k<3;k++) mod+=dr[k]*dr[k];
   mod=1./l_sqrt(mod);
   for(k=0;k<3;k++) u0[k]=dr[k]*mod;

   res=center+1;
   for(i=0;i<objptr->idunion.pers.distmax;i++) 
   {
    atom1=prm->Ipres[res-1];
    atom2=prm->Ipres[res];
    atom1+=objptr->idunion.pers.offset;
    atom2+=objptr->idunion.pers.offset;
    r1 = r+3*(atom1-1);
    r2 = r+3*(atom2-1);
    dr[0]=r1[0]-r2[0];
    dr[1]=r1[1]-r2[1];
    dr[2]=r1[2]-r2[2];
    images(dr,serptr);
    for(k=0,mod=0.;k<3;k++) mod+=dr[k]*dr[k];
    mod=1./l_sqrt(mod);
    for(k=0;k<3;k++) u[k]=dr[k]*mod;
    cos_theta = u[0]*u0[0]+u[1]*u0[1]+u[2]*u0[2];
    objptr->idunion.pers.distr[i]+=cos_theta;
    res++;
   }

    r+=3*prm->Natom;
   }
   r-=3*(prm->Natom*prm->Nmol);

}
/************************************************************************/

void
do_pers_r(idstruct *objptr,idstruct *serptr)
{
   parmstruct   *prm = objptr->idunion.pers.prm;
   _REAL   *r = serptr->idunion.ser.set->crd;
   _REAL   *box = serptr->idunion.ser.set->boxc;
   int center=objptr->idunion.pers.center;
   int i,im,ir,nn,k,atom1,atom2,atom1_2,atom2_2,kk;
   _REAL *r1,*r2,dr[3],u0[3],u[3],mod,cos_theta;
   void images(_REAL *dist, idstruct *serptr);

   for(im=0;im<prm->Nmol;im++)
   {

   objptr->idunion.pers.cases++;
   atom1=prm->Ipres[center-1];
   atom2=prm->Ipres[center];
   atom1+=objptr->idunion.pers.offset;
   atom2+=objptr->idunion.pers.offset;
   nn=prm->Nres-2*objptr->idunion.pers.center;

   for(ir=0;ir<nn;ir++)
   {

   r1 = r+3*(atom1-1);
   r2 = r+3*(atom2-1);
   dr[0]=r1[0]-r2[0];
   dr[1]=r1[1]-r2[1];
   dr[2]=r1[2]-r2[2];
   images(dr,serptr);
   for(k=0,mod=0.;k<3;k++) mod+=dr[k]*dr[k];
   mod=1./l_sqrt(mod);
   for(k=0;k<3;k++) u0[k]=dr[k]*mod;

   for(i=ir,kk=0;i<nn;i++,kk++) 
   {
    atom1_2=prm->Ipres[center+ir+kk-1];
    atom2_2=prm->Ipres[center+ir+kk];
    atom1_2+=objptr->idunion.pers.offset;
    atom2_2+=objptr->idunion.pers.offset;
    r1 = r+3*(atom1_2-1);
    r2 = r+3*(atom2_2-1);
    dr[0]=r1[0]-r2[0];
    dr[1]=r1[1]-r2[1];
    dr[2]=r1[2]-r2[2];
    images(dr,serptr);
    for(k=0,mod=0.;k<3;k++) mod+=dr[k]*dr[k];
    mod=1./l_sqrt(mod);
    /*printf ("mod= %e\n",mod);*/
    for(k=0;k<3;k++) u[k]=dr[k]*mod;
    cos_theta = u[0]*u0[0]+u[1]*u0[1]+u[2]*u0[2];
    objptr->idunion.pers.distr[kk]+=cos_theta;
    objptr->idunion.pers.count[kk]++;
   }

   atom1=prm->Ipres[center+ir];
   atom2=prm->Ipres[center+ir+1];
   atom1+=objptr->idunion.pers.offset;
   atom2+=objptr->idunion.pers.offset;

   }

    r+=3*prm->Natom;
   }
   r-=3*(prm->Natom*prm->Nmol);

}


/************************************************************************/

void
writepers(dptr)
idstruct	*dptr;
{
	int	i,j;
	int	states = dptr->idunion.pers.cases;
	FILE 	*file = dptr->idunion.pers.file;
	_REAL	fac,distrib;

	fac = 1./ ((_REAL)states);

        fprintf(file, "%d %e\n", 0,1.);
        for (i=0; i<dptr->idunion.pers.distmax; i++) {
         distrib = dptr->idunion.pers.distr[i]*fac;
         fprintf(file, "%d %e\n", i+1,distrib);
        }
        fclose(file);
}

/************************************************************************/

void
writepers_r(dptr)
idstruct	*dptr;
{
	int	i,j;
	int	states = dptr->idunion.pers.cases;
	FILE 	*file = dptr->idunion.pers.file;
	_REAL	fac,distrib;

        for (i=0; i<dptr->idunion.pers.nn; i++) {
         if(dptr->idunion.pers.count[i]!=0)
         {
          fac=1./(_REAL)dptr->idunion.pers.count[i];
          distrib = dptr->idunion.pers.distr[i]*fac;
          fprintf(file, "%d %e %e\n", i+1,fac,distrib);
         }
        }
        fclose(file);
}

/************************************************************************/
void
write_pers_flory(idstruct *dptr)
{
 parmstruct   *prm = dptr->idunion.pers.prm;
 int	states = dptr->idunion.pers.cases,i;
 FILE 	*file = dptr->idunion.pers.file;
 _REAL	fac,pers;

 fac=1./(_REAL)states;

 for(i=0;i<prm->Nres-1;i++)
 {
  pers=fac*dptr->idunion.pers.distr[i];
  fprintf (file,"%d %e\n",i,pers);
 }

}
