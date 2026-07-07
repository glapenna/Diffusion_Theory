#include "head.h"

// extern char		tok[];
// char		tok[];
 
idstruct	*getstruct(), *startid();

void writebond();

// FILE *ftmp;

/************************************************************************/

void
initbond()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(BOND, id);

	dptr->idunion.work.type = 0;
	dptr->idunion.work.nboxes = 0;
	dptr->idunion.work.nitems = 0;
	dptr->idunion.work.cases = 0;
	dptr->idunion.work.casesw = 0.;
	dptr->idunion.work.scale = 0.;
	dptr->idunion.work.min = 0.;
	dptr->idunion.work.max = 5.;
        dptr->idunion.work.resunion.scalar.value = 0.;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("BOND: expected file name","");
	strcpy(dptr->idunion.work.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.work.fname, 0, &exists, "BOND"))
		exit(1);
	if ((dptr->idunion.work.file = 
				  fopen(dptr->idunion.work.fname, "w")) == NULL){
		perror(dptr->idunion.work.fname);
		exit(1);
	}
	gettok();
        /* ftmp=fopen("bond.time","wb"); */

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************/

void
setupbond(int npoints)
{
        extern parmstruct *prm;
	idstruct	*objptr, *idptr,*serptr;
	char		errmsg[TOKLEN];
	char		id[TOKLEN];
	int i, j, k, ib,jb, nitems, points = 0,ifound;
//	char atnam1[2*ATOMTOKLEN],atnam2[2*ATOMTOKLEN],*atnamt;
	char atnam1[ATOMTOKLEN],atnam2[ATOMTOKLEN],*atnamt;
        int resnum1,resnum2,resnumt;
 
        gettok();
        if ((objptr = getstruct(tok)) == NULL)
                inerr("BOND: no such id: ", tok);
        if (objptr->id_type != BOND)
                inerr("BOND: wrong type of id: ", tok);

        gettok();
        if ((serptr = getstruct(tok)) == NULL)
                inerr("BOND: no such id: ", tok);
        if (serptr->id_type != STREAM)
                inerr("BOND: wrong type of id: ", tok);


	/* get elements; atoms */
 
	for (nitems=0;;) {
  
		gettok();
  
		if (nitems > npoints)
			inerr(errmsg,"");
 
		if (tok[0] == ';') {
			if (points  != npoints ) 
				inerr(errmsg,"");
			break;
		} else if (!strcmp(tok, "NINT")) {
			gettok();
	                if (!tokint(&objptr->idunion.work.nboxes))
                         inerr("BOND: expected integer # of boxes", "");
		} else if (!strcmp(tok, "MIN")) {
			gettok();
	                if (!tokREAL(&objptr->idunion.work.min))
                         inerr("BOND: expected min", "");
		} else if (!strcmp(tok, "MAX")) {
			gettok();
	                if (!tokREAL(&objptr->idunion.work.max))
                         inerr("BOND: expected max", "");

		} else if ((idptr = getstruct(tok)) == NULL) {

          // sprintf(atnam1,"%s    ",tok);
	  strcpy(atnam1,tok);
          gettok();
          tokint(&resnum1,tok);
          printf("atnam1: %s %d\n",atnam1,resnum1);
			points++;
			nitems++;
          atnamt = serptr->idunion.ser.prm->AtomNames;
          for(i=0;i<prm->Nsit;i++) {
           // atnamt = &prm->AtomNames[i*ATOMTOKLEN];
           resnumt=prm->ResNum[i];
           if (!strncmp(atnam1,atnamt,ATOMTOKLEN)&&resnum1==resnumt) {
            objptr->idunion.work.r1 = serptr->idunion.ser.set->crd;
            printf("atom found:\n");
            ifound++;
            printf("atom1: %s %d\n",atnam1,resnum1);
           }
	   atnamt+=ATOMTOKLEN;
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*prm->Nsit;

          gettok();
          // sprintf(atnam2,"%s    ",tok);
	  strcpy(atnam2,tok);
          gettok();
          tokint(&resnum2,tok);
          printf("atnam2: %s %d\n",atnam2,resnum2);
			points++;
			nitems++;
          atnamt = serptr->idunion.ser.prm->AtomNames;
          for(i=0;i<prm->Nsit;i++) {
           atnamt = &prm->AtomNames[i*ATOMTOKLEN];
           resnumt=prm->ResNum[i];
           // printf("atnam2t: %s %d\n",atnamt,resnumt);
           if (!strncmp(atnam2,atnamt,ATOMTOKLEN)&&resnum2==resnumt) {
            objptr->idunion.work.r2 = serptr->idunion.ser.set->crd;
            printf("atom found:\n");
            ifound++;
            printf("atom2: %s %d\n",atnam2,resnum2);
           }
	   atnamt+=ATOMTOKLEN;
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*prm->Nsit;

          continue;
         }
        }
/* set memory for boxes */
        		objptr->idunion.work.boxes =
                         (_REAL *) get(sizeof(_REAL)*(objptr->idunion.work.nboxes));
			for(k=0;k<objptr->idunion.work.nboxes;k++) 
				objptr->idunion.work.boxes[k] = 0.;
			objptr->idunion.work.scale = 
			 ( objptr->idunion.work.max - objptr->idunion.work.min) 
				/ objptr->idunion.work.nboxes;

	objptr->idunion.work.nitems = nitems;
	objptr->idunion.work.type = BOND;
        objptr->idunion.work.jump = 3*prm->Natom;
        initstat(&objptr->idunion.work.resunion.scalar.stat);
        if(serptr->idunion.ser.iweight) objptr->idunion.work.resunion.scalar.stat.weight=1;

}

/************************************************************************/

void
distribbond(idstruct *objptr,idstruct *serptr)
{
   parmstruct   *prm = serptr->idunion.ser.prm;
   int  nmol=prm->Nmol;
   _REAL distancebox();
   void images(_REAL *dist, idstruct *serptr);
   void images_to(_REAL *dist, _REAL *box, _REAL *boxm1);
   _REAL distance();

   _REAL   *box = serptr->idunion.ser.set->boxc;

   _REAL *r1 = objptr->idunion.work.r1;
   _REAL *r2 = objptr->idunion.work.r2;
   int  jump = objptr->idunion.work.jump;
   _REAL   scale = objptr->idunion.work.scale,
           min = objptr->idunion.work.min,
           max = objptr->idunion.work.max;
   _REAL   *boxes = objptr->idunion.work.boxes;

   int im,k,ibox;
   _REAL boxm1[3],dist12[3],dist,*rt1,*rt2,weight,dist2;

   if(serptr->idunion.ser.iweight)
   {
    weight=serptr->idunion.ser.w->idunion.serw.w;
   } else {
    weight=1.;
   }
   for(k=0;k<3;k++) boxm1[k]=1./box[k];
   rt1=r1; rt2=r2;

   for (im=0;im<nmol;im++) {
    ++objptr->idunion.work.cases;
    objptr->idunion.work.casesw+=weight;
/*
    dist = distance(rt1,rt2,box,boxm1);
    dist = distancebox(rt1,rt2,box,boxm1);
*/
    
    dist12[0]=rt1[0]-rt2[0];
    dist12[1]=rt1[1]-rt2[1];
    dist12[2]=rt1[2]-rt2[2];
#ifdef PBC
    images(dist12,serptr);
#endif
#ifdef PBCTO
    images_to(dist12,box,boxm1);
#endif
    dist2 = (dist12[0]*dist12[0]+dist12[1]*dist12[1]+dist12[2]*dist12[2]);
    dist = sqrt(dist2);
/*
    dist = dist*dist*dist;
    dist = 1./dist;
*/

    if(dist>=min && dist<=max)
    {
     ibox = (int)((dist-min)/scale);
     boxes[ibox]+=weight;
    }
    /* fwrite(&dist,sizeof(_REAL),1,ftmp); */
    dist *= weight;
    objptr->idunion.work.resunion.scalar.value = dist;
    updatestat(&objptr->idunion.work.resunion.scalar.stat, dist);
    rt1 += jump;
    rt2 += jump;
   } /* enddo im */
}

/************************************************************************
							WRITEbond()
*************************************************************************/

void
writebond(idstruct *dptr)
{
	int	i, nboxes = dptr->idunion.work.nboxes;
	_REAL   *boxes = dptr->idunion.work.boxes;
	int	states = dptr->idunion.work.cases;
	FILE 	*file = dptr->idunion.work.file;
	_REAL	scale = dptr->idunion.work.scale, min=dptr->idunion.work.min;

	_REAL	ri,distrib,fac;

	fac = 1./ dptr->idunion.work.casesw;

	for (i=0; i<nboxes; i++) {
          ri = min+(_REAL)i*scale+0.5*scale;
          distrib = boxes[i]*fac;
 	  fprintf(file, "%e %e\n", ri,distrib); 
	}
	fclose(file);
	/* fclose(ftmp); */
}
