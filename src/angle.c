#include "head.h"

// extern char		tok[];
// char		tok[];
 
idstruct	*getstruct(), *startid();

void writeang();

/************************************************************************/

void
initang()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(ANGLE, id);

	dptr->idunion.work.type = 0;
	dptr->idunion.work.nboxes = 0;
	dptr->idunion.work.nitems = 0;
	dptr->idunion.work.cases = 0;
	dptr->idunion.work.scale = 0.;
	dptr->idunion.work.min = 0.;
	dptr->idunion.work.max = 5.;
        dptr->idunion.work.resunion.scalar.value = 0.;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("ANGLE: expected file name","");
	strcpy(dptr->idunion.work.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.work.fname, 0, &exists, "ANGLE"))
		exit(1);
	if ((dptr->idunion.work.file = 
				  fopen(dptr->idunion.work.fname, "w")) == NULL){
		perror(dptr->idunion.work.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/*************************************************************************/

void
setupang(int npoints)
{
        extern parmstruct *prm;
	idstruct	*objptr, *idptr, *serptr;
	char		errmsg[TOKLEN];
	char		id[TOKLEN];
	int i, j, k, ifound, nitems, points = 0;
        char atnam1[2*ATOMTOKLEN],atnam2[2*ATOMTOKLEN],atnam3[2*ATOMTOKLEN],*atnamt;
        int resnum1,resnum2,resnum3,resnumt;
 
        gettok();
        if ((objptr = getstruct(tok)) == NULL)
                inerr("ANGLE: no such id: ", tok);
        if (objptr->id_type != ANGLE)
                inerr("ANGLE: wrong type of id: ", tok);

        gettok();
        if ((serptr = getstruct(tok)) == NULL)
                inerr("ANGLE: no such id: ", tok);
        if (serptr->id_type != STREAM)
                inerr("ANGLE: wrong type of id: ", tok);


	/* get elements; atoms */
 
	for (nitems=0;;) {
 
  
		gettok();
  
		if (nitems > npoints)
			inerr(errmsg,"nitems");
 
		if (tok[0] == ';') {
			if (points  != npoints ) 
				inerr(errmsg,"points");
			break;
		} else if (!strcmp(tok, "NINT")) {
			gettok();
	                if (!tokint(&objptr->idunion.work.nboxes))
                         inerr("ANGLE: expected integer # of boxes", "");
		} else if (!strcmp(tok, "MIN")) {
			gettok();
	                if (!tokREAL(&objptr->idunion.work.min))
                         inerr("ANGLE: expected min", "");
		} else if (!strcmp(tok, "MAX")) {
			gettok();
	                if (!tokREAL(&objptr->idunion.work.max))
                         inerr("ANGLE: expected max", "");

		} else if ((idptr = getstruct(tok)) == NULL) {
 
          sprintf(atnam1,"%s    ",tok);
          gettok();
          tokint(&resnum1,tok);
                        points++;
                        nitems++;
          for(i=0;i<prm->Nsit;i++) {
           atnamt = &prm->AtomNames[i*ATOMTOKLEN];
           resnumt=prm->ResNum[i];
           if (!strncmp(atnam1,atnamt,ATOMTOKLEN)&&resnum1==resnumt) {
            objptr->idunion.work.r1 = serptr->idunion.ser.set->crd;
            printf("atom found:\n");
            ifound++;
            printf("atom1: %s %d\n",atnam1,resnum1);
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*prm->Nsit;

          gettok();
          sprintf(atnam2,"%s    ",tok);
          gettok();
          tokint(&resnum2,tok);
                        points++;
                        nitems++;
          for(i=0;i<prm->Nsit;i++) {
           atnamt = &prm->AtomNames[i*ATOMTOKLEN];
           resnumt=prm->ResNum[i];
           if (!strncmp(atnam2,atnamt,ATOMTOKLEN)&&resnum2==resnumt) {
            objptr->idunion.work.r2 = serptr->idunion.ser.set->crd;
            printf("atom found:\n");
            ifound++;
            printf("atom2: %s %d\n",atnam2,resnum2);
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*prm->Nsit;

          gettok();
          sprintf(atnam3,"%s    ",tok);
          gettok();
          tokint(&resnum3,tok);
                        points++;
                        nitems++;
          for(i=0;i<prm->Nsit;i++) {
           atnamt = &prm->AtomNames[i*ATOMTOKLEN];
           resnumt=prm->ResNum[i];
           if (!strncmp(atnam3,atnamt,ATOMTOKLEN)&&resnum3==resnumt) {
            objptr->idunion.work.r3 = serptr->idunion.ser.set->crd;
            printf("atom found:\n");
            ifound++;
            printf("atom3: %s %d\n",atnam3,resnum3);
           }
           serptr->idunion.ser.set->crd+=3;
          }
          serptr->idunion.ser.set->crd-=3*prm->Nsit;
	}
	}
/* set memory for boxes */
        		objptr->idunion.work.boxes =
                         (_REAL *) get(sizeof(_REAL)*(objptr->idunion.work.nboxes));
			for(k=0;k<objptr->idunion.work.nboxes;k++) 
				objptr->idunion.work.boxes[k] = 0.;
/*
*/
	 objptr->idunion.work.max *= TORAD;
	 objptr->idunion.work.min *= TORAD;
			objptr->idunion.work.scale = 
	 ( objptr->idunion.work.max - objptr->idunion.work.min) 
				/ objptr->idunion.work.nboxes;

	objptr->idunion.work.nitems = nitems;
        objptr->idunion.work.jump = prm->Natom;

}

/************************************************************************/

void
distribang(idstruct *objptr,idstruct *serptr)
{
   parmstruct   *prm = serptr->idunion.ser.prm;
   int  nmol=prm->Nmol;
   _REAL anglebox(),cos_theta;

   _REAL   *box = serptr->idunion.ser.set->boxc;

   _REAL *r1 = objptr->idunion.work.r1;
   _REAL *r2 = objptr->idunion.work.r2;
   _REAL *r3 = objptr->idunion.work.r3;
   int  jump = objptr->idunion.work.jump;
   _REAL   scale = objptr->idunion.work.scale,
           min = objptr->idunion.work.min;
   _REAL   *boxes = objptr->idunion.work.boxes;

   int im,k,ibox;
   _REAL boxm1[3],dist12,*rt1,*rt2,*rt3;

   for(k=0;k<3;k++) boxm1[k]=1./box[k];
   rt1=r1; rt2=r2; rt3=r3;

   for (im=0;im<nmol;im++) {
    ++objptr->idunion.work.cases;
    dist12 = anglebox(rt1,rt2,rt3,box,boxm1);
/*
    cos_theta = l_cos(dist12);
    if(cos_theta<-1.) cos_theta=-1.;
    if(cos_theta> 1.) cos_theta= 1.;
    objptr->idunion.work.resunion.scalar.value = cos_theta;
*/
    objptr->idunion.work.resunion.scalar.value = dist12;
    updatestat(&objptr->idunion.work.resunion.scalar.stat, dist12);
    ibox = (int)((dist12-min)/scale);
    boxes[ibox]+=1.;
    rt1 += jump;
    rt2 += jump;
    rt3 += jump;
   } /* enddo im */
}

/************************************************************************/

void
writeang(idstruct *dptr)
{
	int	i, nboxes = dptr->idunion.work.nboxes;
	_REAL   *boxes = dptr->idunion.work.boxes;
	int	states = dptr->idunion.work.cases;
	FILE 	*file = dptr->idunion.work.file;
	_REAL	scale = dptr->idunion.work.scale, min=dptr->idunion.work.min;

	_REAL	ri,distrib,fac;

	fac = 1./ ((_REAL)states);

	for (i=0; i<nboxes; i++) {
          ri = min+(_REAL)i*scale+0.5*scale;
          ri *= TODEG;
          distrib = boxes[i]*fac;
 	  fprintf(file, "%e %e\n", ri,distrib); 
	}
	fclose(file);
}
