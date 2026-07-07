#include "head.h"

// extern char		tok[];
// char		tok[];
 
idstruct	*getstruct(), *startid();

void writedensity1();

/************************************************************************/

void
initdensity1()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(DENSITY1, id);

	dptr->idunion.work.type = 0;
	dptr->idunion.work.nboxes = 0;
	dptr->idunion.work.nitems = 0;
	dptr->idunion.work.cases = 0;
	dptr->idunion.work.scale = 0.;
	dptr->idunion.work.min = 0.;
	dptr->idunion.work.max = 5.;
	dptr->idunion.work.nmax = 0;
	dptr->idunion.work.nelem = 0;
	dptr->idunion.work.nlist = 0;
	dptr->idunion.work.glrate = 0;
        dptr->idunion.work.resunion.scalar.value = 0.;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("DENSITY1: expected file name","");
	strcpy(dptr->idunion.work.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.work.fname, 0, &exists, "DENSITY1"))
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
setupdensity1()
{
        extern parmstruct *prm;
	idstruct	*objptr, *idptr, *serptr;
	char		errmsg[TOKLEN];
	char		id[TOKLEN];
	int i, j, k;
	char atnamt[ATOMTOKLEN], atnam[ATOMTOKLEN];
 
        gettok();
        if ((objptr = getstruct(tok)) == NULL)
                inerr("DENSITY1: no such id: ", tok);
        if (objptr->id_type != DENSITY1)
                inerr("DENSITY1: wrong type of id: ", tok);
        gettok();
        if ((serptr = getstruct(tok)) == NULL) 
                inerr("setupdensity1: no such id\n","");
        if (serptr->id_type != STREAM)
                inerr("setupdensity1: stream id expected\n","");

	for (i=0; i<ATOMTOKLEN; i++)
		atnam[i] = NULL;

	/* get elements; atoms */
 
	for (;;) {
  
		gettok();
  
		if (tok[0] == ';') {
			break;
		} else if (!strcmp(tok, "NINT")) {
			gettok();
	                if (!tokint(&objptr->idunion.work.nboxes))
                         inerr("DENSITY1: expected integer # of boxes", "");
		} else if (!strcmp(tok, "MIN")) {
			gettok();
	                if (!tokREAL(&objptr->idunion.work.min))
                         inerr("DENSITY1: expected min", "");
		} else if (!strcmp(tok, "MAX")) {
			gettok();
	                if (!tokREAL(&objptr->idunion.work.max))
                         inerr("DENSITY1: expected max", "");
                } else if (!strcmp(tok, "NAT")) {
                        objptr->idunion.work.type=0;
                        gettok();
                        if (!tokint(&objptr->idunion.work.nelem))
                         inerr("DENSITY1: expected integer # of following atoms", "");
                } else if (!strcmp(tok, "ATOMNAME")) {
                        gettok();
			strcat(atnam,tok);
                        for(k=strlen(tok);k<ATOMTOKLEN;k++)
                         strcat(atnam," ");
/*
*/
        		printf("atnam  %s\n",atnam);
		}
	}
        objptr->idunion.work.jump = 3*prm->Natom;
        objptr->idunion.work.nlist = prm->Nmol;
/* set memory for boxes */
       		objptr->idunion.work.boxes =
                        (_REAL *) get(sizeof(_REAL)*(objptr->idunion.work.nboxes));
		for(k=0;k<objptr->idunion.work.nboxes;k++) 
			objptr->idunion.work.boxes[k] = 0.;
		objptr->idunion.work.scale = 
	 ( objptr->idunion.work.max - objptr->idunion.work.min) 
				/ objptr->idunion.work.nboxes;

/* look for the atnam1 string in angle array */
       for (i=0; i< prm->Natom; i++) {
        for (j=0; j<ATOMTOKLEN;j++) atnamt[j]= prm->AtomNames[ATOMTOKLEN*i+j];
        if (!strncmp(atnam,atnamt,ATOMTOKLEN)) {
         objptr->idunion.work.r1 = serptr->idunion.ser.set->crd;
        }
        serptr->idunion.ser.set->crd+=3;
       }
       serptr->idunion.ser.set->crd-=3*prm->Natom;
       initstat(&objptr->idunion.work.resunion.scalar.stat);
}

/************************************************************************/

void
distribdensity1(idstruct *objptr,idstruct *serptr)
{
   extern int it;
   parmstruct   *prm = serptr->idunion.ser.prm;
   int	ngroup=objptr->idunion.work.nelem, nmol=prm->Nmol;

   _REAL   *dir = serptr->idunion.ser.set->dir;
   _REAL   *box = serptr->idunion.ser.set->boxc;

   _REAL *r1 = objptr->idunion.work.r1;
   _REAL  *rt,*rtt;
   int	jump = objptr->idunion.work.jump;
   int  nboxes = objptr->idunion.work.nboxes;
   _REAL   scale = objptr->idunion.work.scale, 
           min = objptr->idunion.work.min;
   _REAL   *boxes = objptr->idunion.work.boxes;

   int im,i,k,ibox;
   _REAL x,y,z,scp;

   rt = r1;
   for (im=0;im<nmol;im++) {
    rtt=rt;
    for (i=0;i<ngroup;i++) {
     ++objptr->idunion.work.cases;
     x=*rtt++; y=*rtt++; z=*rtt++;
     
     scp = dir[0]*x+dir[1]*y+dir[2]*z; 
   /*  printf("scp= %f \n",scp);*/
     updatestat(&objptr->idunion.work.resunion.scalar.stat, Nabs(scp)); 
     ibox = (int)((scp-min)/scale);
     boxes[ibox]+=1.;
    } /* enddo i */
    rt+=jump;
   } /* enddo im */
}

/************************************************************************/

void
writedensity1(idstruct *dptr)
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
          distrib = boxes[i]*fac;
 	  fprintf(file, "%e %e\n", ri,distrib); 
	}
	fclose(file);
}
