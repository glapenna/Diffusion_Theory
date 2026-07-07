#include "head.h"

// extern char		tok[];
// char		tok[];
 
idstruct	*getstruct(), *startid();

void writedihe(), writetorch();

/************************************************************************/

void
initdihe()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(DIHE, id);

	dptr->idunion.work.type = 0;
	dptr->idunion.work.nboxes = 0;
	dptr->idunion.work.nitems = 0;
	dptr->idunion.work.cases = 0;
	dptr->idunion.work.casesw = 0.;
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
		inerr("DIHE: expected file name","");
	strcpy(dptr->idunion.work.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.work.fname, 0, &exists, "DIHE"))
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
setupdihe()
{
        extern parmstruct *prm;
	idstruct	*objptr, *serptr;
	char		errmsg[TOKLEN];
	char		id[TOKLEN];
	int i,j,k,npoints,points,nitems;
        int resnum1,resnum2,resnum3,resnum4,resnumt;
	char atnam1[2*ATOMTOKLEN],atnam2[2*ATOMTOKLEN],atnam3[2*ATOMTOKLEN],
             atnam4[2*ATOMTOKLEN], *atnamt;
        int **imatrix();
 
        gettok();
        if ((objptr = getstruct(tok)) == NULL)
                inerr("DIHE: no such id: ", tok);
        if (objptr->id_type != DIHE)
                inerr("DIHE: wrong type of id: ", tok);
        gettok();
        if ((serptr = getstruct(tok)) == NULL)
                inerr("DIHE: no such id: ", tok);
        if (serptr->id_type != STREAM)
                inerr("DIHE: wrong type of id: ", tok);

	/* get elements; atoms */
 
        npoints=4;
        points=0;
	for (nitems=0;;) {
		gettok();
		if (nitems > npoints)
			inerr(errmsg,"nitems too large");
		if (tok[0] == ';') {
			if (points  != npoints ) 
				inerr(errmsg,"points ne npoints");
			break;
		} else if (!strcmp(tok, "NINT")) {
			gettok();
	                if (!tokint(&objptr->idunion.work.nboxes))
                         inerr("DIHE: expected integer # of boxes", "");
		} else if (!strcmp(tok, "MIN")) {
			gettok();
	                if (!tokREAL(&objptr->idunion.work.min))
                         inerr("DIHE: expected min", "");
		} else if (!strcmp(tok, "MAX")) {
			gettok();
	                if (!tokREAL(&objptr->idunion.work.max))
                         inerr("DIHE: expected max", "");
                } else if (!strcmp(tok, "NTOR")) {
                        objptr->idunion.work.type=1;
                        gettok();
                        if (!tokint(&objptr->idunion.work.nelem))
                         inerr("DIHE: expected integer # of following torsions", "");
                } else if (!strcmp(tok, "ATOMS")) {
                     gettok();
                     sprintf(atnam1,"%s    ",tok);
                     gettok();
                     tokint(&resnum1,tok);
		     points++;
		     nitems++;
	for (i=0; i< prm->Nsit; i++) {
         atnamt=&prm->AtomNames[i*ATOMTOKLEN];
         resnumt=prm->ResNum[i];
           if (!strncmp(atnam1,atnamt,ATOMTOKLEN)&&resnum1==resnumt) {
            objptr->idunion.work.r1 = serptr->idunion.ser.set->crd;

/*
             printf("atom found:\n"); 
             printf("atom1: %s %d\n",atnam1,resnum1); 
             printf("atomt: %s %d\n",atnamt,resnumt); 
             printf("ATOMI %d ",i+1); 
*/
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
	for (i=0; i< prm->Nsit; i++) {
         atnamt=&prm->AtomNames[i*ATOMTOKLEN];
         resnumt=prm->ResNum[i];
           if (!strncmp(atnam2,atnamt,ATOMTOKLEN)&&resnum2==resnumt) {
            objptr->idunion.work.r2 = serptr->idunion.ser.set->crd;
/*
             printf("atom found:\n"); 
             printf("atom2: %s %d\n",atnam2,resnum2); 
             printf("atomt: %s %d\n",atnamt,resnumt); 
             printf("%d ",i+1); 
*/
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
	for (i=0; i< prm->Nsit; i++) {
         atnamt=&prm->AtomNames[i*ATOMTOKLEN];
         resnumt=prm->ResNum[i];
           if (!strncmp(atnam3,atnamt,ATOMTOKLEN)&&resnum3==resnumt) {
            objptr->idunion.work.r3 = serptr->idunion.ser.set->crd;
/*
             printf("atom found:\n");
             printf("atom3: %s %d\n",atnam3,resnum3);
             printf("atomt: %s %d\n",atnamt,resnumt);
             printf("%d 180.\n",i+1); 
*/
           }
           serptr->idunion.ser.set->crd+=3;
        }
        serptr->idunion.ser.set->crd-=3*prm->Nsit;

                     gettok();
                     sprintf(atnam4,"%s    ",tok);
                     gettok();
                     tokint(&resnum4,tok);
		     points++;
		     nitems++;
	for (i=0; i< prm->Nsit; i++) {
         atnamt=&prm->AtomNames[i*ATOMTOKLEN];
         resnumt=prm->ResNum[i];
           if (!strncmp(atnam4,atnamt,ATOMTOKLEN)&&resnum4==resnumt) {
            objptr->idunion.work.r4 = serptr->idunion.ser.set->crd;
/*
             printf("atom found:\n"); 
             printf("atom4: %s %d\n",atnam4,resnum4); 
             printf("atomt: %s %d\n",atnamt,resnumt); 
*/
           }
           serptr->idunion.ser.set->crd+=3;
        }
        serptr->idunion.ser.set->crd-=3*prm->Nsit;

		     continue;
		} else {
                 printf("dihe: case not found.\n");
                 exit(1);
		}
	}
/* set memory for boxes */
        if(objptr->idunion.work.type==1) {
         objptr->idunion.work.nmax =
           (int)pow((_REAL)objptr->idunion.work.nboxes,
                  (_REAL)objptr->idunion.work.nelem);
        		objptr->idunion.work.boxes =
                         (_REAL *) get(sizeof(_REAL)*(objptr->idunion.work.nboxes));
			for(k=0;k<objptr->idunion.work.nmax;k++) 
				objptr->idunion.work.boxes[k] = 0.;
        } else {
        		objptr->idunion.work.boxes =
                         (_REAL *) get(sizeof(_REAL)*(objptr->idunion.work.nboxes));
			for(k=0;k<objptr->idunion.work.nboxes;k++) 
				objptr->idunion.work.boxes[k] = 0.;
        }
/*
*/
	 objptr->idunion.work.max *= TORAD;
	 objptr->idunion.work.min *= TORAD;

			objptr->idunion.work.scale = 
	 ( objptr->idunion.work.max - objptr->idunion.work.min) 
				/ objptr->idunion.work.nboxes;

	objptr->idunion.work.nitems = nitems;
        objptr->idunion.work.jump = 3*prm->Natom;
        objptr->idunion.work.nlist = prm->Nmol;

/* settings for chain torsional distribution */
/*
*/
       if(objptr->idunion.work.type ==1) {
        objptr->idunion.work.value =
             (_REAL *) get(sizeof(_REAL)*(objptr->idunion.work.nelem));
        objptr->idunion.work.gauche =
             (int *) get(sizeof(int)*(objptr->idunion.work.nelem+1));
        objptr->idunion.work.statesp =
             (int *) get(sizeof(int)*(prm->Nmol));
        objptr->idunion.work.ratemat = imatrix(0,objptr->idunion.work.nmax-1,
                              0,objptr->idunion.work.nmax-1);

        for(k=0;k<objptr->idunion.work.nelem;k++)
         objptr->idunion.work.value[k] = 0.;
        for(k=0;k<objptr->idunion.work.nelem+1;k++)
         objptr->idunion.work.gauche[k] = 0;
        for(k=0;k<objptr->idunion.work.nmax;k++)
         objptr->idunion.work.statesp[k] = 9999;
        for(i=0;i<objptr->idunion.work.nmax;i++)
         for(k=0;k<objptr->idunion.work.nmax;k++)
          objptr->idunion.work.ratemat[i][k] = 0;

       }
       initstat(&objptr->idunion.work.resunion.scalar.stat);

}

/************************************************************************/

void
distribdihe(idstruct *objptr,idstruct *serptr)
{
   extern int it;
   parmstruct   *prm = serptr->idunion.ser.prm;
   int	nbond=prm->Ndihe, nmol=prm->Nmol;
   _REAL torsionbox();

   _REAL   *r = serptr->idunion.ser.set->crd;
   _REAL   *box = serptr->idunion.ser.set->boxc;
   _REAL weight,cos_phi,xx;

   int	jump = objptr->idunion.work.jump;
   int  nboxes = objptr->idunion.work.nboxes;
   _REAL   scale = objptr->idunion.work.scale, 
           min = objptr->idunion.work.min;
   _REAL   *boxes = objptr->idunion.work.boxes;
   int     ntor = objptr->idunion.work.nelem;
   int     nmax = objptr->idunion.work.nmax;
   _REAL   *tors = objptr->idunion.work.value;
   int     *gauche = objptr->idunion.work.gauche;
   int     *statesp = objptr->idunion.work.statesp;
   int     state;
   _REAL *r1=objptr->idunion.work.r1;
   _REAL *r2=objptr->idunion.work.r2;
   _REAL *r3=objptr->idunion.work.r3;
   _REAL *r4=objptr->idunion.work.r4;
   _REAL *rt1,*rt2,*rt3,*rt4,*rt1t,*rt2t,*rt3t,*rt4t;

   int im,k,ibox,igauche,ip;
   _REAL boxm1[3],phi;

   if(serptr->idunion.ser.iweight)
   {
    weight=serptr->idunion.ser.w->idunion.serw.w;
   } else {
    weight=1.;
   }

   for(k=0;k<3;k++) boxm1[k]=1./box[k];
   rt1=r1; rt2=r2; rt3=r3; rt4=r4;

   for (im=0;im<nmol;im++) {
    ++objptr->idunion.work.cases;
    objptr->idunion.work.casesw+=weight;
    phi = torsionbox(rt1,rt2,rt3,rt4,box,boxm1);
    cos_phi=l_cos(phi)*weight;
    xx= 1.9-1.4*cos_phi+6.4*cos_phi*cos_phi;
    objptr->idunion.work.resunion.scalar.value = phi;
    updatestat(&objptr->idunion.work.resunion.scalar.stat, phi); 

    if(objptr->idunion.work.type==0) {
     ibox = (int)((phi-min)/scale);
     boxes[ibox]+=weight;
     /* boxes[ibox]+=1.; */
    } else if(objptr->idunion.work.type==1) {
     tors[0]=phi;
     rt1t = rt1+3; rt2t = rt2+3; rt3t = rt3+3; rt4t = rt4+3;
     for(k=1; k<ntor; k++) {
      phi = torsionbox(rt1t,rt2t,rt3t,rt3t,box,boxm1);
      tors[k] = phi;
      rt1t += 3; rt2t += 3; rt3t += 3; rt4t += 3;
     }

     for(k=0,ibox=0,igauche=0; k<ntor; k++) {
         igauche+=abs((int)(tors[k]/scale)-1);
/*
         ibox+=(int)(pow((_REAL)nboxes,(_REAL)k)*(int)((tors[k]-min)/scale));
*/
/* reverse the order */
         ibox+=(int)(pow((_REAL)nboxes,(_REAL)(ntor-1-k))*(int)((tors[k]-min)/scale));
      }
           state=ibox;
           gauche[igauche]+=1;
           if (ibox >= nmax) {
            printf("DISTRIBTORCH: ibox= %d exceeds nmax= %d\n",ibox,nmax);
            printf("DISTRIBTORCH: phi= %e\n",tors[0]);
            exit(1);
           }
           boxes[ibox]+=1.; 
           if(it>0) {
            if(statesp[im] != state) {
             ++objptr->idunion.work.glrate; 
/*
           printf("it,state: %d %d %d\n",it,state,objptr->idunion.work.glrate); 
*/
/* fill ratematrix */
             for (k=0;k<nmax;k++) {
              for (ip=0;ip<nmax;ip++) {
               if(statesp[im]==k && state==ip) {
                ++objptr->idunion.work.ratemat[k][ip]; 
               }
              }
             }
            }
           }
           statesp[im]=state;
    } /* endif type */
    rt1+=jump; rt2+=jump; rt3+=jump; rt4+=jump;
   } /* enddo im */
}

/************************************************************************/

void
writedihe(idstruct *dptr)
{
	int	i, nboxes = dptr->idunion.work.nboxes;
	_REAL   *boxes = dptr->idunion.work.boxes;
	int	states = dptr->idunion.work.cases;
	FILE 	*file = dptr->idunion.work.file;
	_REAL	scale = dptr->idunion.work.scale, min=dptr->idunion.work.min;

	_REAL	ri,distrib,fac;

	/* fac = 1./ ((_REAL)states); */
	fac = 1./ dptr->idunion.work.casesw;

	for (i=0; i<nboxes; i++) {
          ri = min+(_REAL)i*scale+0.5*scale;
          ri *= TODEG;
          distrib = boxes[i]*fac;
 	  fprintf(file, "%e %e\n", ri,distrib); 
	}
	fclose(file);
}

/************************************************************************/

void
writetorch(idstruct *dptr)
{
	int	nboxes = dptr->idunion.work.nboxes;
	int	ntorch = dptr->idunion.work.nelem;
	int	nmax = dptr->idunion.work.nmax;
	_REAL   *boxes = dptr->idunion.work.boxes;
	int	states = dptr->idunion.work.cases;
	int	ntor = dptr->idunion.work.nlist;
	int	*gauche = dptr->idunion.work.gauche;
	int	count = dptr->idunion.work.glrate;
	FILE 	*file = dptr->idunion.work.file;
	_REAL	scale = dptr->idunion.work.scale;
/*
*/
	int	**ratemat = dptr->idunion.work.ratemat;
        FILE *filetemp;

	_REAL	distrib,fac,fact;
        int i,k,quoz,rem;

	fac = 1./ ((_REAL)states * (_REAL)ntor);
	fact = 1./ ((_REAL)states);
        filetemp = fopen ("ratematrix","w");

        printf("number of transitions: %d\n",count);
        printf("number of transitions: %12.6e\n",fact*(_REAL)count);
        printf("rate (W) matrix:\n");
	for (i=0,distrib=0.; i<nmax; i++) {
	 for (k=0; k<nmax; k++) {
          fprintf(filetemp,"%d %d %12.6e\n",i,k,fact*(_REAL)
           dptr->idunion.work.ratemat[i][k]);
          distrib+=fact*(_REAL)dptr->idunion.work.ratemat[i][k];
         }
        }
	printf("\n");
	printf("somma: %12.6e\n",distrib);
/*
        filetemp = fopen ("ratematrix","w");
        fwrite(ratemat,sizeof(int),243*243,filetemp);
*/

        printf("gauche population:\n");
	for (i=0; i<ntorch+1; i++) {
         printf("number of gauche: %d population: %d\n",i,gauche[i]);
	 gauche[i]=-1;
        }

	for (i=0; i<nmax; i++) {
         quoz=i;
         for(k=0;k<ntorch;k++) {
          rem = (int)fmod((_REAL)quoz,(_REAL)nboxes);
          quoz = (int)(quoz/nboxes);
          *gauche = rem-1;
          gauche++;
         }
/*
         gauche-=ntorch;
*/
         distrib = boxes[i]*fac;
/*
         for(k=0;k<ntorch;k++) fprintf(file, "%+d ", *gauche++); 
*/
         for(k=0;k<ntorch;k++) fprintf(file, "%+d ", *(--gauche)); 
	 fprintf(file, "%e\n", distrib); 
	}
        gauche-=ntorch;
	fclose(file);
}
