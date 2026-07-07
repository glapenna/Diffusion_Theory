#include "head.h"


// extern char		tok[];
// char		tok[];
// void *get2(int,int);
extern int		it;

extern idstruct		*startid(), *getstruct();

void	writegr();
/************************************************************************
							INITGR()
*************************************************************************/

void
initgr()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		min = 0, exists;

	/*
	 *  SYNTAX:
	 *
	 *		distrib_id file_name 
	 */

	/* 'distrib_id' */

	dptr = startid(GR, id);

	/* INITIALIZE */

	dptr->idunion.db.type = 0;
	dptr->idunion.db.max = 0.0;
	dptr->idunion.db.min = 0.0;
	dptr->idunion.db.maxz = 0.0;
	dptr->idunion.db.minz = 0.0;
	dptr->idunion.db.nboxes = 0;
	dptr->idunion.db.win = 0;
	dptr->idunion.db.winct = 0;
	dptr->idunion.db.cases = 0;
	dptr->idunion.db.casesw = 0.;
	dptr->idunion.db.rho = 0.;
	dptr->idunion.db.grid = 0;
	dptr->idunion.db.smallside = 10000.;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("initgr: expected file name","");
	strcpy(dptr->idunion.db.fname1, tok);
/*
*/
	gettok();
	if (iskeywd())
		inerr("initgr: expected file name","");
	strcpy(dptr->idunion.db.fname2, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.db.fname1, 0, &exists, "GR"))
		exit(1);
	if ((dptr->idunion.db.file1 = 
				  fopen(dptr->idunion.db.fname1, "w")) == NULL){
		perror(dptr->idunion.db.fname1);
		exit(1);
	}
/*
*/
	if (fileprob(dptr->idunion.db.fname2, 0, &exists, "GR"))
		exit(1);
	if ((dptr->idunion.db.file2 = 
				  fopen(dptr->idunion.db.fname2, "w")) == NULL){
		perror(dptr->idunion.db.fname2);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************
							SETUPGR()
*************************************************************************/

void
setupgr()
{
        extern parmstruct *prm;
	idstruct	*dbptr,*serptr;
	int		i,j,k,kk;
        char atnam1[2*ATOMTOKLEN],atnam2[2*ATOMTOKLEN],*atnamt;

	/*
	 *  SYNTAX:
	 *  		distrib_id
	 *		min max nboxes [WIN sets]
	 */

	/*
	 *  get db struct from previous section
	 *	(FILES_OUT DISTRIBUTION tok_dbid ..)
	 */


	/* 'distrib_id' */
	gettok();
	if ((dbptr = getstruct(tok)) == NULL)
		inerr("setupgr: no such id: ", tok);
	if (dbptr->id_type != GR)
		inerr("setupgr: wrong type of id: ", tok);

/* get stream */
	gettok();
	if ((serptr = getstruct(tok)) == NULL)
		inerr("setupgr: no such id: ", tok);
	if (serptr->id_type != STREAM)
		inerr("setupgr: wrong type of id: ", tok);

	/* CHECK */
	/*
	 *  make sure it has not already been referenced in OUTPUT section
	 */

        dbptr->idunion.db.type = 0;

	if (dbptr->id_setup)
		inerr("setupgr: id multiply referenced in OUTPUT", "");

               for (;;) {
	        gettok();
		if (tok[0] == ';') {
                 break;
		} else if (!strcmp(tok, "NINT")) {
		 gettok();
	         if (!tokint(&dbptr->idunion.db.nboxes))
                  inerr("setupgr: expected integer # of boxes", "");
		} else if (!strcmp(tok, "MAX")) {
		 gettok();
	         if (!tokREAL(&dbptr->idunion.db.max))
                  inerr("setupgr: expected max", "");
		} else if (!strcmp(tok, "MIN")) {
		 gettok();
	         if (!tokREAL(&dbptr->idunion.db.min))
                  inerr("setupgr: expected min", "");
		} else if (!strcmp(tok, "NINTZ")) {
		 gettok();
	         if (!tokint(&dbptr->idunion.db.nboxesz))
                  inerr("setupgr: expected integer # of boxes", "");
		} else if (!strcmp(tok, "MAXZ")) {
		 gettok();
	         if (!tokREAL(&dbptr->idunion.db.maxz))
                  inerr("setupgr: expected max", "");
		} else if (!strcmp(tok, "MINZ")) {
		 gettok();
	         if (!tokREAL(&dbptr->idunion.db.minz))
                  inerr("setupgr: expected min", "");
		} else if (!strcmp(tok, "GRID")) {
	         dbptr->idunion.db.grid=1;
		 gettok();
	         if (!tokint(&dbptr->idunion.db.nboxes))
                  inerr("setupgr: expected integer # of boxes", "");
		 dbptr->idunion.db.rmin = 
                  (_REAL *) get2(dbptr->idunion.db.nboxes,sizeof(_REAL));
		 dbptr->idunion.db.rmax = 
                  (_REAL *) get2(dbptr->idunion.db.nboxes,sizeof(_REAL));
                 for(j=0;j<dbptr->idunion.db.nboxes;j++)
                 {
		  gettok();
	          if (!tokREAL(&dbptr->idunion.db.rmin[j]))
                   inerr("setupgr: expected rmin", "");
		  gettok();
	          if (!tokREAL(&dbptr->idunion.db.rmax[j]))
                   inerr("setupgr: expected rmax", "");
                 }
		} else if (!strcmp(tok, "WIN")) {
			gettok();
			if (!tokint(&dbptr->idunion.db.win))
		    		inerr("setupgr/WIN: expected",
					" # of sets per curve");
			if (dbptr->idunion.db.win < 1)
		    		inerr("setupgr/WIN: need > 0 sets",
					" for window");
		} else if (!strcmp(tok, "ATOMNAMES")) {
                  dbptr->idunion.db.type = 1;
/* read atoms */
                        gettok();
                        sprintf(atnam1,"%s    ",tok);
                        gettok();
                        sprintf(atnam2,"%s    ",tok);
		/*
                        gettok();
                        strcat(atnam1,tok);
                        for (k=strlen(tok);k<ATOMTOKLEN;k++)
                           strcat(atnam1," ");
                        gettok();
                        strcat(atnam2,tok);
                        for (k=strlen(tok);k<ATOMTOKLEN;k++)
                           strcat(atnam2," ");
                */
                } else {
                  inerr("setupgr: unexpected token", "");
                  break;
                }
               }

	       if (tok[0] != ';') {
                 printf("setupgr: ; expected\n");
                 exit(1);
	       }

       if(dbptr->idunion.db.type) {
/* look for atom */
          for(i=0;i<prm->Natom;i++) {
           atnamt = &prm->AtomNames[i*ATOMTOKLEN];
           if (!strncmp(atnam1,atnamt,ATOMTOKLEN)) {
            dbptr->idunion.db.listi = serptr->idunion.ser.set->crd;
            printf("atom %d found.\n",i);
           }
           if (!strncmp(atnam2,atnamt,ATOMTOKLEN)) {
            dbptr->idunion.db.listj = serptr->idunion.ser.set->crd;
            printf("atom %d found.\n",i);
           }
           serptr->idunion.ser.set->crd+=3;
          }
/*
        printf("gr of pair %s %s\n",atnam1,atnam2);
        printf("pointers %d %d\n",
                dbptr->idunion.db.listi,
                dbptr->idunion.db.listj);
*/
        serptr->idunion.ser.set->crd-=3*prm->Natom;
        dbptr->idunion.db.jump = 3*(prm->Natom);
       }

		dbptr->idunion.db.boxes = 
			(_REAL *) get(dbptr->idunion.db.nboxes *
						sizeof(_REAL));
		dbptr->idunion.db.boxespa = 
			(_REAL *) get(dbptr->idunion.db.nboxesz *
						sizeof(_REAL));
/*
*/
		dbptr->idunion.db.boxespe = 
			(_REAL *) get(dbptr->idunion.db.nboxes *
						sizeof(_REAL));
		for (i=0; i<dbptr->idunion.db.nboxes; i++) {
			dbptr->idunion.db.boxes[i] = 0.;
/*
*/
			dbptr->idunion.db.boxespe[i] = 0.;
		}
/*
*/
		for (i=0; i<dbptr->idunion.db.nboxesz; i++) {
			dbptr->idunion.db.boxespa[i] = 0.;
		}

/* average box sides */
        dbptr->idunion.db.boxav = (_REAL *)get(sizeof(_REAL)*3);
		for (i=0; i<3; i++) dbptr->idunion.db.boxav[i] = 0.;

		/*
		 *  derive the scale factor (size of the bins)
		 */

		dbptr->idunion.db.scale = 
			(dbptr->idunion.db.max - dbptr->idunion.db.min) 
			/ (_REAL) dbptr->idunion.db.nboxes;
/*
*/
		dbptr->idunion.db.scalez = 
			(dbptr->idunion.db.maxz - dbptr->idunion.db.minz) 
			/ (_REAL) dbptr->idunion.db.nboxesz;

	dbptr->id_setup++;
}

/************************************************************************
							gr()
*************************************************************************/

void
gr(idstruct *dptr,idstruct *serptr)
{
        extern parmstruct *prm;
        int nmol=prm->Nmol;
        int nsit=prm->Nsit;
	_REAL	max = dptr->idunion.db.max;
	_REAL	min = dptr->idunion.db.min;
	_REAL	maxz = dptr->idunion.db.maxz;
	_REAL	minz = dptr->idunion.db.minz;
	_REAL	scale = dptr->idunion.db.scale;
	_REAL	scalez = dptr->idunion.db.scalez;
	_REAL	small = dptr->idunion.db.smallside;
	_REAL	*listi = dptr->idunion.db.listi;
	_REAL	*listj = dptr->idunion.db.listj;
	int	jump = dptr->idunion.db.jump;
	_REAL   *boxav = dptr->idunion.db.boxav;
	_REAL	*boxes = dptr->idunion.db.boxes;
	_REAL	*boxespa = dptr->idunion.db.boxespa;
	_REAL	*boxespe = dptr->idunion.db.boxespe;

	_REAL   *crdi = serptr->idunion.ser.set->crd;
	_REAL   *crdj = serptr->idunion.ser.set->crd;
        int shift= serptr->idunion.ser.shift;
	_REAL   *dir = serptr->idunion.ser.set->dir;
        _REAL   *boxc = serptr->idunion.ser.set->boxc;

	int	k,im,jm,ismall,iz,i,j;
        _REAL boxh[3],boxm1[3],dr[3],drp[3],mdrp,pspa,
              pspe;
	_REAL	d,distancebox();
        void images(_REAL *dist, idstruct *serptr);
        _REAL determinant(_REAL **a),**aa=serptr->idunion.ser.set->cell_to_cart;
        _REAL *crdj0,weight;

        int *excltemp = prm->excltemp, *nexcl = prm->Nexcl,
         *iexcl = prm->Iexcl, ntotexcl = prm->Natexcl,idum,jx,nx;

   if(serptr->idunion.ser.iweight)
   {
    weight=serptr->idunion.ser.w->idunion.serw.w;
   } else {
    weight=1.;
   }

        for (k=0; k<3; k++) {
         boxm1[k] = 1./boxc[k];
         boxh[k] = 0.5*boxc[k];
         boxav[k] += weight*boxc[k];
        }
/* sort the box sides */
        for (k=0; k<3; k++) {
         if (boxh[k]<=small) {
          small=boxh[k];
          ismall=k;
         }
        }
        dptr->idunion.db.smallside = small;

        dptr->idunion.db.cases++;
        dptr->idunion.db.casesw+=weight;
        dptr->idunion.db.rho +=
	  weight*boxm1[0]*boxm1[1]*boxm1[2];
//        dptr->idunion.db.rho += weight/determinant(aa);
/*
         printf("dir: %e %e %e\n",dir[0],dir[1],dir[2]);
         printf("box: %e %e %e\n",boxc[0],boxc[1],boxc[2]);
         printf("small: %e\n",small);
*/
        if(dptr->idunion.db.type==1) {
/* one atom per mol distribution */
         crdi = listi;
         // crdj0 = listj+jump;
         crdj0 = serptr->idunion.ser.set->crd;
         // for(im=0,i=0; im<nmol-1; im++,crdi+=jump,i+=jump) 
         for(im=0,i=0; im<nmol; im++,crdi+=jump,i+=jump) 
         {
            crdj = crdj0;
	 // for (jm=im+1,j=i+jump; jm<nmol; jm++,crdj+=jump,j+=jump) 
	    for (j=0; j<nsit; j++,crdj+=3) 
            {
             dr[0]=crdi[0]-crdj[0];
             dr[1]=crdi[1]-crdj[1];
             dr[2]=crdi[2]-crdj[2];
#ifdef PBC
/*
             dr[0] -= boxc[0]*rint(dr[0]*boxm1[0]);
             dr[1] -= boxc[1]*rint(dr[1]*boxm1[1]);
             dr[2] -= boxc[2]*rint(dr[2]*boxm1[2]);
*/
             images(dr,serptr);
#endif
             pspa=Nabs(dr[0]*dir[0]+dr[1]*dir[1]+dr[2]*dir[2]);
             drp[0]=dr[0]-pspa*dir[0];
             drp[1]=dr[1]-pspa*dir[1];
             drp[2]=dr[2]-pspa*dir[2];
             mdrp=drp[0]*drp[0]+drp[1]*drp[1]+drp[2]*drp[2];
             pspe=l_sqrt(mdrp);
     	     if(pspa<= maxz && pspa>=minz) {
              boxespa[(int)((pspa-minz)/scalez)]+=weight;
             }
             if(pspe<= max && pspe>=min)
               if(pspa<=scalez)
                 boxespe[(int)((pspe-min)/scale)]+=2.;
             d=l_sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
             /* d = distancebox(crdi,crdj,boxc,boxm1); */
      	     if(d <= max && d >=min) 
//	      boxes[(int) ((d-min) / scale)]+=2.*weight;
	      boxes[(int) ((d-min) / scale)]+=
	        prm->Charge[j]/0.3727337818133816e+02;
            }
//            crdj0+=jump;
	 }
        } else if(dptr->idunion.db.type==0) {
/* all atoms distribution */
         for(im=0,i=0; im<nsit-1; im++,crdi+=3,i++) {

          for (k=0;k<prm->Nsit;k++) *prm->excltemp++ = -1;
          prm->excltemp -= prm->Nsit;
          nx = prm->Nexcl[i];
          for(jx=0;jx<nx;jx++) {
           idum = *iexcl++;
           if (idum ==0) idum = 1;
           excltemp[idum-1] = i;
          }

            crdj = crdi+3;
	    for (jm=im+1,j=jm; jm<nsit; jm++,crdj+=3,j++) {
             if(prm->excltemp[j] != i)
             {

#ifdef PBC
               d = distancebox(crdi,crdj,boxc,boxm1);
#else
               d = distance(crdi,crdj,boxc,boxm1);
#endif
      	       if(d <= max && d >=min) boxes[(int) ((d-min) / scale)]+=2.*weight;
             }
            }
	 }
         crdi-=3*(nsit-1);
         crdj-=3*nsit;
        } else {
         printf("gr: type not found\n");
         exit(1);
        }

	/*
	 *  if this is a window distribution & window done, write it
	 */

	if (dptr->idunion.db.win) {
		if (++dptr->idunion.db.winct == dptr->idunion.db.win) {

			int	i, nboxes = dptr->idunion.db.nboxes;

			/*
			 *  write out distribution curve and reinit
			 */

			dptr->idunion.db.winct = 0;
			writegr(dptr, 0);
			for (i=0; i<nboxes; i++) {
				dptr->idunion.db.boxes[i] = 0.;
				dptr->idunion.db.boxespa[i] = 0.;
				dptr->idunion.db.boxespe[i] = 0.;
			}
			dptr->idunion.db.winct = 0;

		}
	}
}

void
gr_grid(idstruct *dptr,idstruct *serptr)
{
        extern parmstruct *prm;
        int nmol=prm->Nmol;
        int nsit=prm->Nsit;
	_REAL	small = dptr->idunion.db.smallside;
	_REAL	*listi = dptr->idunion.db.listi;
	_REAL	*listj = dptr->idunion.db.listj;
	int	jump = dptr->idunion.db.jump;
	_REAL   *boxav = dptr->idunion.db.boxav;
	_REAL	*boxes = dptr->idunion.db.boxes;
	_REAL	*boxespa = dptr->idunion.db.boxespa;
	_REAL	*boxespe = dptr->idunion.db.boxespe;

	_REAL   *crdi = serptr->idunion.ser.set->crd;
	_REAL   *crdj = serptr->idunion.ser.set->crd;
        int shift= serptr->idunion.ser.shift;
	_REAL   *dir = serptr->idunion.ser.set->dir;
        _REAL   *boxc = serptr->idunion.ser.set->boxc;

	int	k,im,jm,ismall,iz,i,j;
        _REAL boxh[3],boxm1[3];
	_REAL	d,distancebox();
        void images(_REAL *dist, idstruct *serptr);
        _REAL determinant(_REAL **a),**aa=serptr->idunion.ser.set->cell_to_cart;
        _REAL *crdj0,weight;

        int *excltemp = prm->excltemp, *nexcl = prm->Nexcl,
         *iexcl = prm->Iexcl, ntotexcl = prm->Natexcl,idum,jx,nx;

   if(serptr->idunion.ser.iweight)
   {
    weight=serptr->idunion.ser.w->idunion.serw.w;
   } else {
    weight=1.;
   }

        for (k=0; k<3; k++) {
         boxm1[k] = 1./boxc[k];
         boxh[k] = 0.5*boxc[k];
         boxav[k] += weight*boxc[k];
        }
/* sort the box sides */
        for (k=0; k<3; k++) {
         if (boxh[k]<=small) {
          small=boxh[k];
          ismall=k;
         }
        }
        dptr->idunion.db.smallside = small;

        dptr->idunion.db.cases++;
        dptr->idunion.db.casesw+=weight;
        dptr->idunion.db.rho +=
	  boxm1[0]*boxm1[1]*boxm1[2];
//      dptr->idunion.db.rho += weight/determinant(aa);
/*
         printf("dir: %e %e %e\n",dir[0],dir[1],dir[2]);
         printf("box: %e %e %e\n",boxc[0],boxc[1],boxc[2]);
         printf("small: %e\n",small);
*/
        printf("weight: %e\n",weight);

        if(dptr->idunion.db.type==0) {
/* all atoms distribution */
         for(im=0,i=0; im<nsit-1; im++,crdi+=3,i++) {

          for (k=0;k<prm->Nsit;k++) *prm->excltemp++ = -1;
          prm->excltemp -= prm->Nsit;
          nx = prm->Nexcl[i];
          for(jx=0;jx<nx;jx++) {
           idum = *iexcl++;
           if (idum ==0) idum = 1;
           excltemp[idum-1] = i;
          }

            crdj = crdi+3;
	    for (jm=im+1,j=jm; jm<nsit; jm++,crdj+=3,j++) {

             if(prm->excltemp[j] != i)
             {
#ifdef PBC
              d = distancebox(crdi,crdj,boxc,boxm1);
#else
              d = distance(crdi,crdj,boxc,boxm1);
#endif
              for(k=0;k<dptr->idunion.db.nboxes;k++)
              {
      	       if(d >= dptr->idunion.db.rmin[k] && 
                  d  < dptr->idunion.db.rmax[k])
               {
                boxes[k]+=2.*weight;
                boxespe[k]+=2.*weight;
               }
              }
             }
            }
	 }
         crdi-=3*(nsit-1);
         crdj-=3*nsit;
        }

        for(k=0;k<dptr->idunion.db.nboxes;k++)
        {
         boxes[k]/=(_REAL)prm->Nsit;
         fwrite(&boxes[k],1,sizeof(_REAL),dptr->idunion.db.file1);
         boxes[k]=0.;
        }

}


/************************************************************************
							WRITEgr()
*************************************************************************/

/*  normalize and print g(r) */

void
writegr(idstruct *dptr, int done)
{
	int	i, nboxes = dptr->idunion.db.nboxes;
	int	nboxesz = dptr->idunion.db.nboxesz;
	_REAL	*boxes = dptr->idunion.db.boxes;
	_REAL	*boxespa = dptr->idunion.db.boxespa;
	_REAL	*boxespe = dptr->idunion.db.boxespe;
        extern parmstruct *prm;
	int	npart;
	FILE 	*file1 = dptr->idunion.db.file1;
	FILE 	*file2 = dptr->idunion.db.file2;
	_REAL	scale = dptr->idunion.db.scale;
	_REAL	scalez = dptr->idunion.db.scalez;
	_REAL	min = dptr->idunion.db.min;
	_REAL	minz = dptr->idunion.db.minz;
	_REAL	rho = dptr->idunion.db.rho,rhoc;
	_REAL   *boxav = dptr->idunion.db.boxav;
	_REAL   small = dptr->idunion.db.smallside;
        _REAL cases,ri,ri3,fact,fac,facpa,facpe,gr,grpe,grpa,sidexy,
              fact_om;

        if (dptr->idunion.db.type==0) {
         npart=prm->Nsit;
        } else if (dptr->idunion.db.type==1) {
         npart=prm->Nmol;
        } else {
         printf("writedb: type not found\n");
         exit(1);
        }

	if (dptr->idunion.db.winct) {
		/* incomplete window */
		cases = (_REAL)(dptr->idunion.db.winct * npart);
		printf( "Warning: last distrib %s has only %d cases, %d %s",
		    dptr->id_name, dptr->idunion.db.winct, dptr->idunion.db.win,
		    "specified for window\n");
	} else if (dptr->idunion.db.win) {
		/* complete window */
		cases = (_REAL)(dptr->idunion.db.win * npart);
	} else
		cases = dptr->idunion.db.casesw;

        boxav[0]/=cases;
        boxav[1]/=cases;
        boxav[2]/=cases;

/* constants */
        rho = rho * (_REAL)npart / cases;
        rhoc = rho*boxav[2]/40.;
        // fac = 4.*_PI*rho*cases*(_REAL)npart/3.;
        fac = 4.*_PI*rho*cases*(_REAL)(prm->Nsit-1)/3.;
        fact_om = cases*(_REAL)npart;
        sidexy=boxav[0]*boxav[1];
        facpa = sidexy*rho*cases*(_REAL)npart;
        printf("npart=%d\n",npart);
        small=boxav[2];
        facpe = small*_PI*rho*cases*(_REAL)npart;
        facpe *= 0.5;
/*
        facpe = 2.*scalez*_PI*rho*cases*(_REAL)npart;
        facpe = facpe/((_REAL)nboxesz/(_REAL)nboxes);
*/
        printf("cases,rho,smallside %e %e %e\n",cases,rho,small);

	/*
	 *  unless finally done with a window already printed, print curve
	 */

	if (!(done && dptr->idunion.db.win && !dptr->idunion.db.winct)) {

         if(dptr->idunion.db.grid==0)
         {
			for (i=0; i<nboxes; i++) {
                          ri = min + (_REAL)i*scale;
                          ri3 = pow((ri + scale),3.) - pow(ri,3.);
                          fact = fac*ri3;
                           gr = boxes[i] / fact;
                          /* gr = boxes[i] / fact_om; */
                          ri3 = pow((ri + scale),2.) - pow(ri,2.);
                          fact = facpe*ri3;
                          grpe = boxespe[i] / fact;
 			  fprintf(file1, "%e %e %e\n", (ri+0.5*scale)*0.1,gr,
                                       grpe); 
			}

			for (i=0; i<nboxesz; i++) {
                          ri = minz + (_REAL)i*scalez;
                          fact = facpa*scalez;
                          grpa = boxespa[i] / fact;
 			  fprintf(file2, "%e %e\n", (ri+0.5*scalez)*0.1,
                                       grpa); 
			}
         }
         else
         {
			for (i=0; i<nboxes; i++) {
                          scale = dptr->idunion.db.rmax[i]
                                 -dptr->idunion.db.rmin[i];
                          ri = dptr->idunion.db.rmin[i];
                          ri3 = pow((ri + scale),3.) - pow(ri,3.);
                          fact = fac*ri3;
                          /* gr = boxespe[i] / fact; */
                          gr = boxespe[i] / fact_om;
 			  fprintf(file2, "%e %e\n", (ri+0.5*scale)*0.1,gr);
			}
         }
	}
	if (done) {
		fclose(file1);
		fclose(file2);
        }
}
