#include "head.h"

// extern char		tok[];
// char		tok[];

/**** extern ****/
extern idstruct		*startid(),*getstruct();

// FILE *ftmp;

/************************************************************************/

void
initorient()
{
	idstruct	*dptr;
	char		id[TOKLEN];
	int		exists;

	dptr = startid(ORIENT, id);

	dptr->idunion.orient.type = 0;
	dptr->idunion.orient.nboxesa = 0;
	dptr->idunion.orient.nboxesb = 0;
	dptr->idunion.orient.cases = 0;
	dptr->idunion.orient.scalea = 0.;
	dptr->idunion.orient.scaleb = 0.;
	dptr->idunion.orient.mina = 0.;
	dptr->idunion.orient.minb = 1.;
	dptr->idunion.orient.maxa = TWOPI;
	dptr->idunion.orient.maxb = -1.;

	/* 'file_name' */

	gettok();
	if (iskeywd())
		inerr("ORIENT: expected file name","");
	strcpy(dptr->idunion.orient.fname, tok);

	/*
	 *  open the file(s)
	 */

	if (fileprob(dptr->idunion.orient.fname, 0, &exists, "ORIENT"))
		exit(1);
	if ((dptr->idunion.orient.file = 
				  fopen(dptr->idunion.orient.fname, "w")) == NULL){
		perror(dptr->idunion.orient.fname);
		exit(1);
	}
	gettok();

/*
printf("dptr %d id %s\n", dptr, dptr->id_name);
*/
}

/************************************************************************
							setuporient()
*************************************************************************/

void
setuporient()
{
	idstruct	*dptr,*grp1,*grp2;
        extern parmstruct *prm;
	extern FILE *ftmp;
	int i,j,k,grptype;

/* get id */
        gettok();
        if ((dptr = getstruct(tok)) != NULL) {
                if (dptr->id_type != ORIENT)
                        inerr("setuporient: expected id or keyword\n","");
        } else {
                printf("setuporient: id %s not found in idlist\n", tok);
        }
        gettok();
        if ((grp1 = getstruct(tok)) == NULL) 
                inerr("setuporient: no such id\n","");
        if (grp1->id_type != GROUP)
                inerr("setuporient: group1 id expected\n","");
        grptype = grouptype(grp1);

        /* gettok(); */
        if ((grp2 = getstruct(tok)) == NULL) 
                inerr("setuporient: no such id\n","");
        if (grp2->id_type != GROUP)
                inerr("setuporient: group2 id expected\n","");
        grptype = grouptype(grp2);
        dptr->idunion.orient.grp1=grp1;
        dptr->idunion.orient.grp2=grp2;

        dptr->idunion.orient.type = 0;

               for(;;) {
	        /* gettok(); */
        	if (tok[0] == ';') {
                 break;
                } else if (!strcmp(tok, "X")) {
                 dptr->idunion.orient.type = 1;
                } else if (!strcmp(tok, "NINTA")) {
                        gettok();
                        if (!tokint(&dptr->idunion.orient.nboxesa))
                         inerr("ORIENT: expected integer # of boxes", "");
                } else if (!strcmp(tok, "MINA")) {
                        gettok();
                        if (!tokREAL(&dptr->idunion.orient.mina))
                         inerr("ORIENT: expected mina", "");
                } else if (!strcmp(tok, "MAXA")) {
                        gettok();
                        if (!tokREAL(&dptr->idunion.orient.maxa))
                         inerr("ORIENT: expected maxa", "");
                } else if (!strcmp(tok, "NINTB")) {
                        gettok();
                        if (!tokint(&dptr->idunion.orient.nboxesb))
                         inerr("ORIENT: expected integer # of boxes", "");
                } else if (!strcmp(tok, "MINB")) {
                        gettok();
                        if (!tokREAL(&dptr->idunion.orient.minb))
                         inerr("ORIENT: expected minb", "");
                } else if (!strcmp(tok, "MAXB")) {
                        gettok();
                        if (!tokREAL(&dptr->idunion.orient.maxb))
                         inerr("ORIENT: expected maxb", "");
		} else inerr("setuporient: option not found", "");
                gettok();
	       }

/* allocate and initialize */
         dptr->idunion.orient.mina *= TORAD;
         dptr->idunion.orient.maxa *= TORAD;
         dptr->idunion.orient.scalea =
          ( dptr->idunion.orient.maxa - dptr->idunion.orient.mina)
                                / dptr->idunion.orient.nboxesa;
/*
         dptr->idunion.orient.minb *= TORAD;
         dptr->idunion.orient.maxb *= TORAD;
*/
         dptr->idunion.orient.scaleb =
          ( dptr->idunion.orient.maxb - dptr->idunion.orient.minb)
                                / dptr->idunion.orient.nboxesb;
         dptr->idunion.orient.nboxes = dptr->idunion.orient.nboxesa*
                                      dptr->idunion.orient.nboxesb;

		dptr->idunion.orient.boxes = 
                 (int *) get(dptr->idunion.orient.nboxes*sizeof(int));
		for (i=0; i<dptr->idunion.orient.nboxes; i++) {
		 dptr->idunion.orient.boxes[i] = 0;
                }
        ftmp = fopen("bonds.time","w");
}


/************************************************************************/

void
distriborient(idstruct *objptr)
{
   idstruct *grp1=objptr->idunion.orient.grp1,*grp2=objptr->idunion.orient.grp2;
   void qeuler();
   extern FILE *ftmp;

   int     nboxesa = objptr->idunion.orient.nboxesa;
   int     nboxesb = objptr->idunion.orient.nboxesb;
   _REAL   scalea = objptr->idunion.orient.scalea, 
           mina = objptr->idunion.orient.mina,
           scaleb = objptr->idunion.orient.scaleb, 
           minb = objptr->idunion.orient.minb;
   int     *boxes = objptr->idunion.orient.boxes;

   int k,ial,ibe,ibox;
   _REAL q1[4],q2t[4],dq[4],q1t[4],q2[4],cbh,cb,ca,sa,cg,sg,beta,alpha;
   _REAL ux0[]={0.,1.,0.,0.},uz0[]={0.,0.,0.,1.},utmp[4],qtmp[4],phase(),dqt[4];
   _REAL r1tmp[4],r2tmp[4],scp;

/* alpha spans (-pi,pi) if atan2 is used */
/* alpha spans (0,2*pi) if phase is used */

    ++objptr->idunion.orient.cases;
    dogrp(grp1);
    dogrp(grp2);

    q1[0]= grp1->idunion.group.q[0];
    q1[1]= grp1->idunion.group.q[1];
    q1[2]= grp1->idunion.group.q[2];
    q1[3]= grp1->idunion.group.q[3];
    q1t[0]= q1[0];
    q1t[1]=-q1[1];
    q1t[2]=-q1[2];
    q1t[3]=-q1[3];
    q2t[0]= grp2->idunion.group.q[0];
    q2t[1]=-grp2->idunion.group.q[1];
    q2t[2]=-grp2->idunion.group.q[2];
    q2t[3]=-grp2->idunion.group.q[3];
    q2[0]= q2t[0];
    q2[1]=-q2t[1];
    q2[2]=-q2t[2];
    q2[3]=-q2t[3];
    /* proquat(q1,q2t,dq); */
    proquat(q2,q1t,dq);
    /* proquat(q1t,q2,dq); */
    /* proquat(q2t,q1,dq); */

    dqt[0]= dq[0];
    dqt[1]=-dq[1];
    dqt[2]=-dq[2];
    dqt[3]=-dq[3];

    if(objptr->idunion.orient.type==0) {
     proquat(dqt,uz0,qtmp);
     proquat(qtmp,dq,utmp);
    } else {
     proquat(dqt,ux0,qtmp);
     proquat(qtmp,dq,utmp);
    }

/*
     proquat(q1t,uz0,qtmp);
     proquat(qtmp,q1,r1tmp);
     proquat(q2t,uz0,qtmp);
     proquat(qtmp,q2,r2tmp);
     scp=r1tmp[1]*r2tmp[1]+r1tmp[2]*r2tmp[2]+r1tmp[3]*r2tmp[3];
     printf("bb= %f\n",TODEG*acos(scp));
*/

    cb=utmp[3];
    alpha=phase(utmp[1],utmp[2]);
/*
    qeuler(dq,&cbh,&ca,&sa,&cg,&sg);
    cb=2.*cbh*cbh-1.;
    alpha=atan2(sa,ca);
*/
 printf("beta,alpha= %f %f\n",TODEG*acos(cb),TODEG*alpha);
 fprintf(ftmp,"%e %e %e\n",utmp[1],utmp[2],utmp[3]);

    ial = (int)((alpha-mina)/scalea);
    ibe = (int)((cb-minb)/scaleb);
    ibox = nboxesa*ibe+ial;
    boxes[ibox]++;

}

/************************************************************************/

void
writeorient(idstruct *dptr)
{
	int	i,j,ibox;
        int     nboxesa = dptr->idunion.orient.nboxesa,
                nboxesb = dptr->idunion.orient.nboxesb;
	int	*boxes = dptr->idunion.orient.boxes;
	int	states = dptr->idunion.orient.cases;
	FILE 	*file = dptr->idunion.orient.file;
	_REAL	scalea = dptr->idunion.orient.scalea,
                mina=dptr->idunion.orient.mina,
	        scaleb = dptr->idunion.orient.scaleb,
                minb=dptr->idunion.orient.minb;

	_REAL	alpha,beta,distrib,fac;

	fac = 1./ ((_REAL)states);

        for (i=0,ibox=0; i<nboxesb; i++) {
         beta = minb+(_REAL)i*scaleb;
         for (j=0; j<nboxesa; j++,ibox++) {
           alpha = TODEG*(_REAL)j*scalea;
           distrib = (_REAL)boxes[ibox]*fac;
           fprintf(file, "%e %e %e\n", beta,alpha,distrib);
         }
//         fprintf(file, "\n");
        }
        fclose(file);
}
