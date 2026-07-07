#include "head.h"

// extern char     prmtop[],tok[];
char     prmtop[];
extern idstruct *startid();
extern FILE     *genopen();
extern char     *get();

parmstruct      *prm = NULL;
parmstruct      *defprm = NULL;

void
rdprm()
{
   FILE *fileptr;
   char *name,*Temp;
   int i,j,k,ic,natomtot,ntype2;
#ifdef DOUBLE
   static char scanform[] = "%le";
#else
   static char scanform[] = "%e";
#endif
   _REAL Xdum[2];
   char Cdum[1];
   void preadln();
   idstruct        *prmptr;
   char            id[TOKLEN];
   int swap=0;

   prmptr = startid(PARM, id);

   gettok();
   if (!strcmp(tok, "SWAP")) {
     swap=1;
     gettok();
   }

   /* allow room for .Z */
   prmptr->idunion.parm.name = get(strlen(tok) + 3);
   strcpy(prmptr->idunion.parm.name, tok);
   if((fileptr = genopen(tok,"PARM")) == NULL) exit(1);

   prmptr->idunion.parm.prm = (parmstruct *) get(sizeof(parmstruct));
   prm = prmptr->idunion.parm.prm;
   if(defprm==NULL) {
    defprm=prm;
    strncpy(prmtop,tok,TOKLEN-1);
   }

   gettok();
   if (strcmp(tok,";")) inerr("; expected","");

/* read headers */
        fscanf(fileptr, "%d%d%d%d%d%d%d%d%d%d",
                &prm->Ntypes, &prm->Natom, &prm->Nmol,
                &prm->Nbond, &prm->Nangle, &prm->Ndihe, &prm->Ndihei,
                &prm->Natexcl, &prm->Ncharge, &prm->Nres);
        while(getc(fileptr) != 10);
/*
        printf("%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
                prm->Ntypes, prm->Nmol, prm->Natom,
                prm->Nbond, prm->Nangle, prm->Ndihe, prm->Ndihei,
                prm->Natexcl, prm->Ncharge, prm->Nres);
*/

        fscanf(fileptr, "%d%d%d%d%d%d",
                &prm->Nsolu, &prm->Natom_solu,&prm->Nions,&prm->Natom_ions,
                &prm->Nsolv,&prm->Natom_solv);
        while(getc(fileptr) != 10);
        prm->Nsolu_tot=prm->Nsolu*prm->Natom_solu;
        prm->Nions_tot=prm->Nions*prm->Natom_ions;
        prm->Nsolv_tot=prm->Nsolv*prm->Natom_solv;
/*
        printf("%6d%6d%6d%6d%6d%6d\n",
                prm->Nsolu, prm->Natom_solu,prm->Nions,prm->Natom_ions,
                prm->Nsolv,prm->Natom_solv);
*/

        fscanf(fileptr, "%d%d%d%d",
                &prm->Nbondp, &prm->Nanglep, &prm->Ndihep, &prm->Ndiheip);
        while(getc(fileptr) != 10);
/*
        printf("%6d%6d%6d%6d\n",
                prm->Nbondp, prm->Nanglep, prm->Ndihep, prm->Ndiheip);
*/

        fscanf(fileptr, scanform,
                Xdum,Xdum);
        while(getc(fileptr) != 10);

/* set memory */

   natomtot = prm->Natom * prm->Nmol;
   prm->Nsit = natomtot;
   if(prm->Nsit!=(prm->Nsolu_tot+prm->Nions_tot+prm->Nsolv_tot))
   {
    printf("rdprm: Nsit and Natom_tot are different.\n");
    exit(1);
   }
   ntype2 = prm->Ntypes * prm->Ntypes;
   prm->atomtype = (int *) get(sizeof(int)*natomtot);
   prm->Masses = (_REAL *) get(sizeof(_REAL)*(prm->Natom));
   prm->Charge = (_REAL *) get(sizeof(_REAL)*(prm->Ncharge));
   prm->Ibond = (int *) get(sizeof(int)*3*(prm->Nbond));
   prm->Kbond = (_REAL *) get(sizeof(_REAL)*2*(prm->Nbondp));
   prm->Iangle = (int *) get(sizeof(int)*4*(prm->Nangle));
   prm->Kangle = (_REAL *) get(sizeof(_REAL)*2*(prm->Nanglep));
   prm->Idihe = (int *) get(sizeof(int)*5*(prm->Ndihe));
   prm->Idihei = (int *) get(sizeof(int)*5*(prm->Ndihei));
   prm->Kdihe = (_REAL *) get(sizeof(_REAL)*MAXDIHK*(prm->Ndihep));
   prm->Kdihei = (_REAL *) get(sizeof(_REAL)*MAXDIHKI*(prm->Ndihep));
   prm->Nexcl = (int *) get(sizeof(int)*(prm->Natom));
   prm->excltemp = (int *) get(sizeof(int)*(prm->Natom));
    for (i=0;i<prm->Natom;i++) prm->excltemp[i] = -1;
   prm->Iexcl = (int *) get(sizeof(int)*(prm->Natexcl));
   prm->ljAB = (_REAL *) get(sizeof(_REAL)*2*ntype2);
/* labels */
        prm->Swapat = (int *) get(sizeof(int)*(prm->Natom));
	prm->AtomNames = (char *) get(ATOMTOKLEN*sizeof(char)*prm->Natom+81);
	prm->ResNames = (char *) get(ATOMTOKLEN*sizeof(char)*prm->Nres+81);
	prm->ResNamesa = (char *) get(ATOMTOKLEN*sizeof(char)*prm->Natom+81);
	prm->Ipres = (int *) get(sizeof(int)*(prm->Nres+1));
	prm->ResNum = (int *) get(sizeof(int)*(prm->Natom));

/* read parameters */
        for (i=0; i<prm->Natom; i++)
                fscanf(fileptr, "%d", &prm->atomtype[i]);
        for (i=0; i<prm->Natom; i++) prm->atomtype[i] -= 1;
        while(getc(fileptr) != 10);
        for (i=0; i<prm->Natom; i++)
                fscanf(fileptr, scanform, &prm->Masses[i]);

/*
    for (i=0; i<prm->Natom; i++) {
     printf("%e ", prm->Masses[i]);
     if( (i+1)%12 == 0) printf("\n");
    }
    printf("\n");
*/

        while(getc(fileptr) != 10);
        for (i=0; i<prm->Ncharge; i++)
                fscanf(fileptr, scanform, &prm->Charge[i]);
        while(getc(fileptr) != 10);
        for (i=0; i<3*(prm->Nbond); i++)
                fscanf(fileptr, "%d", &prm->Ibond[i]);
        while(getc(fileptr) != 10);
        for (i=0,ic=0; i<prm->Nbond; i++,ic+=3) {
           prm->Ibond[ic]   -= 1;
           prm->Ibond[ic+1] -= 1;
           prm->Ibond[ic+2] -= 1;
        }
/*
        for (i=0,ic=0; i<prm->Nbond; i++,ic+=3) 
         printf("%d %d %d\n",
           prm->Ibond[ic],
           prm->Ibond[ic+1],
           prm->Ibond[ic+2]);
*/

        for (i=0; i<4*(prm->Nangle); i++)
                fscanf(fileptr, "%d", &prm->Iangle[i]);
        for (i=0,ic=0; i<prm->Nangle; i++,ic+=4) {
           prm->Iangle[ic]   -= 1;
           prm->Iangle[ic+1] -= 1;
           prm->Iangle[ic+2] -= 1;
           prm->Iangle[ic+3] -= 1;
        }
        while(getc(fileptr) != 10);
        for (i=0; i<5*(prm->Ndihe); i++)
                fscanf(fileptr, "%d", &prm->Idihe[i]);
        for (i=0,ic=0; i<prm->Ndihe; i++,ic+=5) {
           prm->Idihe[ic]   -= 1;
           prm->Idihe[ic+1] -= 1;
           prm->Idihe[ic+2] -= 1;
           prm->Idihe[ic+3] -= 1;
           prm->Idihe[ic+4] -= 1;
        }
        while(getc(fileptr) != 10);
        for (i=0; i<5*(prm->Ndihei); i++)
                fscanf(fileptr, "%d", &prm->Idihei[i]);
        for (i=0,ic=0; i<prm->Ndihei; i++,ic+=5) {
           prm->Idihei[ic]   -= 1;
           prm->Idihei[ic+1] -= 1;
           prm->Idihei[ic+2] -= 1;
           prm->Idihei[ic+3] -= 1;
           prm->Idihei[ic+4] -= 1;
        }
        while(getc(fileptr) != 10);
        for (i=0; i<2*(prm->Nbondp); i++)
                fscanf(fileptr, scanform, &prm->Kbond[i]);
        while(getc(fileptr) != 10);
        for (i=0; i<2*(prm->Nanglep); i++)
                fscanf(fileptr, scanform, &prm->Kangle[i]);
        while(getc(fileptr) != 10);
        for (i=0; i<MAXDIHK*(prm->Ndihep); i++)
                fscanf(fileptr, scanform, &prm->Kdihe[i]);
        while(getc(fileptr) != 10);
        for (i=0; i<MAXDIHKI*(prm->Ndiheip); i++)
                fscanf(fileptr, scanform, &prm->Kdihei[i]);
        while(getc(fileptr) != 10);
        for (i=0; i<prm->Natom; i++)
                fscanf(fileptr, "%d", &prm->Nexcl[i]);
        while(getc(fileptr) != 10);
        for (i=0; i<prm->Natexcl; i++)
                fscanf(fileptr, "%d", &prm->Iexcl[i]);
        while(getc(fileptr) != 10);
/*
        for (i=0; i<prm->Natexcl; i++) {
           prm->Iexcl[i]   -= 1;
        }
*/
        for (i=0; i<2*(prm->Ntypes)*(prm->Ntypes); i++)
                fscanf(fileptr, scanform, prm->ljAB++);
        prm->ljAB -= 2*(prm->Ntypes)*(prm->Ntypes);
        while(getc(fileptr) != 10);
                fscanf(fileptr, scanform, &prm->sig0r);
                fscanf(fileptr, scanform, &prm->eps0r);
                fscanf(fileptr, scanform, &prm->mass0r);
        while(getc(fileptr) != 10);
        printf("sig0r = %e ", prm->sig0r);
/*
        for (i=0; i<(prm->Ntypes)*(prm->Ntypes); i++) {
                printf("%d %e ", i,*prm->ljAB++);
                printf("%e\n", *prm->ljAB++);
        }
        prm->ljAB -= 2*(prm->Ntypes)*(prm->Ntypes);
        fflush(stdout);
*/

/* READ ATOM NAMES */
        j=1;
	for (i=0; i<(j/20 + (j%20 ? 1 : 0)); i++)
		preadln(fileptr, "", Cdum);
        printf("Nres = %d\n", prm->Nres); fflush(stdout);
        for (i=0; i<prm->Nres+1; i++)
        {
                fscanf(fileptr, "%d", &prm->Ipres[i]);
    
        }
        while(getc(fileptr) != 10);
/*
        for (i=0; i<prm->Nres+1; i++)
                printf("%d\n", prm->Ipres[i]);
*/
/*
        for (i=0; i<prm->Nres+1; i++) {
           prm->Ipres[i]   -= 1;
        }
*/

	for (i=0; i<(prm->Nres/20 + (prm->Nres%20 ? 1 : 0)); i++)
		preadln(fileptr, "", &prm->ResNames[i*80]);
/*
	for (i=0; i<(prm->Nres/20 + (prm->Nres%20 ? 1 : 0)); i++) 
		printf("%.80s\n", &prm->ResNames[i*80]);
*/
/*
        for (i=0; i<prm->Natom; i++)
                fscanf(fileptr, "%d", &prm->Swapat[i]);
        while(getc(fileptr) != 10);
        for (i=0; i<prm->Natom; i++) prm->Swapat[i] -= 1;
*/

	for (i=0; i<(prm->Natom/20 + (prm->Natom%20 ? 1 : 0)); i++)
		preadln(fileptr, "", &prm->AtomNames[i*80]);
/*
 { char atomname[ATOMTOKLEN];
   for(i=0;i<prm->Natom;i++) {
    strncpy(atomname,&prm->AtomNames[i*ATOMTOKLEN],ATOMTOKLEN);
    printf("atomname(%d)= %s\n",i,atomname);
   }
 }
*/

/*
        printf("Atom names old\n");
	for (i=0; i<(prm->Natom/20 + (prm->Natom%20 ? 1 : 0)); i++) 
		printf("%.80s\n", &prm->AtomNames[i*80]);
*/
/* set AtomNames according to the site numbering */
        if (swap==1) {
	 Temp = (char *) get(ATOMTOKLEN*sizeof(char)*prm->Natom+81);
         for (i=0; i<prm->Natom; i++) {
		j = prm->Swapat[i];
		for (k=0;k<ATOMTOKLEN;k++)
		 Temp[j*ATOMTOKLEN+k] = prm->AtomNames[i*ATOMTOKLEN+k];
 	 }
         for (i=0; i<prm->Natom; i++) {
		for (k=0; k<ATOMTOKLEN; k++)
		 prm->AtomNames[i*ATOMTOKLEN+k] = Temp[i*ATOMTOKLEN+k];
	 }
         free((char *) Temp);
/*
         printf("Atom names new\n");
	 for (i=0; i<(prm->Natom/20 + (prm->Natom%20 ? 1 : 0)); i++) 
		printf("%.80s\n", &prm->AtomNames[i*80]);
*/
        }

/* end reading */
/* set residue numbuer for each atom */
       for(i=0,k=0;i<prm->Nres;i++) {
        for(j=prm->Ipres[i];j<prm->Ipres[i+1];j++,k++) {
         prm->ResNum[k]=i+1;
         strncpy(&prm->ResNamesa[k*ATOMTOKLEN],
            &prm->ResNames[i*ATOMTOKLEN],ATOMTOKLEN);
        }
       }
       prm->Nmov=prm->Nsit;

}
 
void
preadln(FILE *file, char *name, char *string)
{
	int 	i, j;

	for (i=0; i<81; i++) {
		if ((j = getc(file)) == EOF) {
			printf("Error: unexpected EOF in %s\n", name);
			exit(1);
		}
		string[i] = (char) j;
		if (string[i] == '\n') {
			break;
		}
	}
	if (i == 80  &&  string[i] != '\n') {
		printf("Error: line too long in %s:\n%.80s", name, string);
		exit(1);
	}
}

int
firstWAT(parmstruct *prm)
{
        char    *restr = prm->ResNames;
        char    *lastres = prm->ResNames + prm->Nres * 4 + 1;
        int     res;

        /*
         *  find 1st water residue
         */

        for (res=0; restr<lastres; restr+=4, res++) {
                if (!strncmp(restr, "WAT ", 4)) {
                  printf("first water: res = %d, atom = %d (%.4s)\n",
                                res+1, prm->Ipres[res]-1,
                                &prm->AtomNames[4 * (prm->Ipres[res]-1)]);
                  fflush(stdout);
                  return(prm->Ipres[res]-1);
                }
        }
        return(0);
}
