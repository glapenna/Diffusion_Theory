#include "head.h"

// extern char     tok[];
// char     tok[];
extern idstruct *startid();
extern FILE     *genopen();
extern char     *get();

void
rdprm_his2()
{
   extern parmstruct *prm;
   FILE *fileptr;
   int i,j,k,ir,bytes,bytest,nitems,itest;
   int Idum[1],sizeAtName;
   char Cdum[1];
   double Ddum[1];
   float Fdum[1];
   int swap;
   idstruct        *prmptr;
   char            id[TOKLEN];
   char atomname[ATOMTOKLEN];
   char atomnamer[ATOMTOKLEN];

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

   gettok();
   if (strcmp(tok,";")) inerr("; expected","");

   prmptr->idunion.parm.prm = (parmstruct *) get(sizeof(parmstruct));
   prm = prmptr->idunion.parm.prm;

   printf("rdprm:reading headers...\n");
/* read headers */
/* control for first frame */
   nitems=1;
   bytest=sizeof(int);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   fread(&bytes,sizeof(int),1,fileptr);

/* Discover version */
   nitems=20;
   bytest=4*sizeof(char);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Cdum,bytest,1,fileptr);
   nitems=1;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   fread(&bytes,sizeof(int),1,fileptr);

/* title */
   nitems=20;
   bytest=4*sizeof(char);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Cdum,bytest,1,fileptr);
   fread(&bytes,sizeof(int),1,fileptr);

/* title */
   nitems=20;
   bytest=4*sizeof(char);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Cdum,bytest,1,fileptr);
   fread(&bytes,sizeof(int),1,fileptr);

/* Forcefield types and masses */
   nitems=1;
   bytest=sizeof(int);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(&prm->Ntypes,bytest,1,fileptr);
   printf("rdprm:Ntypes %6d\n", prm->Ntypes);
   nitems=prm->Ntypes;
   bytest=4*sizeof(char);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Cdum,bytest,1,fileptr);
   nitems=prm->Ntypes;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   fread(&bytes,sizeof(int),1,fileptr);

/* Residue names */
   nitems=1;
   bytest=sizeof(int);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(&prm->Nres,bytest,1,fileptr);
   printf("rdprm:Nres %6d\n", prm->Nres);
/*  prm->ResNames = (char *) get(sizeof(char)*4*prm->Nres); */
   nitems=prm->Nres;
   bytest=4*sizeof(char);
   for(i=0,ir=0;i<nitems;i++) {
/*
    ir+=fread(prm->ResNames,bytest,1,fileptr);
    prm->ResNames+=4;
*/
    ir+=fread(Cdum,bytest,1,fileptr);
   }
/*   prm->ResNames-=4*prm->Nres; */
   fread(&bytes,sizeof(int),1,fileptr);

/* index to forcefield and atom names */
   nitems=1;
   bytest=sizeof(int);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(&prm->Nsit,bytest,1,fileptr);
   printf("rdprm:Nsit %6d\n", prm->Nsit);
   prm->Natom = prm->Nsit;
   nitems=prm->Nsit;
   bytest=sizeof(int);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
/* sizeAtName = 5 if Discover 2.x, = 4 if Discover 3.0 or later */
   sizeAtName = 5;
   if(sizeAtName>ATOMTOKLEN) {
    printf("rdprm: sizeAtName too large.\n");
    exit(1);
   }
   prm->AtomNames = (char *) get(sizeof(char)*ATOMTOKLEN*prm->Nsit);
   nitems=prm->Nsit;
   bytest=sizeAtName*sizeof(char);
   for(i=0,ir=0;i<nitems;i++) {
    ir+=fread(atomname,bytest,1,fileptr);
    for(j=strlen(atomname)+1;j<ATOMTOKLEN;j++) strcat(atomname," ");
/*    printf("atomname(%d) %s\n",i,atomname); */
    strncpy(&prm->AtomNames[i*ATOMTOKLEN],atomname,ATOMTOKLEN);
   }
   fread(&bytes,sizeof(int),1,fileptr);

/* moveable atoms */
   nitems=1;
   bytest=sizeof(int);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   nitems=1;
   bytest=sizeof(int);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(&prm->Nmov,bytest,1,fileptr);
   printf("rdprm:Nmov %6d\n", prm->Nmov);
   prm->Imov_at = (int *) get(sizeof(int)*prm->Nmov);
   prm->Iat_mov = (int *) get(sizeof(int)*prm->Nsit);
   nitems=prm->Nmov;
   bytest=sizeof(int);
   for(i=0,ir=0;i<nitems;i++) {
    ir+=fread(prm->Imov_at,bytest,1,fileptr);
    prm->Imov_at++;
   }
   prm->Imov_at-=prm->Nsit;
   fread(&bytes,sizeof(int),1,fileptr);

   for(i=0;i<prm->Nsit;i++) {
    *prm->Iat_mov++=-1;
   }
   prm->Iat_mov-=prm->Nsit;
   for(i=0;i<prm->Nmov;i++) {
    ir=*prm->Imov_at++;
    prm->Iat_mov[ir-1]=i;
   }
   prm->Imov_at-=prm->Nmov;

/* molecules */
   nitems=1;
   bytest=sizeof(int);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(&prm->Nmolt,bytest,1,fileptr);
   printf("rdprm:Nmolt %6d\n", prm->Nmolt);
   prm->Nmol=1;
   nitems=prm->Nmolt;
   bytest=sizeof(int);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   nitems=prm->Nmolt;
   bytest=sizeof(int);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   fread(&bytes,sizeof(int),1,fileptr);

/* residues */
   nitems=1;
   bytest=sizeof(int);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(&prm->Nrest,bytest,1,fileptr);
   printf("rdprm:Nrest %6d\n", prm->Nrest);
   nitems=prm->Nrest;
   bytest=sizeof(int);
   prm->Ipres = (int *) get(sizeof(int)*prm->Nrest);
   prm->iend = (int *) get(sizeof(int)*prm->Nrest);
   for(i=0,ir=0;i<nitems;i++) {
    ir+=fread(&prm->Ipres[i],bytest,1,fileptr);
    ir+=fread(&prm->iend[i],bytest,1,fileptr);
   }
   nitems=prm->Nrest;
   bytest=sizeof(int);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   fread(&bytes,sizeof(int),1,fileptr);
/* set ResNum */
   prm->ResNum = (int *) get(sizeof(int)*prm->Nsit);
   for(i=0;i<prm->Nrest;i++) {
    for(j=prm->Ipres[i]-1;j<prm->iend[i];j++) {
     prm->ResNum[j]=i+1;
    }
   }

/* bonds */
   nitems=1;
   bytest=sizeof(int);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(&prm->Nbond,bytest,1,fileptr);
   printf("rdprm:Nbond %6d\n", prm->Nbond);
   fread(&bytes,sizeof(int),1,fileptr);

   if(prm->Nbond!=0) {
   nitems=2*prm->Nbond;
   bytest=sizeof(int);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   fread(&bytes,sizeof(int),1,fileptr);
   }

/* cell */
   nitems=6;
   bytest=sizeof(double);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=3*3*196;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=3*196;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=3*3*196;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=1;
   bytest=sizeof(int);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   nitems=4;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
/*
   nitems=1;
   bytest=sizeof(int);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
*/
   nitems=6;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=6;
   bytest=sizeof(int);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   fread(&bytes,sizeof(int),1,fileptr);

/* steps */
   nitems=1;
   bytest=sizeof(int);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(&prm->Nener,bytest,1,fileptr);
   printf("rdprm:Nener: %d\n",prm->Nener);
   prm->Ener = (double *) get(sizeof(double)*prm->Nener);
   nitems=1;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   printf("rdprm:dt: %e\n",*Ddum);
   nitems=1;
   bytest=sizeof(int);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   printf("rdprm:nsavex: %d\n",*Idum);
   nitems=1;
   bytest=sizeof(int);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Idum,bytest,1,fileptr);
   fread(&bytes,sizeof(int),1,fileptr);

/* energies */
   nitems=1;
   bytest=sizeof(double);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   printf("rdprm:Etot: %e\n",*Ddum);
   nitems=1;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   printf("rdprm:Epot: %e\n",*Ddum);
   nitems=1;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   printf("rdprm:Ekin: %e\n",*Ddum);
   nitems=prm->Nener;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=prm->Nmolt;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=prm->Nener*prm->Nmolt;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=prm->Nmolt;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=prm->Nmolt;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=prm->Nmolt;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=prm->Nmolt;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=1;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   printf("rdprm:Press: %e\n",*Ddum);
   nitems=1;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   nitems=9;
   bytest=sizeof(double);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Ddum,bytest,1,fileptr);
   fread(&bytes,sizeof(int),1,fileptr);

/* coordinates and velocities */
   nitems=3*prm->Nsit;
   bytest=sizeof(float);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) {
    ir+=fread(Fdum,bytest,1,fileptr);
/*  if(i<12) printf("%e\n",*Fdum); */
   }
   fread(&bytes,sizeof(int),1,fileptr);
   nitems=3*prm->Nsit;
   bytest=sizeof(float);
   fread(&bytes,sizeof(int),1,fileptr);
   for(i=0,ir=0;i<nitems;i++) ir+=fread(Fdum,bytest,1,fileptr);
   fread(&bytes,sizeof(int),1,fileptr);

   rewind(fileptr);
   genclose(fileptr, prmptr->idunion.parm.name);

   printf("end reading topology.\n");
/*
   for(i=0;i<prm->Nsit;i++) {
    strncpy(atomname,&prm->AtomNames[i*ATOMTOKLEN],ATOMTOKLEN);
    printf("%4s\n",atomname);
   }
*/
   prm->sig0r=1./3.8;
   prm->eps0r=500.;
   prm->mass0r=0.012;

/* end reading */
}
