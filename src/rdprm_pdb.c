#include "head.h"

// extern char     prmtop[],tok[];
// char     prmtop[],tok[];
extern idstruct *startid();
extern FILE     *genopen();
extern char     *get();

void
rdprm_pdb()
{
   extern parmstruct *prm;
   extern parmstruct *defprm;
   FILE *fileptr;
   int swap;
   idstruct        *prmptr;
   char            id[TOKLEN];

   int i,j,k,start,line_count,read_count,line_good;
   int resnum;
   char lab[]="MOL";
   char atomname[ATOMTOKLEN],resname[ATOMTOKLEN];
   char atomnamer[ATOMTOKLEN];
   // char resnamer[ATOMTOKLEN];
   _REAL x,y,z,charge,rlj;
   char buffer[512],first_word[512];
   int NAMD_read_line(FILE *, char *), NAMD_blank_string(char *);
   void NAMD_find_first_word(char *, char *);
#ifdef DOUBLE
   static char scanform[]="%*s %*s %s %s %d %lf %lf %lf";
#else
   static char scanform[]="%*s %*s %s %s %d %f %f %f";
#endif

   prmptr = startid(PARM, id);

   /* allow room for .Z */
   gettok();
   prmptr->idunion.parm.name = get(strlen(tok) + 3);
   strcpy(prmptr->idunion.parm.name, tok);
   if((fileptr = genopen(tok,"PARM")) == NULL) exit(1);

   prmptr->idunion.parm.prm = (parmstruct *) get(sizeof(parmstruct));
   prm = prmptr->idunion.parm.prm;
   if(defprm==NULL) {
    defprm=prm;
    strncpy(prmtop,tok,TOKLEN-1);
   }

   prm->sig0r = 10.;
   prm->eps0r = 500.;
   prm->mass0r= 1.000;

   for (;;) {
    gettok();
    if (tok[0]==';') {
     break;
    } else if (!strcmp(tok, "SWAP")) {
     swap=1;
    } else if (!strcmp(tok, "NORED")) {
     prm->sig0r = 1.;
    }
   }

   printf("rdprm:reading PDB file...\n");
/* figure out how many atoms are there */
   line_count = 0;
   line_good = 0;
   start=1;
   prm->Nres=0;
   while (NAMD_read_line(fileptr, buffer) != -1)
   {
           line_count++;
           NAMD_find_first_word(buffer, first_word);
           if ((buffer[0] != '!') &&
               !NAMD_blank_string(buffer) &&
               (strncasecmp(first_word, "REMARK", 6) != 0) &&
               (strncasecmp(first_word, "CRYST", 5) != 0) &&
               (strncasecmp(first_word, "HEADER", 6) != 0) &&
               (strncasecmp(first_word, "COMPND", 6) != 0) &&
               (strncasecmp(first_word, "AUTHOR", 6) != 0) &&
               (strncasecmp(first_word, "MASTER", 6) != 0) &&
               (strncasecmp(first_word, "TER", 3) != 0) &&
               (strncasecmp(first_word, "END", 3) != 0) )
           {
            line_good++;
/*
            read_count=sscanf(buffer,scanform,
             atomname,resname,&resnum,&x,&y,&z,&rlj,&charge);
*/
            read_count=sscanf(buffer,scanform,
             atomname,resname,&resnum,&x,&y,&z);
	    // printf("resnum= %d\n",resnum); fflush(stdout);
            if (read_count != 6)
            {
             char err_msg[512];
             sprintf(err_msg, "BAD LINE IN PDB FILE\nLINE=*%s*\n",
                buffer);
             printf("%s\n",err_msg);
             exit(1);
            }
            if(resnum!=start) {
             prm->Nres++;
             start=resnum;
            }
           }
   }
   prm->Nsit = line_good;
   printf("rdprm:Nsit: %d\n",prm->Nsit);
   printf("rdprm:Nres: %d\n",prm->Nres);
   fflush(stdout);
   rewind(fileptr);
   prm->AtomNames = (char *)get(prm->Nsit*sizeof(char)*ATOMTOKLEN);
   prm->ResNum = (int *)get(prm->Nsit*sizeof(int));
   prm->Ipres = (int *)get((prm->Nres+1)*sizeof(int));
   prm->ResNames = (char *)get((prm->Nres+1)*sizeof(char)*ATOMTOKLEN);
   prm->ResNamesa = (char *)get(prm->Nsit*sizeof(char)*ATOMTOKLEN);
   prm->Masses = (_REAL *)get(prm->Nsit*sizeof(_REAL));
   prm->Charge = (_REAL *)get(prm->Nsit*sizeof(_REAL));
   prm->ljAB = (_REAL *)get(prm->Nsit*sizeof(_REAL));
   for(i=0;i<prm->Nsit;i++) {
    strncpy(&prm->ResNamesa[ATOMTOKLEN*i],lab,3);
   }

   /*  Keep reading in lines until we hit the EOF                  */
   line_count = 0;
   line_good = 0;
   start=1;
   k=0;
   prm->Ipres[0]=1;
   while (NAMD_read_line(fileptr, buffer) != -1)
    {
           /*  Get the first word of the line                      */
           line_count++;
   /* printf("reading line %d.\n",line_count); */
           NAMD_find_first_word(buffer, first_word);

           /*  First, screen out things that we ignore, such as    */
           /*  blank lines, lines that start with '!', ...         */
           if ((buffer[0] != '!') &&
               !NAMD_blank_string(buffer) &&
               (strncasecmp(first_word, "REMARK", 6) != 0) &&
               (strncasecmp(first_word, "CRYST", 5) != 0) &&
               (strncasecmp(first_word, "HEADER", 6) != 0) &&
               (strncasecmp(first_word, "COMPND", 6) != 0) &&
               (strncasecmp(first_word, "AUTHOR", 6) != 0) &&
               (strncasecmp(first_word, "MASTER", 6) != 0) &&
               (strncasecmp(first_word, "TER", 3) != 0) &&
               (strncasecmp(first_word, "END", 3) != 0) )
           {
            /*  Now, get the required things from buffer    */
            /*  Use sscanf to parse up the input line       */
            line_good++;
/*
            read_count=sscanf(buffer,scanform,
             atomname,resname,&resnum,&x,&y,&z,&rlj,&charge);
*/
            read_count=sscanf(buffer,scanform,
             atomname,resname,&resnum,&x,&y,&z);
            if (read_count != 6)
            {
             char err_msg[512];
             sprintf(err_msg, "BAD LINE IN PDB FILE\nLINE=*%s*\n",
                buffer);
             printf("%s\n",err_msg);
             exit(1);
            }
            // sprintf(resnamer,"%s    ",resname);
            if(line_good==1) {
             strncpy(&prm->ResNames[ATOMTOKLEN*k],resname,ATOMTOKLEN);
            }
            if(resnum!=start) {
             k++;
             strncpy(&prm->ResNames[ATOMTOKLEN*k],resname,ATOMTOKLEN);
             prm->Ipres[k]=line_good;
             start=resnum;
            }
            // sprintf(atomnamer,"%s    ",atomname);
            strncpy(&prm->AtomNames[ATOMTOKLEN*(line_good-1)],atomname,
              ATOMTOKLEN);
/*
            prm->ljAB[line_good-1]=rlj;
            prm->Charge[line_good-1]=charge;
*/
/*
            strncpy(&prm->ResNamesa[ATOMTOKLEN*(line_good-1)],resname,
              ATOMTOKLEN);
*/
            prm->ResNum[line_good-1]=resnum;
/*
            printf("atomname %d %d\n", line_good,k);
            printf("atomname %s %s %d %s %d\n",
             &prm->AtomNames[ATOMTOKLEN*(line_good-1)],
             &prm->ResNames[ATOMTOKLEN*k],
             prm->Ipres[k],
             &prm->ResNamesa[ATOMTOKLEN*(line_good-1)],
             prm->ResNum[line_good-1] );
*/
           }
   }
   fclose(fileptr);
/* assumes some defaults... */
   prm->Ipres[prm->Nres] = prm->Nsit+1;
   prm->Ntypes = 1;
   prm->Nmol = 1;
   prm->Natom = prm->Nsit;
   prm->Nsolu = 1;
   prm->Natom_solu = prm->Nsit;
//   prm->Natom_solu = 123;
   prm->Nsolu_tot = prm->Nsolu*prm->Natom_solu;
   prm->atomtype = (int *)get(prm->Nsit*sizeof(int));

/* assign masses */
   for(i=0;i<prm->Nsit;i++) {
    prm->atomtype[i]=1;
    strncpy(atomname,&prm->AtomNames[i*ATOMTOKLEN],ATOMTOKLEN);
    if(atomname[0]=='H') {
     prm->Masses[i]=1.008/prm->mass0r;
    } else if(atomname[0]=='C') {
     prm->Masses[i]=12.01/prm->mass0r;
     if(atomname[1]=='L') {
      prm->Masses[i]=35.45/prm->mass0r;
     }
    } else if(atomname[0]=='O') {
     prm->Masses[i]=16.00/prm->mass0r;
    } else if(atomname[0]=='N') {
     prm->Masses[i]=14.01/prm->mass0r;
    } else if(atomname[0]=='S') {
     prm->Masses[i]=32.06/prm->mass0r;
    } else if(atomname[0]=='P') {
     prm->Masses[i]=30.97/prm->mass0r;
    } else if(atomname[0]=='F') {
     prm->Masses[i]=55.85/prm->mass0r;
    } else if(atomname[0]=='R') {
     prm->Masses[i]=101.07/prm->mass0r;
    } else if(atomname[0]=='Z') {
     prm->Masses[i]=65.38/prm->mass0r;
    } else {
     printf("Mass not found for atom %d.\n",i);
     printf("Assuming hydrogen.\n");
     prm->Masses[i]=1.008/prm->mass0r;
    }
   }

   printf("end reading topology.\n");
/* end reading */
/*
   for(i=0;i<prm->Nsit;i++) {
     strncpy(atomname,&prm->AtomNames[i*ATOMTOKLEN],ATOMTOKLEN);
     printf("%4s ",atomname);
     strncpy(resname,&prm->ResNamesa[i*ATOMTOKLEN],ATOMTOKLEN);
     printf("%4s ",resname);
     resnum=prm->ResNum[i];
     printf("%d\n",resnum);
   }
*/
}
