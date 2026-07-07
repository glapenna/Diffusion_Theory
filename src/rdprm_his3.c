#include "head.h"

// extern char     tok[];
extern idstruct *startid();
extern FILE     *genopen();
extern char     *get();

void
rdprm_his3()
{
   extern parmstruct *prm;
   FILE *fileptr;
   int swap;
   idstruct        *prmptr;
   char            id[TOKLEN];

   int i,j,line_count,read_count,line_good;
   int resnum;
   char atomname[ATOMTOKLEN],resname[ATOMTOKLEN],atomz[3],atomtype[9];
   char atomnamer[ATOMTOKLEN];
   char resnamer[ATOMTOKLEN];
   _REAL x,y,z,charge;
   char buffer[512],first_word[512];
   int NAMD_read_line(FILE *, char *), NAMD_blank_string(char *);
   void NAMD_find_first_word(char *, char *);
#ifdef DOUBLE
   static char scanform[]="%s %lf %lf %lf %s %d %s %s %lf";
#else
   static char scanform[]="%s %f %f %f %s %s %d %s %f";
#endif

   prmptr = startid(PARM, id);

   /* allow room for .Z */
   gettok();
   prmptr->idunion.parm.name = get(strlen(tok) + 3);
   strcpy(prmptr->idunion.parm.name, tok);
   if((fileptr = genopen(tok,"PARM")) == NULL) exit(1);

   prmptr->idunion.parm.prm = (parmstruct *) get(sizeof(parmstruct));
   prm = prmptr->idunion.parm.prm;

   prm->sig0r = 3.8;
   prm->eps0r = 500.;
   prm->mass0r= 0.012;

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

   printf("rdprm:reading MSI car file...\n");
/* figure out how many atoms are there */
   line_count = 0;
   line_good = 0;
   while (NAMD_read_line(fileptr, buffer) != -1)
   {
           line_count++;
           NAMD_find_first_word(buffer, first_word);
           if ((buffer[0] != '!') &&
               !NAMD_blank_string(buffer) &&
               (strncasecmp(first_word, "!BIOSYM", 7) != 0) &&
               (strncasecmp(first_word, "!DATE", 5) != 0) &&
               (strncasecmp(first_word, "PBC", 3) != 0) &&
               (strncasecmp(first_word, "end", 3) != 0) &&
               (line_count != 3) )
           {
            line_good++;
           }
   }
   prm->Nsit = line_good;
   printf("rdprm:Nsit: %d\n",prm->Nsit);
   rewind(fileptr);
   prm->AtomNames = (char *)get(prm->Nsit*sizeof(char)*ATOMTOKLEN);
   prm->ResNum = (int *)get(prm->Nsit*sizeof(int));
   prm->ResNamesa = (char *)get(prm->Nsit*sizeof(char)*ATOMTOKLEN);

   /*  Keep reading in lines until we hit the EOF                  */
   line_count = 0;
   line_good = 0;
   while (NAMD_read_line(fileptr, buffer) != -1)
    {
           /*  Get the first word of the line                      */
           line_count++;
           NAMD_find_first_word(buffer, first_word);

           /*  First, screen out things that we ignore, such as    */
           /*  blank lines, lines that start with '!', ...         */
           if ((buffer[0] != '!') &&
               !NAMD_blank_string(buffer) &&
               (strncasecmp(first_word, "!BIOSYM", 7) != 0) &&
               (strncasecmp(first_word, "!DATE", 5) != 0) &&
               (strncasecmp(first_word, "PBC", 3) != 0) &&
               (strncasecmp(first_word, "end", 3) != 0) &&
               (line_count != 3) )
           {
            /*  Now, get the required things from buffer    */
            /*  Use sscanf to parse up the input line       */
            line_good++;
            read_count=sscanf(buffer,scanform,
             atomname,&x,&y,&z,resname,&resnum,atomtype,atomz,&charge);
            sprintf(atomnamer,"%-4s",atomname);
            sprintf(resnamer,"%-5s",resname);
            if (read_count != 9)
            {
             char err_msg[512];
             sprintf(err_msg, "BAD LINE IN CAR FILE\nLINE=*%s*\n",
                buffer);
             printf("%s\n",err_msg);
             exit(1);
            }
            strncpy(&prm->AtomNames[ATOMTOKLEN*(line_good-1)],atomnamer,
              ATOMTOKLEN);
            strncpy(&prm->ResNamesa[ATOMTOKLEN*(line_good-1)],resname,
              ATOMTOKLEN);
            prm->ResNum[line_good-1]=resnum;
/*
            printf("atomname %d %s %s %d\n",
             line_good,
             &prm->AtomNames[ATOMTOKLEN*(line_good-1)],
             &prm->ResNamesa[ATOMTOKLEN*(line_good-1)],
             prm->ResNum[line_good-1] );
*/
           }
   }
   fclose(fileptr);
/* assumes all the atoms move */
   prm->Nmov = prm->Nsit;
/* assumes 18 energy values */
   prm->Nener = 18;
/* assumes some defaults... */
   prm->Ntypes = 1;
   prm->Nmolt = 1;
   prm->Nmol = 1;
   prm->Natom = prm->Nsit;
   prm->Nres = 1;
   prm->Nrest = 1;
   prm->Nbond = 0;

   printf("end reading topology.\n");
   rewind(fileptr);
/* end reading */
/*
   for(i=0;i<prm->Nsit;i++) {
     strncpy(atomname,&prm->AtomNames[i*ATOMTOKLEN],ATOMTOKLEN);
     printf("%4s\n",atomname);
   }
*/
}
