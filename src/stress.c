#include "head.h"

extern int	errno;

/***********************************
	global vars
************************************/

/****	default crd set for measurements   ****/

idstruct	*firstserw = NULL;

/***********************************
	external vars
************************************/

/****	util.c	****/
// extern char	tok[];
// char	tok[];

/****	doit.c	****/
extern int	itW;

/****   misc routines	****/

extern idstruct	*startid(), *getstruct();
extern FILE	*genopen();
// extern char	*get2(),*get();
char	*get();
void	*get2(int,int);
extern void reverse4ByteOrder();

void setupstreamW()
{
        extern parmstruct *prm;
	idstruct	*serptrW, *idptr;
	namestruct	*nameptr, *lastname = NULL;
	char		id[TOKLEN];
	int		i,gotfmt=0;
        void startserW();

	serptr = startid(STRESS, id);
 
	if (firstserW == NULL)
		firstserW = serptr;

        serptr->idunion.serW.ifile = 0;
        serptr->idunion.serW.conv = 0;
        serptr->idunion.serW.fmt = STDF;
        serptr->idunion.serW.nfiles = 0;
        serptr->idunion.serW.names = NULL;
	serptr->idunion.serW.fileptr = NULL;
        serptr->idunion.serW.skip = 0;
        serptr->idunion.serW.win = 0;
        serptr->idunion.serW.winct = 0;

	for (;;) {
                gettok();
		if (!strcmp(tok, "REVERSE_FLOAT")) {
                        serptr->idunion.serW.conv=1;
                } else if (!strcmp(tok, "FORM")) {
                        gettok();
                        if (!strcmp(tok, "STDF")) {
                          serptr->idunion.serW.fmt = STDF;
                        } else if (!strcmp(tok, "STDD")) {
                          serptr->idunion.serW.fmt = STDD;
                        } else if (!strcmp(tok, "XYZ")) {
                          serptr->idunion.serW.fmt = XYZ;
                        } else {
                         inerr("STRESS: fmt unknown", "");
                         exit(1);
                        }
                        gotfmt++;
		} else if (!strcmp(tok, "WIN")) {

			/*
			 *  window mode = {skip m, read n} forever
			 */

			if (serptr->idunion.serW.skip)
			  inerr("STRESS: 2 WIN statements?", "");
			gettok();
			if (!tokint(&serptr->idunion.serW.skip)) 
			  inerr("STREAMW/WIN: expected # of sets to skip", "");
			if (serptr->idunion.serw.skip < 1)
			  inerr("STREAMW/WIN: skip # needs to be > 0","");
			serptr->idunion.serw.skip = -serptr->idunion.serw.skip;
			serptr->idunion.serw.winct = serptr->idunion.serw.skip;
			
			gettok();
			if (!tokint(&serptr->idunion.serw.win)) 
			  inerr("STREAMW/WIN: expected # of sets in window", "");
			if (serptr->idunion.serw.win < 1)
			  inerr("STREAMW/WIN: window # needs to be > 0","");
		} else
			break;
	}
        if (gotfmt > 1)
                inerr("STREAMW: more than 1 format given","");

	/*
	 * get file names, starting w/ current token
	 */

	for (;;) {

		/*
		 * check for done 
		 */

		if (tok[0] == ';') {
			if (!serptr->idunion.serw.nfiles)
				inerr("No files given in STREAMW stmt, id=", id);
			break;
		} 

		/*
		 *  assume token is file name: check for kwd, increment count 
		 * 	& add file name to list
		 */

		if (iskeywd())
			inerr("STREAMW: expected file name, got keyword: ", tok);
		serptr->idunion.serw.nfiles++;
		nameptr = (namestruct *) get2(1,sizeof(namestruct));
		/* allow for .Z extension to name */
		nameptr->name = get(strlen(tok) + 3);
		strcpy(nameptr->name, tok);
		if (lastname == NULL)
			serptr->idunion.serw.names = nameptr;
		else
			lastname->next = nameptr;
		lastname = nameptr;
		nameptr->next = NULL;

		gettok();
	}

	/*
	 * check files: can we open them & determine format?
	 */

	for (i=1; i<serptr->idunion.serw.nfiles; i++) {
		startserw(serptr, i);
		genclose(serptr->idunion.serw.fileptr, serptr->idunion.serw.name);
	}

	/*
	 *  leave the 1st one started  
	 *	TODO: leave this to a second pass over the idlist
	 */

	startserw(serptr, 0);

	serptr->id_setup++;
}


void startserw(idstruct *serptr, int i)
{
	namestruct 	*nameptr;
	int		j = 0;
        char title[TOKLEN];

	for (nameptr=serptr->idunion.serw.names; 
	     nameptr!=NULL; 
	     nameptr=nameptr->next) 
		if (j++ == i)
			break;

	if (nameptr == NULL) {
	  printf("Internal error: stream file name %d out of range\n", i);
	  exit(1);
	}

	/*
	 *  if something is already open, close it
	 */

	if (serptr->idunion.serw.fileptr != NULL) 
		genclose(serptr->idunion.serw.fileptr,serptr->idunion.serw.name);

	/*
	 *  now open the file
	 */

	serptr->idunion.serw.name = nameptr->name;
	if ((serptr->idunion.serw.fileptr = genopen(nameptr->name, "streamw"))
								 == NULL)
		exit(1);

/* read only binary or ascii files without skipping any header */
	serptr->idunion.serw.fileset = 0;
}

int
nextserw(idstruct *serptr)
{
	FILE	*fileptr = serptr->idunion.serw.fileptr;
	float	w_s = serptr->idunion.serw.w_s;
        int fmt=serptr->idunion.serw.fmt;
        int i5,ifin=0,itest=0;
        void startserw();

	printf("nextserw itw %d fileset %d\n", 
		itw, serptr->idunion.serw.fileset);
        fflush(stdout);

	/*
	 *  try reading set / open next file if necc;
	 *	keep going as req'd if WINdow mode has been specified
	 */

	for (;;) {

		/*
		 *  try to get another set from current file
		 */

		serptr->idunion.serw.fileset++;

      i5=0;
      switch(fmt)
      {

      case(STDF):

/* read C binary file, single precision then convert into d */
         i5+=fread(&w_s,sizeof(float),1,fileptr);
         if(serptr->idunion.serw.conv) reverse4ByteOrder(1,(char *)&w_s);
         serptr->idunion.serw.w=(_REAL)w_s;
         printf("w= %f\n",serptr->idunion.serw.w); fflush(stdout);
        if (i5!=1) {
         printf("something wrong reading w_s: %d\n",i5);
        }

        itest=1;
        ifin=i5;
        break;

      case(STDD):

/* read C binary file, double precision */
         i5+=fread(&serptr->idunion.serw.w,sizeof(_REAL),1,fileptr);
         printf("w= %f\n",serptr->idunion.serw.w); fflush(stdout);
        if (i5!=1) {
         printf("something wrong reading w: %d\n",i5);
        }

        itest=1;
        ifin=i5;
        break;

       case(XYZ):
        i5+=fscanf(fileptr, "%e ", &w_s);
        serptr->idunion.serw.w=(_REAL)w_s;
        printf("w= %f\n",serptr->idunion.serw.w); fflush(stdout);
        if (i5!=1) {
         printf("something wrong reading w_s: %d\n",i5);
        }
        itest=1;
        ifin=i5;
        break;

       default:
        printf("STREAM: fmt not found\n");
        exit(1);
        break;

       }

                if (ifin==itest) {

			int	wct;

		        printf("coord set ok:\n");

			/*
			 *  success: done unless in WINdow mode
			 */

			if (!serptr->idunion.serw.skip)
				return(1);

			/*
			 *  this stream is in WINdow mode - consider 
			 *	whether to skip this set
			 */

			wct = ++serptr->idunion.serw.winct;
			if (wct < 0) 
				continue;	/* skip another set */
			if (wct ==  serptr->idunion.serw.win) {
				serptr->idunion.serw.winct = 
						serptr->idunion.serw.skip;
				continue;
			}
			return(1);
		} else {
			printf("incomplete coord set: %s\n",
					serptr->idunion.serw.name); 
			printf("\tstarting another file if any\n");
		}
			
		/*
		 *  exhausted file: start next file
		 */

                rewind(fileptr);
		if (++serptr->idunion.serw.ifile 
                      == serptr->idunion.serw.nfiles) {
		    printf("STREAMW %s: no more files/crds (iteration %d)\n", 
						serptr->id_name, itw);
	            serptr->idunion.serw.fileset = 0;
	            serptr->idunion.serw.ifile = 0;
	            serptr->idunion.serw.winct = 
	             serptr->idunion.serw.skip;
                    startserw(serptr, 0);
		    return(0);
		}
                startserw(serptr, serptr->idunion.serw.ifile);
		fileptr = serptr->idunion.serw.fileptr;

		printf("nextserw next file (%d) started\n", 
					serptr->idunion.serw.ifile);

		/*
		 *  continue loop: read coordinates
		 */
	}

#ifdef lint
	return(0);
#endif
}
