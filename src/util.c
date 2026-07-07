#include "head.h"
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/errno.h>
extern int      errno;
/*
extern int      sys_nerr;
extern char     *sys_errlist[];
*/

extern char            analin[], analout[];
extern int debug;
extern idstruct         *firstser;

idstruct        *idlist = NULL, *lastid = NULL;
// char            tok[TOKLEN];
char            line[TOKLEN];
int             linenum = 0;
static int      lineindex = 0;
static int      tab = 0;

void *get2(int,int);
char     *get();
/*
void     gettok(), inerr();
*/
char     *syserr();
setstruct       *getstatset();


char *tokarray[] = {
        "FILES_IN",     "FILES_OUT",    "DECLARE",  
        "STOP",         "END",      
        "PARM",         "STREAM",       "CRD",          "TABLE",
        "BOND", "ANGLE", "DIHE", "GR",  "ORDERB", "DENSITY1",
        "BEAD", "MCARIG", "RIG", "ACFF", 
        ";",            ",",            "&",            "|",
        "!",            "(",            ")",            "-",
        ".",            "?" };

/* NTOK is the count of items in tokarray */
#define NTOK 29


/***********************************************************************
 							FILEARGS()
************************************************************************/

/*
 * fileargs() - Unix style get args for filenames
 */

void
fileargs(int argc, char *argv[])
{
	extern  char    *optarg;
	char		*arg0 = argv[0];
	int		cc, in = 0, out = 0, zcat = 0;

	analin[0] = 0;
        while ((cc = getopt(argc, argv, "i:o:")) != EOF) {
		switch (cc) {
			case 'i':
				strcpy (analin, optarg);
				if (stdin != freopen(analin, "r", stdin)) {
					perror(analin);
					exit(1);
				}
				in = 1;
				break;
			case 'o':
				strcpy (analout, optarg);
				if (stdout != freopen(analout, "w", stdout)) {
					perror(analout);
					exit(1);
				}
				out = 1;
				break;
			default:
				usage(arg0);
				break;
		}
	}

	/*
	 * fill in filenames where necc & open if necc
	 */

	if (!in) {
		if (!isatty(fileno(stdin)) && !zcat)
			/* shell has redirected something to stdin */
			strcpy(analin, "stdin");
		else {
			/* use default filename */
			strcpy(analin, "analin");
			if (stdin != freopen(analin, "r", stdin)) {
				perror(analin);
				exit(1);
			}
		}
	}
	if (!out)
		strcpy(analout, "stdout");
}

/***********************************************************************
 							READLN()
************************************************************************/

int
readln(FILE *file, char *name, char *string)
{
	int 	i, j;

	for (i=0; i<81; i++) {
		if ((j = fgetc(file)) == EOF) {
		    printf("Error: unexpected EOF in %s (%d char in line)\n", 
						name, i);
		    exit(1);
		}
		string[i] = (char) j;
		if (string[i] == '\n') {
			string[i+1] = '\0';
			break;
		}
	}
	if (i == 80  &&  string[i] != '\n') {
		printf("Error: line too long in %s:\n%.80s", name, string);
		exit(1);
	}

	return(i);
}

/***********************************************************************
                                                        GETTOK()
************************************************************************/

void gettok()
{
	static int	linelen;
	int		i;

	for (i=0; ; ) {
		switch (line[lineindex]) {
			case ';':
			case ',':
			case '|':
			case '&':
			case ')':
			case '(':
/*	   	        case '-':    */
			case '!':
			case '%':
				if (!i) {
					/* first char */
					tok[0] = line[lineindex++];
					tok[1] = '\0';
					return;
				} else
					/* abutting prev tok */
					tok[i] = '\0';
					return;
			case '\t':
				tab++;
			case ' ':
			case '\0':
				lineindex++;
				if (i) {
					tok[i] = '\0';
					return;
				}
				break;
			case '\n':
				if (!i) {
					/* still waiting for 1st char -
				   	 get another line */
					for (;;) {
						linelen = 
						  readln(stdin, analin, line);
						linenum++;
						printf("> %.*s\n",linelen,line);
						if (line[0] !='#')
							break;
					}
					lineindex = 0;
					tab = 0;
					continue;
				} else {
					tok[i] = '\0';
					return;
				}
/*
  	   	        case '-':  
				tok[i++] = line[lineindex];
				lineindex++;
*/
			default:
				tok[i++] = line[lineindex];
				lineindex++;
				break;
		}
	}
}

/***********************************************************************
 							TOKINT()
************************************************************************/
/*
 *  tokint() - see if token is an int
 */

int
tokint(int *num)
{
	int	i, j = strlen(tok);

	for (i=0; i<j; i++)
		if (!isdigit(tok[i]))
			return(0);
	*num = atoi(tok);
	return(1);
}

/***********************************************************************
 							TOKREAL()
************************************************************************/

int
tokREAL(_REAL *num)
{
	int	i, j = strlen(tok), dot = 0;

	for (i=0; i<j; i++) {
                if (tok[i] == '-') i++;
		if (!isdigit(tok[i])) {
			if (tok[i] == '.') {
				if (dot++)
					return(0);
			} else
				return(0);
		}
	}
#ifdef DOUBLE
	sscanf(tok, "%lf", num);
#else
	sscanf(tok, "%f", num);
#endif
	/* "*num = (_REAL)atof(tok)" : wrong result */
	return(1);
}

/***********************************************************************
                                                        INERR()
************************************************************************/

/*
 *  inerr() - report error in input
 */

void
inerr(char *str1, char *str2)
{
        int     i;

        printf("%s", line);
        for (i=1; i<lineindex+7*tab; i++)
                printf(".");
        printf("^\n");
        printf("Error in %s line %d: %s %s\n    Last token was %s\n",
                analin, linenum, str1, str2, tok);
        exit(1);
}

/***********************************************************************
                                                        FILEPROB()
************************************************************************/

/*
 * fileprob() - make sure file doesn't exist unless APPENDing
 */

int
fileprob(char *path, int append, int *exists, char *caller)
{
        struct stat     buf;

        if (stat(path, &buf) == -1) {
                *exists = 0;
                switch (errno) {
                    case ENOENT:
                        if (append)
                                printf("%s warning: creating APPEND file %s\n",
                                                caller, path);
                        return(0);
                    case ENOTDIR:
                        printf("%s: missing directory in path\n", path);
                        return(1);
                    default:
                        printf("%s: sys err %s\n", path, syserr());
                        return(1);
                }
        }
        *exists = 1;
        if (!append) {
                printf("%s: file %s already exists\n", caller, path);
                return(1);
        }
        return(0);
}

/***********************************************************************
                                                        ISKEYWD()
************************************************************************/

int
iskeywd()
{
        int             i;

        for (i=0; i<NTOK; i++)
                if (!strcmp(tok, tokarray[i]))
                        return(1);
        return(0);
}

#define FILES_OUT_TYPE(TYPE) \
		( TYPE == ACFF || \
		  TYPE == TABLE)

/***********************************************************************
 							STARTID()
************************************************************************/

idstruct *
startid(int type, char *id)
{
	idstruct	*idptr, *lastx = idlist;

	/*
	 * get id & make sure it's not a reserved word or a duplicate
	 *		& set lastx to last non-CRD/TABLE
	 */

	gettok();
	if (iskeywd())
		inerr("Expected id, got keyword: ", tok);
/*
	printf("startid %s\n", tok);
*/
	for (idptr=idlist; idptr!=NULL; idptr=idptr->next) {

/*
		printf("idptr %d  id  %s type %d\n", 
					idptr, idptr->id_name, idptr->id_type);
*/
		if (!strcmp(tok, idptr->id_name))
			inerr("duplicate id:", tok);
		if (!FILES_OUT_TYPE(idptr->id_type))
			lastx = idptr;
	}

	/*
	 *  Get space for id record, enque & set its type.
	 *
	 *	Queueing order is based on runtime considerations:
	 *	the queue of idstructs amounts to a compiled list
	 *	of commands.
	 *
	 *	The FILES_OUT items (CRD, DIST, TABLE,..) go last so 
	 *	that at runtime the things that are directed at them 
	 *	are already calcd when they are reached.
	 *
	 *	The other ids go in order up to the FILES_OUT.
	 */

	idptr = (idstruct *) get2(1,sizeof(idstruct));
	idptr->next = NULL;
/*
	if (lastx != NULL)
		printf("lastx %s\n", lastx->id_name);
	printf("lastid %d, this id %d\n", lastid, idptr);
*/

	if (idlist == NULL) {

		/*
		 *  First id (a FILES_IN id)
		 */

		idlist = idptr;
		lastid = idptr;
	} else if (FILES_OUT_TYPE(type)) { 

		/* 
		 *  FILES_OUT: put at end 
		 */
/*
		printf("put at end\n");
*/
		lastid->next = idptr;
		lastid = idptr;
	} else {	

		/* 
		 *  before FILES_OUT ids if any
		 */
/*
		printf("put b4 crd/tbl\n");
*/
		idptr->next = lastx->next;
		lastx->next = idptr;

		/*
		 *  if no FILES_OUT type yet, move lastid to end
		 */

		if (!FILES_OUT_TYPE(lastid->id_type))
			lastid = idptr;
	}

	/*
	 *  set up struct
	 */
		
	idptr->id_type = type;
	idptr->id_subtype = 0;
	strcpy(idptr->id_name, tok);
	strcpy(id, tok);
	idptr->id_setup = 0;

	return(idptr);
}

/***********************************************************************
 							GETSTRUCT()
************************************************************************/

	/* getstruct() - get a structure with given id */

idstruct *
getstruct(char *id)
{
	idstruct	*idptr;

	for (idptr=idlist; ; idptr=idptr->next) {
		if (idptr == NULL)
			return(NULL);
		if (!strcmp (idptr->id_name, id)) {
			return(idptr);
		}
	}
}

/***********************************************************************
 							GENOPEN()
************************************************************************/

/*
 *  genopen() - fopen regular or popen compressed file for reading.
 *	Filename assumed to have room for .Z extension if no .Z.
 */

FILE *
genopen(char *name, char *caller)
{
	struct stat	buf;
	char		*cbuf;
	int		length, compressed;
	FILE		*fp;

	printf("%s: opening %s\n", caller, name);
	length = strlen(name);
	compressed = iscompressed(name);
	cbuf = get(length+3);	/* allow for .Z */
	strcpy(cbuf, name);

	/*
	 *  if file doesn't exist, maybe it has been compressed/decompressed
	 */

	if (stat(cbuf, &buf) == -1) {
		switch (errno) {
		  case ENOENT:	{
			if (!compressed) {
				/*
				 *  see if it's been compressed
				 */
				strcat(cbuf, ".Z");
				if (stat(cbuf, &buf) == -1) {
					printf("%s, %s: does not exist\n", 
						name, cbuf);
					free(cbuf);
					return(NULL);
				}
				printf("\tusing %s\n", cbuf);
				strcat(name, ".Z");
				compressed++;
			} else {
				/*
				 *  see if it's been uncompressed
				 */
				cbuf[length-2] = '\0';
				if (stat(cbuf, &buf) == -1) {
					printf("%s, %s: does not exist\n", 
							name, cbuf);
					free(cbuf);
					return(NULL);
				}
				printf("\tusing %s\n", cbuf);
				name[length-2] = '\0';
				compressed = 0;
			}
			break;
		}
		default:
			printf("%s: sys err %s\n", name, syserr());
			free(cbuf);
			return(NULL);
		}
	}

	/*
	 *  open the file
	 */

	if (compressed) {
		char *pcmd = get(strlen("zcat ") + strlen(cbuf) + 1);

		sprintf(pcmd, "zcat %s", cbuf);
		if ((fp = popen(pcmd, "r")) == NULL) {
			perror(pcmd);
			exit(1);
		}
		free(pcmd);
	} else {
		if ((fp = fopen(cbuf, "r")) == NULL) {
			perror(cbuf);
			exit(1);
		}
	}
	free(cbuf);
	return(fp);
}

/***********************************************************************
 							GENCLOSE()
************************************************************************/

/*
 *  genclose() - close fopened or popened file
 */

void genclose(FILE *fileptr, char *name)
{
 int i;
	if (iscompressed(name)) {
		if (pclose(fileptr) == -1)
			perror("pclose");
	} else {
		if (i=fclose(fileptr) == -1) {
			perror("fclose");
                }
                printf("i=%d\n",i);
	}
}

/***********************************************************************
                                                        ISCOMPRESSED()
************************************************************************/

/* 
 *  iscompressed() - look for .Z at end of name
 */

int
iscompressed(char *name)
{
        int     i = strlen(name) - 1;   /* last char in name */

        if (i < 0) {
                fprintf(stderr, "programming error: name w/ length %d\n", i);
                exit(1);
        }
        if (i < 3)
                return(0);
        if (name[i] == 'Z'  &&  name[i-1] == '.')
                return(1);
        return(0);
}

/***********************************************************************
                                                        GET()
************************************************************************/

char *
get(int size)
{
        char    *ptr;

        if (size ==0)
                return((char *) NULL);

        if ((ptr = (char *) malloc((unsigned)size)) == NULL) {
                printf("malloc %d", size);
                fflush(stdout);
                perror("malloc err:");
                exit(1);
        }
        return(ptr);
}

void *
get2(int n,int size)
{
        char    *ptr;

        if (size ==0)
                return((void *) NULL);
        if (n ==0)
                return((void *) NULL);

        if ((ptr = (void *) calloc(n,(unsigned)size)) == NULL) {
                printf("calloc %d", n*size);
                fflush(stdout);
                perror("calloc err:");
                exit(1);
        }
        return(ptr);
}



/***********************************************************************
                                                        SYSERR()
************************************************************************/

/*
 * syserr() - return string describing current errno
 */

char *
syserr()
{
        static char     num[8];
// old version
/*

        if (errno < sys_nerr)
                return(sys_errlist[errno]);
*/
// new version
        if (errno < ERANGE)
                return(strerror(ERANGE));
        sprintf(num, "%d", errno);
        return(num);
}


/***********************************************************************
                                                        PARR()
************************************************************************/

/*
 *  parr() - print 3x3 array
 */

void parr(_REAL *a)
{
        int i, i3, j;

        for (i=i3=0; i<3; i++,i3+=3) {
                for (j=0; j<3; j++)
                        printf("%e    ",a[i3 + j]);
                printf("\n");
        }
}

/**********************************************************************/
/*
static void
*/
void
reverse4ByteOrder(int nele, char *dtaptr)
/* nele specifies number of elements to convert */
/* dtaptr specify and return conversion data    */
/*----------------------------------------------------------------------
Reverse byte order in 4 byte word (inplace).
----------------------------------------------------------------------*/
{
  register int *ip;
  register char *ca;
  register int m;
  union {
    char    ch[4];
    int     i;
  }       bs;

  ca = (char *) dtaptr;
  ip = (int *) dtaptr;
  for (m = 0; m < nele; m++) {
    bs.ch[3] = *ca++;
    bs.ch[2] = *ca++;
    bs.ch[1] = *ca++;
    bs.ch[0] = *ca++;
    *ip++ = bs.i;
  }

}

/************************************************************************                                                        GETSETADDR()
*************************************************************************/

/*
 *  getsetaddr() - starting w/ current tok, see if a set is specd, otherwise
 *      return default
 */

setstruct *
getsetaddr()
{
        idstruct        *idptr;

        /*
         * tok is possibly a stream or stat id; if so figure out set addr         *      leave fresh tok for next stage
         */

        gettok();
        if ((idptr = getstruct(tok)) != NULL) {
                        gettok();
                        return(idptr->idunion.ser.set);
        }
        /*
         *  default
         */

        return(firstser->idunion.ser.set);
}

/***********************************************************************
                                                        SETATOM()
************************************************************************/
/*
 *  setatom() - given an atom number, figure out crd set and return
 *      pointer to X crd of atom in set
 */

 _REAL *
setatom(int atom)
{
        setstruct       *set;

        if (atom < 1)
                inerr("Atom number less than 1","");

        /*
         * bind set
         */

        set = getsetaddr();
        if (atom > set->max)
                inerr("Atom number greater than crd set","");
#ifdef DBG
        printf("$$atom # got set, addr %d \n", set);
#endif
        atom--;
        return(&set->crd[3*atom]);
}

/***********************************************************************
 							GETSETATOM()
************************************************************************/
/*
 *  getsetatom() - starting w/ current tok, figure out an atom # & return
 *	address of X crd of atom, leaving a fresh token. 
 *	Assumes current tok is atom name (e.g. not id)
 */
 
_REAL *
getsetatom()
{
	int		num = 0, res, i, a0, a1;
	char		atomname[2*TOKLEN], *AtomName;
	idstruct	*idptr;
	setstruct	*set;
 
	/*
	 *  maybe atom name
	 */
 
#ifdef DBG
	printf("$$let me see.. %.4s\n", tok);
#endif
	if (strlen(tok) > ATOMTOKLEN)
		inerr("Not an atom name (or anything else)", "");
	sprintf(atomname, "%s    ", tok);
 
	/*
	 * assuming it's an atom name for now, look for a residue number
	 */
 
	gettok();
	if (!tokint(&res))
		inerr("Expected /atomname resnumber/ at least", "");
 
#ifdef DBG
	printf("$$res is %d\n", res);
#endif
	/*
	 *  assuming atom name & res number, see if there's a stream/parm ref
	 */
 
	gettok();
	if ((idptr = getstruct(tok)) == NULL) {
 
		/*
		 *  no stream specd
		 */
 
		set = firstser->idunion.ser.set;
#ifdef DBG
		printf("$$no stream specd\n");
#endif
	} else if (idptr->id_type == STREAM) {
		set = idptr->idunion.ser.set;
		gettok();
	} else if (idptr->id_type == STATIC) {
		gettok();
		set = getstatset(idptr);
	} else
		inerr("Expected STREAM or STATIC id to bind atom name/res", "");
 
#ifdef DBG
	printf("$$look up atom name in residue:\n");
#endif
	/*
	 *  look up the atom name in the residue
	 */
 
	a0 = set->prm->Ipres[res-1] - 1;	/* user's res counts from 1 */
	a1 = set->prm->Ipres[res] - 1;
	AtomName = &set->prm->AtomNames[a0 * ATOMTOKLEN];
	num = -1;
#ifdef DBG
	printf("atomname %.4s  res %d   a0 %d  a1 %d\n", atomname, res, a0, a1);
#endif
	for (i=a0; i<a1; i++) {
#ifdef DBG
		printf("%s  AtomName (%.4s)(%c)\n", 
				atomname, AtomName, AtomName[0]);
#endif
		if(!strncmp(AtomName, atomname, ATOMTOKLEN)) {
			num = i;
			break;
		}
		AtomName += ATOMTOKLEN;
	}
	if (num < 0) {
		printf("No atom %.*s in residue %d\n",ATOMTOKLEN,atomname, res);
		return(NULL);
	}
 
	/*
	 *  bind set
	 */
 
	if (num >= set->max)
		inerr("Atom number greater than crd set max","");
#ifdef DBG
	printf("$$user atom # %d got set\n", num+1);
#endif
	return(&set->crd[3*num]);
}

/***********************************************************************
                                                        GETSTATSET()
************************************************************************/

setstruct *
getstatset(idstruct *idptr)
{
        setstruct       *setptr;
        int             setnum, i;

        if (tok[0] != '%')
                return(idptr->idunion.stat.sets);

        /*
         *  else set is itemized: get item # look up addr & advance tok
         */

        gettok();
        if (!tokint(&setnum))
                inerr("expected STATIC set #","");
        if (setnum < 1)
                inerr("stat set # < 1","");
        if (setnum > idptr->idunion.stat.nsets)
                inerr("set # > max STATIC set","");
        gettok();
        setptr = idptr->idunion.stat.sets;
        setnum--;
        for (i=0; i<setnum; setptr = setptr->next);
        return(setptr);
}

/***********************************************************************
 							CHECKATOMS()
************************************************************************/

int
checkatoms(parmstruct *prm, int max)
{
	int	atoms;

	/*
	 * if tok == ATOM, read # & check it & advance to next token
	 * 	else return default max or prm->Natom
	 *	max is > 0 if the set is already defined
	 */

	if (strcmp(tok, "ATOM")) {	/* i.e. tok != ATOM */
		if (max)
			return(max);
		if (prm == NULL)
			inerr(
			  "Need 'ATOM #' or 'PARM file' or '-p prmtop' arg\n",
			  "\t\t\tto determine # atoms for array size");
		return(prm->Natom);
	}
	gettok();
	if (!tokint(&atoms))
		inerr("ATOM: Expected # of atoms", "");
	if (atoms < 1)
		inerr("ATOM: need >= 1 atom","");

	if (prm != NULL) {
		if (max  &&  atoms > max)
			inerr("More ATOMs than in set","");
		if (atoms > prm->Natom  &&  !max) 
			printf("Warning: more atoms than in prmtop file");
	}
	gettok();
	return(atoms);
}

/***********************************************************************
 							FMTNSTART()
************************************************************************/

/*
 * fmtnstart() - figure out file format 
 *	& position FILE for reading coords
 */

int
fmtnstart(FILE **file, char *name, int atoms)
{
	char	line1[81], line2[81], line3[81];
	int	nat, ifmt = 0;
	int	popn = iscompressed(name);

	/* assuming all files have >= 3 lines */

	readln(*file, name, line1);
	readln(*file, name, line2);
	readln(*file, name, line3);

	/*
	 * figure out format by inspecting 1st 3 lines for "ATOM" or 
	 *	decimal points in appropriate places
	 */
	/* might be a pdb */
	if (!strncmp(line1, "ATOM ", 5)) {
		ifmt = PDB;
	} else if (line2[4] == '.'  &&  line2[12] == '.'  &&  
	    		    line2[20] == '.'  && line2[28] == '.'  &&  
			    line2[36] == '.'  && line2[44] == '.'  &&
	    		    line2[52] == '.'  && line2[60] == '.'  &&  
			    line2[68] == '.'  && line2[76] == '.') {
		ifmt = CRD;
	} else {
		if (sscanf(line2, "%d", &nat)  &&  
			(line3[4] == '.'   &&  line3[16] == '.'  &&  
			 line3[28] == '.'  &&  line3[40] == '.'  &&  
			 line3[52] == '.'  &&  line3[64] == '.') ||
			(line3[5] == '.'   &&  line3[17] == '.'  &&  
			 line3[29] == '.'  &&  line3[41] == '.'  &&  
			 line3[53] == '.'  &&  line3[65] == '.') ) {
				ifmt = RST;
				if (atoms != nat) {
					if (atoms > nat) {
				printf("Error in %s: \n %d %s %d %s %s\n",
				    analin, atoms, "atoms specified, more than",
				    nat, "found in restrt file", name);
						exit(1);
					}
			printf("Warning: restrt file %s has %d %s %d %s\n",
					name, nat, "atoms, only using", atoms,
					"according to ATOM stmt or prmtop");
		}}}
	if (!ifmt) {
		/* default: it wasn't anything! */
		printf("Error: file: %s: format is neither pdb, crd or rst\n",
			name);
		exit(1);
	}

	/*
	 * now rewind & skip any header/noncoordinate lines
	 */

	if (popn) {
		char	pcmd[120];

		sprintf(pcmd, "zcat %s", name);
		if (pclose(*file) == -1)
			perror(name);
		if ((*file = popen(pcmd, "r")) == NULL) {
			perror(pcmd);
			exit(1);
		}
	} else {
		freopen(name, "r", *file);
	}
	if (ifmt == CRD) {
		readln (*file, name, line1);
		if (debug)
			printf("Format is CRD\n");
        } else if (ifmt == RST) {
		readln (*file, name, line1);
		readln (*file, name, line1);
		if (debug)
			printf("Format is RST\n");
	}
	else if (debug)
		printf("Format is PDB\n");
	return(ifmt);
}


/***********************************************************************
 							SCANFMT()
************************************************************************/

/* 
 *  scanfmt() - read xyz for 1 atom from pdb or other format file 
 */

int
scanfmt(int ifmt, FILE *file, _REAL xyz[])
{
	int	i,c;

	if (ifmt != PDB) {
/* TODO - read in a line at atime, & check for only 3 #s on line at beginning
   of set in case user forgot to use BOX in CRD format  */
   
#ifdef DOUBLE
		i = fscanf(file, " %lf %lf %lf", 
#else
		i = fscanf(file, " %f %f %f", 
#endif
				&xyz[0], &xyz[1], &xyz[2]); 
	} else {
		i = fscanf(file, 
#ifdef DOUBLE
			"%*s %*s %*s %*s %*s %lf %lf %lf", 
#else
			"%*s %*s %*s %*s %*s %f %f %f", 
#endif
				&xyz[0], &xyz[1], &xyz[2]);
		for(;;) {	/* skip to eoln */
			c = getc(file); 
			if ( c == 10  ||  c == EOF)
				break;
		}
	}

	if (i == -1  &&  errno) {
		printf("read error (%s) - ignoring if assuming eof works ok\n", 
							syserr());
	}
	return(i-3);
}

