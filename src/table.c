#include "head.h"

/*
 * Table operations - 'wrap' table has a single row starting with an 'A'
 *	that continues on the next line with a 'B' & so on; to get a
 *	page width of columns, grep the appropriate letter from the table.
 *	Non-wrap is intended to write a row across multiple files
 *	numbered off the table id - tabid0, tabid1...
 *
 */

/* #define DBG */

// extern	char		tok[];
// char		tok[];
extern	idstruct	*idlist;

static idstruct		*tabptr;
static	printstruct	*lastprint;
idstruct		*startid(), *getstruct();
char			*get();

void	enqueitem(), REALtabwrite(), tablereset();
/***********************************************************************
 							INITTABLE()
************************************************************************/

/*
 * inittable() - set up table idstruct in idqueue
 */

void
inittable()
{
	char		id[TOKLEN];
	int		exists;

	/*
	 * get table record in idlist & start filling in
	 */

	tabptr = startid(TABLE, id);

	gettok();
	if (iskeywd())
		inerr("TABLE: expected file name, got keyword: ", tok);
	strcpy(tabptr->idunion.tab.filename, tok);

	tabptr->idunion.tab.nitems = 0;

	/*
	 * for now only 'wrap' format w/ fixed # columns 
	 */

	tabptr->idunion.tab.nfiles = 1;
	tabptr->idunion.tab.append = 0;
	tabptr->idunion.tab.files = (filestruct *) get(sizeof(filestruct)); 
	tabptr->idunion.tab.curfile = tabptr->idunion.tab.files;
	tabptr->idunion.tab.files->next = NULL;
	strcpy(tabptr->idunion.tab.files->f_name, tok);

	gettok();
	if (!strcmp(tok, "APPEND")) {
		tabptr->idunion.tab.append = 1;
		gettok();
	}
	if (fileprob(tabptr->idunion.tab.filename, tabptr->idunion.tab.append, 
							&exists, "TABLE"))
		exit(1);

	if ((tabptr->idunion.tab.files->f_file 
			= fopen(tabptr->idunion.tab.filename, "a+")) == NULL) {
		perror(id);
		exit(1);
	}
	tabptr->idunion.tab.fptr = tabptr->idunion.tab.files->f_file;
	tabptr->idunion.tab.wrap = 1;
	strcpy(tabptr->idunion.tab.format, "%4e ");
	tabptr->idunion.tab.itemsline = 10;
	tabptr->idunion.tab.sofar = 0;
	tabptr->idunion.tab.linedelim = 'A';
	tabptr->idunion.tab.print = NULL;
	
	if (tok[0] != ';') 
		inerr("TABLE: unexpected token","");
}

/************************************************************************
							SETUPTAB()
*************************************************************************/

/*
 * setuptab() - look up table id in idqueue and build its print queue
 */

void
setuptab()
{
	idstruct	*idptr;
	int		items = 0;

	/*
	 * get table struct
	 */

	gettok();
	if ((tabptr = getstruct(tok)) == NULL)
		inerr("TABLE: no such id:", tok);
	if (tabptr->id_type != TABLE)
		inerr("TABLE: id isn't table:", tok);

	if (tabptr->id_setup)
		printf(
		  "Warning: TABLE %s multiply referenced in OUTPUT (adding)\n", 
									tok);
	lastprint = tabptr->idunion.tab.print;

	/*
	 * now put the items in the table's print queue 
	 */

	gettok();
	for (;;) {
		if (tok[0] == ';') {
			if (!items) {
				printf("OUTPUT TABLE: no items specified");
				exit(1);
			}
			tabptr->idunion.tab.nitems = items;
			tabptr->id_setup++;
			return;
		} else if (!strcmp(tok, "BOND")) {
			items += typenque(BOND);
		} else if (!strcmp(tok, "ANGLE")) {
			items += typenque(ANGLE);
		} else if (!strcmp(tok, "DIHE")) {
			items += typenque(DIHE);
		} else if (!strcmp(tok, "ORDERB")) {
			items += typenque(ORDERB);
		} else if (!strcmp(tok, "GYR")) {
			items += typenque(GYR);
		} else if (!strcmp(tok, "STREAM")) {
			items += typenque(STREAM);
                } else if (!strcmp(tok, "RMS")) {
                        items += typenque(RMS);
		} else if ((idptr = getstruct(tok)) != NULL) {
			enqueitem(idptr, 1);
			items++;
		} else
			inerr("TABLE: unexpected token","");
	}
}

/************************************************************************
							TYPENQUE()
*************************************************************************/

/*
 * typenque() - put all calcns of a given type on table's print queue
 */

int
typenque(int type)
{
	idstruct	*ptr;
	int		count = 0;

	for (ptr=idlist; ptr!=NULL; ptr=ptr->next) {
/*
		printf("--%d %d %d\n", ptr->id_type, ptr->id_subtype, type);
*/
/*
		if (ptr->id_type == STREAM  ||
			 ptr->id_type == TABLE)
			continue;
*/
		if (ptr->id_type == TABLE)
			continue;
		if (ptr->id_type == type) {
			enqueitem(ptr, 0);
			count++;
		}
	}
	if (!count)
		printf("Warning: no %s items specified\n", tok);
/*
	printf("%s: %d item(s) added\n", tok, count);
*/
	gettok();
	return(count);
}

/************************************************************************
							ENQUEITEM()
*************************************************************************/

/*
 * enqueitem() - put a particular idstruct item on table's print queue
 */

void
enqueitem(idstruct *ptr, int tget)
{	
	printstruct	*printptr;

	printptr = (printstruct *) get(sizeof(printstruct));

	printptr->next = NULL;
	switch (ptr->id_type) {
		case BOND:
		case ANGLE:
		case DIHE:
			printptr->p_value = (char *)
				&ptr->idunion.work.resunion.scalar.value;
			printptr->p_type = ptr->id_type;
			if (tget)
				gettok();
			break;
		case ORDERB:
			printptr->p_value = (char *)
				&ptr->idunion.work.resunion.scalar.stat.avg;
			printptr->p_type = ptr->id_type;
			if (tget)
				gettok();
			break;
		case GYR:
			printptr->p_value = (char *)
				&ptr->idunion.work3.gyr;
			printptr->p_type = ptr->id_type;
			if (tget)
				gettok();
			break;
		case STREAM:
			printptr->p_value1 = (char *)
				&ptr->idunion.ser.set->Etot;
			printptr->p_value2 = (char *)
				&ptr->idunion.ser.set->Epot;
			printptr->p_value3 = (char *)
				&ptr->idunion.ser.set->Ekin;
			printptr->p_value4 = (char *)
				&ptr->idunion.ser.set->Press;
			printptr->p_type = ptr->id_type;
			if (tget)
				gettok();
			break;
                case RMS:
                        printptr->p_type = ptr->id_type;
                        printptr->p_value = (char *)
                                &ptr->idunion.rms.rms;
                        if (tget)
                                gettok();
                        break;
		default:
			inerr("TABLE: wrong type of item","");
	}
	if (tabptr->idunion.tab.print == NULL) {
		tabptr->idunion.tab.print = printptr;
	} else {
		lastprint->next = printptr;
	}
	lastprint = printptr;
}

/************************************************************************
							WRITETAB()
*************************************************************************/

/*
 * writetab() - write out table entries
 */

void
writetab(idstruct *tabpt)
{
	printstruct	*print = tabpt->idunion.tab.print;
	int		i, nitems;
	_REAL		*result;

	tablereset(tabpt);
	for (; print != NULL; print = print->next) {
/*
		printf("ADDR= %d\n", print->p_value);
*/
		switch (print->p_type) {
		case BOND:
			REALtabwrite(tabpt, 
				* ((_REAL *) print->p_value) * 0.1 );
/*
			REALtabwrite(tabpt, 
				* ((_REAL *) print->p_value) * 0.01 );
*/
			break;
		case ANGLE:
		case DIHE:
			REALtabwrite(tabpt, 
				* ((_REAL *) print->p_value) * TODEG);
/*
			REALtabwrite(tabpt, 
				* ((_REAL *) print->p_value) );
*/
			break;
		case ORDERB:
			REALtabwrite(tabpt, 
				* ((_REAL *) print->p_value) );
			break;
		case GYR:
			REALtabwrite(tabpt, 
				* ((_REAL *) print->p_value) * 0.1 );
			break;
		case STREAM:
			REALtabwrite(tabpt, 
				* ((_REAL *) print->p_value1) );
			REALtabwrite(tabpt, 
				* ((_REAL *) print->p_value2) );
			REALtabwrite(tabpt, 
				* ((_REAL *) print->p_value3) );
			REALtabwrite(tabpt, 
				* ((_REAL *) print->p_value4) );
			break;
                case RMS:
                        REALtabwrite(tabpt, * ((_REAL *) print->p_value));
                        break;
		default:
			printf("Programming error: p_type %d", print->p_type);
			exit(1);
		}
	}
}

/************************************************************************
							REALTABWRITE()
*************************************************************************/
/*
 * REALtabwrite() - write a single floating point value to table, wrap line
 *	if necc
 */

void
REALtabwrite (idstruct *tabpt, _REAL value)
{
/*
	printf("X= %f\n",value);
*/

	/*
	 * start new line if necc, writing 1st char if wrap,
	 *	go on to next file if not wrap
	 */

	if (!(tabpt->idunion.tab.sofar % tabpt->idunion.tab.itemsline)) {
		if (tabpt->idunion.tab.wrap) {
/*
			fprintf(tabpt->idunion.tab.fptr, "\n%c ",
						tabpt->idunion.tab.linedelim);
*/
			fprintf(tabpt->idunion.tab.fptr, "\n");
			tabpt->idunion.tab.linedelim++; /* check range on setup */
		} else 
                {
			fprintf(tabpt->idunion.tab.fptr, "\n");
			tabpt->idunion.tab.curfile = 
					tabpt->idunion.tab.curfile->next;
			tabpt->idunion.tab.fptr = 
					tabpt->idunion.tab.curfile->f_file;
		}
	}

	/*
	 * write value & increment items sofar in the line
	 */

	fprintf(tabpt->idunion.tab.fptr, "%13.6e ", value);
	tabpt->idunion.tab.sofar++;
}

/************************************************************************
							TABLERESET()
*************************************************************************/

/*
 * tablereset() - get ready to start another logical line (row) in table
 */

void
tablereset(idstruct *tabpt)
{
	if (tabpt->idunion.tab.wrap) {
		tabpt->idunion.tab.linedelim = 'A';
		tabpt->idunion.tab.sofar = 0;
	} else {
		tabpt->idunion.tab.curfile = tabpt->idunion.tab.files;
		tabpt->idunion.tab.fptr = tabpt->idunion.tab.files->f_file;
	}
}

/************************************************************************
							CLOSETAB()
*************************************************************************/

/*
 * closetab() - print newline & close file
 */

void
closetab(idstruct *tabpt)
{
	fprintf(tabpt->idunion.tab.fptr, "\n");
	fclose(tabpt->idunion.tab.fptr);
}
