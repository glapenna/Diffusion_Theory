#include "head.h"

/************************************************************************/
/*                                                                      */
/*                      FUNCTION NAMD_read_line                         */
/*                                                                      */
/*   INPUTS:                                                            */
/*      fd - FILE to read line from                                     */
/*      buf - buffer to put line into                                   */
/*                                                                      */
/*   OUTPUTS:                                                           */
/*      this function returns a 0 if the line is read successfully or   */
/*   a -1 if an EOF is encountered                                      */
/*                                                                      */
/*      NAMD_read_line reads in a line from an open file.  It places    */
/*   the line that is read in into the buffer passed in via buf.  If    */
/*   an EOF is encountered, a -1 is returned.  If an EOF is encountered */
/*   in the middle of a line, the program will terminate abnormally.    */
/*   Also, any comments of the form {blah blah blah} will automatically */
/*   be skipped by NAMD_read_line, even if usch comments span lines.    */
/*   Also, the string will be left justified.  That is, any spaces at   */
/*   the begining of the line will be removed.                          */
/*                                                                      */
/************************************************************************/

int NAMD_read_line(FILE *fd, char *buf)

{
        int i=0;        /*  Current position in buf */
        int c;          /*  Character read in from file */

        /*  Loop and read characters until we get either an EOF or a    */
        /*  newline                                                     */
        while ( ((c=fgetc(fd)) != EOF) && (c != '\n') )
        {
                /*  If we encounter a bracketed comment, skip it.  This */
                /*  basically means read EVERYTHING until the next } and*/
                /*  throw it into the big bit bucket                    */
                if (c == '{')
                {
                        while ( ((c=fgetc(fd)) != EOF) && (c!='}') )
                        {
                        }

                        if (c==EOF)
                        {
                                char err_msg[512];
                                sprintf
                                 (err_msg, "ABNORMAL EOF FOUND - buffer=*%s*\n", buf);
                                printf("%s\n",err_msg);
                                exit(1);
                        }
                        continue;
                }

                /*  Also, use this little piece of logic to remove      */
                /*  any leading spaces from the line                    */
                if ((i>0) || !isspace(c))
                {
                        buf[i] = c;
                        i++;
                }
        }

        /*  NULL terminate the string                                   */
        buf[i]=(char )NULL;

        /*  Check for an EOF in the middle of a line                    */
        if ((c==EOF) && (i!=0))
        {
                char err_msg[512];
                sprintf(err_msg, "ABNORMAL EOF FOUND - buffer=*%s*\n",
                   buf);
                printf("%s\n",err_msg);
                exit(1);
        }

        /*  Return the appropriate value                                */
        if (c==EOF)
                return(-1);
        else
                return(0);
}
/*                      END OF FUNCTION NAMD_read_line                  */

/************************************************************************/
/*                                                                      */
/*                      FUNCTION NAMD_find_first_word                   */
/*                                                                      */
/*   INPUTS:                                                            */
/*      source - string to obtain word from                             */
/*      word - buffer to place word into                                */
/*                                                                      */
/*   OUTPUTS:                                                           */
/*      word is returned containing the first word of source            */
/*                                                                      */
/*      This function finds the first word in a string.  The first word */
/*   is defined to be the first set of continuous non-space charaters   */
/*   in a string.  So in the string "   AB14^  FDGFD GFDG"  the first   */
/*   word would be "AB14^".  The word is returned in the string pointed */
/*   to by word.                                                        */
/*                                                                      */
/************************************************************************/

void NAMD_find_first_word(char *source, char *word)

{
        int i=0;        /*  Position within source */
        int j=0;        /*  Position within word */

        /*  Skip leading white space                                    */
        while ( (source[i] != (char )NULL) && isspace(source[i]))
                i++;

        /*  Copy the word                                               */
        while ( (source[i] != (char )NULL) && !isspace(source[i]))
        {
                word[j]=source[i];
                i++;
                j++;
        }

        word[j]=(char )NULL;

        return;
}
/*                      END OF FUNCTION NAMD_find_first_word            */

/************************************************************************/
/*                                                                      */
/*                      FUNCTION NAMD_blank_str                         */
/*                                                                      */
/*   INPUTS:                                                            */
/*      str - string to test                                            */
/*                                                                      */
/*   OUTPUTS:                                                           */
/*      a 1 is returned if the string is blank, otherwise a 0 is        */
/*   returned                                                           */
/*                                                                      */
/*      NAMD_blank_str tests to see if a string is blank.  That is,     */
/*   contains only characters where isspace() is true                   */
/*                                                                      */
/************************************************************************/

int NAMD_blank_string(char *str)
{
        int i;          /* Current position in str */
        int blank=1;    /* Flag 1-> string is blank */
        int len;        /* length of the string str */

        len=strlen(str);

        for (i=0; i<len && blank; i++)
        {
                blank = isspace(str[i]);
        }

        return(blank);
}
/*                      END OF FUNCTION NAMD_blank_string               */
