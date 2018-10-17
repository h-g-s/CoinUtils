#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <ctype.h>
#include "macros.h"

char charIsInvisible( const char c )
{
    return ( (c==' ') || (c=='\t') || (c=='\r') || (c=='\n') );
}

char* applyInversion(char* str)
{
    int head, tail;
    for ( head=0,tail=(strlen(str)-1); head<tail; head++,tail-- ) {
        char t;
        t = str[head];
        str[head] = str[tail];
        str[tail] = t;
    }
    return str;
}

char *getFileName(char *destiny, const char *fileWithPath)
{
    int i,pos;
    /* Returning till found a slash */
    for ( i=(strlen(fileWithPath)-1),pos=0; (i>=0); i--,pos++ ) {
        if (fileWithPath[i]=='/')
            break;
        destiny[pos] = fileWithPath[i];
    }
    destiny[pos]='\0';
    applyInversion(destiny);
    /* Now, removing the . */
    const int endIdx = strlen(destiny)-1;
    for ( i=endIdx; (i>=0); i-- ) {
        if (destiny[i]=='.')
            break;
    }

    return destiny;
}

char *strRemoveSpsSol( char *dest, const char *str )
{
    char *startDest = dest;
    const char *send = str + strlen( str );
    const char *s = str;
    while ( (charIsInvisible(*s)) && (s<send) )
        s++;
    while ( (s<send) ) {
        *dest = *s;
        ++dest;
        ++s;
    }
    *dest = '\0';

    return startDest;
}


char* getParamName(char* target, const char* str)
{
    unsigned int i;
    unsigned int size;
    unsigned int pos;
    pos = 0;
    size = strlen(str);
    for ( i=0; i<size; i++ ) {
        if (str[i] != '-') {
            if (str[i] != '=')
                target[pos++] = str[i];
            else
                break;
        }
    }
    target[pos] = '\0';
    return target;
}

int digitsInLine( const char *str, int str_s )
{
    int result = 0;

    {
        int i;
        for ( i=0 ; (i<str_s) ; ++i ) {
            switch (str[i]) {
                case '\0' :
                case '\n' :
                    return result;
                    break;
                default:
                    if ( isdigit( str[i] ) )
                        result++;
                    break;
            };
        }
    }

    return result;
}

char* getParamValue(char* target, const char* str)
{
    unsigned int i;
    unsigned int size;
    size = strlen(str);
    unsigned int destSize = 0;
    for ( i=0; i<size; i++) {
        if (str[i]=='=')
            break;
    }
    i++;
    for (; i<size; i++) {
        target[destSize] = str[i];
        destSize++;
    }
    target[destSize] = 0;
    return target;
}

char *strRemoveSpsEol( char *str )
{
    char *p;

    p = str + strlen(str)-1;
    while ( p>=str ) {
        if (charIsInvisible(*p))
            *p = '\0';
        else
            break;

        --p;
    }

    return str;
}

void strAllToUpper( char *str )
{
    int l = strlen( str );
    {
        int i;
        for ( i=0; (i<l); ++i )
            str[i] = toupper(str[i]);
    }

}

void strFillSpacesLeft( char* dest, const char* str, int n )
{
    int len = strlen( str );
    int toFill = n-len;
    toFill = MAX( 0, toFill );
    if (toFill)
        memset( dest,' ', sizeof(char)*toFill );
    strcpy( dest+toFill, str );
}

void strFillSpacesRight(char* dest, const char* str, int n)
{
    int len = strlen( str );
    int toFill = n-len;
    toFill = MAX( 0, toFill );
    strcpy( dest, str );
    if (toFill)
        memset( dest+len,' ', sizeof(char)*toFill );
    *(dest+len+toFill) = '\0';
}

void strFillSpacesBoth( char *dest, const char *str, int n )
{
    strFillSpacesLeft( dest, str, n/2);
    strFillSpacesRight( dest, str, n/2);
}

char *strClearSpaces(char *str)
{
    char *c;
    char *cNoSpaces;
    c = str;
    cNoSpaces = str;
    while (*c != '\0') {
        if ( (*c != ' ') && (*cNoSpaces!=*c) ) {
            *cNoSpaces = *c;
            cNoSpaces++;
        }
        c++;
    }
    *cNoSpaces = '\0';
    return str;
}

char *strRemoveDblSpaces( char *str )
{
    /* removing spaces at extremes */
    char dest[LINE_SIZE];
    strRemoveSpsSol( dest, str );
    strRemoveSpsEol( dest );

    strcpy( str, dest );
    if (strlen(str)<=3)   /* no double spaces possible  A A */
        return str;

    char *s, *v, *prev;

    /* transforming all intermediate
    invisible chars in spaces */
    for ( s=str+1; (*s!='\0'); ++s )
        if (charIsInvisible(*s))
            *s = ' ';

    v    = str+2,
    s    = str+2,
    prev = str+1;

SKIP_CHAR:
    /* jump just one if valid is a repeated space or tab */
    if ( *prev==' ' )
        while ( *v==' ' )
            ++v;

    if ( *v == '\0' )
        goto CLOSE_STR;

    prev = v;

    if ( v!=s )
        *s = *v;

    ++v;
    ++s;

    if ( *v != '\0' )
        goto SKIP_CHAR;
CLOSE_STR:
    *s = '\0';

    return str;
}

int countChar( const char *str, const char c )
{
    int res = 0;

    const char *s = str;

    while ( *s!='\0' ) {
        res += (int) ( *s == c );
        ++s;
    }

    return res;
}

int firstOccurrence( const char *str, const char c )
{
    const char *s = str;
    while ( *s != '\0' ) {
        if ( *s == c )
            return s - str;

        ++s;
    }

    return INT_MAX;
}

int lastOccurrence( const char *str, const char c )
{
    int last = INT_MAX;

    const char *s = str;
    while ( *s != '\0' ) {
        if ( *s == c )
            last = s-str;

        ++s;
    }

    return last;
}


int splitString( char **columns, const char *str, const char delimiter,
                 const int maxColumns, const int columnSize, const char multDel )
{
    int sizeColumn;
    int ncolumn = 0;
    const char *send = str + strlen(str);
    const char *s = str;
    if (str[0] == '\0')
        return 0;
    const char *ns = s;
PROCESS_COLUMN:
    if ( ncolumn+1 == maxColumns )
        return ncolumn;

    /* finds the next delimiter */
FIND_DELIMITER:
    if ( ns == send )
        goto FOUND_COLUMN;
    if ( *ns != delimiter ) {
        ns++;
        goto FIND_DELIMITER;
    }
FOUND_COLUMN:
    sizeColumn = ns - s;
    if ((!multDel)||(sizeColumn>0)) {
        if (sizeColumn)
            memcpy( columns[ncolumn], s, sizeColumn );
        columns[ncolumn][sizeColumn] = '\0';
        ncolumn++;
    }
    if ( ns == send )
        return ncolumn;
    ++ns;
    s = ns;
    if ( ns != send )
        goto PROCESS_COLUMN;

    return ncolumn;
}
