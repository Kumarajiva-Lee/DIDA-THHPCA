//---------------------------
//by Shishupeng
//shupengshi@outlook.com
//---------------------------
#include <stdio.h>
#include <fcntl.h>

#define STANDARD_ERROR 2

#define STANDARD_OUTPUT 1

#include "mpi.h"

# ifdef NOUNDERSCORE
#      define RSL_LITE_ERROR_DUP1 rsl_error_dup1
#   else
#      define RSL_LITE_ERROR_DUP1 rsl_error_dup1_
# endif

int RSL_LITE_ERROR_DUP1 ( int *me )
{
    int newfd ;
    char filename[256] ;
    //char hostname[256] ;

    //gethostname( hostname, 256 );

/* redirect standard out*/
    sprintf(filename,"std.out.%04d",*me) ;
    if ((newfd = open( filename, O_CREAT | O_WRONLY, 0666 )) < 0 ) {
        perror("error_dup: cannot open std.out.nnnn") ;
        fprintf(stderr,"...sending output to standard output and continuing.\n") ;
        return 0;
    }
    if( dup2( newfd, STANDARD_OUTPUT ) < 0 ) {
        perror("error_dup: dup2 fails to change output descriptor") ;
        fprintf(stderr,"...sending output to standard output and continuing.\n") ;
        close(newfd) ;
        return 0;
    }

/* redirect standard error */
    sprintf(filename,"std.error.%04d",*me) ;
    if ((newfd = open( filename, O_CREAT | O_WRONLY, 0666 )) < 0 ) {
        perror("error_dup: cannot open std.error.log") ;
        fprintf(stderr,"...sending error to standard error and continuing.\n") ;
        return 0;
    }
    if( dup2( newfd, STANDARD_ERROR ) < 0 ) {
        perror("error_dup: dup2 fails to change error descriptor") ;
        fprintf(stderr,"...sending error to standard error and continuing.\n") ;
        close(newfd) ;
        return 0;
    }
    fprintf( stdout, "taskid: %d \n",*me) ;
    fprintf( stderr, "taskid: %d \n",*me) ;
    //fprintf( stdout, "taskid: %d hostname: %s\n",*me,hostname) ;
    //fprintf( stderr, "taskid: %d hostname: %s\n",*me,hostname) ;
}

