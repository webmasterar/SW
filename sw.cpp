
/**
    Smith-Waterman
    Copyright (C) 2015 Solon P. Pissis, Ahmad Retha, Fatima Vayani 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include "sw.h"

unsigned int EDNA[90];
unsigned int BLOSUM[91];

/**
 * Initialises the substitution scoring tables (EDNA and BLOSUM)
 */
void init_substitution_score_tables ( void )
{
    unsigned int i;
    for ( i = 0; i < strlen ( ALPHABET_DNA ) - 1; i ++ ) {
	    EDNA[(int)ALPHABET_DNA[i]] = i;
    }
    EDNA[(int)'N'] = 14; //Setting R/DNA N char
    EDNA[(int)'U'] = 1;  //Setting RNA U = T
    for ( i = 0; i < strlen ( ALPHABET_PROT ); i ++ ) {
	    BLOSUM[(int)ALPHABET_PROT[i]] = i;
    }
}

/**
 * Smith-Waterman algorithm with Traceback
 * @param a score matrix to use: 0 = DNA or RNA, 1 = PROT
 * @param p doubled-up p cstring
 * @param m length of x
 * @param t second cstring
 * @param n length of t
 * @param o open gap penalty
 * @param e extend gap penalty
 * @param score the score of the match
 * @param ii the i position to rotate at
 * @param jj the j position to rotate at
 */
unsigned int sw_tb ( int a, unsigned char * p, unsigned int m, unsigned char * t, unsigned int n, double o, double e, double * score, unsigned int * ii, unsigned int * jj )
{
	init_substitution_score_tables();
	int i;
	int j;
	double g = o;
	double h = e;
	double max_score = -DBL_MAX;
	double u, v, w;
	double matching_score = 0;
	const char UP   = 'U';
	const char DIAG = 'D';
	const char LEFT = 'L';
	const char ZERO = 'Z';

	//init matrices - T = Scores, I = Insertion, D = Deletion, TB = Traceback
	double ** T;
	double ** I;
	double ** D;
	char   ** TB;
	if ( ( T = ( double ** ) calloc ( ( m + 1 ) , sizeof( double * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: T could not be allocated!\n");
                return ( 0 );
        }
        for ( i = 0; i < m + 1; i ++ )
        {
                if ( ( T[i] = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: T could not be allocated!\n");
                        return ( 0 );
                }
        }
	if ( ( I = ( double ** ) calloc ( ( m + 1 ) , sizeof( double * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: I could not be allocated!\n");
                return ( 0 );
        }
        for ( i = 0; i < m + 1; i ++ )
        {
                if ( ( I[i] = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: I could not be allocated!\n");
                        return ( 0 );
                }
        }
	if ( ( D = ( double ** ) calloc ( ( m + 1 ) , sizeof( double * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: D could not be allocated!\n");
                return ( 0 );
        }
        for ( i = 0; i < m + 1; i ++ )
        {
                if ( ( D[i] = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: D could not be allocated!\n");
                        return ( 0 );
                }
        }
        if ( ( TB = ( char ** ) calloc ( ( m + 1 ) , sizeof( char * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: TB could not be allocated!\n");
                return ( 0 );
        }
        for ( i = 0; i < m + 1; i ++ )
        {
                if ( ( TB[i] = ( char * ) calloc ( ( n + 1 ) , sizeof( char ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: TB could not be allocated!\n");
                        return ( 0 );
                }
        }
        for ( i = 0; i < m + 1; i++ ) TB[i][0] = ZERO;
	for ( j = 0; j < n + 1; j++ ) TB[0][j] = ZERO;

        //Do DP
	for( i = 1; i < m + 1; i++ )
	{
		//max_score = 0; //@todo figure out why he zeros this?!!!!!

        	for( j = 1; j < n + 1; j++ )
        	{
			D[i][j] = max ( D[i - 1][j] + h, T[i - 1][j] + g );
			u = D[i][j];

			I[i][j] = max ( I[i][j - 1] + h, T[i][j - 1] + g );
			v = I[i][j];

			matching_score = T[i - 1][j - 1] + ( (a == 0) ? nuc_delta( t[j - 1], p[i - 1] ) : prot_delta( t[j - 1], p[i - 1] ) );
			w = matching_score;

			T[i][j] = max ( 0.0, max ( w, max ( u, v ) ) );
			
			if ( T[i][j] == 0.0 )
			{
			    TB[i][j] = ZERO;
			}
			else if ( T[i][j] == w )
			{
			    TB[i][j] = DIAG;
			}
			else if ( T[i][j] == u )
			{
			    TB[i][j] = UP;
			}
			else
			{
			    TB[i][j] = LEFT;
			}

			if ( T[i][j] > max_score )
			{
				max_score = T[i][j];
				( * score ) = max_score;
				( * ii ) = i;
				( * jj ) = j;
			}
        	}
    	}

    	//Do Traceback
    	char d = 1;
    	i = ( * ii );
	j = ( * jj );
	while ( d != ZERO )
	{
	    d = TB[i][j];
	    if ( d == DIAG )
	    {
		i--;
		j--;
	    }
	    else if ( d == LEFT )
	    {
		j--;
	    }
	    else if ( d == UP )
	    {
		i--;
	    }
	}
	( * ii ) = i;
	( * jj ) = j;

        for ( i = 0; i < m + 1; i ++ )
	{
		free ( D[i] );
		free ( I[i] );
		free ( T[i] );
		free ( TB[i] );
	}
	free ( I );
	free ( D );
	free ( T );
	free ( TB );
	
	return EXIT_SUCCESS;
}

/**
 * Smith-Waterman algorithm (uses linear space)
 * @param a score matrix to use: 0 = DNA or RNA, 1 = PROT
 * @param p doubled-up p cstring
 * @param m length of x
 * @param t second cstring
 * @param n length of t
 * @param o open gap penalty
 * @param e extend gap penalty
 * @param score the score of the match
 */
unsigned int sw_ls ( int a, unsigned char * p, unsigned int m, unsigned  char * t, unsigned int n, double o, double e, double * score )
{
	init_substitution_score_tables();
	int i, j;
	double g = o;
	double h = e;
	double max_score = -DBL_MAX;
	double u, v, w;

        double * d0;
        double * d1;
        double * t0;
        double * t1;
        double * in;
        if ( ( d0 = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL )
        {
            fprintf( stderr, " Error: 'd0' could not be allocated!\n");
            return ( 0 );
        }
        if ( ( d1 = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL  )
        {
            fprintf( stderr, " Error: 'd1' could not be allocated!\n");
            return ( 0 );
        }
        if ( ( t0 = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL )
        {
            fprintf( stderr, " Error: 't0' could not be allocated!\n");
            return ( 0 );
        }
        if ( ( t1 = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL )
        {
            fprintf( stderr, " Error: 't1' could not be allocated!\n");
            return ( 0 );
        }
        if ( ( in = ( double * ) calloc ( ( n + 1 ) , sizeof ( double ) ) ) == NULL )
        {
            fprintf( stderr, " Error: 'in' could not be allocated!\n");
            return ( 0 );
        }

    	for( i = 1; i < m + 1; i++ )
    	{
        	for( j = 1; j < n + 1; j++ )
        	{

                    switch ( i % 2 ) {

                        case 0:

                            u = d0[j] = max ( d1[j] + h, t1[j] + g );

                            v = in[j] = max ( in[j - 1] + h, t0[j - 1] + g ); //i0

                            w = t1[j - 1] + ((a == 0) ? nuc_delta( t[j - 1], p[i - 1] ) : prot_delta( t[j - 1], p[i - 1] ));

                            t0[j] = max ( 0.0, max ( w, max ( u, v ) ) );

                            if ( t0[j] > max_score )
                            {
                                    max_score = t0[j];
                            }

                            break;

                        case 1:

                            u = d1[j] = max ( d0[j] + h, t0[j] + g );

                            v = in[j] = max ( in[j - 1] + h, t1[j - 1] + g ); //i1

                            w = t0[j - 1] + ((a == 0) ? nuc_delta( t[j - 1], p[i - 1] ) : prot_delta( t[j - 1], p[i - 1] ));

                            t1[j] = max ( 0.0, max ( w, max ( u, v ) ) );

                            if ( t1[j] > max_score )
                            {
                                    max_score = t1[j];
                            }                   

                            break;

                    }
	      }
    	}

        free( d0 );
        free( d1 );
        free( t0 );
        free( t1 );
        free( in );

	( * score ) = max_score;

        return EXIT_SUCCESS;
}

/**
 * Decode the input switches 
 */
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
{
      int          oi;
      int          opt;
      double       val;
      int          args;

      /* initialisation */
      sw -> alphabet                       = NULL;
      sw -> method                         = NULL;
      sw -> input_filename                 = NULL;
      sw -> open_gap_penalty               = -10.0;
      sw -> extend_gap_penalty             = -0.5;
      args = 0;

      while ( ( opt = getopt_long ( argc, argv, "a:m:i:O:E:h", long_options, &oi ) ) != - 1 )
      {
	switch ( opt )
	  {
	    case 'a':
	      sw -> alphabet = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
	      strcpy ( sw -> alphabet, optarg );
	      args ++;
	      break;

	    case 'm':
	      sw -> method = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
	      strcpy ( sw -> method, optarg );
	      break;
	      
	    case 'i':
	      sw -> input_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
	      strcpy ( sw -> input_filename, optarg );
	      args ++;
	      break;

	    case 'O':
	      val = atof ( optarg );
	      sw -> open_gap_penalty = -abs(val);
	      args ++;
	      break;

	    case 'E':
	      val = atof ( optarg );
	      sw -> extend_gap_penalty = -abs(val);
	      args ++;
	      break;

	    case 'h':
	      return ( 0 );
	  }
      }
      
      if ( sw -> method == NULL )
      {
	  sw -> method = (char * ) METHOD_LS;
      }

      if ( args < 4 )
	{
	  usage ();
	  exit ( 1 );
	}
      else
	return ( optind );
}

/** 
 * Usage of the tool 
 */
void usage ( void )
{
	fprintf ( stdout, " Usage: Smith-Waterman <options>\n" );
	fprintf ( stdout, " Standard (Mandatory):\n" );
	fprintf ( stdout, "  -a, --alphabet            <str>     Alphabet to use (DNA, RNA or PROT).\n" );
	fprintf ( stdout, "  -m, --method              <str>     Method to use - LS=linear space, TB=traceback. (Default: LS)\n" );
	fprintf ( stdout, "  -i, --input-file          <str>     (Multi)FASTA input filename.\n" );
	fprintf ( stdout, "  -O, --open-gap-penalty    <float>   The cost of opening a gap (e.g. -10.0).\n");
	fprintf ( stdout, "  -E, --extend-gap-penalty  <float>   The cost of extending a gap (e.g. -0.5).\n");
}

/**
 * Time
 */
double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
};

/**
 * Main
 */
int main(int argc, char **argv)
{

	struct TSwitch  sw;

	FILE *           in_fd;                  // the input file descriptor
	FILE *           out_fd;                 // the input file descriptor
	char *           alphabet;               // alphabet (DNA, RNA or PROT)
	char *           alphabet_letters;       // alphabet letters
	char *           method;                 // method, either LS or TB
        char *           input_filename;         // the input file name
        char *           output_filename;        // the output file name
        unsigned char ** seq    = NULL;          // the sequence in memory
        unsigned char ** seq_id = NULL;          // the sequence id in memory
        double           open_gap_penalty;       // open gap penalty
	double           extend_gap_penalty;     // extend gap penalty
	unsigned int     num_args;               // num args
	unsigned int     a;                      // matrix to use, 0 = R/DNA, 1 = PROT

	num_args = decode_switches ( argc, argv, &sw );

	if ( num_args < 4 ) {
	    usage();
	    return EXIT_FAILURE;
	} else {
		alphabet                = sw . alphabet;
		method                  = sw . method;
		open_gap_penalty        = sw . open_gap_penalty;
		extend_gap_penalty      = sw . extend_gap_penalty;
                input_filename          = sw . input_filename;
	}

	if ( strcmp ( alphabet, "DNA" ) == 0 ) {
		a = 0;
		alphabet_letters = (char *) ALPHABET_DNA;
	} else if ( strcmp ( alphabet, "RNA" ) == 0 ) {
		a = 0;
		alphabet_letters = (char *) ALPHABET_RNA;
	} else if ( strcmp ( alphabet, "PROT" ) == 0 ) {
		a = 1;
		alphabet_letters = (char *) ALPHABET_PROT;
	} else {
		usage();
		return EXIT_FAILURE;
	}

        /* Read the (Multi)FASTA file in memory */
        fprintf ( stderr, " Reading the (Multi)FASTA input file: %s\n", input_filename );
        if ( ! ( in_fd = fopen ( input_filename, "r") ) )
        {
                fprintf ( stderr, " Error: Cannot open file %s!\n", input_filename );
                return ( 1 );
        }

        char c;
        unsigned int num_seqs = 0;           	// the total number of sequences considered
        unsigned int total_length = 0;          // the total number of sequences considered
        unsigned int max_alloc_seq_id = 0;
        unsigned int max_alloc_seq = 0;
        c = fgetc( in_fd );
        do
        {
                if ( c != '>' )
                {
                        fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", input_filename );
                        return ( 1 );
                }
                else
                {
                        if ( num_seqs >= max_alloc_seq_id )
                        {
                                seq_id = ( unsigned char ** ) realloc ( seq_id,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
                                max_alloc_seq_id += ALLOC_SIZE;
                        }

                        unsigned int max_alloc_seq_id_len = 0;
                        unsigned int seq_id_len = 0;

                        seq_id[ num_seqs ] = NULL;

                        while ( ( c = fgetc( in_fd ) ) != EOF && c != '\n' )
                        {
                                if ( seq_id_len >= max_alloc_seq_id_len )
                                {
                                        seq_id[ num_seqs ] = ( unsigned char * ) realloc ( seq_id[ num_seqs ],   ( max_alloc_seq_id_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
                                        max_alloc_seq_id_len += ALLOC_SIZE;
                                }
                                seq_id[ num_seqs ][ seq_id_len++ ] = c;
                        }
                        seq_id[ num_seqs ][ seq_id_len ] = '\0';

                }
		if ( num_seqs >= max_alloc_seq )
                {
                        seq = ( unsigned char ** ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
                        max_alloc_seq += ALLOC_SIZE;
                }

                unsigned int seq_len = 0;
                unsigned int max_alloc_seq_len = 0;

                seq[ num_seqs ] = NULL;

                while ( ( c = fgetc( in_fd ) ) != EOF && c != '>' )
                {
                        if( seq_len == 0 && c == '\n' )
                        {
                                fprintf ( stderr, " Omitting empty sequence in file %s!\n", input_filename );
                                c = fgetc( in_fd );
                                break;
                        }
                        if( c == '\n' || c == ' ' ) continue;

                        c = toupper( c );

                        if ( seq_len >= max_alloc_seq_len )
                        {
                                seq[ num_seqs ] = ( unsigned char * ) realloc ( seq[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
                                max_alloc_seq_len += ALLOC_SIZE;
                        }

                        if( strchr ( alphabet_letters, c ) )
                        {
                                seq[ num_seqs ][ seq_len++ ] = c;
                        }
                        else
                        {
                                fprintf ( stderr, " Error: input file %s contains an unexpected character %c!\n", input_filename, c );
                                return ( 1 );
                        }
                }

                if( seq_len != 0 )
                {
                        if ( seq_len >= max_alloc_seq_len )
                        {
                                seq[ num_seqs ] = ( unsigned char * ) realloc ( seq[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
                                max_alloc_seq_len += ALLOC_SIZE;
                        }
                        seq[ num_seqs ][ seq_len ] = '\0';
                        total_length += seq_len;
                        num_seqs++;
                }

        } while( c != EOF );

	if ( fclose ( in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	if ( num_seqs > 2 )
	{
        	fprintf( stderr, " Warning: %d sequences were read from file %s.\n", num_seqs, input_filename );
        	fprintf( stderr, " Warning: Only the first two (%s, %s) will be processed!\n", seq_id[0], seq_id[1] );
	}

	unsigned int m = strlen ( ( char * ) seq[0] );
	unsigned int n = strlen ( ( char * ) seq[1] );
	double distance = 0.0;
	unsigned int ii;
	unsigned int jj;
	unsigned int rot;

	/* Run the algorithm */
	double start = gettime();
	if ( strcmp ( method, METHOD_LS ) == 0 )
	{
		sw_ls ( a, seq[0], m, seq[1], n, open_gap_penalty, extend_gap_penalty, &distance );
	}
	else
	{
		sw_tb ( a, seq[0], m, seq[1], n, open_gap_penalty, extend_gap_penalty, &distance, &ii, &jj );
		/* If ii >= jj then it is easy to find the rotation: ii - jj */
		if ( ii >= jj )
		{
			rot = ii - jj;
		}
		/* Otherwise we would need to shift jj - ii or as many letters as we have on the right of position i - sw. w + 1 */
		else
		{
			int a = jj - ii;
			int b = m - ii - 1;
			rot = m - min ( a, b );
		}
	}
	double end = gettime();

        fprintf( stderr, " Max similarity score: %1.2f\n", distance );
	if ( strcmp ( method, METHOD_TB ) == 0 ) {
		fprintf ( stderr, " Rotation: %u\n", rot );
	}
        fprintf( stderr, " Elapsed time for comparing sequences: %lf secs\n", ( end - start ) );

	/* De-allocate */
	free ( sw . method );
	free ( sw . alphabet );
	free ( sw . input_filename );
	unsigned int i;
        for ( i = 0; i < num_seqs; i ++ )
        {
                free ( seq[i] );
                free ( seq_id[i] );
        }
        free ( seq );
        free ( seq_id );

	return EXIT_SUCCESS;
}
