/*----------------------------------------------------------------------*
 * File:    radplex.c
 * Author:  Richard Leggett (richard.leggett@tgac.ac.uk)
 * Purpose: Demultiplex RADSeq
 * Created: 27 Nov 2012
 *----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h> 
#include <ctype.h>

/*----------------------------------------------------------------------*
 * Constants
 *----------------------------------------------------------------------*/
#define MAX_READ_LENGTH 10000

/*----------------------------------------------------------------------*
 * Structures
 *----------------------------------------------------------------------*/
typedef struct {
    char sequence_header[MAX_READ_LENGTH];
    char sequence[MAX_READ_LENGTH];
    char qualities_header[MAX_READ_LENGTH];
    char qualities[MAX_READ_LENGTH];
} FastqRead;

typedef struct {
    char* input_filename[2];
    FILE* input_fp[2];
    FastqRead read[2];
    int pairs_of_reads;
} FastqReadPair;

/*----------------------------------------------------------------------*
 * Globals
 *----------------------------------------------------------------------*/
int allowed_mismatches = 1;
int verbose = 0;

// Need to make this a more generic solution...
FILE* xmai_r1_fp = 0;
FILE* xmai_r2_fp = 0;
FILE* psti_r1_fp = 0;
FILE* psti_r2_fp = 0;
FILE* not_r1_fp = 0;
FILE* not_r2_fp = 0;
FILE* psti_p1_fp = 0;


/*----------------------------------------------------------------------*
 * Function:   chomp
 * Purpose:    Remove hidden characters from end of line
 * Parameters: str -> String to change
 * Returns:    None
 *----------------------------------------------------------------------*/
void chomp(char* str)
{
    int i = strlen(str) - 1;
    
    while ((i > 0) && (str[i] < ' ')) {
        str[i--] = 0;
    }    
}

/*----------------------------------------------------------------------*
 * Function:   usage
 * Purpose:    Report program usage.
 * Parameters: None
 * Returns:    None
 *----------------------------------------------------------------------*/
void usage(void)
{
    printf("\nradplex v0.1\n" \
           "richard.leggett@tgac.ac.uk\n" \
           "\nDemultiplex RADSeq runs.\n" \
           "\nOptions:\n" \
           "    [-h | --help] This help screen.\n" \
           "    [-a | --one] First FASTQ file.\n" \
           "    [-b | --two] Second (pair) FASTQ file.\n" \
           "    [-m | --mismatches] Number of allowed mismatches (default 1).\n" \
           "    [-v | --verbose] Verbose output.\n" \
           "\n");
}

/*----------------------------------------------------------------------*
 * Function:   parse_command_line
 * Purpose:    Parse command line options
 * Parameters: argc = number of arguments
 *             argv -> array of arguments
 * Returns:    None
 *----------------------------------------------------------------------*/
void parse_command_line(int argc, char* argv[], FastqReadPair* read_pair)
{
    static struct option long_options[] = {
        {"one", required_argument, NULL, 'a'},
        {"two", required_argument, NULL, 'b'},
        {"help", no_argument, NULL, 'h'},
        {"mismatches", required_argument, NULL, 'm'},
        {"verbose", no_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };
    int opt;
    int longopt_index;

    while ((opt = getopt_long(argc, argv, "a:b:hi:j:m:n:p:q:tv", long_options, &longopt_index)) > 0)
    {
        switch(opt) {
            case 'h':
                usage();
                exit(0);
                break;
            case 'a':
                if (optarg==NULL) {
                    printf("Error: Option requires an argument.\n");
                    exit(1);
                }
                read_pair->input_filename[0] = malloc(strlen(optarg)+1);
                strcpy(read_pair->input_filename[0], optarg);
                break;
            case 'b':
                if (optarg==NULL) {
                    printf("Error: Option requires an argument.\n");
                    exit(1);
                }
                read_pair->input_filename[1] = malloc(strlen(optarg)+1);
                strcpy(read_pair->input_filename[1], optarg);
                break;
            case 'm':
                if (optarg==NULL) {
                    printf("Error: Option requires an argument.\n");
                    exit(1);
                }
                allowed_mismatches=atoi(optarg);
                break;
            case 'v':
                verbose = 1;
                break;
        }
    }
    
    if ((read_pair->input_filename[0] == 0) || (read_pair->input_filename[1] == 0)) {
        printf("Error: you must specify both reads.\n");
        exit(2);
    }
}


/*----------------------------------------------------------------------*
 * Function:   
 * Purpose:    
 * Parameters: 
 * Returns:    
 *----------------------------------------------------------------------*/
int get_next_pair(FastqReadPair* read_pair)
{
    int i;
    
    for (i=0; i<2; i++) {
        if (!fgets(read_pair->read[i].sequence_header, MAX_READ_LENGTH, read_pair->input_fp[i])) {
            printf("Error reading input file\n");
            return 1;
        }
        if (!fgets(read_pair->read[i].sequence, MAX_READ_LENGTH, read_pair->input_fp[i])) {
            printf("Error reading input file\n");
            return 2;
        }
        if (!fgets(read_pair->read[i].qualities_header, 1024, read_pair->input_fp[i])) {
            printf("Error reading input file\n");
            return 3;
        }
        if (!fgets(read_pair->read[i].qualities, MAX_READ_LENGTH, read_pair->input_fp[i])) {
            printf("Error reading input file\n");
            return 4;                   
        }
        
        chomp(read_pair->read[i].sequence_header);
        chomp(read_pair->read[i].sequence);
        chomp(read_pair->read[i].qualities_header);
        chomp(read_pair->read[i].qualities);
    }

    read_pair->pairs_of_reads++;
    
    return 0;    
}

/*----------------------------------------------------------------------*
 * Function:   
 * Purpose:    
 * Parameters: 
 * Returns:    
 *----------------------------------------------------------------------*/
int compare_sequence(char* a, char* b, int l)
{
    int i=0;
    int differences=0;
    
    for (i=0; i<l; i++) {
        if (tolower(a[i]) != tolower(b[i])) {
            differences++;
        }
    }
    
    return differences;
}

/*----------------------------------------------------------------------*
 * Function:   
 * Purpose:    
 * Parameters: 
 * Returns:    
 *----------------------------------------------------------------------*/
void write_read(FastqRead* read, FILE* fp)
{
    fprintf(fp, "%s\n", read->sequence_header);
    fprintf(fp, "%s\n", read->sequence);
    fprintf(fp, "%s\n", read->qualities_header);
    fprintf(fp, "%s\n", read->qualities);
}

/*----------------------------------------------------------------------*
 * Function:   
 * Purpose:    
 * Parameters: 
 * Returns:    
 *----------------------------------------------------------------------*/
void check_current_read_for_adaptors(FastqReadPair* read_pair)
{
    int i;
    char p1[16];
    char p2[16];
    char temp[128];
    int m;
    int matched = 0;
    
    strncpy(p2, read_pair->read[1].sequence, 7);
    p2[7] = 0;
    
    m = compare_sequence(read_pair->read[0].sequence + 6, "CCGGG", 5);
    if (m <= allowed_mismatches) {
        strncpy(p1, read_pair->read[0].sequence, 6);
        p1[6]=0;
        printf("    Detected XmaI from base 7 with %d mismatches: %s-CCGGG p2 is %s\n", m, p1, p2);
        write_read(&read_pair->read[0], xmai_r1_fp);
        write_read(&read_pair->read[1], xmai_r2_fp);
        matched = 1;
    } else {
        int o;
        for (o=4; o<=7; o++) {
            m = compare_sequence(read_pair->read[0].sequence + o, "TGCAG", 5);
            if (m <= allowed_mismatches) {
                strncpy(p1, read_pair->read[0].sequence, o);
                printf("    Detected PstI from base %d with %d mismatches: %s-TGCAG p2 is %s\n", o+1, m, p1, p2);   
                fprintf(psti_p1_fp, "%s\n", p1);           
                write_read(&read_pair->read[0], psti_r1_fp);
                write_read(&read_pair->read[1], psti_r2_fp);
                matched = 1;
            }
        }
    }
    
    if (!matched) {
        write_read(&read_pair->read[0], not_r1_fp);
        write_read(&read_pair->read[1], not_r2_fp);
    }
}

/*----------------------------------------------------------------------*
 * Function:   
 * Purpose:    
 * Parameters: 
 * Returns:    
 *----------------------------------------------------------------------*/
void read_files(FastqReadPair* read_pair)
{
    int i;
    int rc = 0;

    for (i=0; i<2; i++) {
        read_pair->input_fp[i] = fopen(read_pair->input_filename[i], "r");
        if (!read_pair->input_fp[i]) {
            printf("Error: can't open %s\n", read_pair->input_filename[i]);
            exit(2);
        }
    }
    
    xmai_r1_fp = fopen("XmaI_R1.fastq", "w");
    xmai_r2_fp = fopen("XmaI_R2.fastq", "w");
    psti_r1_fp = fopen("PstI_R1.fastq", "w");
    psti_r2_fp = fopen("PstI_R2.fastq", "w");
    psti_p1_fp = fopen("PstI_p1.fastq", "w");
    not_r1_fp = fopen("NotRecognised_R1.fastq", "w");
    not_r2_fp = fopen("NotRecognised_R2.fastq", "w");

    if (!xmai_r1_fp || !xmai_r2_fp || !psti_r1_fp || !psti_r2_fp || !not_r1_fp || !not_r2_fp) {
        printf("Error: can't open output files\n");
        exit(2);
    }
    
            
    while ((!feof(read_pair->input_fp[0])) && (rc==0)) {
        rc - get_next_pair(read_pair);
        if (verbose) {
            printf("\nPair %d: %s\n", read_pair->pairs_of_reads, read_pair->read[0].sequence_header);
            printf("    Read 1: %s\n", read_pair->read[0].sequence);
            printf("    Read 2: %s\n", read_pair->read[1].sequence);
        }
        check_current_read_for_adaptors(read_pair);
    }
    
    for (i=0; i<2; i++) {
        fclose(read_pair->input_fp[i]);
    }
    
    fclose(xmai_r1_fp);
    fclose(xmai_r2_fp);
    fclose(psti_r1_fp);
    fclose(psti_r2_fp);
    fclose(psti_p1_fp);
    fclose(not_r1_fp);
    fclose(not_r2_fp);
}
                   
/*----------------------------------------------------------------------*
* Function:   
* Purpose:    
* Parameters: 
* Returns:    
*----------------------------------------------------------------------*/
void initialise_read_pair_struct(FastqReadPair* r)
{
    int i;
    r->pairs_of_reads = 0;
    for (i=0; i<2; i++) {
        r->input_filename[i] = 0;
        r->input_fp[i] = 0;
    }
}
            
/*----------------------------------------------------------------------*
 * Function:   main
 *----------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    FastqReadPair read_pair;
    
    initialise_read_pair_struct(&read_pair);
    parse_command_line(argc, argv, &read_pair);
    read_files(&read_pair);
    
    return 0;
}
