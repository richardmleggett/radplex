/*----------------------------------------------------------------------*
 * File:    radplex.c
 * Author:  Richard Leggett (richard.leggett@tgac.ac.uk)
 * Purpose: Demultiplex RADSeq
 * Created: 27 Nov 2012
 * Updated:  7 Oct 2013
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
#define MAX_ADAPTORS 100
#define MAX_PATH_LENGTH 1024

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
    char* input_filename[3];
    FILE* input_fp[3];
    FastqRead read[3];
    int pairs_of_reads;
} FastqReadPair;

/*----------------------------------------------------------------------*
 * Globals
 *----------------------------------------------------------------------*/
int allowed_mismatches = 1;
int verbose = 0;
char adaptor_filename[2][MAX_PATH_LENGTH];
char output_prefix[MAX_PATH_LENGTH];
char* adaptors[2][MAX_ADAPTORS];
int n_adaptors[2];
FILE* undetermined_fp[2];
FILE* out_fp[MAX_ADAPTORS][MAX_ADAPTORS][2];

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
    printf("Demultiplex RADSeq runs.\n" \
           "\nOptions:\n" \
           "    [-h | --help] This help screen.\n" \
           "    [-a | --one] FASTQ R1.\n" \
           "    [-b | --two] FASTQ R2.\n" \
           "    [-c | --index] FASTQ index read.\n" \
           "    [-m | --mismatches] Number of allowed mismatches (default 1).\n" \
           "    [-p | --output_prefix] Output filename prefix.\n" \
           "    [-v | --verbose] Verbose output.\n" \
           "    [-1 | --p1] p1 Adaptor file.\n" \
           "    [-2 | --p2] p2 Adaptor file.\n" \
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
        {"index", required_argument, NULL, 'c'},
        {"help", no_argument, NULL, 'h'},
        {"mismatches", required_argument, NULL, 'm'},
        {"output_prefix", required_argument, NULL, 'p'},
        {"verbose", no_argument, NULL, 'v'},
        {"p1", required_argument, NULL, '1'},
        {"p2", required_argument, NULL, '2'},
        {0, 0, 0, 0}
    };
    int opt;
    int longopt_index;

    while ((opt = getopt_long(argc, argv, "a:b:c:hm:p:v1:2:", long_options, &longopt_index)) > 0)
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
            case 'c':
                if (optarg==NULL) {
                    printf("Error: Option requires an argument.\n");
                    exit(1);
                }
                read_pair->input_filename[2] = malloc(strlen(optarg)+1);
                strcpy(read_pair->input_filename[2], optarg);
                break;
            case 'm':
                if (optarg==NULL) {
                    printf("Error: Option requires an argument.\n");
                    exit(1);
                }
                allowed_mismatches=atoi(optarg);
                break;
            case 'p':
                if (optarg==NULL) {
                    printf("Error: Option requires an argument.\n");
                    exit(1);
                }
                strcpy(output_prefix, optarg);
                break;
            case 'v':
                verbose = 1;
                break;
            case '1':
                if (optarg==NULL) {
                    printf("Error: Option requires an argument.\n");
                    exit(1);
                }
                strcpy(adaptor_filename[0], optarg);
                break;
            case '2':
                if (optarg==NULL) {
                    printf("Error: Option requires an argument.\n");
                    exit(1);
                }
                strcpy(adaptor_filename[1], optarg);
                break;
        }
    }
    
    if ((read_pair->input_filename[0] == 0) || (read_pair->input_filename[1] == 0) || (read_pair->input_filename[2] == 0)) {
        printf("Error: you must specify both reads.\n");
        exit(2);
    }
    
    if ((adaptor_filename[0][0] == 0) || (adaptor_filename[1][0] == 0)) {
        printf("Error: you must specify adaptor files.\n");
        exit(3);
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
    
    for (i=0; i<3; i++) {
        if (!fgets(read_pair->read[i].sequence_header, MAX_READ_LENGTH, read_pair->input_fp[i])) {
            printf("End of file\n");
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
int match_adaptor(char*seq, int n)
{
    int i;
    int index = -1;
    
    for (i=0; i<n_adaptors[n]; i++) {
        if (compare_sequence(seq, adaptors[n][i], strlen(adaptors[n][i])) <=allowed_mismatches) {
            index = i;
            break;
        }
    }
    
    return index;
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
    int o;
    int matched = 0;
    int p1_index = -1;
    int p2_index = -1;
    FILE* out_r1 = undetermined_fp[0];
    FILE* out_r2 = undetermined_fp[1];
    
    // Get p2 from index read
    strncpy(p2, read_pair->read[2].sequence, 7);
    p2[7] = 0;
    
    p2_index = match_adaptor(p2, 1);
    
    // Deprecated XmaI detection
    //m = compare_sequence(read_pair->read[0].sequence + 6, "CCGGG", 5);
    //if (m <= allowed_mismatches) {
    //    strncpy(p1, read_pair->read[0].sequence, 6);
    //    p1[6]=0;
    //    printf("    Detected XmaI from base 7 with %d mismatches: %s-CCGGG p2 is %s\n", m, p1, p2);
    //    write_read(&read_pair->read[0], xmai_r1_fp);
    //    write_read(&read_pair->read[1], xmai_r2_fp);
    //    matched = 1;
    //} else {

    for (o=4; o<=7; o++) {
        m = compare_sequence(read_pair->read[0].sequence + o, "TGCAG", 5);
        if (m <= allowed_mismatches) {
            strncpy(p1, read_pair->read[0].sequence, o);
            p1[o] = 0;
            p1_index = match_adaptor(p1, 0);
            //printf("    Detected PstI from base %d with %d mismatches: %s-TGCAG p2 is %s\n", o+1, m, p1, p2);
            matched = 1;
        }
    }
    
    if ((matched) && (p1_index >=0) && (p2_index >=0)) {
        //printf("p1=%s (%d)\tp2=%s (%d)\n", p1, p1_index, p2, p2_index);
        if (out_fp[p1_index][p2_index][0] == 0) {
            int i;
            for (i=0; i<2; i++) {
                char filename[MAX_PATH_LENGTH];
                sprintf(filename, "%s_%c%d_R%d.fastq", output_prefix, p2_index+'A', p1_index+1, i);
                out_fp[p1_index][p2_index][i] = fopen(filename, "w");
                if (!out_fp[p1_index][p2_index][i]) {
                    printf("Can't open %s\n", filename);
                    exit(6);
                } else {
                    printf("Created %s\n", filename);
                }
            }
        }
        out_r1 = out_fp[p1_index][p2_index][0];
        out_r2 = out_fp[p1_index][p2_index][1];
    } else {
        //printf("No match\n");
        out_r1 = undetermined_fp[0];
        out_r2 = undetermined_fp[1];
    }

    write_read(&read_pair->read[0], out_r1);
    write_read(&read_pair->read[1], out_r2);
}

/*----------------------------------------------------------------------*
 * Function:   
 * Purpose:    
 * Parameters: 
 * Returns:    
 *----------------------------------------------------------------------*/
void read_files(FastqReadPair* read_pair)
{
    int i, j, k;
    int rc = 0;
    char filename[MAX_PATH_LENGTH];

    // Clear output file handles
    for (i=0; i<MAX_ADAPTORS; i++) {
        for (j=0; j<MAX_ADAPTORS; j++) {
            for (k=0; k<2; k++) {
                out_fp[i][j][k] = 0;
            }
        }
    }
    
    for (i=0; i<2; i++) {
        sprintf(filename, "%s_undetermined_R%d.fastq", output_prefix, i+1);
        undetermined_fp[i] = fopen(filename, "w");
        if (!undetermined_fp[i]) {
            printf("Error: Can't open %s\n", filename);
            exit(5);
        }
    }
    
    for (i=0; i<3; i++) {
        read_pair->input_fp[i] = fopen(read_pair->input_filename[i], "r");
        if (!read_pair->input_fp[i]) {
            printf("Error: can't open %s\n", read_pair->input_filename[i]);
            exit(2);
        }
    }
    
    while ((!feof(read_pair->input_fp[0])) && (rc==0)) {
        rc = get_next_pair(read_pair);
        if (verbose) {
            printf("\nPair %d: %s\n", read_pair->pairs_of_reads, read_pair->read[0].sequence_header);
            printf("    Read 1: %s\n", read_pair->read[0].sequence);
            printf("    Read 2: %s\n", read_pair->read[1].sequence);
        }
        check_current_read_for_adaptors(read_pair);
    }
    
    for (i=0; i<3; i++) {
        fclose(read_pair->input_fp[i]);
    }
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
void load_adaptor_files(void)
{
    int i;
    
    for (i=0; i<2; i++) {
        FILE* fp = fopen(adaptor_filename[i], "r");
        char string[1024];
        
        n_adaptors[i] = 0;
        if (fp) {
            printf("Reading P%d adaptors...\n", i+1);
            while (!feof(fp)) {
                if (fgets(string, 1024, fp)) {
                    chomp(string);
                    if (strlen(string) > 1) {
                        adaptors[i][n_adaptors[i]] = calloc(strlen(string)+1, sizeof(char));
                        if (adaptors[i][n_adaptors[i]] == 0) {
                            printf("Error: can't allocate memory.\n");
                            exit(5);
                        }
                        strcpy(adaptors[i][n_adaptors[i]], string);
                        if (i == 0) {
                            printf("%d. %s\n", n_adaptors[i], adaptors[i][n_adaptors[i]]);
                        } else {
                            printf("%c. %s\n", n_adaptors[i]+'A', adaptors[i][n_adaptors[i]]);
                        }
                        n_adaptors[i]++;
                    }
                }
            }
        } else {
            printf("Error: Can't open %s\n", adaptor_filename[i]);
            exit(4);
        }
    }
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
    for (i=0; i<3; i++) {
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
    
    printf("\nRADplex v0.2\n\n");
    
    adaptor_filename[0][0] = 0;
    adaptor_filename[1][0] = 0;
    strcpy(output_prefix, "RADplex_output");
    
    initialise_read_pair_struct(&read_pair);
    parse_command_line(argc, argv, &read_pair);
    load_adaptor_files();
    read_files(&read_pair);
    
    printf("\nDone.\n");
    
    return 0;
}
