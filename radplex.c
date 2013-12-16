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
#include <math.h>

/*----------------------------------------------------------------------*
 * Constants
 *----------------------------------------------------------------------*/
#define MAX_READ_LENGTH 10000
#define MAX_ADAPTORS 100
#define MAX_PATH_LENGTH 1024
#define MAX_UNDETERMINED 1000000
#define MAX_HASH 7
#define UNDETERMINED_HASH_SIZE 279936

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
int adaptor_counts[MAX_ADAPTORS][MAX_ADAPTORS];
int undetermined_read_count = 0;
char undetermined_indices[2][UNDETERMINED_HASH_SIZE];
int total_read_count = 0;
int clip_psti = 0;

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
           "    [-z | --clip_psti] Clip PstI sequence too.\n" \
           "    [-1 | --p1] p1 Adaptor file.\n" \
           "    [-2 | --p2] p2 Adaptor file.\n" \
           "\n");
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
void initialise_main(void)
{
    int i,j;
    
    adaptor_filename[0][0] = 0;
    adaptor_filename[1][0] = 0;
    strcpy(output_prefix, "RADplex_output");

    for (i=0; i<2; i++) {
        for (j=0; j<UNDETERMINED_HASH_SIZE; j++) {
            undetermined_indices[i][j] = 0;
        }
    }
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
int base_to_n(char base)
{
    int n = -1;
    
    switch(base) {
        case 'A':
        case 'a':
            n=1;
            break;
        case 'C':
        case 'c':
            n=2;
            break;
        case 'G':
        case 'g':
            n=3;
            break;
        case 'T':
        case 't':
            n=4;
            break;
        case 'N':
        case 'n':
            n=5;
            break;
        default:
            printf("Error: Unknown base %c\n", base);
            exit(3);
            break;
    }
    
    return n;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
char n_to_base(int n)
{
    char base;
    
    switch (n) {
        case 1: base = 'A'; break;
        case 2: base = 'C'; break;
        case 3: base = 'G'; break;
        case 4: base = 'T'; break;
        case 5: base = 'N'; break;
        default: printf("Error: unknown n %d\n", n); exit(3); break;
    }
    
    return base;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
int generate_hash(char* sequence)
{
    int multiplier = 0;
    int hash = 0;
    int i, m;
    
    printf("Sequence: %s\n", sequence);
    
    if (!sequence) {
        printf("Error: Bad sequence\n");
        exit(1);
    }
    
    if ((strlen(sequence) == 0) || (strlen(sequence) > 7)) {
        printf("Error: Bad length for hash sequence\n");
        exit(1);
    }
    
    for (i=0; i<strlen(sequence); i++) {
        m = (int)pow(6.0, (double)(MAX_HASH-1-i));
        hash = hash + (base_to_n(sequence[i]) * m);
    }
    
    return hash;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
char* hash_to_string(int hash, char* hash_string)
{
    int i, n, m;
    int running = hash;
    int length = 0;
    
    for (i=0; i<MAX_HASH; i++) {
        m = (int)pow(6.0, (double)(MAX_HASH-1-i));
        if (running >= m) {
            n = (running / m);
            running -= (n * m);
            hash_string[i] = n_to_base(n);
        } else {
            hash_string[i] = 0;
        }
    }

    return hash_string;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
void store_undetermined(int p, char* index)
{
    if (index[0] != 0) {
        undetermined_indices[p][generate_hash(index)]++;
    }
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
void output_undetermined_indices(void)
{
    int i,j;
    char hash_string[8];
    FILE* fp;
    char filename[1024];
    
    
    for (i=0; i<2; i++) {
        sprintf(filename, "%s_p%d_undetermined_counts.txt", output_prefix, i+1);
        fp = fopen(filename, "w");
        if (fp) {
            for (j=0; j<UNDETERMINED_HASH_SIZE; j++) {
                if (undetermined_indices[i][j] > 0) {
                    hash_to_string(j, hash_string);
                    fprintf(fp, "%s\t%d\n", hash_string, undetermined_indices[i][j]);
                }
            }
        } else {
            printf("ERROR: Can't open %s\n", filename);
        }
    }
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
char* assign_string(char* s)
{
    char* r = malloc(strlen(s)+1);
    if (!r) {
        printf("Error: can't assign string %s\n",s);
        exit(12);
    }
    
    strcpy(r, s);
    return r;
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
void setup_default_adaptors(void)
{
    n_adaptors[0] = 12;
    adaptors[0][0]  = assign_string("TGAGTGCAG");
    adaptors[0][1]  = assign_string("ACGTATGCAG");
    adaptors[0][2]  = assign_string("CTCCGATGCAG");
    adaptors[0][3]  = assign_string("GATACCATGCAG");
    adaptors[0][4]  = assign_string("GGCATGCAG");
    adaptors[0][5]  = assign_string("CTAGGTGCAG");
    adaptors[0][6]  = assign_string("ACGCACTGCAG");
    adaptors[0][7]  = assign_string("TATTCAATGCAG");
    adaptors[0][8]  = assign_string("GTATTGCAG");
    adaptors[0][9]  = assign_string("TACGTTGCAG");
    adaptors[0][10] = assign_string("CCGCACTGCAG");
    adaptors[0][11] = assign_string("AGTAGAATGCAG");

    n_adaptors[1] = 8;
    adaptors[1][0] = assign_string("AATAGTT");
    adaptors[1][1] = assign_string("ACCGACC");
    adaptors[1][2] = assign_string("ATGGCAA");
    adaptors[1][3] = assign_string("CCGGTCG");
    adaptors[1][4] = assign_string("GACCTGG");
    adaptors[1][5] = assign_string("GTTCGGT");
    adaptors[1][6] = assign_string("TGAACTA");
    adaptors[1][7] = assign_string("TGATAAC");
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
void write_read(FastqRead* read, int trim_start, FILE* fp)
{
    fprintf(fp, "%s\n", read->sequence_header);
    fprintf(fp, "%s\n", (read->sequence) + trim_start);
    fprintf(fp, "%s\n", read->qualities_header);
    fprintf(fp, "%s\n", (read->qualities) + trim_start);
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
    int clip_size = 0;
    FILE* out_r1 = undetermined_fp[0];
    FILE* out_r2 = undetermined_fp[1];
    
    total_read_count++;
    
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
    //    write_read(&read_pair->read[0], 0, xmai_r1_fp);
    //    write_read(&read_pair->read[1], 0, xmai_r2_fp);
    //    matched = 1;
    //} else {

    for (i=0; i<n_adaptors[0]; i++) {
        if (compare_sequence(read_pair->read[0].sequence, adaptors[0][i], strlen(adaptors[0][i])) <= allowed_mismatches) {
            p1_index = i;
            strncpy(p1, read_pair->read[0].sequence, strlen(adaptors[0][i]));
            p1[strlen(adaptors[0][i]) - 5] = 0;
            matched = 1;
            break;
        }
    }
        
    
    if (p1_index < 0) {
        p1[0] = 0;
        for (o=4; o<=7; o++) {
            m = compare_sequence(read_pair->read[0].sequence + o, "TGCAG", 5);
            if (m <= allowed_mismatches) {
               strncpy(p1, read_pair->read[0].sequence, o);
               p1[o] = 0;
            }
        }
    }
    
    /*
    for (o=4; o<=7; o++) {
        m = compare_sequence(read_pair->read[0].sequence + o, "TGCAG", 5);
        if (m <= allowed_mismatches) {
            strncpy(p1, read_pair->read[0].sequence, o);
            p1[o] = 0;
            p1_index = match_adaptor(p1, 0);
            //printf("    Detected PstI from base %d with %d mismatches: %s-TGCAG p2 is %s\n", o+1, m, p1, p2);
            matched = 1;
        }
    }*/
    
    if ((matched) && (p1_index >=0) && (p2_index >=0)) {
        //printf("p1=%s (%d)\tp2=%s (%d)\n", p1, p1_index, p2, p2_index);
        if (out_fp[p1_index][p2_index][0] == 0) {
            int i;
            for (i=0; i<2; i++) {
                char filename[MAX_PATH_LENGTH];
                sprintf(filename, "%s_%c%d_R%d.fastq", output_prefix, p2_index+'A', p1_index+1, i+1);
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
        clip_size = strlen(adaptors[0][p1_index]);
        if (clip_psti == 0) {
            clip_size -= 5;
        }
        adaptor_counts[p1_index][p2_index]++;
    } else {
        //printf("No match\n");
        
        store_undetermined(0, p1);
        store_undetermined(1, p2);
        
        out_r1 = undetermined_fp[0];
        out_r2 = undetermined_fp[1];
        clip_size = 0;
        undetermined_read_count++;
    }

    strcat(read_pair->read[0].sequence_header, " ");
    strcat(read_pair->read[0].sequence_header, p1);
    strcat(read_pair->read[0].sequence_header, "-");
    strcat(read_pair->read[0].sequence_header, p2);
    
    write_read(&read_pair->read[0], clip_size, out_r1);
    write_read(&read_pair->read[1], 0, out_r2);
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
void display_adaptors(void)
{
    int i,j;
    
    for (i=0; i<2; i++) {
        printf("P%d adaptors:\n", i);
        for (j=0; j<n_adaptors[i]; j++) {
            if (i == 0) {
                printf("  %d. %s\n", j, adaptors[i][j]);
            } else {
                printf("  %c. %s\n", j+'A', adaptors[i][j]);
            }
            adaptor_counts[i][j] = 0;
        }
    }
    
    printf("\n");
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
void load_adaptor_files(void)
{
    int i,j;
    
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
                            strcat(adaptors[i][n_adaptors[i]], "TGCAG");
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
    
    printf("\n");
}

/*----------------------------------------------------------------------*
 * Function:
 * Purpose:
 * Parameters:
 * Returns:
 *----------------------------------------------------------------------*/
void display_counts(void)
{
    double percent = 0.0;
    int i, j;
    
    printf("\nCat\tP1\tP2\tCount\tPercent\n");
    
    for (j=0; j<n_adaptors[1]; j++) {
        for (i=0; i<n_adaptors[0]; i++) {
            percent = 0.0;
            if (adaptor_counts[i][j] > 0) {
                percent = (100.0 * adaptor_counts[i][j]) / total_read_count;
            }
            printf("%c%d\t%s\t%s\t%d\t%.2f\n", j+'A', i+1, adaptors[0][i], adaptors[1][j], adaptor_counts[i][j], percent);
        }
    }
    
    if (undetermined_read_count > 0) {
        percent = (100.0 * undetermined_read_count) / total_read_count;
    }
    printf("Und\t\t\t%d\t%.2f\n", undetermined_read_count, percent);
    printf("Total\t\t\t%d\t100\n", total_read_count);
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
        {"clip_psti", no_argument, NULL, 'z'},
        {"p1", required_argument, NULL, '1'},
        {"p2", required_argument, NULL, '2'},
        {0, 0, 0, 0}
    };
    int opt;
    int longopt_index;
    
    while ((opt = getopt_long(argc, argv, "a:b:c:hm:p:vz1:2:", long_options, &longopt_index)) > 0)
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
            case 'z':
                clip_psti = 1;
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
        printf("Using default adaptors.\n");
        setup_default_adaptors();
    } else {
        load_adaptor_files();
    }
    
    printf("Allowed mismatches: %d\n\n", allowed_mismatches);
}

/*----------------------------------------------------------------------*
 * Function:   main
 *----------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    FastqReadPair read_pair;
    
    printf("\nRADplex v0.4\n\n");
    
    initialise_main();

    //store_undetermined(0, "AATAGTT");
    //store_undetermined(0, "AATAGTT");
    //store_undetermined(0, "AAT");

    initialise_read_pair_struct(&read_pair);
    parse_command_line(argc, argv, &read_pair);
    display_adaptors();
    read_files(&read_pair);
    display_counts();

    output_undetermined_indices();
    
    printf("\nDone.\n");
    
    return 0;
}
