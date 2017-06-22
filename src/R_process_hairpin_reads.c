#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#define BLOCKSIZE 10000000
#define ALLOC(n,t) (t *) calloc(n, sizeof(t))
#define VALID_MATCH(seq1, seq2, len) valid_match(seq1, seq2, len, g_barcode_n_mismatch)

typedef struct {
  char   *seq1;
  char   *seq2;
  char   *seq_rev;
  int    original_pos;
} a_barcode;

typedef struct {
   char   *seq;
   int    original_pos;
} a_hairpin;

a_barcode **barcodes;
a_hairpin **hairpins;

//* Globals *//
int g_is_paired_reads;
int g_is_dual_indexing_reads;
int g_num_barcode;
int g_num_hairpin;
long g_num_read;
long **g_summary;
int g_barcode_start;
int g_barcode_end;
int g_barcode2_start;
int g_barcode2_end;
int g_barcode_start_rev;
int g_barcode_end_rev;
int g_barcode_length;
int g_barcode2_length;
int g_barcode_length_rev;
int g_hairpin_start;
int g_hairpin_end;
int g_hairpin_length;
int g_allow_shifting;
int g_shifting_n_base;
int g_allow_mismatch;
int g_barcode_n_mismatch;
int g_hairpin_n_mismatch;
int g_allow_shifted_mismatch;
int g_is_verbose;

long g_barcode_count;
long g_hairpin_count;
long g_bc_hp_count;
//** End Globals **//

// get number of lines in file
int get_lines_in_file(FILE* fin) {
  int n_lines = 0, ch, last_ch = '\n';
  while (1) {
    ch = fgetc(fin);
    if (ch == '\n') {
      n_lines++;
    } else if (ch == EOF) {
      if (last_ch != '\n') {
        n_lines++;
      } // Capture non-newline-terminated last line.
      break;
    }
    last_ch = ch;
  }
  rewind(fin);
  return n_lines;
}


void read_in_barcodes(char* filename) {
  FILE *fin;
  char *line = NULL;
  size_t len = 1000;
  char *readline;

  fin = fopen(filename,"r");

  // get number of lines in the file.
  g_num_barcode = get_lines_in_file(fin);
  barcodes = (a_barcode**)R_alloc(g_num_barcode+1, sizeof(a_barcode*));

  line = ALLOC(len+1, char);
  a_barcode *new_barcode;
  int count = 0;
  char * token;

  while ((readline = fgets(line, len, fin)) != NULL) {
    count++;
    new_barcode = ALLOC(1, a_barcode);
    new_barcode->seq1 = ALLOC(g_barcode_length, char);
    strncpy(new_barcode->seq1, line, g_barcode_length);
    new_barcode->original_pos = count;
    if (g_is_paired_reads > 0) {
      token = strtok(line, "\t");
      token = strtok(NULL, "\t"); // next token
      new_barcode->seq_rev = ALLOC(g_barcode_length_rev, char);
      strncpy(new_barcode->seq_rev, token, g_barcode_length_rev);
    } else if (g_is_dual_indexing_reads > 0) {
      token = strtok(line, "\t");
      token = strtok(NULL, "\t"); // next token
      new_barcode->seq2 = ALLOC(g_barcode_length_rev, char);
      strncpy(new_barcode->seq2, token, g_barcode2_length);
    } else {
      new_barcode->seq_rev = NULL;
    };
    barcodes[count] = new_barcode;
  }
  fclose(fin);
  free(line);

  Rprintf(" -- Number of Barcodes : %d\n", g_num_barcode);
}


bool valid_match(char *sequence1, char *seq2, int length, int threshold) {
  int mismatch_base_count = 0;
  for (int i = 0; i < length; i++) {
    if (sequence1[i] != seq2[i]) {
      mismatch_base_count++;
    }
  }

  if (mismatch_base_count <= threshold) {
    return true;
  } else {
    return false;
  }
}


int locate_barcode(char *a_barcode) {
  int imin = 1;
  int imax = g_num_barcode;

  while (imax >= imin) {
    int imid = (imax + imin) / 2;

    if (strncmp(barcodes[imid]->seq1, a_barcode, g_barcode_length) < 0) {
      imin = imid + 1;
    } else if (strncmp(barcodes[imid]->seq1, a_barcode, g_barcode_length) > 0) {
      imax = imid - 1;
    } else {
      return barcodes[imid]->original_pos;
    }
  }

  if (g_allow_mismatch > 0) {
    for (int i = 1; i <= g_num_barcode; i++) {
      if (VALID_MATCH(a_barcode, barcodes[i]->seq1, g_barcode_length)) {
        return barcodes[i]->original_pos;
      }
    }
  }

  return -1;
}


int locate_barcode_paired(char *a_barcode, char *a_barcode_rev) {
  int imin = 1;
  int imax = g_num_barcode;

  while (imax >= imin) {
    int imid = (imax + imin) / 2;
    // compare forward barcode sequence
    if (strncmp(barcodes[imid]->seq1, a_barcode, g_barcode_length) < 0) {
      imin = imid + 1;
    } else if (strncmp(barcodes[imid]->seq1, a_barcode, g_barcode_length) > 0) {
      imax = imid - 1;
    } else {
      // same forward sequence, compare reverse barcode sequence
      if (strncmp(barcodes[imid]->seq_rev, a_barcode_rev, g_barcode_length_rev) < 0) {
        imin = imid + 1;
      } else if (strncmp(barcodes[imid]->seq_rev, a_barcode_rev, g_barcode_length_rev) > 0) {
        imax = imid - 1;
      } else {
        return barcodes[imid]->original_pos;
      }
    }
  }

  if (g_allow_mismatch > 0) {
    for (int i = 1; i <= g_num_barcode; i++) {
      if (VALID_MATCH(a_barcode, barcodes[i]->seq1, g_barcode_length) &&
          VALID_MATCH(a_barcode_rev, barcodes[i]->seq_rev, g_barcode_length_rev)) {
        return barcodes[i]->original_pos;
      }
    }
  }

  return -1;
}


int locate_barcode_dual_indexing(char *a_barcode, char *a_barcode2) {
  int imin = 1;
  int imax = g_num_barcode;

  while (imax >= imin) {
    int imid = (imax + imin) / 2;
    // compare forward barcode sequence 1
    if (strncmp(barcodes[imid]->seq1, a_barcode, g_barcode_length) < 0) {
      imin = imid + 1;
    } else if (strncmp(barcodes[imid]->seq1, a_barcode, g_barcode_length) > 0) {
      imax = imid - 1;
    } else {
      // same forward sequence 1, compare forward barcode sequence 2
      if (strncmp(barcodes[imid]->seq2, a_barcode2, g_barcode2_length) < 0) {
        imin = imid + 1;
      } else if (strncmp(barcodes[imid]->seq2, a_barcode2, g_barcode2_length) > 0) {
        imax = imid - 1;
      } else {
        return barcodes[imid]->original_pos;
      }
    }
  }

  if (g_allow_mismatch > 0) {
    for (int i = 1; i <= g_num_barcode; i++) {
      if (VALID_MATCH(a_barcode, barcodes[i]->seq1, g_barcode_length) &&
          VALID_MATCH(a_barcode2, barcodes[i]->seq2, g_barcode2_length)) {
        return barcodes[i]->original_pos;
      }
    }
  }

  return -1;
}


int locate_exact_match_hairpin(char *a_hairpin) {
  int imin = 1;
  int imax = g_num_hairpin;

  while (imax >= imin) {
    int imid = (imax + imin) / 2;

    if (strncmp(hairpins[imid]->seq, a_hairpin, g_hairpin_length) < 0) {
      imin = imid + 1;
    } else if (strncmp(hairpins[imid]->seq, a_hairpin, g_hairpin_length) > 0) {
      imax = imid - 1;
    } else {
      return hairpins[imid]->original_pos;
    }
  }
  return -1;
}


int locate_mismatch_hairpin(char *a_hairpin) {
  for (int i = 1; i <= g_num_hairpin; i++) {
    if (VALID_MATCH(a_hairpin, hairpins[i]->seq, g_hairpin_length)) {
      return hairpins[i]->original_pos;
    }
  }
  return -1;
}


int locate_hairpin(char *a_hairpin, char *read) {
  // check if a perfect match exists
  int hairpin_index = locate_exact_match_hairpin(a_hairpin);
  if (hairpin_index > 0) {
    return hairpin_index;
  }

  // if a perfect match doesn't exist, check if a mismatch exists
  if (g_allow_mismatch > 0) {
    hairpin_index = locate_mismatch_hairpin(a_hairpin);
    if (hairpin_index > 0) {
      return hairpin_index;
    }
  }

  if (g_allow_shifting > 0) {
    // Check if given hairpin can be mapped to a shifted location.
    char *shifted_hairpin_str;
    shifted_hairpin_str = ALLOC(g_hairpin_length, char);

    // check shifting leftwards
    for (int index = 1; index <= g_shifting_n_base; index++) {
      strncpy(shifted_hairpin_str, read + g_hairpin_start - 1 - index, g_hairpin_length);

      hairpin_index = locate_exact_match_hairpin(shifted_hairpin_str);

      // check if mismatch on a shifted position is allowed and try to match
      if ((hairpin_index <= 0) && (g_allow_shifted_mismatch)) {
        hairpin_index = locate_mismatch_hairpin(shifted_hairpin_str);
      }

      if (hairpin_index > 0) {
        return hairpin_index;
      }
    }

    // check shifting rightwards
    for (int index = 1; index <= g_shifting_n_base; index++) {
      strncpy(shifted_hairpin_str, read + g_hairpin_start - 1 + index, g_hairpin_length);

      hairpin_index = locate_exact_match_hairpin(shifted_hairpin_str);

      if ((hairpin_index <= 0) && (g_allow_shifted_mismatch)) {
        hairpin_index = locate_mismatch_hairpin(shifted_hairpin_str);
      }

      if (hairpin_index > 0) {
        return hairpin_index;
      }
    }
  }

  return -1;
}


int barcode_compare(a_barcode *barcode1, a_barcode *barcode2) {
  int ans = strncmp(barcode1->seq1, barcode2->seq1, g_barcode_length);
  if (ans == 0) {
    if (g_is_paired_reads > 0) {
      ans = strncmp(barcode1->seq_rev, barcode2->seq_rev, g_barcode_length_rev);
    } else if (g_is_dual_indexing_reads > 0) {
      ans = strncmp(barcode1->seq2, barcode2->seq2, g_barcode2_length);
    }
  }

  return ans;
}


void sort_barcodes(void) {
  a_barcode *temp;
  for (int i = 1; i < g_num_barcode; i++) {
    for (int j = i + 1; j <= g_num_barcode; j++) {
      if (barcode_compare(barcodes[i], barcodes[j]) > 0) {
        temp = barcodes[i];
        barcodes[i] = barcodes[j];
        barcodes[j] = temp;
      }
    }
  }
}


void read_in_hairpins(char *filename) {
  FILE *fin;
  char *line = NULL;
  size_t len = 1000;
  char *readline;

  fin = fopen(filename, "r");

  // Getting number of lines in the file.
  g_num_hairpin = get_lines_in_file(fin);
  hairpins = (a_hairpin **)R_alloc(g_num_hairpin + 1, sizeof(a_hairpin *));

  line = ALLOC(len + 1, char);
  a_hairpin *new_hairpin;
  int count = 0;

  while ((readline = fgets(line, len, fin)) != NULL) {
    count++;
    new_hairpin = ALLOC(1, a_hairpin);
    new_hairpin->seq = ALLOC(g_hairpin_length, char);
    new_hairpin->original_pos = count;
    strncpy(new_hairpin->seq, line, g_hairpin_length);
    hairpins[count] = new_hairpin;
  }
  fclose(fin);
  free(line);

  Rprintf(" -- Number of Hairpins : %d\n", g_num_hairpin);
}


void sort_hairpins(void) {
  a_hairpin *temp;
  for (int i = 1; i < g_num_hairpin; i++) {
    for (int j = i + 1; j <= g_num_hairpin; j++) {
      if (strncmp(hairpins[i]->seq, hairpins[j]->seq, g_hairpin_length) > 0) {
        temp = hairpins[i];
        hairpins[i] = hairpins[j];
        hairpins[j] = temp;
      }
    }
  }
}


void process_hairpin_reads(char *filename, char *filename2) {
  FILE *fin = NULL;
  FILE *fin_rev = NULL;
  char *line = NULL;
  char *line2 = NULL;
  size_t len = 1000;
  char *readline;
  char *readline2;
  long num_read_this_file = 0;

  char *this_barcode_for = NULL;
  char *this_barcode_2 = NULL;
  char *this_barcode_rev = NULL;
  char *this_hairpin = NULL;

  line = ALLOC(len+1, char);
  fin = fopen(filename,"r");
  if (g_is_paired_reads > 0) {
    fin_rev = fopen(filename2, "r");
    line2 = ALLOC(len + 1, char);
  }

  if (g_is_verbose > 0) {
    if (g_is_paired_reads > 0) {
      Rprintf("Processing reads in %s and %s.\n", filename, filename2);
    } else {
      Rprintf("Processing reads in %s.\n", filename);
    }
  }

  this_barcode_for = ALLOC(g_barcode_length, char);
  if (g_is_dual_indexing_reads > 0) {
    this_barcode_2 = ALLOC(g_barcode2_length, char);
  }
  if (g_is_paired_reads > 0) {
    this_barcode_rev = ALLOC(g_barcode_length_rev, char);
  }
  this_hairpin = ALLOC(g_hairpin_length, char);
  long line_count = 0;

  while ((readline = fgets(line, len, fin)) != NULL) {
    if (g_is_paired_reads > 0) {
      readline2 = fgets(line2, len, fin_rev);
      if (readline2 == NULL) {
        break;
      }
    }
    line_count++;
    if ((line_count % 4) != 2) {
      continue;
    }

    if ((g_is_verbose > 0) && (num_read_this_file % BLOCKSIZE == 0)) {
      Rprintf(" -- Processing %d million reads\n", (num_read_this_file / BLOCKSIZE + 1) * 10);
    }
    g_num_read++;
    num_read_this_file++;

    strncpy(this_barcode_for, line + g_barcode_start - 1, g_barcode_length);

    int barcode_index;
    if (g_is_paired_reads > 0) {
      strncpy(this_barcode_rev, line2 + g_barcode_start_rev - 1, g_barcode_length_rev);
      barcode_index = locate_barcode_paired(this_barcode_for, this_barcode_rev);
    } else if (g_is_dual_indexing_reads > 0) {
      strncpy(this_barcode_2, line + g_barcode2_start - 1, g_barcode2_length);
      barcode_index = locate_barcode_dual_indexing(this_barcode_for, this_barcode_2);
    } else {
      barcode_index = locate_barcode(this_barcode_for);
    }

    strncpy(this_hairpin, line + g_hairpin_start - 1, g_hairpin_length);
    int hairpin_index = locate_hairpin(this_hairpin, line);

    if (barcode_index > 0) {
      g_barcode_count++;
    }

    if (hairpin_index > 0) {
      g_hairpin_count++;
    }

    if ((barcode_index > 0) && (hairpin_index > 0)) {
      g_summary[hairpin_index][barcode_index]++;
      g_bc_hp_count++;
    }
  }

  if (g_is_verbose > 0) {
    if (g_is_paired_reads > 0) {
      Rprintf("Number of reads in file %s and %s: %ld\n", filename, filename2, num_read_this_file);
    } else {
      Rprintf("Number of reads in file %s : %ld\n", filename, num_read_this_file);
    }
  }

  fclose(fin);
  free(line);
  free(this_barcode_for);
  free(this_hairpin);

  if (g_is_paired_reads > 0) {
    fclose(fin_rev);
    free(line2);
    free(this_barcode_rev);
  }
}


void initialise(int is_paired, int is_dual_indexing,
                int bar_code_start, int barcode_end,
                int barcode2_start, int barcode2_end,
                int bar_code_start_rev, int barcode_endrev,
                int hairpin_start, int hairpin_end,
                int allow_shifting, int shifting_n_base,
                int allow_mismatch, int barcode_mismatch, int hairpin_mismatch,
                int allow_shifted_mismatch, int verbose) {

  g_num_barcode = 0;
  g_num_hairpin = 0;

  g_is_paired_reads = is_paired;
  g_is_dual_indexing_reads = is_dual_indexing;
  g_barcode_start = bar_code_start;
  g_barcode_end = barcode_end;
  g_barcode2_start = barcode2_start;
  g_barcode2_end = barcode2_end;
  g_barcode_start_rev = bar_code_start_rev;
  g_barcode_end_rev = barcode_endrev;
  g_hairpin_start = hairpin_start;
  g_hairpin_end = hairpin_end;
  g_barcode_length = barcode_end - bar_code_start + 1;
  g_barcode2_length = barcode2_end - barcode2_start + 1;
  g_barcode_length_rev = barcode_endrev - bar_code_start_rev + 1;
  g_hairpin_length = hairpin_end - hairpin_start + 1;

  g_allow_shifting = allow_shifting;
  g_shifting_n_base = shifting_n_base;
  g_allow_mismatch = allow_mismatch;
  g_barcode_n_mismatch = barcode_mismatch;
  g_hairpin_n_mismatch = hairpin_mismatch;
  g_allow_shifted_mismatch = allow_shifted_mismatch;
  g_is_verbose = verbose;

  g_num_read = 0;
  g_barcode_count = 0;
  g_hairpin_count = 0;
  g_bc_hp_count = 0;
}


void output_summary_table(char *output) {
  FILE *fout;
  fout = fopen(output, "w");
  for (int i = 1; i <= g_num_hairpin; i++) {
    fprintf(fout, "%ld", g_summary[i][1]);
    for (int j = 2; j <= g_num_barcode; j++) {
      fprintf(fout, "\t%ld", g_summary[i][j]);
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
}


void check_hairpins(void) {
  for (int p = 1; p <= g_num_hairpin; p++) {
    for (int q = 0; q < g_hairpin_length; q++) {
      char base = hairpins[p]->seq[q];
      if ((base != 'A') && (base != 'T') && (base != 'G') && (base != 'C')) {
        Rprintf("Hairpin no.%d: %s contains invalid base %c\n", p, hairpins[p]->seq, base);
      }
    }
  }
}


void clean_up(void) {
  for (int index = 1; index <= g_num_barcode; index++) {
    free(barcodes[index]->seq1);

    if (g_is_paired_reads > 0) free(barcodes[index]->seq_rev);

    if (g_is_dual_indexing_reads > 0) free(barcodes[index]->seq2);

    free(barcodes[index]);
  }

  for (int index = 1; index <= g_num_hairpin; index++) {
    free(hairpins[index]->seq);
    free(hairpins[index]);
  }

  for (int index = 0; index <= g_num_hairpin; index++) {
    free(g_summary[index]);
  }

  free(g_summary);
}


void allocate_summary_table(void) {
  g_summary = ALLOC(g_num_hairpin + 1, long *);
  for (int i = 0; i <= g_num_hairpin; i++) {
    g_summary[i] = ALLOC(g_num_barcode + 1, long);
  }

  for (int i = 0; i <= g_num_hairpin; i++) {
    for (int j = 0; j <= g_num_barcode; j++) {
      g_summary[i][j] = 0;
    }
  }
}


void processHairpinReads(int *g_is_paired_reads, int *g_is_dual_indexing_reads,
                         char **file, char **file2, int *filecount,
                         char **barcode_seqs, char **hairpin_seqs,
                         int *g_barcode_start, int *g_barcode_end,
                         int *g_barcode2_start, int *g_barcode2_end,
                         int *g_barcode_start_rev, int *g_barcode_end_rev,
                         int *g_hairpin_start, int *g_hairpin_end,
                         int *g_allow_shifting, int *g_shifting_n_base,
                         int *g_allow_mismatch, int *barcode_mismatch, int *hairpin_mismatch,
                         int *g_allow_shifted_mismatch,
                         char **output, int *verbose) {
  // initialise global variables
  initialise(*g_is_paired_reads, *g_is_dual_indexing_reads,
             *g_barcode_start, *g_barcode_end,
             *g_barcode2_start, *g_barcode2_end,
             *g_barcode_start_rev, *g_barcode_end_rev,
             *g_hairpin_start, *g_hairpin_end,
             *g_allow_shifting, *g_shifting_n_base,
             *g_allow_mismatch, *barcode_mismatch, *hairpin_mismatch,
             *g_allow_shifted_mismatch, *verbose);

  read_in_barcodes(*barcode_seqs);
  sort_barcodes(); // () is necessary to invoke function without parameters

  read_in_hairpins(*hairpin_seqs);

  check_hairpins();
  sort_hairpins();

  allocate_summary_table();

  int fin;

  for (fin = 0; fin < *filecount; fin++) {
    process_hairpin_reads(file[fin], file2[fin]);
  }

  Rprintf("\nThe input run parameters are: \n");
  Rprintf(" -- Barcode: start position %d\t end position %d\t length %d\n", g_barcode_start, g_barcode_end, g_barcode_length);
  if (g_is_dual_indexing_reads) {
    Rprintf(" -- Second Barcode in forward read: start position %d\t end position %d\t length %d\n", g_barcode2_start, g_barcode2_end, g_barcode2_length);
  }
  if (g_is_paired_reads) {
    Rprintf(" -- Barcode in reverse read: start position %d\t end position %d\t length %d\n", g_barcode_start_rev, g_barcode_end_rev, g_barcode_length_rev);
  }
  Rprintf(" -- Hairpin: start position %d\t end position %d\t length %d\n", g_hairpin_start, g_hairpin_end, g_hairpin_length);
  if (g_allow_shifting) {
    Rprintf(" -- Allow hairpin sequences to be matched to a shifted position, <= %d base left or right of the specified positions. \n", g_shifting_n_base);
  } else {
    Rprintf(" -- Hairpin sequences need to match at specified positions. \n");
  }
  if (g_allow_mismatch) {
    Rprintf(" -- Allow sequence mismatch, <= %d base in barcode sequence and <= %d base in hairpin sequence. \n", g_barcode_n_mismatch, g_hairpin_n_mismatch );
  } else {
    Rprintf(" -- Mismatch in barcode/hairpin sequences not allowed. \n");
  }

  Rprintf("\nTotal number of read is %ld \n", g_num_read);
  Rprintf("There are %ld reads (%.4f percent) with barcode matches\n", g_barcode_count, 100.0*g_barcode_count/g_num_read);
  Rprintf("There are %ld reads (%.4f percent) with hairpin matches\n", g_hairpin_count, 100.0*g_hairpin_count/g_num_read);
  Rprintf("There are %ld reads (%.4f percent) with both barcode and hairpin matches\n", g_bc_hp_count, 100.0*g_bc_hp_count/g_num_read);

  output_summary_table(*output);

  clean_up();
}
