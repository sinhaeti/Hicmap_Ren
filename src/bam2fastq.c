#include <zlib.h>
#include <stdio.h>
#include <math.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "utils.h"
#include <unistd.h>
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "1.0.0"
#endif

typedef enum {true, false} bool;
/*
 * reverse a string
 */
char* strrev(char *s){
	if(s == NULL) return NULL;
	int l = strlen(s);
	char *ss = strdup(s);
	free(s);
	s = mycalloc(l, char);
	int i; for(i=0; i<l; i++){
		s[i] = ss[l-i-1];
	}
	s[l] = '\0';
	return s;
}

const int get_read_idx(const bam1_t *b) {
    return !(b->core.flag & BAM_FREAD1);
}

char* get_pair_name(const bam1_t *b) {
    return bam_get_qname(b);
}

const char* get_read_name(const bam1_t *b) {
    char* suffix[2] = {"/1", "/2"};
    char* name = bam_get_qname(b);
    if (b->core.flag & BAM_FPAIRED)
		strcpy(name, suffix[get_read_idx(b)]);
    return name;
}

const char* get_sequence(const bam1_t *b){
	if(b == NULL) die("get_sequence: parameter error\n");
	const uint8_t *seq = bam_get_seq(b);
	size_t len = b->core.l_qseq;
	char* sequence;
	sequence = malloc(len*sizeof(char));
	uint8_t offset = (b->core.flag & BAM_FREVERSE) ? 16 : 0;
    size_t i; for (i=0; i<len; i++) {
		switch(bam_seqi(seq, i) + offset)
		{
			case 1:
				sequence[i] = 'A';
				break;
			case 2:
				sequence[i] = 'C';
				break;
			case 4:
				sequence[i] = 'G';
				break;
			case 8:
				sequence[i] = 'T';
				break;
			case 15:
				sequence[i] = 'N';
				break;
				//Complements (original index + 16)
			case 17:
				sequence[i] = 'T';
				break;
			case 18:
				sequence[i] = 'G';
				break;
			case 20:
				sequence[i] = 'C';
				break;
			case 24:
				sequence[i] = 'A';
				break;
			case 31:
				sequence[i] = 'N';
				break;
			default:
				sequence[i] = 'N';
				break;
		}
     }
	 if (offset) sequence = strrev(sequence);
	 return sequence;
}

char* get_qualities(const bam1_t *b) {
    const uint8_t *qual = bam_get_qual(b);
    size_t len = b->core.l_qseq;
    char* quality;
	quality = mycalloc(len, char);
    for (size_t i=0; i<len; i++) {
		quality[i] = (char)((int)qual[i] + 33);
    }
	if (b->core.flag & BAM_FREVERSE)
		quality = strrev(quality);
    return quality;
}

int get_lane_id(const bam1_t *b){
	if(b==NULL) die("get_lane_id: input error");
	kstring_t *name = mycalloc(1, kstring_t);
	name->s = get_pair_name(b);
	name->l = strlen(name->s);
	int *fields, n, i;
	fields = ksplit(name, (int)':', &n);
	if(n==1) return 0;
	int lane_id = 0;
	sscanf(name->s + fields[1], "%d", &lane_id);
	return lane_id;
}

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: bam2fastq (convert bam to fastq file)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Rongxin Fang <r3fang@ucsd.edu>\n\n");
	fprintf(stderr, "Usage:   bam2fastq [options] <in.sam>\n\n");
	fprintf(stderr, "Options: -o FILENAME   Specifies the name of the FASTQ file(s) that will be generated\n");
	fprintf(stderr, "         -f            overwriting existing files if necessary\n");
	fprintf(stderr, "         -a            Reads in the BAM that are aligned will not be extracted\n");
	fprintf(stderr, "         -u            Reads in the BAM that are not aligned will not be extracted\n");
	fprintf(stderr, "\n");
	return 1;
}

// DP matrix
typedef struct {
	bool f;
	bool a;
	bool u;
	char* in_name;	
	char* out_name;
} opt_t;

void init_opt(opt_t *s){
	s->f = false;
	s->a = false;
	s->u = false;
}

void opt_destory(opt_t *s){
	free(s->in_name);
	free(s->out_name);
	free(s);
}

void bam_parser(opt_t *opt){
	samFile *in = sam_open(opt->in_name, "r");
	if(in == NULL) die("bam_parser: fail to open file '%s'", opt->in_name);
	// if output file exists but not force to overwrite
	if(access(opt->out_name, F_OK)!=-1 && opt->f==false) die("bam_parser: %s exists, use opetion -f to overwrite", opt->out_name);
	bam_hdr_t *header = sam_hdr_read(in);
	bam1_t *aln = bam_init1();
	int8_t *p;
	int32_t n;
	int ret;
	while((ret=sam_read1(in, header, aln)) >= 0)
		printf("name=%s\nflag=%d\nseq=%s\nqual=%s\nlane_id=%d\n", get_read_name(aln), aln->core.flag, get_sequence(aln), get_qualities(aln), get_lane_id(aln));

	bam_destroy1(aln);
	sam_close(in);
}

int main(int argc, char *argv[]){
	if (argc < 2) return usage();
	int c;
	opt_t *opt = mycalloc(1, opt_t);
	init_opt(opt);
	while ((c = getopt (argc, argv, "fauo:")) != -1){
		switch(c){
			case 'f': opt->f = true; break;
			case 'a': opt->a = true; break;
			case 'u': opt->u = true; break;
			case 'o': opt->out_name = strdup(optarg); break;
			default: return 1;
		}
	}
	opt->in_name = strdup(argv[argc-1]);
	bam_parser(opt);
	opt_destory(opt);
	return 0;
}
