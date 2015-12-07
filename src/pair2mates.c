#include <zlib.h>
#include <stdio.h>
#include <math.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/bgzf.h"
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


static int strnum_cmp(const char *_a, const char *_b)
{
    const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
    const unsigned char *pa = a, *pb = b;
    while (*pa && *pb) {
        if (isdigit(*pa) && isdigit(*pb)) {
            while (*pa == '0') ++pa;
            while (*pb == '0') ++pb;
            while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
            if (isdigit(*pa) && isdigit(*pb)) {
                int i = 0;
                while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
                return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
            } else if (isdigit(*pa)) return 1;
            else if (isdigit(*pb)) return -1;
            else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
        } else {
            if (*pa != *pb) return (int)*pa - (int)*pb;
            ++pa; ++pb;
        }
    }
    return *pa? 1 : *pb? -1 : 0;
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
	char* sequence="";
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
    char* quality="";
	quality = mycalloc(len, char);
	size_t i;
    for (i=0; i<len; i++) {
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
	fprintf(stderr, "Program: pair2mates (pair up two mates of Hi-C)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Rongxin Fang <r3fang@ucsd.edu>\n\n");
	fprintf(stderr, "Usage:   pair2mates [options] <R1.bam> <R2.bam>\n\n");
	fprintf(stderr, "Options: -o FILENAME   Specifies the name of output file (e.g. name.paired.bam)\n");
	fprintf(stderr, "         -d            Min distance for valid pairs [10000]\n");
	fprintf(stderr, "\n");
	return 1;
}

static char *get_cigar_str(const bam1_t *b){
	uint32_t *cigar = bam_get_cigar(b);
	char* str1 = mycalloc(30, char);
	strcpy(str1, "");
	char buf[10];
	int k;
	for (k = 0; k < b->core.n_cigar; ++k){
	    if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
			sprintf(buf, "%d%c", bam_cigar_oplen(cigar[k]), bam_cigar_opchr(cigar[k]));
			strcat(str1, buf);
	}
	return str1;
}
// DP matrix

typedef struct {
	bool f;
	char* input1;	
	char* input2;
	char* output;
	int dist;	
} opt_t;

void init_opt(opt_t *s){
	s->f = false;
	s->dist = 10000;
}

void opt_destory(opt_t *s){
	free(s->input1);
	free(s->input2);
	free(s->output);
	free(s);
}

void pairUp(opt_t *opt){
	samFile *in1 = sam_open(opt->input1, "r");
	samFile *in2 = sam_open(opt->input2, "r");
	BGZF *out = bgzf_open(opt->output, "w");
	
	if(in1 == NULL) die("pairUp: fail to open file '%s'", opt->input1);
	if(in2 == NULL) die("pairUp: fail to open file '%s'", opt->input2);
	if(out == NULL) die("pairUp: fail to open file '%s'", opt->output);
	
	//if(access(opt->output, F_OK)!=-1 && opt->f==false) die("pairUp: %s exists, use opetion -f to overwrite", opt->output);
	bam_hdr_t *header1 = sam_hdr_read(in1);
	bam_hdr_t *header2 = sam_hdr_read(in2);

	bam1_t *aln1 = bam_init1();
	bam1_t *aln2 = bam_init1();
	
	int8_t *p;
	int32_t n;
	int ret1, ret2;

	ret1=sam_read1(in1, header1, aln1);
	ret2=sam_read1(in2, header2, aln2);
	
	kstring_t *name1 = mycalloc(1, kstring_t); 
	kstring_t *name2 = mycalloc(1, kstring_t); 
	bam_hdr_write(out, header1);
	
	while (ret1 >= 0 && ret2 >= 0){
		name1->s = get_pair_name(aln1); name1->l = strlen(name1->s);
		name2->s = get_pair_name(aln2); name2->l = strlen(name2->s);
		if(strnum_cmp(name1->s, name2->s) < 0){
			ret1=sam_read1(in1, header1, aln1);
		}
		else if(strnum_cmp(name1->s, name2->s) == 0){
			if(strcmp(header2->target_name[aln2->core.tid],  header2->target_name[aln2->core.tid])==0 && abs(aln1->core.pos - aln2->core.pos) > opt->dist){
				aln1->core.flag = 83;
				aln1->core.mpos = aln2->core.pos;
				aln1->core.isize = aln1->core.pos - aln2->core.pos;
				
				aln2->core.flag = 163;				
				aln2->core.mpos = aln1->core.pos;
				aln2->core.isize = aln2->core.pos - aln1->core.pos;

				bam_write1(out, aln1);
				bam_write1(out, aln2);
			}			
			// 83/163
			ret1=sam_read1(in1, header1, aln1);
			ret2=sam_read1(in2, header2, aln2);
		}
		else if(strnum_cmp(name1->s, name2->s) > 0){
			ret2=sam_read1(in2, header2, aln2);
		}
	}	
	bam_destroy1(aln1);
	bam_destroy1(aln2);
	sam_close(in1);
	sam_close(in2);
	if ( bgzf_close(out)<0 ) error("Close failed\n");
	free(name1);
	free(name2);
}

int main(int argc, char *argv[]){
	if (argc < 2) return usage();
	int c;
	opt_t *opt = mycalloc(1, opt_t);
	init_opt(opt);
	while ((c = getopt (argc, argv, "o:m:")) != -1){
		switch(c){
			case 'o': opt->output = strdup(optarg); break;
			case 'm': opt->dist = atoi(optarg); break;
			default: return 1;
		}
	}
	opt->input1 = strdup(argv[argc-2]);
	opt->input2 = strdup(argv[argc-1]);
	pairUp(opt);
	opt_destory(opt);
	return 0;
}
