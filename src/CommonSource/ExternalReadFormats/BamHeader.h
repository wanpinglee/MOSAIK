#ifndef BAMHEADER_H_
#define BAMHEADER_H_

// our zlib constants
#define GZIP_ID1             31
#define GZIP_ID2            139
#define CM_DEFLATE            8
#define FLG_FEXTRA            4
#define OS_UNKNOWN          255
#define BGZF_XLEN             6
#define BGZF_ID1             66
#define BGZF_ID2             67
#define BGZF_LEN              2
#define GZIP_WINDOW_BITS    -15
#define Z_DEFAULT_MEM_LEVEL   8

// our BZGF constants
#define BLOCK_HEADER_LENGTH    18
#define BLOCK_FOOTER_LENGTH     8
#define MAX_BLOCK_SIZE      65536
#define DEFAULT_BLOCK_SIZE  65536

// our BAM constants
#define BAM_CORE_SIZE  32
#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6
#define BAM_CIGAR_SHIFT 4

#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)

// define some SAM/BAM FLAGS
#define BAM_SEQUENCED_AS_PAIRS         1
#define BAM_PROPER_PAIR                2
#define BAM_QUERY_UNMAPPED             4
#define BAM_MATE_UNMAPPED              8
#define BAM_QUERY_REVERSE_COMPLEMENT  16
#define BAM_MATE_REVERSE_COMPLEMENT   32
#define BAM_QUERY_FIRST_MATE          64
#define BAM_QUERY_SECOND_MATE        128
#define BAM_SECONDARY_ALIGNMENT      256
#define BAM_FAILS_PLATFORM_QUALITY   512
#define BAM_PCR_OPTICAL_DUPLICATE   1024

// define some tag lengths
#define MISMATCH_TAG_LEN 7

#endif
