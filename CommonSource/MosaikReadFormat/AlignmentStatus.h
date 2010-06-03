// ***************************************************************************
// AlignmentStatus.h - stores the flags used in our alignments.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

// define our alignment status flags (relevant before MosaikMerge)
typedef unsigned char AlignmentStatus;
#define AS_UNKNOWN                      0   // specifies an unset status flag
#define AS_SINGLE_END_READ              1   // transferred from the read format
#define AS_PAIRED_END_READ              2   // transferred from the read format
#define AS_UNSORTED_READ                4   // expected in MosaikAligner data
#define AS_SORTED_ALIGNMENT             8   // expected in MosaikSort data
#define AS_ALL_MODE                     16  // enables non-unique PE resolution
#define AS_UNIQUE_MODE                  32  // disables non-unique PE resolution
#define AS_RESERVED1                    64  // reserved, not currently used
#define AS_RESERVED2                    128 // reserved, not currently used

// define our read flags
#define RF_UNKNOWN                      0   // specifies an unset read flag
#define RF_FAILED_QUALITY_CHECK         1   // reserved, not currently used
#define RF_HAVE_MATE1                   2   // specifies if mate1 is stored
#define RF_HAVE_MATE2                   4   // specifies if mate2 is stored
#define RF_IS_LONG_READ                 8   // specifies any read that is longer than 255 bases
#define RF_IS_PAIRED_IN_SEQUENCING      16  // originates from a paired-end data set
#define RF_IS_PCR_OR_OPTICAL_DUPLICATE  32  // reserved, not currently used
#define RF_IS_UNALIGNED                 64  // reserved, not currently used
#define RF_RESOLVED_AS_PAIR             128 // seen only after MosaikSort or MosaikMerge

// define our alignment flags
#define AF_UNKNOWN                      0   // specifies an unset alignment flag
#define AF_IS_FIRST_MATE                1   // reserved, not currently used
#define AF_IS_MATE_REVERSE_STRAND       2   // specifies the orientation of the mate
#define AF_IS_NOT_PRIMARY               4   // reserved, not currently used
#define AF_IS_REVERSE_STRAND            8   // specifies the orientation of the current alignment
#define AF_IS_SECOND_MATE               16  // reserved, not currently used
#define AF_WAS_RESCUED                  32  // specifies if the alignment was rescued during local alignment search
#define AF_RESERVED1                    64  // reserved, not currently used
#define AF_RESERVED2                    128 // reserved, not currently used

// define our alignment tags
#define AT_UNKNOWN                      0

// define our header tags
#define HT_UNKNOWN                      0

// define our reference sequence tags
#define RST_UNKNOWN                     0

// define our read group tags
#define RGT_UNKNOWN                     0

// define our tag data structure
typedef unsigned char TagType;

const TagType TT_CHAR   = 'c';
const TagType TT_DOUBLE = 'd';
const TagType TT_FLOAT  = 'f';
const TagType TT_INT16  = 's';
const TagType TT_INT32  = 'i';
const TagType TT_INT64  = 'l';
const TagType TT_STRING = 'z';
const TagType TT_UCHAR  = 'C';
const TagType TT_UINT16 = 'S';
const TagType TT_UINT32 = 'I';
const TagType TT_UINT64 = 'L';

#define TAG_STRING_LEN 512

struct Tag {
	unsigned char ID;
	TagType Type;
	union {
		char Char;
		char String[TAG_STRING_LEN];
		double Double;
		float Float;
		int Int32;
		long long Int64;
		short Int16;
		unsigned char UChar;
		unsigned int UInt32;
		unsigned long long UInt64;
		unsigned short UInt16;
	};
};
