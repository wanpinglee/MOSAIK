// ***************************************************************************
// SequencingTechnologies.h - stores the internal codes relating to different
//                            sequencing technologies.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

typedef unsigned short SequencingTechnologies;

// we have space for 10 additional sequencing technologies
#define ST_UNKNOWN               0
#define ST_454                   1
#define ST_HELICOS               2
#define ST_ILLUMINA              4
#define ST_PACIFIC_BIOSCIENCES   8
#define ST_SOLID                16
#define ST_SANGER               32
