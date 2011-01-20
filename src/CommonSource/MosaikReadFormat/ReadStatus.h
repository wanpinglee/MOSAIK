// ***************************************************************************
// ReadStatus.h - stores the sequencing arrangement: single-end or paire-end.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

typedef unsigned char ReadStatus;

#define RS_UNKNOWN           0 

#define RS_SINGLE_END_READ   1
#define RS_PAIRED_END_READ   2
