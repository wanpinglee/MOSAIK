// ***************************************************************************
// CSRF - imports reads from the SRF file format.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg [Adapted from io_lib: James Bonfield]
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "SRF.h"

// constructor
CSRF::CSRF(void)
: mIsOpen(false)
{
	// note this LUT provides phred base qualities
	for(int i = -128; i < 128; i++)
		mBaseQualityLUT[i + 128] = (int)(-10.0 * log(1.0 / (2.0 + pow(10.0, i / 10.0))) / (log(2.0) + log(5.0)) + 0.5);

	// provides Solexa base qualities
	//mBaseQualityLUT[i + 128] = (int)((10.0 * log(1.0 + pow(10.0, i / 10.0)) / log(10.0) + 0.499));
}

// destructor
CSRF::~CSRF(void) {
	Close();
}

// validates the supplied SRF file
bool CSRF::CheckFile(const string& filename, const bool showError) {

	// read in the first few characters
	char signature[CHECK_BYTES + 1];
	signature[CHECK_BYTES] = 0;
	bool foundError = false;

	// open the SRF file
	FILE* checkStream = NULL;
	errno_t err;

	if((err = fopen_s(&checkStream, filename.c_str(), "rb")) != 0) {
		if(showError) {
			cout << "ERROR: Could not open " << filename << " when validating reference sequence archive." << endl;
			exit(1);
		}

		foundError = true;
	}

	// retrieve the SRF signature
	if(!foundError) {
		fread(signature, CHECK_BYTES, 1, checkStream);

		// check if the magic values match
		if(strncmp(signature, SRF_MAGIC, 4) != 0) {
			if(showError) {
				cout << "ERROR: It seems that the input file (" << filename << ") is not in the SRF format." << endl;
				exit(1);
			}

			foundError = true;
		}
	}

	// check the version
	if(!foundError) {
		
		// check if the versions match
		if(strncmp(signature + 9, SRF_VERSION, 3) != 0) {

			if(showError) {
				cout << "ERROR: It seems that the input file (" << filename << ") was prepared for a different version of the SRF format. Expecting version: " << SRF_VERSION << endl;
				exit(1);
			}

			foundError = true;
		}
	}

	// close the file
	fclose(checkStream);

	// return the appropriate values
	if(foundError) return false;
	return true;
}

// closes the SRF file
void CSRF::Close(void) {
	if(mIsOpen) {
		mIsOpen = false;
		fclose(mSrfData.fp);
	}
}

// gets the next read from the SRF file
bool CSRF::GetRead(Mosaik::Read& mr) {

	if(!mIsOpen) {
		cout << "ERROR: An attempt was made to get reads from an SRF file that hasn't been opened yet." << endl;
		exit(1);
	}

	char name[NAMELEN];
	ZTR_t* ztr = FetchNextTrace(&mSrfData, name);
	if(ztr == NULL) return false;

	// retrieve the bases
	int numChunks;
	ZtrChunk** ztrChunks = FindZtrChunks(ztr, ZTR_TYPE_BASE, &numChunks);

	if(numChunks != 1) {
		cout << "ERROR: Zero or greater than one BASE chunks found." << endl;
		exit(1);
	}

	UncompressChunk(ztr, ztrChunks[0]);
	unsigned int numBases = ztrChunks[0]->dlength - 1;

	mr.Mate1.Bases.Copy(ztrChunks[0]->data + 1, numBases);

	// retrieve the qualities
	ztrChunks = FindZtrChunks(ztr, ZTR_TYPE_CNF4, &numChunks);

	if(numChunks != 1) {
		cout << "ERROR: Zero or greater than one CNF4 chunks found." << endl;
		exit(1);
	}

	UncompressChunk(ztr, ztrChunks[0]);

	mr.Mate1.Qualities.Copy(ztrChunks[0]->data + 1, numBases);

	for(unsigned int i = 0; i < numBases; i++) {
		mr.Mate1.Qualities[i] = mBaseQualityLUT[mr.Mate1.Qualities[i] + 128];
	}

	// replace any periods
	mr.Mate1.Bases.Replace('.', 'N');

	// clean up
	DeleteZtr(ztr);

	// assign our data
	mr.Name = name;

	return true;
}

// opens the SRF file
void CSRF::Open(const string& filename) {

	if(mIsOpen) {
		cout << "ERROR: An attempt was made to open an already open SRF file." << endl;
		exit(1);
	}

	// validate the SRF file
	CheckFile(filename, true);

	FILE* fp = NULL;
	fopen_s(&fp, filename.c_str(), "rb");
	if(fp == NULL) {
		cout << "ERROR: Unable to open the SRF file (" << filename << ") for reading." << endl;
		exit(1);
	}

	mSrfData.fp = fp;
	mIsOpen = true;
}

// Searches for chunks of a specific type.
CSRF::ZtrChunk** CSRF::FindZtrChunks(ZTR_t* ztr, unsigned int type, int* nchunks_p) {
	ZtrChunk **chunks = NULL;
	int nchunks = 0;
	int i;

	for (i = 0; i < ztr->nchunks; i++) {
		if (ztr->chunk[i].type == type) {
			chunks = (ZtrChunk **)realloc(chunks, (nchunks + 1) *
				sizeof(*chunks));
			chunks[nchunks++] = &ztr->chunk[i];
		}
	}
	*nchunks_p = nchunks;
	return chunks;
}

// fetches the next trace from an SRF container as a ZTR object
CSRF::ZTR_t* CSRF::FetchNextTrace(SRF_t* srf, char* name) {

	do {
		int type;
		switch(type = GetBlockType(srf)) {
			case -1:
				/* EOF */
				return NULL;

			case SRFB_CONTAINER:
				if (0 != ReadContainerHeader(srf, &srf->ch)) return NULL;
				break;

			case SRFB_TRACE_HEADER:
				if (0 != ReadTraceHeader(srf, &srf->th)) return NULL;

				srf->mf = mfcreate(NULL, 0);
				if (srf->th.trace_hdr_size)
					mfwrite(srf->th.trace_hdr, 1, srf->th.trace_hdr_size, srf->mf);
				if (srf->ztr)
					DeleteZtr(srf->ztr);
				mfrewind(srf->mf);

				if (NULL != (srf->ztr = DecodePartialZtr(srf, srf->mf, NULL))) {
					srf->mf_pos = mftell(srf->mf);
					mfseek(srf->mf, 0, SEEK_END);
					srf->mf_end = mftell(srf->mf);
				} else {
					/* Maybe not enough to decode or no headerBlob. */
					/* So delay until decoding the body. */
					srf->mf_pos = srf->mf_end = 0;
				}
				break;

			case SRFB_TRACE_BODY: {
				TraceBody tb;
				ZTR_t *ZTR_tmp;

				if(!srf->mf || 0 != ReadTraceBody(srf, &tb, 0)) return NULL;

				if(name) sprintf_s(name, NAMELEN, "%s%s", srf->th.id_prefix, tb.read_id);

				mfseek(srf->mf, srf->mf_end, SEEK_SET);
				if (tb.trace_size) {
					mfwrite(tb.trace, 1, tb.trace_size, srf->mf);
					free(tb.trace);
					tb.trace = NULL;
				}

				mftruncate(srf->mf, mftell(srf->mf));
				mfseek(srf->mf, srf->mf_pos, SEEK_SET);

				if (srf->ztr) {
					ZTR_tmp = DuplicateZtr(srf->ztr); /* inefficient, but simple */
				} else {
					ZTR_tmp = NULL;
				}

				return DecodePartialZtr(srf, srf->mf, ZTR_tmp);
								  }

			default:
				cout << "ERROR: Unknown SRF block type found: " << type << endl;
				exit(1);
				break;
		}
	} while (1);

	return NULL;
}

// uncompresses an individual chunk from all levels of compression
int CSRF::UncompressChunk(ZTR_t* ztr, ZtrChunk* chunk) {
	char *new_data = NULL;
	int new_len;

	while (chunk->dlength > 0 && chunk->data[0] != ZTR_FORM_RAW) {

		switch (chunk->data[0]) {
			case ZTR_FORM_STHUFF:
				new_data = InflateStaticHuffman(ztr, chunk->data, chunk->dlength, &new_len);
				break;
			case ZTR_FORM_XRLE2:
				new_data = ExpandMultiByteRLE(chunk->data, chunk->dlength, &new_len);
				break;
			case ZTR_FORM_QSHIFT:
				new_data = DeinterleaveQualityData(chunk->data, chunk->dlength, &new_len);
				break;
			default:
				cout << "ERROR: Unknown compression algorithm encountered while decompressing ZTR chunk: " << chunk->data[0] << endl;
				exit(1);
				break;
		}

		if(!new_data) return -1;

		chunk->dlength = new_len;
		free(chunk->data);
		chunk->data = new_data;
	}

	return 0;
}

// returns the type of the next block
int CSRF::GetBlockType(SRF_t* srf) {
	int c = fgetc(srf->fp);
	if(c == EOF) return -1;
	ungetc(c, srf->fp);
	return c;
}

// reads a container header and stores the result in 'ch'
int CSRF::ReadContainerHeader(SRF_t* srf, ContainerHeader* ch) {
	char magic[3];
	unsigned int sz;

	if(!ch) return -1;

	/* Check block type */
	if (EOF == (ch->block_type = fgetc(srf->fp))) return -1;
	if (ch->block_type != SRFB_CONTAINER) return -1;

	/* Check magic number && version */
	if (3 != fread(magic, 1, 3, srf->fp)) return -1;
	if (0 != ReadUInt32(srf, &sz))	return -1;
	if (ReadPascalString(srf, ch->version) < 0) return -1;
	if (strncmp(magic, "SRF", 3) || strcmp(ch->version, SRF_VERSION)) return -1;

	/* Containter type, base caller bits */
	if (EOF == (ch->container_type = fgetc(srf->fp)) ||
		ReadPascalString(srf, ch->base_caller) < 0||
		ReadPascalString(srf, ch->base_caller_version) < 0)
		return -1;

	return 0;
}

// reads a data header and stores the result in 'th'
int CSRF::ReadTraceHeader(SRF_t* srf, TraceHeader* th) {
	int z;

	/* Check block type */
	if(EOF == (th->block_type = fgetc(srf->fp))) return -1;

	if(th->block_type != SRFB_TRACE_HEADER) return -1;

	if(0 != ReadUInt32(srf, &th->trace_hdr_size)) return -1;

	th->trace_hdr_size -= 1 + 4 + 1;

	/* Read-id prefix */
	if (EOF == (th->read_prefix_type = fgetc(srf->fp))) return -1;

	if ((z = ReadPascalString(srf, th->id_prefix)) < 0) return -1;
	th->trace_hdr_size -= z+1;

	/* The data header itself */
	if (th->trace_hdr_size) {
		if (th->trace_hdr) free(th->trace_hdr);

		if (NULL == (th->trace_hdr = (unsigned char*)malloc(th->trace_hdr_size))) return -1;

		if (th->trace_hdr_size != fread(th->trace_hdr, 1, th->trace_hdr_size, srf->fp)) {
			free(th->trace_hdr);
			return -1;
		}
	} else {
		th->trace_hdr = NULL;
	}

	return 0;
}

// delete ZTR
void CSRF::DeleteZtr(ZTR_t *ztr) {

	int i;
	if(!ztr) return;

	if(ztr->chunk) {
		for(i = 0; i < ztr->nchunks; i++) {
			if (ztr->chunk[i].data && ztr->chunk[i].ztr_owns) free(ztr->chunk[i].data);
			if (ztr->chunk[i].mdata && ztr->chunk[i].ztr_owns) free(ztr->chunk[i].mdata);
		}
		free(ztr->chunk);
	}

	if(ztr->hcodes) free(ztr->hcodes);
	free(ztr);
}

// decodes a partial ZTR file consisting of data in 'mf'
CSRF::ZTR_t* CSRF::DecodePartialZtr(SRF_t* srf, mFILE* mf, ZTR_t* z) {

	ZTR_t* ztr      = NULL;
	ZtrChunk* chunk = NULL;
	long pos = 0;

	if(z) {

		// Use existing ZTR object => already loaded header
		ztr = z;

	} else {

		// Allocate or use existing ztr
		if((ztr = NewZtr()) == NULL) return NULL;

		// Read the header
		if(ReadZtrHeader(mf, &ztr->header) == -1) {
			if(!z) DeleteZtr(ztr);
			mfrewind(mf);
			return NULL;
		}

		// Check magic number and version
		if(memcmp(ztr->header.magic, ZTR_MAGIC, 8) != 0) {
			if(!z) DeleteZtr(ztr);
			mfrewind(mf);
			return NULL;
		}

		if(ztr->header.version_major != ZTR_VERSION_MAJOR) {
			if(!z) DeleteZtr(ztr);
			mfrewind(mf);
			return NULL;
		}
	}

	// Load chunks
	pos = mftell(mf);
	while((chunk = ReadZtrChunkHeader(mf))) {
		chunk->data = (char *)malloc(chunk->dlength);
		if(chunk->dlength != mfread(chunk->data, 1, chunk->dlength, mf)) break;
		ztr->nchunks++;
		ztr->chunk = (ZtrChunk *)realloc(ztr->chunk, ztr->nchunks * sizeof(ZtrChunk));
		memcpy(&ztr->chunk[ztr->nchunks-1], chunk, sizeof(*chunk));
		free(chunk);
		pos = mftell(mf);
	}

	//int a = ztr->hcodes_checked;
	//int b = ztr->nchunks;
	//int c = ztr->nhcodes;
	//int d = a + b + c;
	//printf("\n");
	//printf("JANE 10 hcodes_checked: %d, nchunks: %d, nhcodes: %d\n", ztr->hcodes_checked, ztr->nchunks, ztr->nhcodes);


	// At this stage we're 'pos' into the mFILE mf with any remainder being
	// a partial block.
	if(ztr->nchunks == 0) {
		if(!z) DeleteZtr(ztr);
		mfrewind(mf);
		return NULL;
	}

	// Ensure we exit at the start of a ztr CHUNK
	mfseek(mf, pos, SEEK_SET);

	// If this is the header part, ensure we uncompress and init. data
	if(!z) {

		// Force caching of huffman code_sets
		SearchCachedHuffmanTables(ztr, CODE_USER);

		// And uncompress the rest
		UncompressZtr(ztr);
	}

	return ztr;
}

// reads a trace header + trace 'blob' and stores the result in 'th'
int CSRF::ReadTraceBody(SRF_t* srf, TraceBody* tb, int no_trace) {
	int z;

	/* Check block type */
	if (EOF == (tb->block_type = fgetc(srf->fp))) return -1;
	if (tb->block_type != SRFB_TRACE_BODY) return -1;

	/* Size */
	if (0 != ReadUInt32(srf, &tb->trace_size)) return -1;
	tb->trace_size -= 6;

	/* Flags */
	if (EOF == (z = fgetc(srf->fp))) return -1;
	tb->flags = z;

	/* Read-id suffix */
	if ((z = ReadPascalString(srf, tb->read_id)) < 0) return -1;
	tb->trace_size -= z+1;

	/* The trace data itself */
	if (!no_trace) {
		if (tb->trace_size) {
			if (NULL == (tb->trace = (unsigned char*)malloc(tb->trace_size)))
				return -1;
			if (tb->trace_size != fread(tb->trace, 1, tb->trace_size,
				srf->fp)) {
					free(tb->trace);
					tb->trace = NULL;
					return -1;
			}
		} else {
			tb->trace = NULL;
		}
	} else {
		/* Skip */
		fseek64(srf->fp, tb->trace_size, SEEK_CUR);
		tb->trace = NULL;
	}

	return 0;
}



// creates a copy of ZTR_t 'src'
CSRF::ZTR_t* CSRF::DuplicateZtr(ZTR_t* src) {

	ZTR_t *dest = NewZtr();
	int i;

	if (!dest)
		return NULL;

	/* Basics */
	*dest = *src;

	/* Mirror chunks */
	dest->chunk = (ZtrChunk*)malloc(src->nchunks * sizeof(ZtrChunk));
	for (i = 0; i < src->nchunks; i++) {
		dest->chunk[i] = src->chunk[i];
		dest->chunk[i].ztr_owns = 0; /* src owns the data/meta_data */
	}

	/* huffman hcodes */
	dest->hcodes = (ztr_hcode_t *)malloc(src->nhcodes * sizeof(ztr_hcode_t));
	for (i = 0; i < src->nhcodes; i++) {
		dest->hcodes[i] = src->hcodes[i];
		dest->hcodes[i].ztr_owns = 0;
	}

	return dest;
}

// Read unsigned 32-bit values in big-endian format
int CSRF::ReadUInt32(SRF_t* srf, unsigned int* val) {
	unsigned char d[4];
	if (1 != fread(d, 4, 1, srf->fp)) return -1;

	*val = (d[0] << 24) | (d[1] << 16) | (d[2] << 8) | (d[3] << 0);
	return 0;
}

// reads a pascal-style string from the srf file
int CSRF::ReadPascalString(SRF_t* srf, char *str) {
	int len;

	if (EOF == (len = fgetc(srf->fp))) return -1;
	if (len != (signed)fread(str, 1, len, srf->fp)) return -1;
	str[len] = '\0';

	return len;
}

// allocates and initialises a ZTR_t structure
CSRF::ZTR_t* CSRF::NewZtr(void) {

	ZTR_t *ztr;

	/* Allocate */
	if (NULL == (ztr = (ZTR_t *)malloc(sizeof(*ztr))))
		return NULL;

	ztr->chunk = NULL;
	ztr->nchunks = 0;

	ztr->nhcodes = 0;
	ztr->hcodes = NULL;
	ztr->hcodes_checked = 0;

	return ztr;
}

// reads a ZTR file header
int CSRF::ReadZtrHeader(mFILE* mf, ZtrHeader* h) {
	if (1 != mfread(h, sizeof(*h), 1, mf)) return -1;
	return 0;
}

// reads a ZTR chunk header and metadata, but not the main data segment
CSRF::ZtrChunk* CSRF::ReadZtrChunkHeader(mFILE* mf) {
	int bei4;
	ZtrChunk *chunk;

	if (NULL == (chunk = (ZtrChunk *)malloc(sizeof(*chunk)))) return NULL;

	/* type */
	if (1 != mfread(&bei4, 4, 1, mf)) {
		free(chunk);
		return NULL;
	}
	chunk->type = be_int4(bei4);

	/* metadata length */
	if (1 != mfread(&bei4, 4, 1, mf)) {
		free(chunk);
		return NULL;
	}
	chunk->mdlength = be_int4(bei4);

	/* metadata */
	chunk->ztr_owns = 1;
	if (chunk->mdlength) {
		if (NULL == (chunk->mdata = (char *)malloc(chunk->mdlength))) {
			free(chunk);
			return NULL;
		}
		if (chunk->mdlength != mfread(chunk->mdata, 1, chunk->mdlength, mf)) {
			free(chunk->mdata);
			free(chunk);
			return NULL;
		}
	} else {
		chunk->mdata = NULL;
	}

	/* data length */
	if (1 != mfread(&bei4, 4, 1, mf)) {
		if (chunk->mdata)
			free(chunk->mdata);
		free(chunk);
		return NULL;
	}
	chunk->dlength = be_int4(bei4);

	return chunk;
}

// Searches through the cached huffman_codeset_t tables looking for a stored huffman code of type 'code_set'
CSRF::ztr_hcode_t* CSRF::SearchCachedHuffmanTables(ZTR_t* ztr, int code_set) {
	int i;

	if (code_set < CODE_USER) return NULL; /* computed on-the-fly or use a hard-coded set */

	/* Check through chunks for undecoded HUFF chunks */
	if (!ztr->hcodes_checked) {
		for (i = 0; i < ztr->nchunks; i++) {
			if (ztr->chunk[i].type == ZTR_TYPE_HUFF) {
				block_t *blk;
				huffman_codeset_t *cs;
				UncompressChunk(ztr, &ztr->chunk[i]);
				blk = CreateBlock((unsigned char *)(ztr->chunk[i].data+2),
					ztr->chunk[i].dlength-2);
				cs = GetHuffmanCodes(blk, NULL);
				if (!cs) {
					DestroyBlock(blk, 1);
					return NULL;
				}
				cs->code_set = (unsigned char)(ztr->chunk[i].data[1]);
				AddUserDefinedHuffmanCodes(ztr, cs, 1);
				DestroyBlock(blk, 1);
			}
		}
		ztr->hcodes_checked = 1;
	}

	/* Check cached copies */
	for (i = 0; i < ztr->nhcodes; i++) {
		if (ztr->hcodes[i].codes->code_set == code_set) {
			return &ztr->hcodes[i];
		}
	}

	return NULL;
}

// uncompresses a ztr (in memory)
int CSRF::UncompressZtr(ZTR_t *ztr) {
	int i;

	for (i = 0; i < ztr->nchunks; i++) 
		UncompressChunk(ztr, &ztr->chunk[i]);

	return 0;
}

// deallocates memory created by CreateBlock()
void CSRF::DestroyBlock(block_t* blk, int keep_data) {
	if (!blk) return;
	if (!keep_data && blk->data) free(blk->data);
	free(blk);
}

// allocates and returns a new block_t struct of a specified default size
CSRF::block_t* CSRF::CreateBlock(unsigned char* data, size_t size) {
	block_t *b = (block_t *)malloc(sizeof(*b));
	if(!b) return NULL;

	b->data = data;
	b->alloc = size;
	b->byte = 0;
	b->bit = 0;

	if(size && !data && NULL == (b->data = (unsigned char*)calloc(size, 1))) {
		free(b);
		return NULL;
	}

	return b;
}

// adds a user-defined huffman_codeset_t code-set to the available code sets used by huffman_decode
CSRF::ztr_hcode_t* CSRF::AddUserDefinedHuffmanCodes(ZTR_t *ztr, huffman_codeset_t* codes, int ztr_owns) {
	if (!codes) return NULL;

	ztr->hcodes = (ztr_hcode_t*)realloc(ztr->hcodes, (ztr->nhcodes+1)*sizeof(*ztr->hcodes));
	ztr->hcodes[ztr->nhcodes].codes = codes;
	ztr->hcodes[ztr->nhcodes].ztr_owns = ztr_owns;

	return &ztr->hcodes[ztr->nhcodes++];
}

// this is the opposite of the store_codes() function
CSRF::huffman_codeset_t* CSRF::GetHuffmanCodes(block_t* block, int* bfinal) {
	int btype;
	huffman_codeset_t *cs;

	/* Header details */
	if (bfinal)
		*bfinal = GetBits(block, 1);
	else
		GetBits(block, 1);
	btype  = GetBits(block, 2);

	cs = (huffman_codeset_t *)malloc(sizeof(*cs));
	cs->code_set = 0;
	cs->blk = NULL;
	cs->bit_num = 0;
	cs->decode_t = NULL;
	cs->decode_J4 = NULL;

	if (btype == 2) {
		/* Standard Deflate algorithm */
		cs->ncodes = 1;
		cs->codes = (huffman_codes_t **)malloc(cs->ncodes*sizeof(*cs->codes));
		cs->codes[0] = GetHuffmanCodesSingle(block);
	} else if (btype == 3) {
		/* Deflate extension - multiple codes */
		int nbits, i;
		nbits  = GetBits(block, 4) + 1;
		cs->ncodes = GetBits(block, nbits) + 1;
		cs->codes = (huffman_codes_t **)malloc(cs->ncodes*sizeof(*cs->codes));
		for (i = 0; i < cs->ncodes; i++) {
			cs->codes[i] = GetHuffmanCodesSingle(block);
		}
	} else {
		fprintf(stderr, "GetHuffmanCodes() only implemented for "
			"BTYPE == DYNAMIC HUFFMAN and INTERLACED HUFFMAN\n");
		return NULL;
	}

	cs->bit_num = block->bit;

	return cs;
}

// reads up to 24-bits worth of data and returns
signed int CSRF::GetBits(block_t* block, int nbits) {
	unsigned int val, bnum = 0;

	if (block->byte*8 + block->bit + nbits > block->alloc * 8) return -1;

	/* Fetch the partial byte of data */
	val = (block->data[block->byte]) >> block->bit;
	bnum = 8 - block->bit;

	/* And additional entire bytes worth as required */
	while (bnum <= (unsigned)nbits) {
		val |= block->data[++block->byte] << bnum;
		bnum += 8;
	}

	block->bit = (block->bit + nbits) % 8;
	return val & ((1 << nbits) - 1);
}

// it restores huffman_codes_t structs from the a serialised data stream
CSRF::huffman_codes_t* CSRF::GetHuffmanCodesSingle(block_t* block) {
	int hlit, hdist, hclen;
	int hclen_order[19] = {
		16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
	};
	int hc_bitlen[19], i;
	huffman_codes_t *bl_cds, *cds;
	int sym, sym_val, last_len;
	int htab[256];

	hlit   = GetBits(block, 5)+257;
	hdist  = GetBits(block, 5)+1;
	hclen  = GetBits(block, 4)+4;

	/*
	fprintf(stderr, "bfinal = %d, btype=%d\n", *bfinal, btype);
	fprintf(stderr, "hlit=0x%x, hdist=0x%x, hclen=0x%x\n",
	hlit, hdist, hclen);
	*/

	/* Read HCLEN code-lengths and construct huffman codes from them */
	for (i = 0; i < hclen; i++)
		hc_bitlen[hclen_order[i]] = GetBits(block, 3);
	for (; i < 19; i++)
		hc_bitlen[hclen_order[i]] = 0;

	bl_cds = (huffman_codes_t *)malloc(sizeof(*bl_cds));
	bl_cds->codes_static = 0;
	bl_cds->ncodes = 0;
	bl_cds->codes = (huffman_code_t *)malloc(19 * sizeof(*bl_cds->codes));
	bl_cds->max_code_len = 7;
	for (i = 0; i < 19; i++) {
		if (hc_bitlen[i]) {
			bl_cds->codes[bl_cds->ncodes].symbol = i;
			bl_cds->codes[bl_cds->ncodes].nbits = hc_bitlen[i];
			bl_cds->ncodes++;
		}
	}

	GenerateCanonicalHuffmanCodes(bl_cds);
	/* output_code_set (stderr, bl_cds); */

	/* Build a lookup table of possible codes to symbols */
	for (i = 0; i < 256; i++)
		htab[i] = -1;
	for (i = 0; i < bl_cds->ncodes; i++) {
		htab[ReverseBitOrder(bl_cds->codes[i].code, bl_cds->codes[i].nbits)
			| (1<<bl_cds->codes[i].nbits)]
		= bl_cds->codes[i].symbol;
	}

	/* Now decode the next HLIT literal codes using bl_cds */
	cds = (huffman_codes_t *)malloc(sizeof(*cds));
	cds->codes_static = 0;
	cds->ncodes = 0;
	cds->codes = (huffman_code_t *)malloc(257 * sizeof(*cds->codes));
	cds->max_code_len = 15;
	sym_val = last_len = 0;
	while ((sym = GetNextSymbol(block, htab)) != -1) {
		int count;
		/* fprintf(stderr, "LIT Sym=%d\n", sym); */

		switch(sym) {
	case 16:
		count = GetBits(block, 2) + 3;
		/* fprintf(stderr, "   +%d\n", count); */
		for (i = 0; i < count; i++) {
			cds->codes[cds->ncodes].symbol = sym_val++;
			cds->codes[cds->ncodes++].nbits = last_len;
		}
		break;

	case 17:
		count = GetBits(block, 3) + 3;
		/* fprintf(stderr, "   +%d\n", count); */
		sym_val += count;
		last_len = 0;
		break;

	case 18:
		count = GetBits(block, 7) + 11;
		/* fprintf(stderr, "   +%d\n", count); */
		sym_val += count;
		last_len = 0;
		break;

	case 0:
		sym_val++;
		last_len = 0;
		break;

	default:
		cds->codes[cds->ncodes].symbol = sym_val++;
		last_len = cds->codes[cds->ncodes++].nbits = sym;
		}

		if (sym_val >= hlit)
			break;
	}
	assert(sym != -1);
	assert(cds->ncodes <= 257);

	/* Skip HDIST codes. Hopefully only 1 of zero length */
	sym_val = 0;
	while ((sym = GetNextSymbol(block, htab)) != -1) {
		/* fprintf(stderr, "DIST Sym=%d\n", sym); */

		switch(sym) {
	case 16:
		sym_val += GetBits(block, 2) + 3;
		break;

	case 17:
		sym_val += GetBits(block, 3) + 3;
		break;

	case 18:
		sym_val += GetBits(block, 7) + 11;
		break;

	default:
		sym_val++;
		}

		if (sym_val >= hdist)
			break;
	}
	assert(sym != -1);

	//huffman_codes_destroy(bl_cds);
	GenerateCanonicalHuffmanCodes(cds);
	/* output_code_set(stderr, cds); */

	return cds;
}

// generates canonical huffman codes given a set of symbol bit lengths
int CSRF::GenerateCanonicalHuffmanCodes(huffman_codes_t* c) {
	int i, j;
	unsigned int code, last_len;
	int clens[33];
	int offs[33];
	huffman_code_t ctmp[258];
	signed int symtab[258];

	/* Sort by bit-length, subfield symbol - much faster than qsort() */
	for (i = 0; i < 258; i++)
		symtab[i] = -1;
	for (i = 0; i < c->ncodes; i++)
		symtab[c->codes[i].symbol] = i;
	for (i = 0; i <= 32; i++)
		offs[i] = clens[i] = 0;
	for (i = 0; i < c->ncodes; i++)
		clens[c->codes[i].nbits]++;
	for (i = 1; i <= 32; i++)
		offs[i] = offs[i-1] + clens[i-1];
	for (i = 0; i < 258; i++) {
		if (symtab[i] != -1)
			ctmp[offs[c->codes[symtab[i]].nbits]++] = c->codes[symtab[i]];
	}
	memcpy(c->codes, ctmp, c->ncodes * sizeof(huffman_code_t));

	/*
	* Force all codes to be <= max_code_len. This is needed due to the
	* 15-bit length limitation of Deflate literal codes and the 7-bit 
	* limit of the code bit-length table.
	*/
	/* Find first point of failure */
	for (i = 0; i < c->ncodes; i++) {
		if (c->codes[i].nbits > c->max_code_len)
			break;
	}

	/*
	* From here on we shrink the length of the current code by increasing
	* the length of an earlier symbol, at last_code.
	*/
	if (i != c->ncodes) {
		int delta = 0;

		/*
		fprintf(stderr, "=== REORDERING %d ===\n", c->code_set);
		output_code_set(stderr, c);
		output_code_set2(stderr, c);
		*/

		for (; i < c->ncodes; i++) {
			int k, cur_len;

			c->codes[i].nbits -= delta;
			if (c->codes[i].nbits <= c->max_code_len)
				continue;

			for (j = i; j >= 0 && c->codes[j].nbits >= c->max_code_len; j--)
				;
			if (j < 0) {
				fprintf(stderr,
					"Too many symbols to fit in bit-length requirements\n");
				fprintf(stderr, "=== FAILING ===\n");
				//output_code_set(stderr, c);
				//output_code_set2(stderr, c);
				abort();
			}

			/*
			fprintf(stderr, "Changing code %d/%d to len %d\n",
			c->codes[i].symbol, c->codes[j].symbol,
			c->codes[j].nbits+1);
			*/
			cur_len = c->codes[i].nbits;
			c->codes[i].nbits = ++c->codes[j].nbits;

			/*
			* Shrink the next code by one, or if none at that bit-length
			* the next 2, and so on
			*/
			delta = 1;
			for (k = i+1; delta && k < c->ncodes; k++) {
				while (c->codes[k].nbits > cur_len) {
					delta *= 2;
					cur_len++;
				}
				c->codes[k].nbits--;
				delta--;
			}
			assert(delta == 0);
		}

		/*
		fprintf(stderr, "=== REORDERED TO %d ===\n", c->code_set);
		output_code_set(stderr, c);
		output_code_set2(stderr, c);
		*/

		/* Ordering is shot - regenerate via brute force way */
		return GenerateCanonicalHuffmanCodes(c);
	}


	/* Generate codes */
	code = last_len = 0; /* stop warning */
	for (i = 0; i < c->ncodes; i++) {
		int nbits = c->codes[i].nbits;

		if (i == 0) {
			code = 0;
			last_len = nbits;
		} else {
			code++;
		}
		if ((unsigned)nbits > last_len) {
			code <<= (nbits - last_len);
			last_len = nbits;
		}
		c->codes[i].code = ReverseBitOrder(code, nbits);
	}

	/* Reindex so the symbol is the primary index into codes */
	for (i = 0; i <= 257; i++) {
		c->lookup[i].nbits = 0;
	}
	for (i = 0; i < c->ncodes; i++) {
		c->lookup[c->codes[i].symbol] = c->codes[i];
	}

	return 0;
}

// reverses the order of bits in the bottom nbits of val
unsigned int CSRF::ReverseBitOrder(unsigned int val, int nbits) {

	unsigned int rev = 0;
	for(int i = 0; i < nbits; i++) {
		rev = (rev << 1) | (val & 1);
		val >>= 1;
	}

	return rev;
}

// A slow version of the above huffman_decode function
int CSRF::GetNextSymbol(block_t* in, int* htab) {
	int b, v = 0, c = 1;
	while ((b = GetBits(in, 1)) != -1) {
		v = (v<<1) | b | (c <<= 1);
		if (htab[v] != -1) return htab[v];
	}
	return -1;
}


// stores nbytes bytes, padding to align on the next byte boundary
void CSRF::StoreBytes(block_t* block, unsigned char* val, int nbytes) {

	// Align
	if (block->bit != 0) {
		block->byte++;
		block->bit = 0;
	}

	// Resize
	ResizeBlock(block, block->byte + nbytes + 1);

	// Store
	memcpy(&block->data[block->byte], val, nbytes);
	block->byte += nbytes;
}

// Decode a huffman stream from 'block' using huffman codes 'c'
CSRF::block_t* CSRF::DecodeHuffmanStream(block_t* in, huffman_codeset_t* cs) {
	block_t *out = NULL;
	int j;
	unsigned int i;
	int node_num;
	unsigned char *cp;
	huffman_codes_t **c;
	int nc;
	h_jump4_t (*J4)[16];
	htree_t *t;

	if (!cs)
		return NULL;
	c = cs->codes;
	nc = cs->ncodes;

	/* Ensure precomputed lookup tables exist */
	if (!cs->decode_t || !cs->decode_J4)
		if (-1 == InitializeDecodeHuffmanTables(cs))
			return NULL;

	t  = cs->decode_t;
	J4 = cs->decode_J4;

	if (NULL == (out = CreateBlock(NULL, 9*(in->alloc+1)))) {
		goto error;
	}

	/*
	* Decoding - part 1
	* We're part way through a byte, so decode bit by bit up to the next
	* whole byte and then we start the fast decoding section.
	*/
	cp = out->data;
	node_num = 0;
	while (in->bit != 0) {
		int b = GetBits(in, 1);
		htree_t *t2 = &t[node_num];
		if (t2->l[b] != -1) {
			if (t2->l[b] != SYM_EOF) {
				*cp++ = (unsigned char)t2->l[b];
			} else {
				out->byte = cp - out->data;
				goto success;
			}
		}
		node_num = t2->c[b];
	}

	/*
	* Decoding - part 2
	*
	* We now handle data nibble by nibble using the nibble to get an
	* h_jump4_t lookup from the J4[] table.
	* If top_bit is clear then we know we have no funny business (SYM_EOF)
	* so we use a fast decoding technique, otherwise we have to do a slower
	* loop with a check.
	*/
	{
		int last_node = node_num;
		unsigned char *last_cp = cp;
		h_jump4_t *x = &J4[node_num][in->data[in->byte] & 0x0f];
		int l = x->nsymbols;
		int b;

		/*
		* This is the tight loop, so we over-optimise here by ignoring EOF
		* and relying on knowing the length of the input data stream.
		* This allows us to ignore the 9-bit data and only operate on
		* the basic 0-255 symbols, glossing over the minor issue that EOF
		* will look like an ordinary symbol.
		*/
		for (i = in->byte; i < in->alloc; i++) {
			last_cp = cp;
			last_node = node_num;

			x = &J4[node_num][in->data[i] & 0x0f];
			l = x->nsymbols;

			/* printf("val=%d\n", in->data[i] & 0x0f); */
			for (j = 0; j < l; j++) {
				*cp++ = x->symbol[j];
			}
			node_num = x->jump;

			if (x->top_bit)
				break;

			x = &J4[node_num][(in->data[i] >> 4) & 0x0f];
			l = x->nsymbols;

			for (j = 0; j < l; j++) {	
				*cp++ = x->symbol[j];
			}
			node_num = x->jump;

			if (x->top_bit)
				break;
		}


		/*
		* Decoding - part 3
		*
		* The above optimisation has unfortunately added EOF to our data
		* along with whatever else was packed in the last byte after the
		* EOF symbol. So we rewind one byte and finish off decoding
		* the slow way - walking the tree.
		*/
		cp = last_cp;
		node_num = last_node;
		in->byte = i;
		in->bit = 0;
		while (-1 != (b = GetBits(in, 1))) {
			htree_t *t2 = &t[node_num];
			if (t2->l[b] != -1) {
				if (t2->l[b] != SYM_EOF) {
					*cp++ = (unsigned char)t2->l[b];
				} else {
					out->byte = cp - out->data;
					goto success;
				}
			}
			node_num = t2->c[b];
		}
	}

success:
	return out;

error:
	if (out)
		DestroyBlock(out, 0);

	return NULL;
}

// ensures a block_t holds at least 'size' bytes
int CSRF::ResizeBlock(block_t* blk, size_t size) {
	unsigned char *newp = NULL;

	if(!blk) return -1;

	/* Grow size to next power of 2, if we're growing */
	if(size > blk->alloc) {
		size--;
		size |= size >> 1;
		size |= size >> 2;
		size |= size >> 4;
		size |= size >> 8;
		size |= size >> 16;
		size++;
	}

	if(NULL == (newp = (unsigned char*)realloc(blk->data, size))) return -1;
	else blk->data = newp;

	if(size > blk->alloc) memset(&blk->data[blk->alloc], 0, size - blk->alloc);
	blk->alloc = size;

	return 0;
}

//
int CSRF::InitializeDecodeHuffmanTables(huffman_codeset_t* cs) {
	int nnodes, i, j, n, nc;
	huffman_codes_t **c;
	int new_node, rec;
	h_jump4_t (*J4)[16] = NULL;
	htree_t *t;

	c = cs->codes;
	nc = cs->ncodes;

	/* Allocate memory for internal nodes (nsyms-1 for each code set) */
	for (nnodes = i = 0; i < nc; i++) {
		nnodes += c[i]->ncodes-1;
	}

	if (NULL == (t = (htree_t *)malloc(nnodes * sizeof(*t))))
		goto error;

	if (NULL == (J4 = (h_jump4_t (*)[16])malloc(nnodes * sizeof(*J4))))
		goto error;

	/*
	* Construct the tree from the codes.
	* We have one tree for all 'nc' huffman codes with each tree pointing
	* to the root of the next one (or first) tree whenever we emit a
	* symbol.
	*
	* This then effectively means the decoding step is identical to the
	* single huffman code function.
	*/
	new_node = 0;
	for (rec = 0; rec < nc; rec++) {
		int root = new_node++;
		int next_root = rec == nc-1
			? 0
			: root + c[rec]->ncodes-1;

		t[root].l[0] = t[root].l[1] = -1;
		t[root].c[0] = t[root].c[1] = next_root;
		for (i = 0; i < c[rec]->ncodes; i++) {
			int n = root;
			unsigned int v = c[rec]->codes[i].code;

			for (j = 0; j < c[rec]->codes[i].nbits-1; j++) {
				int b = v & 1;
				if (t[n].c[b] != next_root) {
					n = t[n].c[b];
				} else {
					n = (t[n].c[b] = new_node++);
					t[n].c[0] = t[n].c[1] = next_root;
					t[n].l[0] = t[n].l[1] = -1;
				}
				v >>= 1;
			}
			/* last bit */
			t[n].l[v & 1] = c[rec]->codes[i].symbol;
		}
	}

	/*
	for (i = 0; i < new_node; i++) {
	printf("t[%d] = {(%d,%d), (%d,%d)}\n",
	i,
	t[i].l[0], t[i].l[1],
	t[i].c[0], t[i].c[1]);
	}
	*/

	/* Build the 16 wide lookup table per node */
	for (n = 0; n < new_node; n++) {
		for (j = 0; j < 16; j++) {
			unsigned int v = j;
			int n2 = n;
			h_jump4_t *hj = &J4[n][j];
			hj->nsymbols = 0;
			hj->top_bit = 0;
			for (i = 0; i < 4; i++) {
				int b = v & 1;
				if (t[n2].l[b] >= 0) {
					hj->symbol[hj->nsymbols++] = (unsigned char)t[n2].l[b];
					if (t[n2].l[b] == SYM_EOF)
						if (!hj->top_bit)
							hj->top_bit |= 1 << (hj->nsymbols-1);
				}
				n2 = t[n2].c[b];
				v >>= 1;
			}
			hj->jump = n2;
			/*
			printf("J4[%d][%d] = {'%.*s', %d}\n",
			n, j, hj->nsymbols, hj->symbol, n2);
			*/
		}
	}

	cs->decode_t = t;
	cs->decode_J4 = J4;

	return 0;

error:
	if (t)
		free(t);

	if (J4)
		free(J4);

	cs->decode_t = NULL;
	cs->decode_J4 = NULL;

	return -1;
}


// =====================
// memory file operators
// =====================

// for creating existing mFILE pointers directly from memory buffers
CSRF::mFILE* CSRF::mfcreate(char* data, int size) {
	mFILE *mf = (mFILE *)malloc(sizeof(*mf));
	mf->fp = NULL;
	mf->data = data;
	mf->alloced = size;
	mf->size = size;
	mf->eof = 0;
	mf->offset = 0;
	mf->flush_pos = 0;
	mf->mode = MF_READ | MF_WRITE;
	return mf;
}

// memory fread
size_t CSRF::mfread(void* ptr, size_t size, size_t nmemb, mFILE* mf) {
	size_t len;
	char *cptr = (char *)ptr;

	if (mf->size <= mf->offset) return 0;
	len = size * nmemb <= mf->size - mf->offset ? size * nmemb : mf->size - mf->offset;
	if (!size) return 0;

	memcpy(cptr, &mf->data[mf->offset], len);
	mf->offset += len;
	cptr += len;

	if (len != size * nmemb) mf->eof = 1;

	return len / size;
}

// memory rewind
void CSRF::mfrewind(mFILE* mf) {
	mf->offset = 0;
	mf->eof    = 0;
}

// memory fseek
int CSRF::mfseek(mFILE* mf, long offset, int whence) {
	switch (whence) {
	case SEEK_SET:
		mf->offset = offset;
		break;
	case SEEK_CUR:
		mf->offset += offset;
		break;
	case SEEK_END:
		mf->offset = mf->size + offset;
		break;
	default:
		errno = EINVAL;
		return -1;
	}

	mf->eof = 0;
	return 0;
}

// memory ftell
long CSRF::mftell(mFILE* mf) {
	return mf->offset;
}

// memory ftruncate
void CSRF::mftruncate(mFILE* mf, long offset) {
	mf->size = offset != -1 ? offset : mf->offset;
	if (mf->offset > mf->size) mf->offset = mf->size;
}

// memory fwrite
size_t CSRF::mfwrite(void* ptr, size_t size, size_t nmemb, mFILE* mf) {

	if(!(mf->mode & MF_WRITE)) return 0;

	/* Append mode => forced all writes to end of file */
	if(mf->mode & MF_APP) mf->offset = mf->size;

	/* Make sure we have enough room */
	while (size * nmemb + mf->offset > mf->alloced) {
		mf->alloced = mf->alloced ? mf->alloced * 2 : 1024;
		mf->data = (char *)realloc(mf->data, mf->alloced);
	}

	/* Record where we need to reflush from */
	if (mf->offset < mf->flush_pos) mf->flush_pos = mf->offset;

	/* Copy the data over */
	memcpy(&mf->data[mf->offset], ptr, size * nmemb);
	mf->offset += size * nmemb;
	if(mf->size < mf->offset) mf->size = mf->offset;

	return nmemb;
}

// ====================
// compression routines
// ====================

// reorders quality data from an interleaved 4-byte aligned format to its RAW format 
char* CSRF::DeinterleaveQualityData(char* qold, int qlen, int* new_len) {

	// Correct input is 4x (nbases+1) bytes
	if(qlen%4 != 0 || *qold != ZTR_FORM_QSHIFT) return NULL;

	int nbases = qlen / 4 - 1;
	char *qnew = (char *)malloc(nbases*4+1);
	qnew[0] = 0; // raw byte

	for(int i = 0, j = 4; i < nbases; i++) {
		qnew[1+i]            = qold[j++];
		qnew[1+nbases+i*3]   = qold[j++];
		qnew[1+nbases+i*3+1] = qold[j++];
		qnew[1+nbases+i*3+2] = qold[j++];
	}

	*new_len = nbases*4+1;
	return qnew;
}

// implements decompression using a set of static huffman codes stored using the Deflate algorithm
char* CSRF::InflateStaticHuffman(ZTR_t* ztr, char* comp, int comp_len, int* uncomp_len) {
	int cset = (unsigned char)(comp[1]);
	huffman_codeset_t *cs = NULL, *cs_free = NULL;
	block_t *blk_in = CreateBlock(NULL, comp_len), *blk_out = CreateBlock(NULL, 1000);
	int bfinal = 1, bit_num = 0;
	char *uncomp = NULL;

	if (cset >= CODE_USER) {
		/* Scans through HUFF chunks */
		ztr_hcode_t* hc = SearchCachedHuffmanTables(ztr, cset);
		if (!hc) return NULL;

		cs = hc->codes;
		bit_num = cs->bit_num;
		blk_in->bit = 0;
	} else if (cset > 0) {
		/* Create some temporary huffman_codes to stringify */
		//cs_free = cs = generate_code_set(cset, 1, NULL, 0, 1, MAX_CODE_LEN, 0);
		if (!cs) return NULL;

		bit_num = cs->bit_num;
		blk_in->bit = 0;
	} /* else inline codes */
	/*
	* We need to know at what bit the huffman codes would have ended on
	* so we can store our huffman encoded symbols immediately following it.
	* For speed though this bit-number is cached.
	*/
	blk_in->data[blk_in->byte++] |= *(comp+2);
	StoreBytes(blk_in, (unsigned char *)comp+3, comp_len-3);


	/* Rewind */
	blk_in->byte = 0;
	blk_in->bit = bit_num;

	do {
		block_t *out = NULL;

		/*
		* We're either at the start of a block with codes to restore
		* (cset == INLINE or the 2nd onwards block) or we've already
		* got some codes in cs and we're at the position where huffman
		* encoded symbols are stored. 
		*/
		if (!cs)
			if (NULL == (cs = cs_free = GetHuffmanCodes(blk_in, &bfinal)))
				return NULL;

		/*
		{int i;
		for (i = 0; i < cs->ncodes; i++) {
		output_code_set(stderr, cs->codes[i]);
		}}
		*/

		if (NULL == (out = DecodeHuffmanStream(blk_in, cs))) {
			//huffman_codeset_destroy(cs);
			return NULL;
		}

		/* Could optimise this for the common case of only 1 block */
		StoreBytes(blk_out, out->data, out->byte);
		DestroyBlock(out, 0);
		//if (cs_free)
		//huffman_codeset_destroy(cs_free);
		cs = cs_free = NULL;
	} while (!bfinal);

	*uncomp_len = blk_out->byte;
	uncomp = (char *)blk_out->data;

	DestroyBlock(blk_in, 0);
	DestroyBlock(blk_out, 1);

	return uncomp;
}

// reverses multi-byte run length encoding
char* CSRF::ExpandMultiByteRLE(char* comp, int comp_len, int* uncomp_len) {
	char *out, *last;
	int out_len, out_alloc, rsz, i, j, run_len;

	out_alloc = comp_len*2; /* just an estimate */
	out_len = 0;
	if (NULL == (out = (char *)malloc(out_alloc)))
		return NULL;

	if (*comp++ != ZTR_FORM_XRLE2)
		return NULL;

	/* Read rsz and swallow padding */
	rsz = *comp++;
	comp_len -= 2;
	for (i = 2; i < rsz; i++) {
		comp++;
		comp_len--;
	}

	/* Uncompress */
	run_len = 0;
	last = comp;
	for (i = 0; i < comp_len;) {
		while (out_len + rsz > out_alloc) {
			out_alloc *= 2;
			if (NULL == (out = (char *)realloc(out, out_alloc)))
				return NULL;
		}
		memcpy(&out[out_len], &comp[i], rsz);

		if (memcmp(&out[out_len], last, rsz) == 0) {
			run_len++;
		} else {
			run_len = 1;
		}

		i += rsz;
		out_len += rsz;

		if (run_len >= 2) {
			/* Count remaining copies */
			run_len = (unsigned char)comp[i];
			i += rsz;

			while (out_len + run_len * rsz > out_alloc) {
				out_alloc *= 2;
				if (NULL == (out = (char *)realloc(out, out_alloc)))
					return NULL;
			}

			for (j = 0; j < run_len; j++) {
				memcpy(&out[out_len], last, rsz);
				out_len += rsz;
			}
			run_len = 0;
		}

		last = &comp[i-rsz];
	}

	/* Shrink back down to avoid excessive memory usage */
	out = (char*)realloc(out, out_len);
	*uncomp_len = out_len;

	return out;
}
