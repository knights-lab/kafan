//#include <windows.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#ifdef WIN32
	#define fseek _fseeki64
	#define ftell _ftelli64
#endif
#pragma GCC push_options
#pragma GCC optimize ("-O2")
static size_t enmuleFb(uint8_t *in, size_t len, uint8_t *out) { // make TCF
	static const uint16_t N = 256, bel = UINT16_MAX;
	typedef struct {uint8_t Child[128];} TN;
	TN *Node = calloc(N+1,sizeof(*Node)); // dict init
	uint16_t loc = 0, nxt = 1, lastN = 0, temp; 
	for (; nxt < 129; ++nxt) Node[lastN].Child[nxt-1] = nxt;
	uint8_t *buffer = malloc(N), *orig = out;
	for (size_t i = 0; i < len; ++i) {
		register uint8_t c = in[i];
		buffer[loc] = c; // consider next byte
		if (c > 127) {
			printf("\nAlright, look. Level '0' is for text. ASCII-7 *text*\n");
			printf(" --> Found non-text at byte %llu (char %u). Bye.\n",i,c);
			exit(1);
		}
		else if (!(temp=Node[lastN].Child[buffer[loc++]])) {
			if (nxt != N) Node[lastN].Child[c]=nxt++;
			*out++ = lastN-1; // store last good
			*buffer = c, loc = 1; // reset heads
			lastN = c+1;
		}
		else lastN = temp; // shortcut to next good
	}
	*out++ = lastN-1; 
	free(Node), free(buffer); 
	return out - orig;
}
static size_t demuleFb(uint8_t *in, size_t inSz, uint8_t *out) {
	static const uint16_t N = 256, dictSz = N+1;
	char **dct = calloc(dictSz,sizeof(*dct));
	uint16_t *lens = calloc(dictSz,sizeof(*lens));
	uint16_t nxt = 0, prv = *in; 
	for ( ; nxt < 128; ++nxt) 
		dct[nxt] = malloc(8), 
		*dct[nxt] = nxt, lens[nxt] = 1;
	*out = prv;
	size_t strIx = 1, j = 1;
	for (; j < inSz && nxt != N; ++j) {
		uint16_t ix = in[j];
		lens[nxt] = lens[prv] + 1;
		dct[nxt] = malloc(lens[nxt]+8);
		dct[nxt] = memcpy(dct[nxt],dct[prv],lens[prv]);
		dct[nxt][lens[prv]] = *dct[dct[ix] ? ix : prv];
		++nxt, prv = ix;
		memcpy(out+strIx,dct[ix],lens[ix]);
		strIx += lens[ix];
	}
	for (; j < inSz; ++j) {
		uint16_t ix = in[j]; 
		uint64_t *D = (uint64_t*)(out+strIx), *S = (uint64_t*)dct[ix];
		for (int x = 0; x < lens[ix]; x+=8) *D++ = *S++;
		strIx += lens[ix];
	}
	for (size_t i = 0; i < dictSz; ++i) free(dct[i]);
	free(dct), free(lens); 
	return strIx; 
}
char NodeQ[(UINT16_MAX+1)*512]; 
uint16_t outBuf[UINT16_MAX+2];
uint8_t buffer[UINT16_MAX+3];
#define ENKFN(b) { \
	static const uint32_t N = (1 << b), \
		endcap = N-1, bel = UINT16_MAX; \
	typedef struct {uint16_t Child[256];} TN; \
	TN *Node = memset(NodeQ,0,(N+1)*sizeof(*Node)); /* calloc(N+1,sizeof(*Node)); */ \
	uint16_t lastN = 0, temp; \
	uint32_t nxt = 1; for (; nxt < 257; ++nxt) \
		Node[lastN].Child[nxt-1] = nxt; \
	uint32_t loc = 0; size_t bix = 0, bix2 = 0; \
	/* uint8_t *buffer = malloc(N+2); */ \
	/* uint16_t *outBuf = malloc((bel+2)*sizeof(*outBuf)); */ \
	int bin = 0; \
	uint64_t mask = 0, i = 0; \
	while (i < len) { \
		if (nxt != N) for (;  i < len && bix != bel; ++i ) { \
			buffer[loc] = in[i]; \
			if ((temp=Node[lastN].Child[buffer[loc++]])==0) { \
				if (nxt != N) Node[lastN].Child[in[i]]=nxt++; \
				/*printf("adding symbol: %u, length = %u\n",nxt-1,loc);*/ \
				outBuf[bix++] = lastN - 1; \
				*buffer = in[i], loc = 1; \
				lastN = in[i] + 1; \
			} \
			else lastN = temp; \
		} \
		else for (;  i < len && bix != bel; ++i ) { \
			buffer[loc] = in[i]; \
			if ((temp=Node[lastN].Child[buffer[loc++]])==0) \
				outBuf[bix++] = lastN - 1, \
				*buffer = in[i], loc = 1, \
				lastN = in[i] + 1; \
			else lastN = temp; \
		} \
		if (bel != bix) break; \
		for (int x = 0; x < bel; ++x) { \
			mask |= (uint64_t)outBuf[x] << bin; \
			bin += b; \
			if (bin >= 64) { \
				out[bix2++] = mask; \
				mask = (uint64_t)outBuf[x] >> (64 + b) - bin; \
				bin -= 64; \
			} \
		} \
		bix = 0; \
	} \
	outBuf[bix++] = lastN - 1; \
	outBuf[bix++] = endcap; \
	for (int x = 0; x < bix; ++x) { \
		mask |= (uint64_t)outBuf[x] << bin; \
		bin += b; \
		if (bin >= 64) { \
			out[bix2++] = mask; \
			mask = (uint64_t)outBuf[x] >> ((64 + b) - bin); \
			bin -= 64; \
		} \
	} \
	if (bin) out[bix2++] = mask; \
	/* free(Node), free(buffer), free(outBuf); */ \
	return bix2*8; \
}

static inline size_t enkfn9(uint8_t *in, size_t len, uint64_t *out) ENKFN(9)
//static inline size_t enkfn10(uint8_t *in, size_t len, uint64_t *out) ENKFN(10)
static inline size_t enkfn11(uint8_t *in, size_t len, uint64_t *out) ENKFN(11)
//static inline size_t enkfn12(uint8_t *in, size_t len, uint64_t *out) ENKFN(12)
static inline size_t enkfn13(uint8_t *in, size_t len, uint64_t *out) ENKFN(13)
static inline size_t enkfn14(uint8_t *in, size_t len, uint64_t *out) ENKFN(14)
//static inline size_t enkfn15(uint8_t *in, size_t len, uint64_t *out) ENKFN(15)
static inline size_t enkfn16(uint8_t *in, size_t len, uint64_t *out) ENKFN(16)

// fast KFN decomp
#define DEKFNFALLC(b) { \
	static const uint32_t N = (1 << b) - 1, \
		bel = (UINT16_MAX/(b*64))*(b*64), b8 = bel*b/64; \
	uint64_t *inB = in, *D, *S, strIx = 1; \
	uint16_t IB = 64-b, code = (*inB << IB) >> IB; \
	if (code == N) return 0; \
	char **dct = calloc(N+1,sizeof(*dct)); \
	uint16_t *lens = calloc(N+1,sizeof(*lens)); \
	uint16_t *buf = malloc(bel*sizeof(*buf)); \
	uint16_t i=1, k=1, j=0, shf, bin=b, prv = code, nxt = 0; \
	for ( ; nxt < 256; ++nxt) \
		dct[nxt] = malloc(8), \
		*dct[nxt] = nxt, lens[nxt] = 1; \
	*out = code; \
	while (code != N) { \
		for (; i < bel && code != N; ++i) { \
			register uint16_t slc = bin > IB; \
			shf = (slc ? IB + 64 : IB) - bin, j += slc; \
			register uint16_t ps = (inB[j] << shf) >> IB; \
			code = !slc | (bin == 64) ? ps : (inB[j-1] >> bin) | ps; \
			bin = slc ? bin - IB : bin + b; \
			buf[i] = code; \
		} \
		for (; k < i; ++k) { \
			if (nxt < N) \
				lens[nxt] = lens[prv] + 1, \
				dct[nxt] = malloc(lens[nxt] + 8), \
				dct[nxt] = memcpy(dct[nxt],dct[prv],lens[prv]), \
				dct[nxt++][lens[prv]] = *dct[dct[buf[k]] ? buf[k] : prv]; \
			prv = buf[k]; \
			D = (uint64_t *)(out+strIx), S = (uint64_t *)dct[buf[k]]; \
			for (int x = 0; x < lens[buf[k]]; x+=8) *D++=*S++; \
			strIx += lens[buf[k]]; \
		} \
		inB += b8; \
		i = j = k = bin = 0; \
	} \
	for (int i = 0; i < N+1; ++i) free(dct[i]); \
	free(dct), free(lens), free(buf); \
	return strIx; \
}
char *dct[UINT16_MAX+1];
uint16_t lens[UINT16_MAX+1], buf[UINT16_MAX+1];

#define DEKFNF(b) { \
	static const uint32_t N = (1 << b) - 1, \
		bel = (UINT16_MAX/(b*64))*(b*64), b8 = bel*b/64; \
	uint64_t *inB = in, *D, *S, strIx = 1; \
	static const uint16_t IB = 64-b; \
	uint16_t code = (*inB << IB) >> IB; \
	uint16_t i=1, k=1, j=0, shf, bin=b, prv = code, nxt = 256; \
	lens[nxt] = 0, *out = code; \
	while (code != N) { \
		for (; i < bel && code != N; ++i) { \
			register uint16_t slc = bin > IB; \
			shf = (slc ? IB + 64 : IB) - bin, j += slc; \
			register uint16_t ps = (inB[j] << shf) >> IB; \
			code = !slc | (bin == 64) ? ps : (inB[j-1] >> bin) | ps; \
			bin = slc ? bin - IB : bin + b; \
			buf[i] = code; \
		} \
		for (; k < i; ++k) { \
			if (nxt < N) \
				lens[nxt] = lens[prv] + 1, \
				dct[nxt] = realloc(dct[nxt],lens[nxt] + 8), \
				dct[nxt] = memcpy(dct[nxt],dct[prv],lens[prv]), \
				dct[nxt++][lens[prv]] = *dct[lens[buf[k]] ? buf[k] : prv]; \
				lens[nxt] = 0; \
			prv = buf[k]; \
			D = (uint64_t *)(out+strIx), S = (uint64_t *)dct[buf[k]]; \
			for (int x = 0; x < lens[buf[k]]; x+=8) *D++=*S++; \
			strIx += lens[buf[k]]; \
		} \
		inB += b8; \
		i = j = k = bin = 0; \
	} \
	return strIx; \
}
#define DEKFNF_PREDEF(b) { \
	static const uint32_t N = (1 << b) - 1, \
		bel = (UINT16_MAX/(b*64))*(b*64), b8 = bel*b/64; \
	uint64_t *inB = in, *D, *S, strIx = 1; \
	static const uint16_t IB = 64-b; \
	uint16_t code = (*inB << IB) >> IB; \
	uint16_t i=1, k=1, j=0, shf, bin=b, prv = code, nxt = 256; \
	lens[nxt] = 0, *out = code; \
	while (code != N) { \
		for (; i < bel && code != N; ++i) { \
			register uint16_t slc = bin > IB; \
			shf = (slc ? IB + 64 : IB) - bin, j += slc; \
			register uint16_t ps = (inB[j] << shf) >> IB; \
			code = !slc | (bin == 64) ? ps : (inB[j-1] >> bin) | ps; \
			bin = slc ? bin - IB : bin + b; \
			buf[i] = code; \
		} \
		for (; k < i; ++k) { \
			if (nxt < N) { \
				lens[nxt] = lens[prv] + 1; \
				uint64_t *D = (uint64_t *)dct[nxt], *S = (uint64_t *)dct[prv]; \
				for (int x = 0; x < lens[prv]; ++x) *D++=*S++; \
				dct[nxt++][lens[prv]] = *dct[lens[buf[k]] ? buf[k] : prv]; \
				lens[nxt] = 0; \
			} \
			prv = buf[k]; \
			D = (uint64_t *)(out+strIx), S = (uint64_t *)dct[buf[k]]; \
			for (int x = 0; x < lens[buf[k]]; x+=8) *D++=*S++; \
			strIx += lens[buf[k]]; \
		} \
		inB += b8; \
		i = j = k = bin = 0; \
	} \
	return strIx; \
}
static inline size_t dekfn9(uint64_t *in, char *out) DEKFNF_PREDEF(9)
//static inline size_t dekfn10(uint64_t *in, char *out) DEKFNF_PREDEF(10)
static inline size_t dekfn11(uint64_t *in, char *out) DEKFNF_PREDEF(11)
//static inline size_t dekfn12(uint64_t *in, char *out) DEKFNF(12)
static inline size_t dekfn13(uint64_t *in, char *out) DEKFNF(13)
static inline size_t dekfn14(uint64_t *in, char *out) DEKFNF(14)
//static inline size_t dekfn15(uint64_t *in, char *out) DEKFNF(15)
static inline size_t dekfn16(uint64_t *in, char *out) DEKFNF(16)

static inline size_t enkfnX(uint8_t *in, size_t len, uint64_t *out, int b, int NP, int FP) { // make TCF
	uint32_t N = (1<<b), bel = UINT16_MAX, loc = 0;
	typedef struct {uint32_t Child[256];} TN;
	TN *Node = calloc(N,sizeof(*Node)); // dict init
	uint32_t *Parent = calloc(N,sizeof(*Parent));
	uint8_t *WhichCh = calloc(N,sizeof(*WhichCh));
	uint32_t *Count = calloc(N,sizeof(*Count)); // DELETE
	uint32_t nxt = 1; for (; nxt < 257; ++nxt) 
		Node[0].Child[nxt-1] = nxt;
	const uint32_t INCBIT = 256, EOFCAP = 257, FIRST = 259;
	nxt = FIRST;
	uint8_t *buffer = malloc(N); 
	uint32_t *outBuf = malloc((bel+32)*sizeof(*outBuf)); 
	uint32_t lastN = 0, temp, crit = 512, bits = 9, bin = 0; 
	uint64_t mask = 0, bix = 0, bix2 = 0, pass = 0; 
	for (size_t i = 0; i < len; ++i) {
		buffer[loc] = in[i]; // consider next byte
		if (!(temp=Node[lastN].Child[buffer[loc++]])) {
			if (nxt != N) {
				Parent[nxt]=lastN;
				WhichCh[nxt] = in[i];
				Node[lastN].Child[in[i]]=nxt;
				Count[nxt] = 1;
				register uint32_t p = nxt; 
				while (p = Parent[p]) 
					if (Count[p] != UINT32_MAX)
						++Count[p];
				p = nxt; // reset p to head of this chain
				pass += N == ++nxt;
				while (pass < NP) { // find next good nxt
					nxt = nxt == N ? FIRST : nxt; // if (nxt == N) {...}
					for (; nxt < N; ++nxt) if (Count[nxt] <= 1) break;
					if (nxt == N) ++pass;
					else if (nxt == p) printf("[%u] lapped at %u\n",pass,nxt), nxt=N, pass=NP;
					else {
						if (!FP) {
							register uint32_t x = nxt; //, cnt = Count[nxt];
							while (x = Parent[x]) --Count[x]; //!=Parent[x]) //for (x = Parent[nxt]; x !=Parent[x]; x = Parent[x])
						}
						Node[nxt] = (TN){0};
						Node[Parent[nxt]].Child[WhichCh[nxt]] = 0;
						break;
					}
				} 
				while (lastN - 1 >= crit) crit+=crit, 
					outBuf[bix++] = INCBIT;
			} 
			outBuf[bix++] = lastN - 1; // store last good
			if (bix >= bel) { 
				for (int x = 0; x < bix; ++x) {
					mask |= (uint64_t)outBuf[x] << bin;
					bin += bits;
					if (bin >= 64) {
						out[bix2++] = mask;
						mask = (uint64_t)outBuf[x] >> ((64 + bits) - bin);
						bin -= 64;
					}
					if (outBuf[x]==INCBIT) ++bits;
				}
				bix = 0;
			}
			*buffer = in[i], loc = 1; // reset heads
			lastN = in[i]+1;
		}
		else lastN = temp; // shortcut to next good
	}
	outBuf[bix++] = lastN - 1; 
	outBuf[bix++] = EOFCAP; 
	for (int x = 0; x < bix; ++x) {
		mask |= (uint64_t)outBuf[x] << bin;
		bin += bits;
		if (bin >= 64) {
			out[bix2++] = mask;
			mask = (uint64_t)outBuf[x] >> ((64 + bits) - bin);
			bin -= 64;
		}
		if (outBuf[x]==INCBIT) ++bits;
	}
	if (bin) out[bix2++] = mask; 
	free(Parent), free(Node), free(buffer), free(outBuf); 
	free(WhichCh), free(Count);
	return bix2*8;
}

// TCF Decompressor
static inline size_t dekfnX(uint64_t *inB, char *out, int b, int NP, int FP) {
	uint32_t dictSz = 1<<b, N = dictSz-1; 
	char **dct = calloc(dictSz,sizeof(*dct));
	uint32_t *lens = calloc(dictSz,sizeof(*lens)), pass=0;
	uint32_t *Parent = calloc(dictSz,sizeof(*Parent)); 
	uint32_t *Count = calloc(dictSz,sizeof(*Count)); 
	uint32_t nxt = 0; for ( ; nxt < 256; ++nxt) 
		dct[nxt] = malloc(8), Parent[nxt] = N,
		*dct[nxt] = nxt, lens[nxt] = 1;
	const uint32_t INCBIT = 256, EOFCAP = 257, FIRST = 258;
	nxt = FIRST;
	size_t strIx = 0, z = 0, b8 = 1024;
	uint32_t *in = malloc(b8*sizeof(*in));
	uint64_t prvChunk = *inB, *D, *S;
	uint32_t prv = N, bits = 9, bin = 0, code = 0;
	while (code != EOFCAP) {
		int i = 0; //if (bin > (64-bits)) --z;
		for (uint32_t shf=0; i < b8 && code != EOFCAP;) { 
			int inv_bits = 64 - bits;
			if (bin > inv_bits) {
				++z;
				code = bin == 64 ? 0 : prvChunk >> bin; // start word with tail of old
				shf = (128 - bits) - bin; 
				code |= (inB[z] << shf) >> inv_bits;
				bin -= inv_bits;
				prvChunk = inB[z];
			}
			else shf = inv_bits - bin,
				code = (inB[z] << shf) >> inv_bits,
				bin += bits;
			if (code == INCBIT) ++bits;
			else in[i++] = code;
		}
		for (int j = 0; j < i; ++j) {
			if (prv < N && nxt < N) {
				
				lens[nxt] = lens[prv] + 1;
				dct[nxt] = malloc(lens[nxt]+8);
				dct[nxt] = memcpy(dct[nxt],dct[prv],lens[prv]);
				//uint64_t *D = (uint64_t *)dct[nxt], *S=(uint64_t *)dct[prv];
				//for (int x = 0; x < lens[prv]; x+=8) *D++=*S++;
				dct[nxt][lens[prv]] = *dct[dct[in[j]] ? in[j] : prv];
				if (NP) {
					Parent[nxt] = prv;
					Count[nxt] = 1;
					register uint32_t p = nxt;
						while ((p = Parent[p])!=N) 
							if (Count[p] != UINT32_MAX)
							++Count[p];
					p = nxt; // reset p
					pass += N == ++nxt;
					if (pass) while (pass < NP) {
						nxt = nxt == N ? FIRST : nxt; 
						for (; nxt < N; ++nxt) if (Count[nxt] <= 1) break;
						if (nxt == N) ++pass;
						else if (nxt == p) printf("[%u] lapped at %u\n",pass,nxt), nxt=N, pass=NP;
						else {
							if (!FP) {
								register uint32_t x = nxt;
								while (x = Parent[x]) --Count[x];
							}
							free(dct[nxt]), dct[nxt]=0; 
							break;
						}
					}
				}
				else ++nxt;
			}
			prv = in[j];
			D = (uint64_t *)(out+strIx), S = (uint64_t *)dct[in[j]];
			for (int x = 0; x < lens[in[j]]; x+=8) *D++=*S++;
			strIx += lens[in[j]];
		}
	}
	for (size_t i = 0; i < dictSz; ++i) free(dct[i]);
	free(dct), free(lens), free(in), free(Parent), free(Count);
	return strIx; 
} 
#pragma GCC pop_options

int main( int argc, char *argv[] ) {
	char *mode, *ref_FN, *out_FN;
	if (argc == 3) mode = "c", ref_FN = argv[1], out_FN = argv[2];
	else if (argc == 4) mode = argv[1], ref_FN = argv[2], out_FN = argv[3];
	else {
		puts("\nNope. You need a mode (cX, pX, or d), input, and output. That order."); 
		puts(" Mode c is \"compress\" and does speculative large-block compression.");
		puts(" Mode p is \"pack\" and does conservative small-block compression.");
		puts(" The 'X' means the compression level, from 0 to 9:");
		puts("   0: Fastest, and for text only. No extended symbols.");
		puts("   1-3: Fast compress. Acts roughly the same for modes c and p.");
		puts("   4-6: Balanced compress. Big difference between c and p.");
		puts("   7-9: Slow compress. Slight difference between c and p.");
		puts(" The fastest modes for decompression are c0 and c2.");
		puts("\nExample: 'kafan c2 orig_file compressed' or 'kafan d compressed orig_file'");
		puts("Guide: try c2 or c3. Too big? c5. Got worse? p5. Need more? c7/p7.");
		puts("If none of this clicks with you, don't touch this program.");
		exit(1);
	}
	int shifts[] = {22,20,20,20,22,24,24,25,26,27};
	int shiftZ[] = {20,16,16,16,16,18,20,23,25,26};
	if (*mode=='c' || *mode=='p') {
		if (mode[1] < '0' || mode[1] > '9') {
			puts("Seriously? No. You need to pick a level from 0 to 9"); exit(1);}
		FILE *testComp = fopen(ref_FN,"rb");
		FILE *compressed = fopen(out_FN,"wb");
		if (!testComp || !compressed) {
			puts("You need actual, accessible paths. Yeah, don't quit your day job."); exit(1);}
		fseek(testComp,0,SEEK_END); size_t sz = ftell(testComp); rewind(testComp);
		printf("File size is %llu\n",sz);
		int sel = mode[1] - 48, bs = 1 << (*mode == 'c' ? shifts[sel] : shiftZ[sel]);
		if (sz < 10) puts("\nUm, you *DO* know such a tiny file won't compress? (No? Shocking.)\n");
		fwrite(&sel,1,1,compressed);
		printf("Compressing using block size %d\n",bs);
		
		clock_t start = clock();
		uint8_t *bufIn = malloc(bs+1);
		uint64_t *bufOut = malloc(2*bs + 2);
		size_t total = 0, totCmp = 1; 
		int32_t read, bad;
		#define ENCLP(COMPLINE) { \
			while (read=fread(bufIn,1,bs,testComp)) { \
				uint32_t comped = COMPLINE; \
				if (comped >= read) \
					comped = read, bad = -read, \
					fwrite(&bad,4,1,compressed), \
					fwrite(bufIn,1,read,compressed); \
				else  \
					fwrite(&comped,4,1,compressed), \
					fwrite(bufOut,1,comped,compressed); \
				total+=read; totCmp+=comped + 4; \
				printf("Block at %llu, size %llu, csize %llu (%f)\n", \
					total,read,comped,(double)comped/read); \
			} \
		}
		switch (mode[1]) {
			case '0': ENCLP(enmuleFb(bufIn,read,(char*)bufOut)); break;
			case '1': ENCLP(enkfn9(bufIn,read,bufOut)); break;
			case '2': ENCLP(enkfn11(bufIn,read,bufOut)); break;
			case '3': ENCLP(enkfn13(bufIn,read,bufOut)); break;
			case '4': ENCLP(enkfn14(bufIn,read,bufOut)); break;
			case '5': ENCLP(enkfn16(bufIn,read,bufOut)); break;
			case '6': ENCLP(enkfnX(bufIn,read,bufOut,18,0,0)); break;
			case '7': ENCLP(enkfnX(bufIn,read,bufOut,20,4,1)); break;
			case '8': ENCLP(enkfnX(bufIn,read,bufOut,21,16,1)); break;
			case '9': ENCLP(enkfnX(bufIn,read,bufOut,21,256,1)); break;
		}
		printf("Total size = %llu, compressed = %llu (%f, %f)\n",
			total,totCmp,(double)totCmp/total, (double)total/totCmp);
		free(bufIn), free(bufOut); //close file
		
		printf("Time: %.3f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC);
		if (totCmp >= total) puts("\nSo how do you feel having wasted your time for nothing?"),
			puts("\"Hey, I'm gonna go and compress near-randomness!\"\nYep, stunning alright.");
		fclose(compressed); 
	}
	
	else if (*mode == 'd') {
		FILE *compressed = fopen(ref_FN,"rb");
		FILE *reflate = fopen(out_FN,"wb");
		if (!compressed || !reflate) {puts("Error in filenames."); exit(1);}
		
		// prep the global structures for the fast decoders
		for (int i = 0; i < 256; ++i) 
			dct[i] = malloc(8), 
			*dct[i] = i, lens[i] = 1;
		for (uint32_t i = 256, b, c; i < (1 << 11); ++i) 
			b = i - 254, c = b >> 3, c = c << 3,
			dct[i] = malloc(c + (c < b ? 8 : 0));

		fseek(compressed,0,SEEK_END); 
		size_t sz = ftell(compressed); rewind(compressed);
		printf("Opening compressed file (%llu bytes)...\n",sz);
		if (!sz) {puts("Nice try genius. File is EMPTY."); exit(1);}
		int32_t read = 0, *readP = malloc(sizeof(*readP));
		fread(readP,1,1,compressed);
		if (*readP < 0 || *readP > 9) {
			puts("Umm, no. You need an *actual KFN file* to decompress."); exit(1);}
		int sel = *readP, bs = 1 << shifts[sel];
		printf("Using block size %d [%d]\n",bs,sel);

		clock_t start = clock();
		uint8_t *bufOut = malloc(bs + 16);
		uint64_t *bufIn = malloc(bs + 16);
		size_t total = 1, totInfl = 0;
		#define DECLP(DECLINE) { \
			while (read=fread(readP,4,1,compressed)) { \
				read = *readP; \
				/*printf("Block at %llu: %d\n",total,read);*/ \
				size_t infl; \
				if (read < 0) \
					read = -read, infl = read, \
					fread(bufIn,1,read,compressed), \
					fwrite(bufIn,1,read,reflate); \
				else  \
					fread(bufIn,1,read,compressed), \
					infl = DECLINE, \
					fwrite(bufOut,1,infl,reflate); \
				total += read + 4; \
				totInfl+=infl; \
			} \
		}
		
		switch (sel) {
			case 0: DECLP(demuleFb((uint8_t*)bufIn,read,bufOut)); break;
			case 1: DECLP(dekfn9(bufIn,bufOut)); break;
			case 2: DECLP(dekfn11(bufIn,bufOut)); break;
			case 3: DECLP(dekfn13(bufIn,bufOut)); break;
			case 4: DECLP(dekfn14(bufIn,bufOut)); break;
			case 5: DECLP(dekfn16(bufIn,bufOut)); break;
			case 6: DECLP(dekfnX(bufIn,bufOut,18,0,0)); break;
			case 7: DECLP(dekfnX(bufIn,bufOut,20,4,1)); break;
			case 8: DECLP(dekfnX(bufIn,bufOut,21,16,1)); break;
			case 9: DECLP(dekfnX(bufIn,bufOut,21,256,1)); break;
		}
		printf("\nTotal size = %llu, comped = %llu (%f, %f)\n",totInfl,total,(double)totInfl/total, (double)total/totInfl);
		printf("\ntime: %.3f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC);
		free(bufOut), free(bufIn);
	}

	else printf("Look, there is no mode '%s'. Read the simple instructions.",mode);
	return 0;
}