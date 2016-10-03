/* gencode.c	program to generate "gene" string and snippets*/
/*  1-Apr-2016	freda creation*/


#include "platform.h"
#include "types.h"
#include "submalloc.h"

typedef struct snippet {	/*snipet of genetic code*/
	int len;		/*length of snippet*/
	struct rope *prope;	/*pt to rope that snippet is in*/
	llong code[2];		/*snippet value*/
} SNIPPET, *PSNIPPET;

#ifdef SNIPLISTPT	/*if using pointers in SNIPLIST, doubles size*/
typedef struct sniplist {	/*list of snippets a unit is in and position*/
	int pos;		/*unit position in snippet or snippet in rope*/
	PSNIPPET psnippet;	/*pt to snippet*/
} SNIPLIST, *PSNIPLIST;
#define GPSNIP(psnipl)	((psnipl)->psnippet)
#define GISNIP(psnipl)	0
#define SPSNIP(psnipl, psnip, snipidx)	((psnipl)->psnippet = (psnip))
#else			/*otherwise using indices instead of pt in sniplist*/
typedef struct sniplist {	/*list of snippets a unit is in and position*/
	int pos;		/*unit position in snippet or snippet in rope*/
	int snipidx;		/*index of snippet in insniplist*/
} SNIPLIST, *PSNIPLIST;
#define GPSNIP(psnipl)	(insniplist[(psnipl)->snipidx])
#define GISNIP(psnipl)	((psnipl)->snipidx)
#define SPSNIP(psnipl, psnip, idx)	((psnipl)->snipidx = (idx))
#endif

#ifdef NOSUBMALLOC	/*if not doing any submallocs*/
#define MALSNIPPET(size)	(malloc(size))
#define FREESNIPPET(ptr, size)	(free(ptr))
#define MALSNIPLIST(size)	(malloc(size))
#define REALSNIPLIST(ptr, size1, size2)	(realloc(ptr, size2))
#define FREESNIPLIST(ptr, size)	(free(ptr))
#define MALROPE(size)		(malloc(size))
#define FREEROPE(ptr, size)	(free(ptr))
#else
#define MXSUBMAL 64
#define MALSNIPPET(size) ((size) <= MXSUBMAL ? submalloc(size) : malloc(size))
#define FREESNIPPET(ptr, size)	((size) <= MXSUBMAL ? subfree((Uint8*)(ptr), size) : free(ptr))
#define MALSNIPLIST(size)	((size) <= MXSUBMAL ? submalloc(size) : malloc(size))
#define REALSNIPLIST(ptr, size1, size2)	(subrealloc((Uint8*)(ptr), size1, size2))
#define FREESNIPLIST(ptr, size)	((size) <= MXSUBMAL ? subfree((Uint8*)(ptr), size) : free(ptr))
#define MALROPE(size)		(submalloc(size))
#define FREEROPE(ptr, size)	(subfree((Uint8*)(ptr), size))
#endif

typedef struct unit {		/*unit length string of codes*/
	llong code;		/*negative or code for this unit*/
	int mxsnip;		/*max number of snippets in current snippet list*/
	int nsnip;		/*number of snippets in current snippet list*/
	PSNIPLIST psniplist;	/*pt to current snippet list*/
} UNIT, *PUNIT;

typedef struct rope {		/*rope made out of many snippets*/
	int mxsnip;		/*max number of snippets in current snippet list*/
	int nsnip;		/*number of snippets in current snippet list*/
	PSNIPLIST psniplist;	/*pt to current snippet list*/
} ROPE, *PROPE;

typedef struct nxtlist {	/*link to possible next overlapping rope*/
	short next;		/*index of next rope*/
	short olap;		/*amount of overap*/
} NXTLIST, *PNXTLIST;
#define MAXTAILS	3	/*max number of possible next ropes to record*/
typedef struct ropelist {	/*list of ropes for endjoinropes*/
	PROPE prope;		/*this rope*/
	PSNIPPET psnips;	/*pt to start snippet*/
	PSNIPPET psnipe;	/*pt to ending snippet*/
	int len;		/*length of rope*/
	short fwdmir;		/*index of forward mirror rope*/
	short rev;		/*index of reverse rope*/
	short revmir;		/*index of reverse mirror rope*/
	short ntails;		/*number of possible next ropes*/
	llong codes;		/*start code*/
	llong codee;		/*end code*/
	NXTLIST nxtlist[MAXTAILS];	/*list of possible next ropes*/
} ROPELIST, *PROPELIST;


int randseed=0x12345678;	/*random seed for generating random genetic code*/
int ncopy=4;		/*number of copies of full code*/
int nsnippet;		/*number of input snippets*/
int nsnipinrope=0;	/*number of snippets in ropes*/
PSNIPPET *insniplist;	/*input snippet list*/
int unitlen=10;		/*number of letters in a reference unit*/
PUNIT unitdb;		/*unit data base*/
int nunitdb;		/*number of unit slots in hash coded unit table*/
PSNIPLIST catsniplist;	/*concatenated snippet lists for unit database*/
int ncatsnip;		/*number of catsniplist entries*/
int nropes=0;		/*number of ropes of code*/
int maxnropes=0;	/*max number of ropes during this solution*/
//Uint8 *minmalloc;	/*first malloc pt*/
//Uint8 *maxmalloc;	/*largest malloc*/
int debugpr=0;		/*if debug printing*/
char *codechar="ACGT";	/*genetic code letters*/
char *bigspace=0;	/*maybe extra space for contig allocation*/
llong bigspacesz=0;	/*size of large contiguous space*/
Uint8 *lowmalst=(Uint8*)0xFFFFFFFFFFFFFFF;	/*low malloc start and end*/
Uint8 *lowmalend=0;
Uint8 *himalst=(Uint8*)0xFFFFFFFFFFFFFFF;	/*high malloc start and end*/
Uint8 *himalend=0;
#define BIGSUBMAL 16384	/*size of large malloc for submalloc allocation*/
#define NSMSIZES 7	/*number of special submalloc sizes*/
int submalsizes[NSMSIZES] = {	/*special sub malloc sizes*/
	16,24,32,40,48,56,64};
int prsubmal=0;	/*if print submalloc statistics*/


void elimsubsets(
/*eliminate snippets that are subsets of other snippets*/
	int *pnsnippet,		/*pt to number of input snippets*/
	PSNIPPET *insniplist,	/*input snippets*/
	int unitlen,		/*length of each unit*/
	PUNIT unitdb);		/*where unit database*/
void makeropes(
/*make ropes out of snippets*/
	int nsnippet,		/*number of input snippets*/
	PSNIPPET *insniplist,	/*input snippets*/
	int unitlen,		/*length of each unit*/
	PUNIT unitdb,		/*where unit database*/
	int minolap);		/*minimum overlap of snippets*/
int chkmatch(
/*check that two snippets match in overlap regions*/
/*return 0 if OK, 1 if mismatch*/
	PSNIPPET psnip1,	/*a snippet*/
	PSNIPPET psnip2,	/*another snippet*/
	int off2);		/*code offset from first to second*/
int add2rope(
/*add snippet to rope*/
/*return 0 for success, 1 for failure*/
	PSNIPPET insnip,	/*snippet already in rope*/
	int inpos,		/*relative position of snippet in rope*/
	PSNIPPET newsnip,	/*snippet to be added to rope*/
	int newpos,		/*relative position of new snippet*/
	int newsnipidx);	/*index of new snippet in insniplist*/
int joinropes(
/*join two overlapping ropes via snippets in each*/
/*return 0 for success, 1 for failure*/
	PSNIPPET psnip1,	/*first snippet*/
	int s1,			/*relative offset of first snippet*/
	PSNIPPET psnip2,	/*second snippet*/
	int s2);		/*relative offset of scond snippet*/
void endjoinropes(
/*join remaining ropes if possible*/
/*this part breaks the minolap rule*/
	int nsnippet,		/*number of input snippets*/
	PSNIPPET *insniplist,	/*input snippets*/
	int ncopy,		/*number of copies of code*/
	int minolap);		/*minimum codes for overlapping snippets*/
void propes(
/*print out ropes*/
	int nsnippet,		/*number of snippets*/
	PSNIPPET *insniplist,	/*input snippet list*/
	int nropes);		/*number of ropes*/
void dounits(
/*create database of code units of length unitlen in the snippets*/
	int nsnippet,		/*number of input snippets*/
	PSNIPPET *insniplist,	/*input snippets*/
	int unitlen,		/*length of each unit*/
	int both,		/*if track both ends of snippets*/
	PUNIT *punitdb);	/*where to return pt to unit database*/
llong extcode(
/*extract section of code*/
	PSNIPPET psnippet,	/*snippet from which to extract*/
	int i,			/*start position*/
	int len);		/*length of section*/
void readsnips(
/*read snippets from file*/
	char *filebase,		/*file base name*/
	int minsniplen,		/*minimum snippet length*/
	int *pnsnippet,		/*where return number of snippets*/
	PSNIPPET **ppsniplist);	/*where to return snippet data*/
static void trackmal(
/*track range of mallocs*/
	void *start,	/*start of malloced range*/
	llong len);	/*length of malloced range in bytes*/
static int getrand(
/*get random value in a specified range*/
	int srange,	/*start of range*/
	int erange);	/*end of range*/
int readfile(
/*open file, read into memory and break into lines*/
/*return 0 or number of lines*/
	char *name,	/*name of file*/
	Uint8 **lines,	/*where to put pointers to lines*/
	Uint8 *bigbuf,	/*where to put file contents*/
	llong maxsize);	/*max size of file*/




main(int ac, char **av)
{
	int mode;	/*split snippets or assemble snippets*/
	char *filebase;	/*base name of genome code files*/
	int i,j,k,n;
	char *ap;	/*point into arg strings*/
	int minsniplen=20;	/*minimum snippet length*/
	int maxsniplen=100;	/*maximim snippet length*/
	int codelen=10000;	/*total length of code*/
	int minolap=10;	/*minimum overlap for allowing match*/
	int noendjoin=0;/*if skip end join*/
	int finddupe=0;	/*if find longest duplicate sequence in code*/
	int size;	/*number of llongs in code*/
	PSNIPPET fullsnip;	/*full code snippet*/
	int usize;	/*size of a unit of codes*/
	int *dutab;	/*duplicate unit table*/
	int idx;	/*index into dutab or occtab*/
	int *occtab;	/*concatenated occurence tables for all units*/
	int maxdup;	/*max duplicate length +1 */
	llong code1;	/*first copy of duplicate code*/
	llong code2;	/*second copy of duplicate code*/
	int dupe1;	/*position of first duplicate*/
	int dupe2;	/*position of second duplicate*/
	char filename[100];	/*name of code and snippet file*/
	llong *codep;	/*pointer to code when generating code*/
	llong code;	/*llong of code*/
	llong lval;	/*byte's worth of code*/
	llong lmask;	/*mask of trailing bits of code*/
	int fd;		/*file descriptor*/
	llong *cp;	/*pt into genetic code*/
	int nib;	/*which nibble*/
	char *obuf;	/*output buffer*/
	char *op;	/*point into output buffer*/
	int dig;	/*digit*/
	char hdig;	/*hex digit*/
	int *snipdef;	/*list of snippet definitions, start, length pairs*/
	int *snipmap;	/*map of which snippet definitions used*/
	int nc;		/*copy number for making snippets*/
	int sniplen;	/*snippet length*/
	int nsnip;	/*number of snippets in this copy*/
	int li;		/*shift index for extracting gencode letter index*/
	int rev;	/*if reverse snippet*/
	int ver;	/*version of code, forward, reverse, & mirror helix*/
	llong ln;	/*long count*/
/*read arguments*/
	if(ac<3)
	{
		printf("%s mode filebase\n",av[0]);
		printf("mode 0=gen code & snippets, mode 1=assemble snippets\n");
		printf("-l CodeLen -i minSnipLen -x maxSniplen -c ncopy -r rand#\n");
		printf("-u unitlen -o minOlap -d debug# -eNoEndJoin\n");
		printf("-? to find dupe -p <#Mbytes to prealloc> -submal_stats\n");
		exit(1);
	}
	mode = atoi(av[1]);
	filebase = av[2];
	for(i=3;i<ac;++i)
	{
		ap = av[i];
		if(*ap++ != '-')
		{
			printf("expecting '-' in argument '%s'\n",ap-1);
			exit(1);
		}
		switch(*ap++)
		{
		case 'c':	/*number of copies*/
			if(!*ap)
				ap = av[++i];
			ncopy = atoi(ap);
			break;
		case 'd':	/*debug dump*/
			if(!*ap)
				ap = av[++i];
			debugpr = atoi(ap);
			break;
		case 'e':	/*no end join*/
			noendjoin = 1;
			break;
		case 'i':	/*minimum snippet length*/
			if(!*ap)
				ap = av[++i];
			minsniplen = atoi(ap);
			break;
		case 'l':	/*total code length*/
			if(!*ap)
				ap = av[++i];
			codelen = atoi(ap);
			break;
		case 'o':	/*minimum overlap length*/
			if(!*ap)
				ap = av[++i];
			minolap = atoi(ap);
			break;
		case 'p':	/*numper of megabytes to preallocate*/
			if(!*ap)
				ap = av[++i];
			bigspacesz = atoi(ap)*1024*1024;
			break;
		case 'r':	/*starting random number*/
			if(!*ap)
				ap = av[++i];
			randseed = atoi(ap);
			break;
		case 's':	/*print submalloc statistics*/
			prsubmal = 1;
			break;
		case 'u':	/*unit length*/
			if(!*ap)
				ap = av[++i];
			unitlen = atoi(ap);
			break;
		case 'x':	/*max snippet length*/
			if(!*ap)
				ap = av[++i];
			maxsniplen = atoi(ap);
			break;
		case '?':	/*find longest duplicate in code*/
			finddupe = 1;
			break;
		default:
			printf("unknown argument, '%s'\n", ap-2);
			exit(1);
			break;
		}
	}
/*check for code and snippet generation*/
	if(mode==0)
	{
	/*randomly generate code and code snippets*/
	/*allocate space for code*/
		size = (codelen+31) / 32;
		fullsnip = (PSNIPPET)malloc(n = sizeof(SNIPPET)+size*sizeof(fullsnip->code[0]));
		fullsnip->len =codelen;
		fullsnip->prope = 0;
		codep = fullsnip->code;
	/*use random number generator to fill code*/
		for(i=0;i<size;++i)
		{
			for(j=56,code=0;j>=0;j-=8)
			{
				lval = (llong) getrand(0,255);
				code |= lval << j;
			}
			*(codep+i) = code;
		}
	/*clear unused tail of code*/
		j = size*32 - codelen;
		lmask = ((llong)1 << (2*j)) - 1;
		*(codep + i -1) &= ~lmask;
	/*optionally find longest duplicate section in code*/
		if(finddupe)
		{
		/*find longest duplicate section in code*/
		/*this is a more complicated, but faster way to search for the max duplicate*/
		/*make a unit table with a count of the number of times each unit appears*/
			if(codelen >= 100000)
				unitlen = 10;
			if(codelen >= 30000000)
				unitlen = 14;
			usize = (llong)1 << (unitlen*2);
			dutab = (int*)malloc(n = usize*sizeof(int));
			memset(dutab, 0, n);;
		/*count the number of times each unit appears in the code*/
			for(i=0; i<codelen-unitlen+1; ++i)
			{
				code1 = extcode(fullsnip, i, unitlen);
				idx = (int)code1;
				dutab[idx]++;
			}
		/*convert that count list into offsets that leave room for */
		/* an int for each occurrence plus and end marker*/
			for(i=0,idx=0;i<usize;++i)
			{
				j = dutab[i];
				dutab[i] = idx;
				idx += j + 1;
			}
		/*allocate space for concatenated occurence list for each unit*/
			occtab = (int*)malloc(ln = idx*sizeof(int));
			memset(occtab, 0xFF, ln);
		/*fill in the occurence lists for each unit*/
			for(i=0; i<codelen-unitlen+1; ++i)
			{
				code1 = extcode(fullsnip, i, unitlen);
				idx = (int)code1;
				idx = dutab[idx];
				for( ;occtab[idx] > 0; ++idx)
					;
				occtab[idx] = i;
			}
		/*using this two level database of units, look for longest duplicate*/
			for(i=0,maxdup=unitlen,dupe1=dupe2=0; i<codelen-maxdup; ++i)
			{
				code1 = extcode(fullsnip, i, unitlen);
				idx = (int)code1;
				idx = dutab[idx];
				code1 = extcode(fullsnip, i, maxdup);
				for( ;occtab[idx] >= 0; ++idx)
				{
					j = occtab[idx];
					if(j <= i)
						continue;
					code2 = extcode(fullsnip, j, maxdup);
					while(code1 == code2)
					{
						++maxdup;
						code1 = extcode(fullsnip, i, maxdup);
						code2 = extcode(fullsnip, j, maxdup);
						dupe1 = i;
						dupe2 = j;
					}
				}
			}
			printf("//longest duplicate is %d units at %d and %d\n",
			maxdup-1, dupe1, dupe2);
		}
	/*write full code*/
	/*form name*/
		strcpy(filename, filebase);
		strcat(filename, ".code");
		fd = creat(filename, 0666);
		if(fd < 0)
		{
			printf("can't create %s\n", filename);
			exit(1);
		}
	/*allocate output line buffer*/
		n = maxsniplen+8 < 200 ? 200 : maxsniplen+8;
		obuf = (char*)malloc(n);
	/*output 4 versions of the genetic code*/
		for(ver=0;ver<4;++ver)
		{
		/*header line*/
			strcpy(obuf,"// ");
			strcat(obuf, ver<2 ? "forward " : "reverse ");
			if(ver&1)
				strcat(obuf, "mirror ");
			strcat(obuf, "code\n");
			write(fd, obuf, strlen(obuf));
		/*loop, writing code, up to 64 letters per line*/
		/* 32 codes / llong,  2 llong / line */
		/* 1 text line = 64 letters = 2 llong */
			for(i=0;i<codelen;i+=64)
			{
				for(j=0, op=obuf;(j<64) && (i+j < codelen); ++j)
				{
					k = ver<2 ? i+j : codelen-i-j-1;
					dig = extcode(fullsnip, k, 1);
					if(ver&1)
						dig = 3 - dig;
					hdig = codechar[dig];
					*op++ = hdig;
				}
				*op++ = '\n';
				if(write(fd, obuf, op - obuf) != op-obuf)
				{
					printf("can't write to %s\n", filename);
					exit(1);
				}
			}
		}
		close(fd);
	/*write snippet file*/
	/*form name*/
		strcpy(filename, filebase);
		strcat(filename, ".snip");
		fd = creat(filename, 0666);
		if(fd < 0)
		{
			printf("can't create %s\n", filename);
			exit(1);
		}
	/*report parameters*/
		sprintf(obuf, "// %d copies of %d letter code in %d to %d length snippets\n",
		ncopy, codelen, minsniplen, maxsniplen);
		write(fd, obuf, strlen(obuf));
	/*allocate space for snippet definitions*/
		snipdef = (int*)malloc((i=codelen/minsniplen+1)*2*sizeof(int));
		snipmap = (int*)malloc((i+31)/32 * sizeof(int));
	/*loop over number of times to split code*/
		for(nc=0;nc<ncopy;++nc)
		{
		/*calculate division of snippets*/
			for(i=0,j=0;i<codelen;i+=sniplen,++j)
			{
			/*calculate snippet length*/
				sniplen = getrand(minsniplen, maxsniplen);
				sniplen = sniplen < codelen - i ? sniplen : codelen - i;
			/*save snippet definition*/
				snipdef[2*j] = i;	//start
				snipdef[2*j+1] = sniplen;
			}
		/*clear map of which snippets written*/
			nsnip = j;
			memset(snipmap, 0, (nsnip+7)/8);
		/*loop writing snippets, one per line*/
			for(n=0;n<nsnip;++n)
			{
			/*pick snippet that has not yet been written*/
			/*first try N random picks*/
				for(k=0;k<10;++k)
				{
					j = getrand(0,nsnip-1);
					if(~snipmap[j>>5] & (1<<(j&31)))
						break;
				}
			/*if have not yet found one, do linear search*/
				if(k>=10)
				{
					for(k=0;1;++k)
					{
						if(snipmap[k] != -1)
							break;
					}
					for(j=k*32;1;++j)
					{
						if(~snipmap[j>>5] & (1<<(j&31)))
							break;
					}
				}
			/*mark this snippet used*/
				snipmap[j>>5] |= (1<<(j&31));
			/*get snippet parameters*/
				i = snipdef[2*j];
				sniplen = snipdef[2*j+1];
			/*print snippet length*/
				sprintf(obuf,"%4d ", sniplen);
			/*maybe reverse snippet, maybe mirror version*/
				rev = getrand(0,3);
			/*loop over letters in snippet*/
				for(j=0,op=obuf+5;j<sniplen;++j)
				{
					dig = extcode(fullsnip, rev>=2 ? i+sniplen-1-j : i+j, 1);
					dig = rev&1 ? 3 - dig : dig;
					hdig = codechar[dig];
					*op++ = hdig;
				}
				*op++ = '\n';
				if(write(fd,obuf,op-obuf) != op-obuf)
				{
					printf("can't write to %s\n", filename);
					exit(1);
				}
			/*maybe tell where code is from*/
				if(debugpr & 2)
				{
					sprintf(obuf, "// snip %d to %d, %s%s\n",
					i, i+sniplen-1, rev&2 ? "rev" : "fwd",
					rev&1 ? "mir" : "");
					write(fd,obuf,strlen(obuf));
				}
			}
		}
		close(fd);
		exit(0);
	}
/*otherwise we are trying to reassemble a code out of snippets*/
/*check reasonable parameters*/
	if(unitlen > minolap)
	{
		printf("can't have unit length, %d, greater than minolap, %d\n",
		unitlen, minolap);
		exit(1);
	}
	if(unitlen < minolap)
		unitlen = minolap;
/*dump parameters*/
	printf("//");
	for(i=0;i<ac;++i)
		printf("%s ", av[i]);
	printf("\n");
/*initialize sub malloc allocation*/
	submalinit(BIGSUBMAL, NSMSIZES, submalsizes);
/*read a list of snippets from input file*/
	readsnips(filebase, minsniplen, &nsnippet, &insniplist);
/*create database of code units of length unitlen in the snippets*/
	dounits(nsnippet, insniplist, unitlen, 0, &unitdb);
/*eliminate complete subset snippets*/
	elimsubsets(&nsnippet, insniplist, unitlen, unitdb);
/*create database of code units of length unitlen in the snippets*/
	dounits(nsnippet, insniplist, unitlen, 1, &unitdb);
/*make ropes out of snippets in multiple passes*/
	nropes = 0;
	for(n=32;n>=minolap;--n)
		makeropes(nsnippet, insniplist, unitlen, unitdb, n);
/*join remaining ropes if possible*/
	if(!noendjoin)
		endjoinropes(nsnippet, insniplist, ncopy, minolap);
/*print out ropes*/
	propes(nsnippet, insniplist, nropes);
/*done*/
	exit(0);
}




PUNIT findunit(
/*find unit in hash table*/
/*return 0 or pt to unit*/
	llong code,	/*code in unit*/
	int doadd)	/*if add entry for missing unit*/
{
	int idx;	/*index into hash table*/
	llong ucode;	/*code in unit table*/
	PUNIT punit;	/*pt to unit in table*/
/*get hash index*/
	idx = code % nunitdb;
/*loop till find entry with proper code or empty entry*/
	while((ucode = (punit = unitdb + idx)->code) >= 0)
	{
		if(ucode == code)
			break;
		if(++idx >= nunitdb)
			idx = 0;
	}
/*maybe return found entry*/
	if(ucode == code)
		return(punit);
/*maybe and entry and return it*/
	if(doadd)
	{
		punit->code = code;
		return(punit);
	}
/*otherwise return 0*/
	return(0);
}




void elimsubsets(
/*eliminate snippets that are subsets of other snippets*/
	int *pnsnippet,		/*pt to number of input snippets*/
	PSNIPPET *insniplist,	/*input snippets*/
	int unitlen,		/*length of each unit*/
	PUNIT unitdb)		/*where unit database*/
{
	int nsnippet;		/*number of snippets*/
	int nsnip;		/*number of snippets processed*/
	PSNIPPET psnippet;	/*current snippet*/
	int len;		/*length of current snippet*/
	int f;			/*if working from first snippet*/
	llong code;		/*code from start or end of snippet*/
	//int unit;		/*Unit ID*/
	PUNIT punit;		/*pt to unit of interest*/
	int s1,e1;		/*range of snippet 1 relative to unit*/
	PSNIPLIST psnipl2;	/*snippet list for this unit*/
	int n2;			/*number of snippets in unit list*/
	PSNIPPET psnip2;	/*one of snippets in this unit list*/
	int s2,e2;		/*range of snippet 2 relative to unit*/
	int s,e;		/*range of overlap of snippets relative to unit*/
	int olap;		/*total potential overlap of 2 snippets*/
	int i, n;
	llong code2;		/*code from second snippet*/
	ROPE rope;		/*non null rope*/
	int nunits;		/*number of units in unit database*/
/*get number of snippets at start*/
	nsnippet = *pnsnippet;
/*loop over each input snippet*/
	for(nsnip=0;nsnip<nsnippet;++nsnip)
	{
	/*point to snippet data*/
		psnippet = insniplist[nsnip];
		len = psnippet->len;
	/*work on first and last unit of each snippet*/
		for(f=2;(--f>0) && !psnippet->prope;)
		{
		/*get first or last unit*/
			code = extcode(psnippet, f ? 0 : len-unitlen, unitlen);
			//unit = (int)code;
		/*pt into unit db for this unit*/
			//punit = unitdb + unit;
			punit = findunit(code, 0);
		/*figure range of possible overlap for first snippet*/
			s1 = f ? 0 : unitlen-len;
			e1 = f ? len : unitlen;
		/*loop over other snippets that this unit is in*/
			for(psnipl2=punit->psniplist, n2=0;
			(n2 < punit->nsnip) && !psnippet->prope;
			++psnipl2, ++n2)
			{
			/*get pointer to another snippet*/
				psnip2 = GPSNIP(psnipl2);
			/*skip if same snippet*/
				if(psnip2 == psnippet)
					continue;
			/*see if snippets match more than just this unit*/
			/*figure range of possible overlap*/
				s2 = -psnipl2->pos;
				e2 = psnip2->len + s2;
				s = s1 < s2 ? s2 : s1;
				e = e1 < e2 ? e1 : e2;
				olap = e - s;
			/*skip if not a full overlap of psnippet or psnip2*/
				if((olap < len) && (olap < psnip2->len))
					continue;
			/*see if codes match around this unit*/
				for(i=0;i<olap;i+=32)
				{
					code = extcode(psnippet, s-s1+i,
						i+32 <= olap ? 32 : olap-i);
					code2 = extcode(psnip2, s-s2+i,
						i+32 <= olap ? 32 : olap-i);
					if(code != code2)
						break;
				}
				if(i < olap)
					continue;
			/*check if one snippet is entirely inside another*/
				if(olap >= psnip2->len)
				{
					psnip2->prope = &rope;
					continue;
				}
				if(olap >= len)
				{
					psnippet->prope = &rope;
					continue;
				}
			}
		}
	}
/*eliminate all the snippets that are subsets of longer, completely overlapping ones*/
/*also eliminate too short snippets*/
	for(nsnip = n2 = 0; nsnip < nsnippet; ++nsnip)
	{
//		if(!(nsnip & 3))
//			printf("%d: %d %lX %lX %lX %lX\n", nsnip, insniplist[nsnip]->len,
//			insniplist[nsnip], insniplist[nsnip+1],
//			insniplist[nsnip+2], insniplist[nsnip+3]);
		psnippet = insniplist[nsnip];
		if(psnippet->prope || (psnippet->len < unitlen))
		{
			//free(psnippet);
			n = ((psnippet->len + 31) >> 5) - 2;	//extra space for code
			n = n < 0 ? 0: n;
			n = sizeof(SNIPPET) + n * 8;
			FREESNIPPET(psnippet, n);
			continue;
		}
	/*copy this good one*/
		insniplist[n2++] = psnippet;
	}
/*reallocate insniplist*/
	insniplist = (PSNIPPET*) realloc( insniplist, n2 * sizeof(psnippet));
/*free up unit database*/
	free(unitdb);
	free(catsniplist);
/*report number of snippets eliminated*/
	printf("// eliminated %d subset snippets out of %d original, leaving %d snippets\n",
	nsnippet-n2, nsnippet, n2);
/*return new number of snippets*/
	*pnsnippet = n2;
}




void makeropes(
/*make ropes out of snippets*/
	int nsnippet,		/*number of input snippets*/
	PSNIPPET *insniplist,	/*input snippets*/
	int unitlen,		/*length of each unit*/
	PUNIT unitdb,		/*where unit database*/
	int minolap)		/*minimum overlap of snippets*/
{
	int nsnip;		/*number of snippets processed*/
	PSNIPPET psnippet;	/*current snippet*/
	int snippetidx;		/*index of current snippet*/
	PSNIPPET psnip2;	/*one of snippets in this unit list*/
	int n,i;
	int nlropes;		/*number of local ropes in hash table*/
	PROPE *propes;		/*hash table of ropes*/
	PROPE prope;		/*current rope*/
	int idx;		/*index into hash table*/
	int nsawropes;		/*number of ropes I've seen in this pass*/
	PSNIPLIST psnipl;	/*snippet list for this rope*/
	int posmax;		/*max position in rope*/
	PSNIPPET psnipmax;	/*longest snippet that goes to end of rope*/
	int snipmaxidx;		/*index of psnipmax*/
/*if no snippets in ropes, link all snippets forward and back*/
	if(!nsnipinrope)
	{
	/*loop over each input snippet*/
		for(nsnip=0;nsnip<nsnippet;++nsnip)
		{
		/*point to snippet data*/
			psnippet = insniplist[nsnip];
		/*try to do best forward and best backward link*/
			extendsnip( psnippet, nsnip, 1, minolap, &psnip2);
			extendsnip( psnippet, nsnip, -1, minolap, &psnip2);
		}
	}
	else
	{
	/*only look for matches at start and end of each rope*/
	/*make a list of ropes seen*/
	/*create a hash table of ropes seen*/
		nlropes = 3*nropes/2;
		propes = (PROPE*)malloc(n = nlropes * sizeof(PROPE));
		memset(propes, 0, n);
	/*loop over each input snippet*/
		for(nsnip=0,nsawropes=0;nsnip<nsnippet;++nsnip)
		{
		/*point to snippet data*/
			psnippet = insniplist[nsnip];
		/*if not in rope, try to add it now*/
			if(!psnippet->prope)
			{
			/*try to do best forward and best backward link*/
				extendsnip( psnippet, nsnip, 1, minolap, &psnip2);
				extendsnip( psnippet, nsnip, -1, minolap, &psnip2);
				continue;
			}
		/*find or insert rope in local list*/
			prope = psnippet->prope;
			idx = ((llong)prope) % nlropes;
			for(;propes[idx] && (propes[idx] != prope); )
			{
				if(++idx >= nlropes)
					idx = 0;
			}
		/*if found in hash table, skip this snippet*/
			if(propes[idx])
				continue;
		/*put entry in hash table*/
			propes[idx] = prope;
			if(++nsawropes >= nropes)
				break;
		/*try to work back from start of rope*/
			psnippet = GPSNIP(prope->psniplist);
			snippetidx = GISNIP(prope->psniplist);
			extendsnip( psnippet, snippetidx, -1, minolap, &psnip2);
		/*try to work forward from end of rope*/
			prope = psnippet->prope;
			n = prope->nsnip - 1;
			psnipl = prope->psniplist + n;
			for(posmax = -100,psnipmax = 0,i=ncopy*6; (--i>=0) && (n>= 0);
			--n, --psnipl)
			{
				if(psnipl->pos + GPSNIP(psnipl)->len >= posmax)
				{
					posmax = psnipl->pos + GPSNIP(psnipl)->len;
					psnipmax = GPSNIP(psnipl);
					snipmaxidx = GISNIP(psnipl);
				}
			}
			extendsnip( psnipmax, snipmaxidx, 1, minolap, &psnip2);
		}
	/*free my list of ropes*/
		free(propes);
	}
}




int extendsnip(
/*routine to find next best extension of snippet*/
/*return 0 if OK, 1 if can't find next*/
	PSNIPPET psnippet,	/*base snippet*/
	int snippetidx,		/*index of base snippet*/
	int dir,		/*desired direction of extension*/
	int minolap,		/*minimum overlap*/
	PSNIPPET *ppnewsnip)	/*where to return new snippet*/
{
	int len;		/*length of current snippet*/
	int f;			/*if working from first snippet*/
	llong code;		/*code from start or end of snippet*/
	int unit;		/*Unit ID*/
	PUNIT punit;		/*pt to unit of interest*/
	int s1,e1;		/*range of snippet 1 relative to unit*/
	PSNIPLIST psnipl2;	/*snippet list for this unit*/
	int n2;			/*number of snippets in unit list*/
	PSNIPPET psnip2;	/*one of snippets in this unit list*/
	int snip2idx;		/*snippet 2 index*/
	int s2,e2;		/*range of snippet 2 relative to unit*/
	int s,e;		/*range of overlap of snippets relative to unit*/
	int olap;		/*total potential overlap of 2 snippets*/
	int i, n;
	llong code2;		/*code from second snippet*/
	int bestolap;		/*best new overlap for this snippet*/
	PSNIPPET bestpsnip2;	/*best other overlapped snippet*/
	int bestsnip2idx;	/*best index of snippet 2*/
	int bests1;		/*best start of snippet 1*/
	int bests2;		/*best start of snippet 2*/
	PROPE prope;		/*rope being created*/
	int mxsnip;		/*initial max snippets in rope*/
	PSNIPLIST propesnipl;	/*snippet list associated with rope*/
/*get snippet length*/
	len = psnippet->len;
/*working from first or last unit in snippet*/
	f = dir > 0 ? 0 : 1;
/*get first or last unit*/
	code = extcode(psnippet, f ? 0 : len-unitlen, unitlen);
	unit = (int)code;
/*pt into unit db for this unit*/
	punit = findunit(code, 0);
/*figure range of possible overlap for first snippet*/
	s1 = f ? 0 : unitlen-len;
	e1 = f ? len : unitlen;
/*loop over other snippets that this unit is in*/
	for(psnipl2=punit->psniplist, n2=0,bestolap=0; n2 < punit->nsnip;
	++psnipl2, ++n2)
	{
	/*get pointer to another snippet*/
		psnip2 = GPSNIP(psnipl2);
		snip2idx = GISNIP(psnipl2);
	/*skip if same snippet*/
		if(psnip2 == psnippet)
			continue;
	/*skip if now in same rope*/
		if(psnippet->prope && (psnippet->prope == psnip2->prope))
			continue;
	/*see if snippets match more than just this unit*/
	/*figure range of possible overlap*/
		s2 = -psnipl2->pos;
		e2 = psnip2->len + s2;
		s = s1 < s2 ? s2 : s1;
		e = e1 < e2 ? e1 : e2;
		olap = e - s;
	/*skip if less than minimum overlap*/
		if(olap < minolap)
			continue;
	/*see if codes match around this unit*/
		for(i=0;i<olap;i+=32)
		{
			code = extcode(psnippet, s-s1+i,
				i+32 <= olap ? 32 : olap-i);
			code2 = extcode(psnip2, s-s2+i,
				i+32 <= olap ? 32 : olap-i);
			if(code != code2)
				break;
		}
		if(i < olap)
			continue;
	/*have matching sets of code*/
	/*save best match for snippet so far*/
		if(olap > bestolap)
		{
			bestolap = olap;
			bestpsnip2 = psnip2;
			bestsnip2idx = snip2idx;
			bests1 = s1;
			bests2 = s2;
		}
	}
/*if I have a best match, process it now*/
	if(bestolap < minolap)
		return(1);
	olap = bestolap;
	psnip2 = bestpsnip2;
	snip2idx = bestsnip2idx;
	s1 = bests1;
	s2 = bests2;
/*form, add to or join ropes*/
/*if neither is in a rope*/
	if(!psnippet->prope && !psnip2->prope)
	{
	/*start a rope with these two*/
		//prope = (PROPE)malloc(n = sizeof(ROPE));
		prope = (PROPE)MALROPE(n = sizeof(ROPE));
		memset(prope, 0, n);
		mxsnip = 2;
		//propesnipl = (PSNIPLIST)malloc(n = mxsnip*sizeof(SNIPLIST));
		n = mxsnip*sizeof(SNIPLIST);
		propesnipl = (PSNIPLIST)MALSNIPLIST(n);
		trackmal((void*)propesnipl, n);
		++nropes;
		if(nropes > maxnropes)
			maxnropes = nropes;
		prope->mxsnip = mxsnip;
		prope->nsnip = 2;
		prope->psniplist = propesnipl;
		if(s1<s2)
		{
			propesnipl->pos = s1;
			SPSNIP(propesnipl, psnippet, snippetidx);
			++propesnipl;
			propesnipl->pos = s2;
			SPSNIP(propesnipl, psnip2, snip2idx);
		}
		else
		{
			propesnipl->pos = s2;
			SPSNIP(propesnipl, psnip2, snip2idx);
			++propesnipl;
			propesnipl->pos = s1;
			SPSNIP(propesnipl, psnippet, snippetidx);
		}
		psnippet->prope = prope;
		psnip2->prope = prope;
		if(debugpr & 1)
			printf("snip%8X + snip%8X => rope%8X\n",
			psnippet, psnip2, prope);
		nsnipinrope += 2;
	}
/*if first is not in a rope*/
	else if(!psnippet->prope)
	{
	/*add first to rope*/
		if(add2rope( psnip2, s2, psnippet, s1, snippetidx))
			return(1);
		if(debugpr & 1)
			printf("snip%8X => rope%8X\n",
			psnippet, psnip2->prope);
		++nsnipinrope;
	}
/*if second is not in rope*/
	else if(!psnip2->prope)
	{
	/*add second to rope*/
		if(add2rope( psnippet, s1, psnip2, s2, snip2idx))
			return(1);
		if(debugpr & 1)
			printf("snip%8X => rope%8X\n",
			psnip2, psnippet->prope);
		++nsnipinrope;
	}
/*else both in ropes*/
	else
	{
	/*join ropes*/
		if(debugpr & 1)
			printf("join ropes%8X &%8X\n",
			psnippet->prope, psnip2->prope);
		if(joinropes( psnippet, s1, psnip2, s2))
			return(1);
	}
/*return new snippet*/
	*ppnewsnip = psnip2;
	return(0);
}




int chkmatch(
/*check that two snippets match in overlap regions*/
/*return 0 if OK, 1 if mismatch*/
	PSNIPPET psnip1,	/*a snippet*/
	PSNIPPET psnip2,	/*another snippet*/
	int off2)		/*code offset from first to second*/
{
	int s1,e1;		/*range of first snippet*/
	int s2,e2;		/*range of second snippet*/
	int s,e;		/*range of overlap*/
	int olap;		/*length of overlap*/
	int i;
	llong code,code2;	/*sections of overlapped code*/
/*figure range of possible overlap for first snippet*/
	s1 = 0;
	e1 = psnip1->len;
/*figure range of possible overlap for second snippet*/
	s2 = off2;
	e2 = off2 + psnip2->len;
/*figure range of possible overlap*/
	s = s1 < s2 ? s2 : s1;
	e = e1 < e2 ? e1 : e2;
	olap = e - s;
/*see if codes match around this unit*/
	for(i=0;i<olap;i+=32)
	{
		code = extcode(psnip1, s-s1+i,
			i+32 <= olap ? 32 : olap-i);
		code2 = extcode(psnip2, s-s2+i,
			i+32 <= olap ? 32 : olap-i);
		if(code != code2)
			break;
	}
	if(i < olap)
		return(1);
	else
		return(0);
}




void prrope(
/*print contents of rope*/
	PROPE prope)	/*rope to print*/
{
	PSNIPLIST psnipl;	/*pt into snippet list for rope*/
	int n,i;
	PSNIPPET psnip;		/*pt to current snippet in rope*/
/*print number of snippets in rope*/
	printf("%d snippets in rope %X\n", prope->nsnip, prope);
	printf("   pos   len  epos    psnip   code\n");
/*print snippets in rope*/
	for(psnipl=prope->psniplist, n=0; n<prope->nsnip; ++psnipl, ++n)
	{
		psnip = GPSNIP(psnipl);
		printf("%6d%5d%6d %8X ",
		psnipl->pos, psnip->len, psnipl->pos+psnip->len-1, psnip);
		for(i=0;i<psnip->len;++i)
		{
			printf("%c", codechar[extcode(psnip, i, 1)]);
		}
		printf("\n");
	}
}




int chkadd2rope(
/*check if OK to add snippet to rope*/
/*return 0 if OK, 1 if error*/
	PSNIPPET insnip,	/*snippet already in rope*/
	int inpos,		/*relative position of snippet in rope*/
	PSNIPPET newsnip,	/*snippet to be added to rope*/
	int newpos)		/*relative position of new snippet*/
{
	PROPE prope;		/*rope to be added to*/
	int mxsnip;		/*new max snippets in rope list*/
	PSNIPLIST psnipl;	/*list of snippets in rope*/
	int n;
	int pos;		/*relative position for newsnip in rope*/
	int epos;		/*end position*/
	int len;		/*length of new snippet*/
	int err;		/*if got mismatch error in rope*/
/*get rope*/
	prope = insnip->prope;
/*find snippet already in rope*/
	for(psnipl=prope->psniplist,n=0;n<prope->nsnip;++n,++psnipl)
	{
		if(insnip == GPSNIP(psnipl))
			break;
	}
	if(n>=prope->nsnip)
		return(1);
/*get relative position of that snippet*/
	pos = psnipl->pos;
/*figure relative position of new snippet*/
	pos += newpos - inpos;
	epos = pos + newsnip->len;
/*check that new snippet matches all currently in rope*/
	for(psnipl=prope->psniplist,n=0,err=0;!err && (n<prope->nsnip);++n,++psnipl)
	{
	/*check if in range*/
		if(psnipl->pos < pos)
		{
			if(psnipl->pos + GPSNIP(psnipl)->len > pos)
				err |= chkmatch( newsnip, GPSNIP(psnipl),
					 psnipl->pos - pos);
		}
		else
		{
			if(psnipl->pos >= epos)
				break;
			err |= chkmatch( newsnip, GPSNIP(psnipl),
				psnipl->pos - pos);
		}
	}
/*check for mismatch*/
	if(err)
		return(1);
	else
		return(0);
}




int add2rope(
/*add snippet to rope*/
/*return 0 for success, 1 for failure*/
	PSNIPPET insnip,	/*snippet already in rope*/
	int inpos,		/*relative position of snippet in rope*/
	PSNIPPET newsnip,	/*snippet to be added to rope*/
	int newpos,		/*relative position of new snippet*/
	int newsnipidx)		/*index of new snippet in insniplist*/
{
	PROPE prope;		/*rope to be added to*/
	int mxsnip;		/*new max snippets in rope list*/
	PSNIPLIST psnipl;	/*list of snippets in rope*/
	int n;
	int pos;		/*relative position for newsnip in rope*/
	int epos;		/*end position*/
	int len;		/*length of new snippet*/
	int err;		/*if got mismatch error in rope*/
/*get rope*/
	prope = insnip->prope;
/*see if I need to expand snippet list for rope*/
	if(prope->nsnip >= prope->mxsnip)
	{
		mxsnip = 2 * prope->mxsnip;
		//prope->psniplist = (PSNIPLIST)realloc(prope->psniplist,
		//			n = mxsnip * sizeof(SNIPLIST));
		n = mxsnip * sizeof(SNIPLIST);
		prope->psniplist = (PSNIPLIST)REALSNIPLIST(prope->psniplist,
				prope->mxsnip * sizeof(SNIPLIST), n);
		trackmal((void*)(prope->psniplist), n);
		prope->mxsnip = mxsnip;
	}
/*find snippet already in rope*/
	for(psnipl=prope->psniplist,n=0;n<prope->nsnip;++n,++psnipl)
	{
		if(insnip == GPSNIP(psnipl))
			break;
	}
	if(n>=prope->nsnip)
	{
		printf("could not find snippet in it's rope\n");
		exit(1);
	}		
/*get relative position of that snippet*/
	pos = psnipl->pos;
/*figure relative position of new snippet*/
	pos += newpos - inpos;
	epos = pos + newsnip->len;
/*check that new snippet matches all currently in rope*/
	for(psnipl=prope->psniplist,n=0,err=0;!err && (n<prope->nsnip);++n,++psnipl)
	{
	/*check if in range*/
		if(psnipl->pos < pos)
		{
			if(psnipl->pos + GPSNIP(psnipl)->len > pos)
				err |= chkmatch( newsnip, GPSNIP(psnipl),
					 psnipl->pos - pos);
		}
		else
		{
			if(psnipl->pos >= epos)
				break;
			err |= chkmatch( newsnip, GPSNIP(psnipl),
				psnipl->pos - pos);
		}
	}
/*check for mismatch*/
	if(err)
	{
		newsnip->prope = 0;
		printf("mismatch, trying to add snippet %X to rope %X\n", newsnip, prope);
		if(debugpr & 2)
			prrope(prope);
		return(1);
	}
/*insert new snippet in list ordered by relative start position*/
	for(n=prope->nsnip, psnipl=prope->psniplist+n-1; (--n>=0) && (pos <= psnipl->pos);
		--psnipl)
	{
		(psnipl+1)->pos = psnipl->pos;
		SPSNIP(psnipl+1, GPSNIP(psnipl), GISNIP(psnipl));
	}
	++psnipl;
	++n;
	psnipl->pos = pos;
	SPSNIP(psnipl, newsnip, newsnipidx);
	prope->nsnip++;
/*mark snippet in rope*/
	newsnip->prope = prope;
	return(0);
}




int joinropes(
/*join two overlapping ropes via snippets in each*/
/*return 0 for success, 1 for failure*/
	PSNIPPET psnip1,	/*first snippet*/
	int s1,			/*relative offset of first snippet*/
	PSNIPPET psnip2,	/*second snippet*/
	int s2)			/*relative offset of second snippet*/
{
	PSNIPPET tsnip;		/*temp snippet*/
	int t;
	PROPE prope2;		/*rope to be copied from*/
	PSNIPLIST psnipl2;	/*list of snippets in second rope*/
	int n2;			/*number of psnip2 in rope 2*/
	int off2;		/*offset for positions in second rope*/
	PROPE prope1;		/*rope to be copied to*/
	PSNIPLIST psnipl1;	/*list of snippets in first rope*/
	int n1;			/*number of psnip1 in rope 1*/
	int off1;		/*offset for positions in first rope*/
	int i1,e1;		/*comparison range in rope 1 */
	int i2,e2;		/*comparison range in rope 2 */
	int i,n;
	int err;		/*if error*/
	int mxsnip;		/*max number of snippets I have room for*/
	PSNIPLIST psniplm;	/*pt into merged snippet list*/
/*find which rope has fewer entries*/
/*make the rope with fewer entries the second rope*/
	if(psnip2->prope->nsnip > psnip1->prope->nsnip)
	{
		tsnip = psnip1;
		psnip1 = psnip2;
		psnip2 = tsnip;
		t = s1;
		s1 = s2;
		s2 = t;
	}
/*find second snippet already in second rope*/
	prope2 = psnip2->prope;
	for(psnipl2=prope2->psniplist,n2=0;n2<prope2->nsnip;++n2,++psnipl2)
	{
		if(psnip2 == GPSNIP(psnipl2))
			break;
	}
	if(n2>=prope2->nsnip)
	{
		printf("couldn't find snippet 2 in it's rope\n");
		exit(1);
	}
/*figure offset from rope two positions to snippet 2 position*/
	off2 = psnipl2->pos - s2;
/*find first snippet already in first rope*/
	prope1 = psnip1->prope;
	for(psnipl1=prope1->psniplist,n1=0;n1<prope1->nsnip;++n1,++psnipl1)
	{
		if(psnip1 == GPSNIP(psnipl1))
			break;
	}
	if(n1>=prope1->nsnip)
	{
		printf("couldn't find snippet 1 in it's rope\n");
		exit(1);
	}
/*figure offset from rope two positions to snippet 2 position*/
	off1 = psnipl1->pos - s1;
/*make sure snippets around psnip2 will work in rope 1*/
/*compare nearby snippets in two ropes*/
	psnipl1=prope1->psniplist;
	psnipl2=prope2->psniplist;
	i2 = n2 - 8;
	e2 = n2 + 8;
	i2 = i2 < 0 ? 0 : i2;
	e2 = e2 < prope2->nsnip ? e2 : prope2->nsnip;
	i1 = n1 - 8;
	e1 = n1 + 8;
	i1 = i1 < 0 ? 0 : i1;
	e1 = e1 < prope1->nsnip ? e1 : prope1->nsnip;
	for(err=0;!err && (i2 < e2); ++i2)
	{
		for(i=i1;!err && (i < e1); ++i)
		{
			err |= chkmatch(GPSNIP(psnipl1+i), GPSNIP(psnipl2+i2),
				(psnipl2+i2)->pos - off2 - ((psnipl1+i)->pos - off1));
		}
	}
	if(err)
	{
		if(debugpr & 1)
			printf("can't join ropes %X and %X\n",
			psnip1->prope, psnip2->prope);
		return(1);
	}
/*add snippets from rope 2 to rope 1*/
/*see if I need to expand snippet list for rope 1*/
	if(prope1->nsnip + prope2->nsnip >= prope1->mxsnip)
	{
		mxsnip = 2 * prope1->mxsnip;
		//prope1->psniplist = (PSNIPLIST)realloc(prope1->psniplist,
		//			n = mxsnip * sizeof(SNIPLIST));
		n = mxsnip * sizeof(SNIPLIST);
		prope1->psniplist = (PSNIPLIST)REALSNIPLIST(prope1->psniplist,
				prope1->mxsnip * sizeof(SNIPLIST), n);
		trackmal((void*)(prope1->psniplist), n);
		prope1->mxsnip = mxsnip;
	}
/*merge snippet lists from 2 ropes in reverse order*/
	n1 = prope1->nsnip - 1;
	psnipl1 = prope1->psniplist + n1;
	n2 = prope2->nsnip - 1;
	psnipl2 = prope2->psniplist + n2;
	n = n1 + n2 + 1;
	psniplm = prope1->psniplist + n;
	for( ;n>=0; --n, --psniplm)
	{
	/*determine which rope to get the next snippet from*/
	/*pick later start after adjusting for differnt bases in 2 ropes*/
		if((n2 >= 0) && ((n1 < 0) || (psnipl2->pos - off2 > psnipl1->pos - off1)))
		{
		/*rope 2 snippet starts later*/
			SPSNIP(psniplm, GPSNIP(psnipl2), GISNIP(psnipl2));
			psniplm->pos = psnipl2->pos - off2 + off1;
			GPSNIP(psniplm)->prope = prope1;
			--psnipl2;
			--n2;
		}
		else
		{
		/*rope 1 snippet starts later*/
			SPSNIP(psniplm, GPSNIP(psnipl1), GISNIP(psnipl1));
			psniplm->pos = psnipl1->pos;
			--psnipl1;
			--n1;
		}
	}
	prope1->nsnip += prope2->nsnip;
/*delete second rope*/
	//free(prope2->psniplist);
	n = prope2->mxsnip * sizeof(SNIPLIST);
	FREESNIPLIST(prope2->psniplist, n);
	//free(prope2);
	FREEROPE(prope2, sizeof(ROPE));
	--nropes;
	return(0);
}




void propes(
/*print out ropes*/
	int nsnippet,		/*number of snippets*/
	PSNIPPET *insniplist,	/*input snippet list*/
	int nropes)		/*number of ropes*/
{
	PROPE *propespr;	/*list of ropes printed*/
	int nropepr;		/*number of ropes printed*/
	llong ln;		/*number of bytes of hi malloc*/
	PSNIPPET *insnipl;	/*pt into input snippet list*/
	int nin;		/*count for input snippets*/
	PSNIPPET snip;		/*current input snippet*/
	PROPE *pprope;		/*pt into list of ropes printed*/
	int i,n;
	PROPE prope;		/*rope being printed*/
	PSNIPLIST rsnipl;	/*rope snippet list*/
	int nsnip;		/*number of snippets in rope*/
	int pos;		/*print position in rope*/
	int maxpos;		/*max position in rope+1 */
	int cct;		/*char on line count*/
	int soff;		/*starting offset within snippet to print*/
	int len;		/*length to print in snippet*/
	int endcode;		/*trailing odd endcode if any*/
	int haveend;		/*if have trailing odd end code*/
/*allocate space for ropes printed*/
	propespr = (PROPE*)malloc(nropes * sizeof(PROPE));
	nropepr = 0;
/*print some overall statixtics*/
	ln = himalst ==(Uint8*)0xFFFFFFFFFFFFFFF ? 0 : himalend - himalst;
	printf("// %d ropes, %d maxnropes, used ~%ld bytes of small mallocs & %ld large\n",
	nropes, maxnropes, lowmalend - lowmalst, ln);
	if(prsubmal)
		prsubmalloc();
/*loop over input snippets*/
	for(insnipl=insniplist,nin=0;nin<nsnippet;++nin,++insnipl)
	{
	/*check for snippet not in rope*/
		snip = *insnipl;
		if(!snip->prope)
		{
			printf("//snippet %d, len %d, %s %s not in rope\n",
			nin/4+1, snip->len, nin&2 ? "rev" : "fwd", nin&1 ? "mir" : "");
			continue;
		}
	/*check for snippet in rope already printed*/
		for(pprope=propespr,i=0;i<nropepr;++i,++pprope)
		{
			if(*pprope == snip->prope)
				break;
		}
		if(i<nropepr)
			continue;
	/*add this rope to list of already printed*/
		*pprope = snip->prope;
		++nropepr;
	/*have an unprinted rope; print it*/
		prope = snip->prope;
		rsnipl = prope->psniplist;
		nsnip = prope->nsnip;
	/*find length of rope*/
		for(n=0,maxpos=rsnipl->pos;n<nsnip;++n,++rsnipl)
		{
			snip = GPSNIP(rsnipl);
			pos = rsnipl->pos + snip->len;
			if(pos > maxpos)
				maxpos = pos;
		}
		rsnipl = prope->psniplist;
		printf("//rope %d, len=%d\n", nropepr-1, maxpos-rsnipl->pos);
	/*loop over snippets in rope*/
		for(n=0,pos=rsnipl->pos,cct=0;n<nsnip;++n,++rsnipl)
		{
		/*figure out available length to print from this snippet*/
			soff = pos - rsnipl->pos;
			snip = GPSNIP(rsnipl);
			len = snip->len - soff;
			if(len <= 0)
				continue;
		/*print out pairs of codes as hex digits*/
			for(;len>0;--len, ++soff, ++pos)
			{
				printf("%c", codechar[(int)extcode(snip, soff, 1)]);
				if(++cct >= 64)
				{
					printf("\n");
					cct = 0;
				}
			}
		}
	/*end line*/
		if(cct)
			printf("\n");
	}
}




#define SETBIT(base,bit)	*((base)+((bit)>>6)) |= (llong)1 << ((bit)&63);
#define CLRBIT(base,bit)	*((base)+((bit)>>6)) &= ~((llong)1 << ((bit)&63));
#define TSTBIT(base,bit)	(*((base)+((bit)>>6)) & (llong)1 << ((bit)&63))
#define MAXBESTSOL	10	/*max best solutions*/
typedef struct solset	{	/*sequenceof ropes solution set*/
	short entry;	/*which rope*/
	short olap;	/*overlap from previous rope*/
} SOLSET, *PSOLSET;

PSOLSET solset;		/*pt to solution set*/
llong *usedropes;	/*where bitmap of used ropes*/
PROPELIST propelist;	/*list of ropes for joining*/
int maxsolropes;	/*max number of ropes in a solution*/
int bestscore[MAXBESTSOL];	/*best scores so far*/
PSOLSET bestsolset;	/*best solutions so far*/
int nbestsol;		/*number of best solutions so far*/
int solutcnt=0;		/*number of solutions found*/



searchseq(
/*search sequence of ropes to find best sequence*/
	int newent,	/*new entry to add*/
	int nent,	/*number of entries before this one*/
	int olap,	/*overlap from previous rope*/
	int minolap,	/*minimum overlap in this sequence*/
	int sumolap,	/*sum of overlaps in this sequence*/
	int seqlen)	/*total code length of this sequence*/
{
	PROPELIST propel;	/*pt to current rope list item*/
	short fwdmir;	/*idx of fwd mirror rope*/
	short rev;	/*idx of reverse rope*/
	short revmir;	/*idx of reverse mirror rope*/
	int score;	/*score for this solution*/
	int n;
/*put entry in list*/
	solset[nent].entry = newent;
	solset[nent].olap = olap;
/*if this is the second entry initialize the overlap data*/
	if(nent == 1)
		minolap = sumolap = olap;
/*update overlap data*/
	else if(nent > 1)
	{
		minolap = olap < minolap ? olap : minolap;
		sumolap += olap;
	}
/*set mask of used ropes and related ropes*/
	SETBIT(usedropes, newent);
	propel = propelist + newent;
	if((fwdmir = propel->fwdmir) >= 0)
		SETBIT(usedropes, fwdmir);
	if((rev = propel->rev) >= 0)
		SETBIT(usedropes, rev);
	if((revmir = propel->revmir) >= 0)
		SETBIT(usedropes, revmir);
/*add to length of combined rope sequence*/
	++nent;
	seqlen += propel->len - olap;
/*if used about 1/4 of total ropes, maybe save a best solution*/
	if(nent >= maxsolropes)
	{
	/*maybe save solution and score in best solutions list*/
		//score = minolap * sumolap;
		score = -seqlen;
		for(n = nbestsol-1;n>=0;--n)
		{
			if(score < bestscore[n])
				break;
		/*slide existing solution down*/
			if(n<MAXBESTSOL-1)
			{
				bestscore[n+1] = bestscore[n];
				memcpy(bestsolset+(n+1)*maxsolropes,
					bestsolset+n*maxsolropes,
					maxsolropes*sizeof(SOLSET));
			}
		}
	/*save new solution*/
		if(n<MAXBESTSOL-1)
		{
			bestscore[n+1] = score;
			memcpy(bestsolset+(n+1)*maxsolropes, solset,
				maxsolropes*sizeof(SOLSET));
		/*increment solution count*/
			if(nbestsol < MAXBESTSOL)
				++nbestsol;
		}
#if 0
	/*maybe print solution count*/
		++solutcnt;
		printf("%d: bscore %d score %d slut: ",
		solutcnt, bestscore[0], score);
		for(n=0;(n<maxsolropes) && (n<10);++n)
		{
			printf("%d ",solset[n].entry);
		}
		printf("\n");
#endif
	}
/*otherwise loop over choices to extend solution*/
	else
	{
		for(n=0;n<propel->ntails;++n)
		{
		/*skip this choice if rope or relative already used*/
			if(TSTBIT(usedropes, propel->nxtlist[n].next))
				continue;
		/*try this continuation*/
			searchseq( propel->nxtlist[n].next, nent, propel->nxtlist[n].olap,
				minolap, sumolap, seqlen);
		}
	}
/*clear used bits*/
	CLRBIT(usedropes, newent);
	if(fwdmir >= 0)
		CLRBIT(usedropes, fwdmir);
	if(rev >= 0)
		CLRBIT(usedropes, rev);
	if(revmir >= 0)
		CLRBIT(usedropes, revmir);
}




void endjoinropes(
/*join remaining ropes if possible*/
/*this part breaks the minolap rule*/
	int nsnippet,		/*number of input snippets*/
	PSNIPPET *insniplist,	/*input snippets*/
	int ncopy,		/*number of copies of code*/
	int minolap)		/*minimum codes for overlapping snippets*/
{
	int joined;		/*if I joined at least two ropes*/
	int n,i;
	int nropel;		/*number of ropes in list*/
	PSNIPPET *insnipl;	/*pt into input snippet list*/
	int nin;		/*count of input snippets*/
	int totlen;		/*total length of ropes*/
	PSNIPLIST psnipl;	/*pt into snippet list for rope*/
	PSNIPPET snip;		/*current snippet pt*/
	PROPE prope;		/*current rope*/
	PROPELIST propel;	/*pt into rope list*/
	int posmax;		/*max relative end position in rope*/
	int len;		/*rope length*/
	llong mirmask;		/*value to calculate mirror of code snippet*/
	llong reve;		/*reverse of end of rope*/
	llong revs;		/*reverse of start of rope*/
	PSNIPPET psnipmax;	/*snippet that includes end of rope*/
	PROPELIST propel1;	/* 1st pt into rope list*/
	int n1;			/* 1st count of rope list*/
	llong code1;		/*tail code*/
	PSNIPPET psnip1;	/*snippet 1*/
	PROPE prope1;		/*rope 1*/
	PROPELIST propel2;	/* 2nd pt into rope list*/
	int n2;			/* 2nd count of rope list*/
	llong code2;		/*tail code*/
	PSNIPPET psnip2;	/*snippet 2*/
	PROPE prope2;		/*rope 2*/
	PROPELIST propel3;	/*another pointer into rope list*/
	int sh2;		/*amount to shift code*/
	llong lmask;		/*mask for tail/head overlap*/
	llong *startropes;	/*where bitmap of used starts for ropes*/
	int olap;		/*amount of overlap*/
	PSOLSET psol;		/*pt into solution set*/
/*report number of ropes at start of end join process*/
	printf("// %d ropes at start of end join process\n",nropes);
/*done if only 4 ropes left*/
	if(nropes <= 4)
		return;
/*give up if over some max ropes*/
	if(nropes > 512)
	{
		printf("// way too many ropes for end join process\n");
		return;
	}
/*make a list of ropes*/
	propelist = (PROPELIST)malloc(n = nropes * sizeof(ROPELIST));
	memset(propelist, 0, n);
	nropel = 0;
/*loop over input snippets*/
	for(insnipl=insniplist,nin=0,totlen=0;nin<nsnippet;++nin,++insnipl)
	{
	/*check for snippet not in rope*/
		snip = *insnipl;
		if(!snip->prope)
			continue;
		prope = snip->prope;
	/*check for snippet in rope already listed*/
		for(propel = propelist, n = 0; n<nropel; ++n, ++propel)
		{
			if(propel->prope == prope)
				break;
		}
		if(n<nropel)
			continue;
	/*add this rope to rope list*/
		propel->prope = prope;
		propel->psnips = GPSNIP(prope->psniplist);
		psnipl = prope->psniplist + prope->nsnip - 1;
		n1 = prope->nsnip - 1;
		for(posmax = -100,psnipmax = 0,i=ncopy*6; (--i>=0) && (n1>= 0);
		--n1, --psnipl)
		{
			if(psnipl->pos + GPSNIP(psnipl)->len >= posmax)
			{
				posmax = psnipl->pos + GPSNIP(psnipl)->len;
				psnipmax = GPSNIP(psnipl);
			}
		}
		propel->psnipe = psnipmax;
		len = posmax - prope->psniplist->pos;
		propel->len = len;
		totlen += len;
		propel->codes = extcode(propel->psnips, 0, minolap);
		propel->codee = extcode(psnipmax, psnipmax->len - minolap, minolap);
		++nropel;
	}
/*make sure I've found all ropes*/
	if(nropel != nropes)
	{
		printf("mismatch rope count in endjoinropes, %d != %d\n",
		nropes, nropel);;
		exit(1);
	}
/*find 3 matches for each rope, reverse, fwd mirror and rev mirror*/
	mirmask = ((llong)1 << (2*minolap)) - 1;
	for(propel1 = propelist, n1 = 0; n1 < nropel; ++propel1, ++n1)
	{
	/*first clear related ropes*/
		propel1->fwdmir = -1;
		propel1->rev = -1;
		propel1->revmir = -1;
		for(i=0,revs=0;i<minolap;++i)
			revs |= ((propel1->codes >> (2*i)) & 3) << (minolap-i-1)*2;
		for(i=0,reve=0;i<minolap;++i)
			reve |= ((propel1->codee >> (2*i)) & 3) << (minolap-i-1)*2;
		for(propel2=propelist,n2=0,n=0;n2<nropel;++propel2,++n2)
		{
		/*skip same*/
			if(n1 == n2)
				continue;
		/*check that lengths match*/
			if(propel1->len != propel2->len)
				continue;
		/*check for forward mirror*/
			if((propel1->codes == (mirmask - propel2->codes)) &&
			   (propel1->codee == (mirmask - propel2->codee)))
			{
				propel1->fwdmir = n2;
				n |= 1;
			}
		/*check for reverse*/
			if((revs == propel2->codee) && (reve == propel2->codes))
			{
				propel1->rev = n2;
				n |= 2;
			}
		/*check for reverse mirror*/
			if((revs == (mirmask - propel2->codee)) &&
			   (reve == (mirmask - propel2->codes)))
			{
				propel1->revmir = n2;
				n |= 4;
			}
		}
	/*report if don't find fwd mir, rev and rev mir*/
		if(n != 7)
		{
			printf("// can't find 3 related ropes for rope %d len %d, %d %d %d\n",
			n1, propel1->len, propel1->fwdmir, propel1->rev, propel1->revmir);
		}
	}
/*make list of best trailing rope from this one in rank of overlap*/
/*nested loops to match tail of this with head of next*/
	for(propel1 = propelist, n1 = 0; n1 < nropel; ++propel1, ++n1)
	{
	/*get tail of this rope*/
		code1 = propel1->codee;
	/*loop over other ropes*/
		for(propel2=propelist,n2=0;n2<nropel;++propel2,++n2)
		{
		/*skip if same rope*/
			if(n1 == n2)
				continue;
		/*get head of this rope*/
			code2 = propel2->codes;
		/*loop over sizes of overlap*/
			for(i=minolap, lmask=((llong)1<<(2*minolap))-1;
			i>2; --i, lmask >>=2)
			{
				sh2 = (minolap-i)*2;
				if((code1 & lmask) == ((code2>>sh2) & lmask))
					break;
			}
		/*record best N tail to head matches for each rope*/
			if(i>2)
			{
			/*insert this entry in list*/
				for(n = propel1->ntails-1;n>=0;--n)
				{
					if(propel1->nxtlist[n].olap > i)
						break;
					if(n < MAXTAILS-1)
					{
						propel1->nxtlist[n+1].olap =
							propel1->nxtlist[n].olap;
						propel1->nxtlist[n+1].next =
							propel1->nxtlist[n].next;
					}
				}
				if(++n < MAXTAILS)
				{
					propel1->nxtlist[n].olap = i;
					propel1->nxtlist[n].next = n2;
					if(propel1->ntails < MAXTAILS)
						propel1->ntails++;
				}
			}
		}
	}
#if 1
/*print data about each rope*/
	for(propel1 = propelist, n1 = 0; n1 < nropel; ++propel1, ++n1)
	{
		printf("//rope%3d: len%9d fm%3d r%3d rm%3d, %d tails: ",
		n1, propel1->len, propel1->fwdmir, propel1->rev, propel1->revmir,
		propel1->ntails);
		for(i=0;i<propel1->ntails;++i)
		{
			printf("%3d %2d, ",propel1->nxtlist[i].next,
				propel1->nxtlist[i].olap);
		}
		printf("\n");
	}
#endif
/*there may be too many ropes to get a good solution in a reasonable search time*/
	if(nropes > 128)
	{
		printf("too many ropes, %d, to get a good solution\n",nropes);
		return;
	}
/*look for best solution set*/
/*initialize for search*/
	n = ((nropes+63) >> 6) * 8;
	startropes = (llong*)malloc(n);
	usedropes = (llong*)malloc(n);
	memset(startropes, 0, n);
	nbestsol = 0;
	maxsolropes = nropes/4;
	solset = (PSOLSET)malloc(maxsolropes*sizeof(SOLSET));
	bestsolset = (PSOLSET)malloc(maxsolropes*MAXBESTSOL*sizeof(SOLSET));
/*loop over ropes as possible start rope*/
	for(n1=0, propel1=propelist;n1<nropel;++n1,++propel1)
	{
	/*see if this rope or a related one already used*/
		if(TSTBIT(startropes, n1))
			continue;
	/*mark this rope and related, used as start*/
		SETBIT(startropes, n1);
		if(propel1->fwdmir >= 0)
			SETBIT(startropes, propel1->fwdmir);
	/*clear ropes used bitmap*/
		memset(usedropes, 0, (nropes+7)>>3);
	/*search for best sequence, starting with this rope*/
		searchseq(n1, 0, 0, 0, 0, 0);
	}
#if 0
/*print best solutions*/
	for(n=0;n<nbestsol;++n)
	{
		printf("//solution %d, score %d\n", n, bestscore[n]);
		for(i=0,psol=bestsolset+n*maxsolropes;i<maxsolropes;++i,++psol)
		{
			printf("//%4d: %4d %4d\n", i, psol->entry, psol->olap);
		}
	}
#endif
/*report solutions found*/
	printf("//found at least %d end join sequences with scores of %d to %d\n",
	nbestsol, bestscore[0], bestscore[nbestsol-1]);
	printf("// using best end join sequence\n");
/*join ropes in best sequence*/
	for(n=1,n1=bestsolset[0].entry;n<maxsolropes;++n,n1=n2)
	{
	/*join forward ropes*/
		olap = bestsolset[n].olap;
		n2 = bestsolset[n].entry;
		propel1 = propelist + n1;
		propel2 = propelist + n2;
		psnip1 = propel1->psnipe;
		psnip2 = propel2->psnips;
		if(debugpr & 1)
			printf("end join ropes%8X &%8X\n",
			propel1->prope, propel2->prope);
		joinropes(psnip1, olap - psnip1->len, psnip2, 0);
	/*maybe join forward mirror ropes*/
		if((propel1->fwdmir >= 0) && (propel2->fwdmir >= 0))
		{
			propel1 = propelist + propel1->fwdmir;
			propel2 = propelist + propel2->fwdmir;
			psnip1 = propel1->psnipe;
			psnip2 = propel2->psnips;
			if(debugpr & 1)
				printf("end join ropes%8X &%8X\n",
				propel1->prope, propel2->prope);
			joinropes(psnip1, olap - psnip1->len, psnip2, 0);
		}
	/*maybe join reverse ropes*/
		propel1 = propelist + n1;
		propel2 = propelist + n2;
		if((propel1->rev >= 0) && (propel2->rev >= 0))
		{
			propel1 = propelist + propel1->rev;
			propel2 = propelist + propel2->rev;
			psnip2 = propel1->psnips;
			psnip1 = propel2->psnipe;
			if(debugpr & 1)
				printf("end join ropes%8X &%8X\n",
				propel1->prope, propel2->prope);
			joinropes(psnip1, olap - psnip1->len, psnip2, 0);
		}
	/*maybe join reverse mirror ropes*/
		propel1 = propelist + n1;
		propel2 = propelist + n2;
		if((propel1->revmir >= 0) && (propel2->revmir >= 0))
		{
			propel1 = propelist + propel1->revmir;
			propel2 = propelist + propel2->revmir;
			psnip2 = propel1->psnips;
			psnip1 = propel2->psnipe;
			if(debugpr & 1)
				printf("end join ropes%8X &%8X\n",
				propel1->prope, propel2->prope);
			joinropes(psnip1, olap - psnip1->len, psnip2, 0);
		}
	}
/*free rope list*/
	free(propelist);
}




void dounits(
/*create database of code units of length unitlen in the snippets*/
	int nsnippet,		/*number of input snippets*/
	PSNIPPET *insniplist,	/*input snippets*/
	int unitlen,		/*length of each unit*/
	int both,		/*if track both ends of snippets*/
	PUNIT *punitdb)		/*where to return pt to unit database*/
{
	PUNIT unitdb;	/*database of where units occur in snippets*/
	int i,j,n;
	int nsnip;	/*which snippet I'm on*/
	PSNIPPET psnippet;	/*pt to snippet I'm on*/
	int len;	/*length of this snippet*/
	llong code;	/*unit length of code*/
	llong idx;	/*which unit*/
	int unit;	/*which of all possible units*/
	PUNIT punit;	/*pt to current unit*/
	int mxsnip;	/*max snippets I currently have room for*/
	PSNIPLIST newsnip;	/*where put snippet data*/
	llong ln;	/*long long n*/
/*check reasonable unit len*/
	if((unitlen<4) || (unitlen>16))
	{
		printf("bad unit length %d, should be between 4 and 16\n",unitlen);
		exit(1);
	}
/*allocate space for unit database*/
/*start, end and spare for each snippet in database*/
/*but don't allocate more an max possible from unitlen*/
	ln = (llong)1 << (2*unitlen);
	nunitdb = 3 * nsnippet;
	if(nunitdb > ln)
		nunitdb = ln;
	unitdb = (PUNIT)malloc(ln = nunitdb * sizeof(UNIT));
	trackmal((void*)unitdb, ln);
	if(!unitdb)
	{
		printf("can't allocate %lX = %ld bytes for %ld unitdb\n",
		ln, ln, nunitdb);
		exit(1);
	}
	if(debugpr & 4)
	{
		printf("//unitdb:        %12lX-%12lX  %12lX  %10d\n",
		unitdb, unitdb+nunitdb, ln, nunitdb);
	}
	//memset(unitdb, 0, ln);
	*punitdb = unitdb;
/*mark entries empty*/
	for(i=0,punit=unitdb;i<nunitdb;++i,++punit)
	{
		punit->code = (llong)(-1);
		punit->nsnip = 0;
	}
/*fill in valid units, start and end of snippets*/
	for(nsnip=0;nsnip<nsnippet;++nsnip)
	{
	/*point to snippet data*/
		psnippet = insniplist[nsnip];
		len = psnippet->len;
	/*enter start and end snippet*/
		code = extcode(psnippet, 0, unitlen);
		punit = findunit(code, 1);
	/*only do front end for first pass to eliminate subset snippets*/
		if(both)
		{
			code = extcode(psnippet, len-unitlen, unitlen);
			punit = findunit(code, 1);
		}
	}
/*count the number of times each unit appears in the snippets*/
/*loop over each input snippet*/
	for(nsnip=0,idx=0;nsnip<nsnippet;++nsnip)
	{
	/*point to snippet data*/
		psnippet = insniplist[nsnip];
		len = psnippet->len;
	/*loop over all units in this snippet*/
		for(i=0;i<len-unitlen+1;++i)
		{
		/*extract unit of code*/
			code = extcode(psnippet, i, unitlen);
			punit = findunit(code, 0);
			if(punit)
			{
				punit->nsnip++;
				++idx;
			}
		}
	}
/*convert that count list into offsets into a concatenated sniplist array*/
	for(i=0,idx=0,punit=unitdb;i<nunitdb;++i,++punit)
	{
		j = punit->nsnip;
		punit->nsnip = idx;
		idx += j;
	}
	ncatsnip = idx;
/*allocate space for concatenated sniplist list for each unit*/
	if(bigspace && (bigspacesz >= (ln = idx*sizeof(SNIPLIST))))
	{
		catsniplist = (PSNIPLIST)realloc(bigspace, ln);
		bigspace = 0;
	}
	else
	{
		if(bigspace)
			free(bigspace);
		catsniplist = (PSNIPLIST)malloc(ln = idx*sizeof(SNIPLIST));
	}
	if(!catsniplist)
	{
		printf("can't allocate %lX = %ld bytes for %ld catsniplist entries\n",
		ln, ln, idx);
		exit(1);
	}
	trackmal((void*)catsniplist, ln);
	if(debugpr & 4)
	{
		printf("//catsniplist:   %12lX-%12lX  %12lX  %10d\n",
		catsniplist, catsniplist+idx, ln, idx);
	}
/*for each unit, fill in pointer in it's snippet list*/
	for(i=0;i<nunitdb;++i)
	{
		punit = unitdb + i;
		punit->psniplist = catsniplist + punit->nsnip;
		punit->nsnip = 0;
	}
/*loop over each input snippet*/
	for(nsnip=0;nsnip<nsnippet;++nsnip)
	{
	/*point to snippet data*/
		psnippet = insniplist[nsnip];
		len = psnippet->len;
	/*loop over all units in this snippet*/
		for(i=0;i<len-unitlen+1;++i)
		{
		/*extract unit of code*/
			code = extcode(psnippet, i, unitlen);
			punit = findunit(code, 0);
		/*if this is a unit of interest, record matching unit in other snippets*/
			if(punit)
			{
			/*add this snippet to the snippet list for this unit*/
				newsnip = punit->psniplist + punit->nsnip;
				newsnip->pos = i;
				SPSNIP(newsnip, psnippet, nsnip);
				punit->nsnip++;
			}
		}
	}
}




llong extcode(
/*extract section of code*/
	PSNIPPET psnippet,	/*snippet from which to extract*/
	int i,			/*start position*/
	int len)		/*length of section*/
{
	int epos;	/*end position of code section*/
	llong code;	/*code to return*/
/*calculate end position*/
	epos = i + len - 1;
/*extract unit of code*/
/*check for crossing code word*/
	if((i>>5) == (epos>>5))
		code = psnippet->code[epos>>5] >> (62-2*(epos&31));
	else
	{
		code = psnippet->code[i>>5] << (((epos&31)+1) * 2);
		code |= (psnippet->code[epos>>5] >> (62-2*(epos&31))) &
			(((llong)1 << (((epos&31)+1) * 2)) - 1);
	}
	if(len < 32)
		code &= ((llong)1 << (2*len)) - 1;
	return(code);
}




#define MXBUF 65536		/*input buffer size*/
Uint8 bigbuf[MXBUF];		/*input buffer*/


void readsnips(
/*read snippets from file*/
	char *filebase,		/*file base name*/
	int minsniplen,		/*minimum snippet length*/
	int *pnsnippet,		/*where return number of snippets*/
	PSNIPPET **ppsniplist)	/*where to return snippet data*/
{
	char filename[100];	/*name of snippet file*/
	struct stat statbuf;	/*file status buffer*/
	llong snipfilesize;	/*size of snippet file*/
	Uint8 *rdbuf;		/*pt to file char buffer*/
	Uint8 **rdlines;	/*pt to list of lines in file buffer*/
	int nlines;		/*number of lines in snippet file*/
	int ifd;		/*input file descriptor*/
	llong fp;		/*index into file*/
	int maxline;		/*max line length*/
	int linelen;		/*length of this line*/
	int rdsize;		/*number of characters to read from file*/
	PSNIPPET *insniplist;	/*list of input snippets and reverse*/
	int line;		/*input line number*/
	int nsnip;		/*number of input snippets and reverse*/
	Uint8 *lp;		/*pointer into line*/
	int len;		/*snippet length*/
	Uint8 c;		/*current character*/
	PSNIPPET pfwdsnip;	/*forward snippet*/
	PSNIPPET prevsnip;	/*reverse snippet*/
	PSNIPPET pfwdsnipmir;	/*forward snippet, mirror of helix*/
	PSNIPPET prevsnipmir;	/*reverse snippet, mirror of helix*/
	int i,n;
	llong *cp;		/*pt into snippet code*/
	int li;			/*where in code word to insert next piece of code*/
	llong code;		/*word of code*/
	int dig;		/*hex digit of code*/
/*form snippet file name*/
	strcpy(filename, filebase);
	strcat(filename, ".snip");
/*get size of snippet file*/
	if(stat(filename, &statbuf) != 0)
	{
		printf("can't get size of %s\n", filename);
		exit(1);
	}
	snipfilesize = statbuf.st_size;
/*maybe allocate a big contiguous block to free up just before I need another block*/
	if(bigspacesz)
	{
		bigspace = (char*)malloc(bigspacesz);
		if(debugpr & 4)
		{
			printf("// bigspace:     %12lX-%12lX  %12lX  %10d\n",
			bigspace, bigspace + bigspacesz, bigspacesz, bigspacesz);
		}
	}
/*count lines in file*/
/*get max line length*/
	if((ifd = open(filename, 0)) < 0)
	{
		printf("can't open %s\n", filename);
		exit(1);
	}
	for(fp = 0, maxline = 0, nlines = 0,linelen = 0; fp < snipfilesize; fp += rdsize)
	{
		if(snipfilesize - fp > MXBUF)
			rdsize = MXBUF;
		else
			rdsize = snipfilesize - fp;
		if(read(ifd, bigbuf, rdsize) != rdsize)
		{
			printf("can't read %d bytes from %s at %ld\n",
			rdsize, filename, fp);
			exit(1);
		}
		for(lp = bigbuf, i=rdsize; --i >= 0; )
		{
			if(((c = *lp++) == '\n') || (c == '\r'))
			{
				if(linelen > 0)
					++nlines;
				if(linelen >= maxline)
					maxline = linelen+1;
				linelen = 0;
			}
			else
				++linelen;
		}
	}
/*check that longest line fits in buffer*/
	if(maxline >= MXBUF)
	{
		printf("max line length of %d char too long for line buffer of %d char\n",
		maxline, MXBUF);
		exit(1);
	}
/*rewind input file*/
	lseek(ifd, 0, 0);
#if 0
/*allocate space to read file*/
	rdbuf = (Uint8*)malloc(snipfilesize+8);
	rdlines = (Uint8**)malloc((n = snipfilesize/(minsniplen+4)+1)*sizeof(Uint8*));
	if(debugpr & 4)
	{
		printf("//file buffer:   %12lX-%12lX  %12lX\n",
		rdbuf, rdbuf+snipfilesize+8, snipfilesize+8);
		printf("//line pointers: %12lX-%12lX  %12lX  %10d\n",
		rdlines, rdlines+n, n*sizeof(Uint8*), n);
	}
	if((nlines = readfile(filename, rdlines, rdbuf, snipfilesize)) < 1)
	{
		printf("can't read all the data in %s\n", filename);
		exit(1);
	}
#endif
/*allocate a list of snippets*/
	insniplist = (PSNIPPET*)malloc(n =4*nlines*sizeof(PSNIPPET));
	if(debugpr & 4)
	{
		printf("// in sniplist:  %12lX-%12lX  %12lX  %10d\n",
		insniplist, insniplist+4*nlines, n, 4*nlines);
	}
	trackmal(insniplist, n);
/*loop over lines in snippet file*/
	for(fp = 0, rdsize = 0, line = 0, lp = bigbuf, nsnip = 0; line < nlines; ++line)
	{
	/*maybe read a new buffer*/
		if((lp - bigbuf + maxline >= rdsize) && (fp < snipfilesize))
		{
		/*copy existing part line down*/
			memcpy(bigbuf, lp, n = rdsize - (lp - bigbuf));
			lp = bigbuf;
			if(snipfilesize - fp > MXBUF - n)
				rdsize = MXBUF - n;
			else
				rdsize = snipfilesize - fp;
			if(read(ifd, bigbuf+n, rdsize) != rdsize)
			{
				printf("can't read %d bytes from %s at %ld\n",
				rdsize, filename, fp);
				exit(1);
			}
			fp += rdsize;
			rdsize += n;
		}
	/*start a new line*/
	/*skip comment lines and blank lines*/
		if(((c = *lp++) == '/') || (c == '\r') || (c == '\n'))
		{
			while(c != '\n')
				c = *lp++;
			continue;
		}
	/*read size of snippet*/
		for(len=0; 1; c = *lp++)
		{
			if((c == ' ') || (c =='\t'))
			{
				if(len==0)
					continue;
				else
					break;
			}
			else if((c>='0') && (c<='9'))
				len = len*10 + c - '0';
			else
			{
				printf("can't read size of snippet on line %d of %s\n",
				line+1, filename);
				exit(1);
			}
		}
		if(!c || !len)
		{
			printf("can't read snippet on line %d of %s\n",
			line+1, filename);
			exit(1);
		}
	/*allocate space for fwd and rev snippets*/
		//pfwdsnip = (PSNIPPET)malloc(n = sizeof(SNIPPET) + (len+3)/4);
		n = ((len + 31) >> 5) - 2;	//extra space for code
		if(n < 0)
			n = 0;
		n = sizeof(SNIPPET) + n * 8;
		//pfwdsnip = (PSNIPPET)malloc(n);
		//pfwdsnip = (PSNIPPET)malloc(n);
		//prevsnip = (PSNIPPET)malloc(n);
		//pfwdsnipmir = (PSNIPPET)malloc(n);
		pfwdsnip = (PSNIPPET)MALSNIPPET(n);
		prevsnip = (PSNIPPET)MALSNIPPET(n);
		pfwdsnipmir = (PSNIPPET)MALSNIPPET(n);
		prevsnipmir = (PSNIPPET)MALSNIPPET(n);
		trackmal(prevsnipmir, n);
		memset(pfwdsnip, 0, n);
		memset(prevsnip, 0, n);
		memset(pfwdsnipmir, 0, n);
		memset(prevsnipmir, 0, n);
	/*link snippets into list*/
		insniplist[nsnip++] = pfwdsnip;
		insniplist[nsnip++] = prevsnip;
		insniplist[nsnip++] = pfwdsnipmir;
		insniplist[nsnip++] = prevsnipmir;
	/*fill in length of snippets*/
		pfwdsnip->len = prevsnip->len = pfwdsnipmir->len = prevsnipmir->len = len;
	/*fill in snippet data from input line*/
		while((c==' ') || (c=='\t'))
			c = *lp++;
		for(cp=pfwdsnip->code,li=62,code=0,n=0; (c!='\n') && (c!= '\r') ; c = *lp++)
		{
			if((c=='A') || (c=='a'))
				dig = 0;
			else if((c=='C') || (c=='c'))
				dig = 1;
			else if((c=='G') || (c=='g'))
				dig = 2;
			else if((c=='T') || (c=='t'))
				dig = 3;
			else
			{
				printf("bad snippet data on line %d of %s\n",
				line+1, filename);
				exit(1);
			}
			code |= (llong)dig << li;
			++n;
			if((li -= 2) < 0)
			{
				*cp++ = code;
				code = 0;
				li = 62;
			}
		}
		if(c == '\r')
			c = *lp++;
		if(li != 62)
			*cp = code;
		if(len != n)
		{
			printf("bad snippet data length on line %d of %s\n",
			line+1, filename);
			exit(1);
		}
	/*convert forward codes to reverse codes*/
		for(n=len;--n>=0;)
		{
			i = len - n - 1;
			prevsnip->code[i>>5] |= extcode(pfwdsnip, n, 1) << ((31-(i&31))*2);
		}
	/*create mirror versions of fwd and reverse codes*/
		for(n=0;n<((len+31)>>5);++n)
		{
			pfwdsnipmir->code[n] = (llong)(-1) - pfwdsnip->code[n];
			prevsnipmir->code[n] = (llong)(-1) - prevsnip->code[n];
		}
	}
/*return snippets*/
	*pnsnippet = nsnip;
	*ppsniplist = insniplist;
}




static void trackmal(
/*track range of mallocs*/
	void *start,	/*start of malloced range*/
	llong len)	/*length of malloced range in bytes*/
{
	Uint8 *s;	/*start of malloced area*/
	Uint8 *e;	/*end of malloced area*/
/*quit if malloc failed*/
	if(!start)
	{
		printf("failed to malloc %ld of space\n", len);
		exit(1);
	}
/*for high and low addresses*/
	s = (Uint8*)start;
	e = s + len;
/*split into two ranges*/
	if((llong)s & ((llong)0xF << (4*11)))
	{
		himalst = s < himalst ? s : himalst;
		himalend = e > himalend ? e : himalend;
	}
	else
	{
		lowmalst = s < lowmalst ? s : lowmalst;
		lowmalend = e > lowmalend ? e : lowmalend;
	}
}




static int lfsrrand(void)
{
/*use a 31 bit Linear Feedback Shift Register (LFSR) with taps at 31 and 28 */
/*process it 8 bits at a time in parallel */
/*this value will not repeat for almost 2^31 cycles*/
/*it distributes the random values evenly over 0-255*/
/*the values have no sequential correlation*/
	int feedback;	/*feedback for LFSR*/
/*do 8 1 bit LFSR steps in parallel*/
	feedback = ((randseed ^ (randseed>>3)) & 0xFF) ^ 0xFF;
	randseed = (randseed>>8) | (feedback<<23);
	return(randseed);
}




static int getrand(
/*get random value in a specified range*/
	int srange,	/*start of range*/
	int erange)	/*end of range*/
{
	int range;	/*range of value to return*/
	int val;	/*value being formed*/
	int mask;	/*mask for power of 2 ranges*/
	range = erange-srange+1;
	val = lfsrrand();
	if(!(range&(mask=range-1)))
		val = val & mask;
	else
		val = val % range;
	val = val+srange;
	return(val);
}




int readfile(
/*open file, read into memory and break into lines*/
/*return 0 or number of lines*/
	char *name,	/*name of file*/
	Uint8 **lines,	/*where to put pointers to lines*/
	Uint8 *bigbuf,	/*where to put file contents*/
	llong maxsize)	/*max size of file*/
{
	int ifd;	/*input file descriptor*/
	llong size;	/*size of file*/
	llong nread;	/*number of bytes to read*/
	llong didread;	/*number of bytes I read*/
	int nlines;	/*number of lines in file*/
	Uint8 *cp;	/*pointer into file*/
	Uint8 **lp;	/*pt into list of line pointers*/
	Uint8 *end;	/*pt to end of file*/
	Uint8 c;	/*char in file*/
/*open file*/
	if((ifd=open(name,0))<0)
	{
		printf("can't open %s\n",name);
		return(0);
	}
/*read whole file into big buffer*/
	for(size=0;size<maxsize;size+=nread)
	{
		nread = maxsize - size;
		nread = nread > 1000000000 ? 1000000000 : nread;
		if((didread = read(ifd, bigbuf + size, nread)) != nread)
		{
			printf("can't read bytes %ld to %ld of %s\n",
			size, size + nread, name);
			exit(1);
		}
	}
	close(ifd);
/*break it into lines*/
/*null terminate the lines*/
/*make a list of pointers to all the lines*/
	for(nlines=0,cp=bigbuf,lp=lines,*lp++=cp,end=bigbuf+size;cp<end;)
	{
	/*check for end of line*/
		if((c= *cp++) == '\n')
		{
			*(cp-1)=0;
			if(*(cp-2) == '\r')
				*(cp-2)=0;
			*lp++=cp;
			++nlines;
		}
	}
/*return number of lines*/
	return(nlines);
}
