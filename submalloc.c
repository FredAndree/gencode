/* famalloc.c	sub malloc allocator to reduce overhead and mem fragmentation*/
/*  5-May-2016	freda creation*/

#include "platform.h"
#include "types.h"

typedef struct malhead	{ /*malloc head for each different size item*/
	Uint16 itemsize;	/*size of low level item*/
	Uint16 itempblock;	/*number of itmes per block*/
	Uint16 bmapsize;	/*bitmap size in bytes*/
	struct malblock *curblock; /*current block to do suballocations from*/
	struct malblock *blocklist; /*linked list of blocks of this size items*/
	int flags;		/*1=didfree */
	int nmalloc;		/*number of big mallocs*/
	int nfree;		/*number of big frees*/
	int nsubmalloc;		/*number of sub mallocs*/
	int nsubfree;		/*number of sub frees*/
	int nblockchk;		/*number of times I checked block for space*/
	int nslotchk;		/*number of times I checked block for slots*/
	int nnewcur;		/*number of times I reset curblock from free*/
} MALHEAD, *PMALHEAD;
#define DIDFREE 1	/*if freed a block of this size*/

typedef struct malblock	{ /*malloc block to do suballocations from*/
	Uint8 *base;	/*base of this block*/
	struct malblock *prevblock;	/*link to previous block*/
	struct malblock *nextblock;	/*link to next block*/
	Uint16 usedslots;	/*number of sub malloc slots used*/
	Uint16 nextslot;	/*next slot to try*/
	Uint16 keycode;		/*key value*/
	Uint16 totslots;	/*total sub malloc slots available*/
	llong bitmap[1];		/*bitmap of used suballocations*/
} MALBLOCK, *PMALBLOCK;

static int powtwo;	/* approx saze of large mallocs*/
static int malsize;	/* exact size of large mallocs*/
static int nsizehash;	/* number of hash entries for different sizes*/
static int firstmalloc=1;	/* if first malloc of run*/
static PMALHEAD pheads;	/* table of headers of different sizes*/
#define KEYCODE	0xAA55	/* key value*/




Uint8* submalloc(
/*try to submalloc a block of memory*/
	int size)	/*size of item*/
{
	int idx;	/*index into hashed header table*/
	PMALHEAD phead;	/*head of submalloc info for this size*/
	int i,j;
	PMALBLOCK pblock; /*pt to current allocation block for this size*/
	PMALBLOCK pbl2;	/*next or previous block*/
	Uint8 *base;	/*base of allocation*/
/*look up items of that size*/
	idx = size % nsizehash;
	phead = pheads + idx;
	if(phead->itemsize != size)
	{
		for(i=nsizehash; --i >= 0; )
		{
			if(++idx >= nsizehash)
				idx = 0;
			phead = pheads + idx;
			if(phead->itemsize == size)
				break;
		}
		if(i < 0)
			return(0);
	}
	phead->nsubmalloc++;
/*if no current block or current block full, allocate and initialize a new block*/
	if(!(pblock = phead->curblock) || (pblock->usedslots == pblock->totslots))
	{
		if(pblock)
		{
			phead->nblockchk++;
			if((pbl2 = pblock->nextblock) && (pbl2->usedslots < pbl2->totslots))
				pblock = pbl2;
			else if((pbl2 = pblock->prevblock) && (pbl2->usedslots < pbl2->totslots))
				pblock = pbl2;
			else if(phead->flags & DIDFREE)
			{
				for(pblock=phead->blocklist;pblock;pblock=pblock->nextblock)
				{
					phead->nblockchk++;
					if(pblock->usedslots < pblock->totslots)
						break;
				}
			}
			if(pblock)
				phead->curblock = pblock;
		}
	/*see if I now have a block with space*/
		if(!pblock || (pblock->usedslots == pblock->totslots))
		{
		/*if this is first malloc of run, allocate enough mem to avoid seg fault*/
			if(firstmalloc)
			{
				Uint8 *base2;
				firstmalloc = 0;
				base = (Uint8*)malloc(malsize);
				pblock = (PMALBLOCK) (((llong)(base) + powtwo - 1) & (-powtwo));
				i = (Uint8*)pblock - base - sizeof(Uint8*);
				i = i < sizeof(Uint8*) ? sizeof(Uint8*) : i;
				base2 = (Uint8*)realloc(base, i);
				//printf("first allocation: base %lX, size %X, base2 %lX\n",
				//	base, i, base2);
			}
		/*clear DIDFREE*/
			phead->flags &= ~DIDFREE;
		/*create and initialize a new malloc block*/
			base = (Uint8*)malloc(malsize);
			if(!base)
				return(0);
			phead->nmalloc++;
		/*point to block data within the malloc, on a power of 2 boundary*/
			pblock = (PMALBLOCK) (((llong)(base) + powtwo - 1) & (-powtwo));
		/*fill in header*/
			pblock->base = base;
			pblock->prevblock = 0;
			pblock->nextblock = phead->blocklist;
			if(phead->blocklist)
				phead->blocklist->prevblock = pblock;
			phead->blocklist = pblock;
			phead->curblock = pblock;
			pblock->totslots = phead->itempblock;
			pblock->usedslots = 0;
			pblock->nextslot = 0;
			pblock->keycode = KEYCODE ^ size;
			memset(pblock->bitmap, 0, phead->bmapsize * sizeof(llong));
		/*take out slots used by block header*/
			i = ((Uint8*)pblock - base) / size;
			j = ((Uint8*)(&(pblock->bitmap[phead->bmapsize])) - 1 - base) / size;
			for(;i <= j; ++i)
			{
				if(i < phead->itempblock)
					pblock->totslots--;
				pblock->bitmap[i >> 6] |= (llong)1 << (i & 63);
			}
		/*check no slots left for user*/
			if(!(pblock->totslots))
			{
				printf("no user submalloc slots of size %d within %d\n",
				size, malsize);
				return(0);
			}
		}
	}
/*find next available subunit*/
	for(i=phead->itempblock; --i >= 0; )
	{
	/*if next unit available*/
		j = pblock->nextslot;
		if(++(pblock->nextslot) >= phead->itempblock)
			pblock->nextslot = 0;
		phead->nslotchk++;
		if(!(pblock->bitmap[j >> 6] & ((llong)1 << (j & 63))))
		{
		/*mark slot used*/
			pblock->bitmap[j >> 6] |= ((llong)1 << (j & 63));
			pblock->usedslots++;
		/*return pointer*/
			return(pblock->base + phead->itemsize * j);
		}
	/*if bitmap very full, advance more quickly*/
		if(pblock->bitmap[j >> 6] == (llong)(-1))
		{
			pblock->nextslot += 63 - (j & 63);
			if(pblock->nextslot >= phead->itempblock)
				pblock->nextslot = 0;
		}
	}
/*should not get here*/
/*count says a slot available, but I could not find it*/
	printf("can't find open sub malloc slot of size %d, %d used of %d in %lX\n",
	size, pblock->usedslots, pblock->totslots, pblock);
	exit(1);
}




#ifdef DEBUG
static int firstfree=1;
#endif
Uint8* subfree(
/*try to submalloc a block of memory*/
	Uint8* addr,	/*address of item to be freed*/
	int size)	/*size of item*/
{
	int idx;	/*index into hashed header table*/
	PMALHEAD phead;	/*head of submalloc info for this size*/
	int i,j;
	PMALBLOCK pblock; /*pt to current allocation block for this size*/
	int item;	/*sub malloc item number in this malloc block*/
	Uint8 *base;	/*base of this malloc block*/
/*look up items of that size*/
	idx = size % nsizehash;
	phead = pheads + idx;
	if(phead->itemsize != size)
	{
		for(i=nsizehash; --i >= 0; )
		{
			if(++idx >= nsizehash)
				idx = 0;
			phead = pheads + idx;
			if(phead->itemsize == size)
				break;
		}
		if(i < 0)
			return(0);
	}
	phead->nsubfree++;
	phead->flags |= DIDFREE;
/*find tentative base for this block of suballocations*/
	pblock = (PMALBLOCK) ((llong)addr & (-powtwo));
#ifdef DEBUG
	if(firstfree)
		printf("subfree: pblock=%lX\n",pblock);
#endif
/*see if this is the right base*/
/*check that I now have the right base*/
	for(i=0;i<2;++i)
	{
		if((addr < pblock->base) ||
		(addr >= (pblock->base + powtwo)) ||
		(pblock->keycode != (KEYCODE ^ size)) ||
		(pblock->totslots > phead->itempblock) ||
		(pblock->usedslots > pblock->totslots) ||
		(pblock->nextslot >= phead->itempblock) ||
		(((addr - pblock->base) % phead->itemsize) != 0))
		{
			if(i == 0)
			{
				pblock = (PMALBLOCK) ((llong)pblock + powtwo);
				continue;
			}
			else
			{
				printf("subfree of non submalloced item, size %d, addr %lX\n",
				size, addr);
				printf("pblock %lX, base %lX, end %lX, keycode %X,\n",
				pblock, pblock->base, pblock->base + powtwo, pblock->keycode);
				printf("uslots %d, tslots %d, nslot %d, itpblk %d\n",
				pblock->usedslots, pblock->totslots,
				pblock->nextslot, phead->itempblock);
				printf("nextblock %lX, prevblock %lX\n",
				pblock->nextblock, pblock->prevblock);
				exit(1);
			}
		}
		break;
	}
#ifdef DEBUG
	if(firstfree)
	{
		printf("subfree: pblock+=%lX\n",pblock);
		firstfree=0;
	}
#endif
/*free the item*/
	item = (addr - pblock->base) / phead->itemsize;
	pblock->bitmap[item >> 6] &= ~((llong)1 << (item & 63));
	pblock->usedslots--;
/*if none used in this block, can free it*/
	if(!(pblock->usedslots))
	{
		base = pblock->base;
		pblock->base = 0;
		pblock->keycode = 0;
		if(pblock->nextblock)
			pblock->nextblock->prevblock = pblock->prevblock;
		if(pblock->prevblock)
			pblock->prevblock->nextblock = pblock->nextblock;
		if(pblock == phead->blocklist)
			phead->blocklist = pblock->nextblock;
		if(pblock == phead->curblock)
			phead->curblock = phead->blocklist;
		free(base);
		phead->nfree++;
	}
/*else if more free slots here than in cur block, make this cur block*/
	else if(pblock->usedslots < phead->curblock->usedslots)
	{
		phead->curblock = pblock;
		phead->nnewcur++;
	}
}




Uint8* subrealloc(
/*try to subrealloc a block of memory*/
	Uint8 *addr,	/*address of existing block*/
	int origsize,	/*original size of item*/
	int newsize)	/*new size of item*/
{
	int idx;	/*index into hashed header table*/
	PMALHEAD orighead;	/*head of submalloc info for orig size*/
	PMALHEAD newhead;	/*head of submalloc info for new size*/
	int i,j;
	int csize;	/*size of data to copy*/
	Uint8 *newaddr;	/*address of new memory*/
/*look up items of original size*/
	idx = origsize % nsizehash;
	orighead = pheads + idx;
	if(orighead->itemsize != origsize)
	{
		for(i=nsizehash; --i >= 0; )
		{
			if(++idx >= nsizehash)
				idx = 0;
			orighead = pheads + idx;
			if(orighead->itemsize == origsize)
				break;
		}
		if(i < 0)
			orighead = 0;
	}
/*look up items of new size*/
	idx = newsize % nsizehash;
	newhead = pheads + idx;
	if(newhead->itemsize != newsize)
	{
		for(i=nsizehash; --i >= 0; )
		{
			if(++idx >= nsizehash)
				idx = 0;
			newhead = pheads + idx;
			if(newhead->itemsize == newsize)
				break;
		}
		if(i < 0)
			newhead = 0;
	}
/*figure copy size*/
	csize = origsize < newsize ? origsize : newsize;
/*handle the various cases of submalloc and regular malloc*/
	if(!newhead & !orighead)
		return(realloc(addr, newsize));
	else if(newhead && !orighead)
	{
		newaddr = submalloc(newsize);
		if(!newaddr)
			return(0);
		memcpy(newaddr, addr, csize);
		free(addr);
		return(newaddr);
	}
	else if(!newhead && orighead)
	{
		newaddr = malloc(newsize);
		if(!newaddr)
			return(0);
		memcpy(newaddr, addr, csize);
		subfree(addr, origsize);
		return(newaddr);
	}
	else
	{
		newaddr = submalloc(newsize);
		if(!newaddr)
			return(0);
		memcpy(newaddr, addr, csize);
		subfree(addr, origsize);
		return(newaddr);
	}
}




int submalinit(
/*init sub malloc mechanism*/
/*return 0 for success, 1 for failure*/
	int bigsize,	/*big malloc size, a power of two*/
	int nsmsize,	/*number of small sizes*/
	int *smsizes)	/*list of small sizes*/
{
	int i,j,n;
	int size;	/*size of current small items*/
	int idx;	/*hash table index*/
	PMALHEAD phead;	/*head data for this size small item*/
	int extrasize;	/*extra space needed in case block head falls at end*/
/*check and save big size*/
	if(bigsize & (bigsize-1))
		return(1);
	powtwo = bigsize;
/*allocate an array of MALHEAD structures to be hash searched*/
	nsizehash = ((nsmsize * 3) / 2) | 1;
	pheads = (PMALHEAD) malloc(n = nsizehash * sizeof(MALHEAD));
	if(!pheads)
		return(1);
	memset(pheads, 0, n);
/*clear malsize*/
	malsize = 0;
/*loop over sizes*/
	for(i=0; i<nsmsize; ++i)
	{
		size = smsizes[i];
	/*set up hash code entry*/
	/*look up items of that size*/
		idx = size % nsizehash;
		phead = pheads + idx;
		if(phead->itemsize)
		{
			for(j=nsizehash; --j >= 0; )
			{
				if(++idx >= nsizehash)
					idx = 0;
				phead = pheads + idx;
				if(!phead->itemsize)
					break;
			}
			if(j < 0)
				return(1);
		}
	/*fill in that header*/
		phead->itemsize = size;
		phead->itempblock = powtwo / size;
		phead->bmapsize = (phead->itempblock + 63) >> 6;
	/*check the amount of extra size I need*/
		extrasize = sizeof(MALBLOCK) + (phead->bmapsize - 1) *
			 sizeof(llong) - sizeof(Uint8*);
		if(powtwo + extrasize > malsize)
			malsize += powtwo + extrasize;
	}
	return(0);
}




void prsubmalloc(void)
/*print submalloc statistics*/
{
	int i;
	PMALHEAD phead;	/*header block for size of intwerest*/
	printf("idx size malloc free submalloc subfree blk_chk slot_chk new_cur\n");
/*loop over sizes*/
	for(i=0; i<nsizehash; ++i)
	{
		phead = pheads + i;
		if(!(phead->itemsize))
			continue;
		printf("%3d:%4d%6d%6d%8d%8d%8d%8d%8d\n",
		i, phead->itemsize, phead->nmalloc, phead->nfree,
		phead->nsubmalloc, phead->nsubfree,
		phead->nblockchk, phead->nslotchk, phead->nnewcur);
	}
}
