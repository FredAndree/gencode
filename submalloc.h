/* famalloc.c	sub malloc allocator to reduce overhead and mem fragmentation*/
/*  5-May-2016	freda creation*/

#ifndef SUBMALLOC_H
#define SUBMALLOC_H

#include "platform.h"
#include "types.h"

Uint8* submalloc(
/*try to submalloc a block of memory*/
	int size);	/*size of item*/

Uint8* subfree(
/*try to submalloc a block of memory*/
	Uint8* addr,	/*address of item to be freed*/
	int size);	/*size of item*/

Uint8* subrealloc(
/*try to subrealloc a block of memory*/
	Uint8 *addr,	/*address of existing block*/
	int origsize,	/*original size of item*/
	int newsize);	/*new size of item*/

int submalinit(
/*init sub malloc mechanism*/
/*return 0 for success, 1 for failure*/
	int bigsize,	/*big malloc size, a power of two*/
	int nsmsize,	/*number of small sizes*/
	int *smsizes);	/*list of small sizes*/

void prsubmalloc(void);
/*print submalloc statistics*/

#endif
