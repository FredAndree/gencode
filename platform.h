/* platform.h	platform dependent manifests*/

#ifndef	PLATFORM_H
#define PLATFORM_H

//#define ONPC		//uncomment this for windows
#define LINUX		//uncomment this for Linux on PC

#ifdef ONPC
#include <io.h>
#include <process.h>
#endif
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdio.h>
#include <memory.h>
#include <string.h>
#include <malloc.h>

#ifdef ONPC
#define open(asdf,qwer)	_open(asdf, _O_BINARY | _O_RDONLY)
#define creat(asdf,qwer) _open(asdf, _O_CREAT|_O_TRUNC|_O_BINARY|_O_RDWR,_S_IWRITE)
#define read	_read
#define write	_write
#define lseek	_lseek
#define close	_close
#define unlink	_unlink
#endif

#endif //PLATFORM_H
