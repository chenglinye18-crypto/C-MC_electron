#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <syslog.h>
#include <errno.h>

#define TO3D(i,j,k) ((j) * numx * numz + ((i) * numz + (k)))
#define MAXLINE 80
#define MAX_ENTRY_NUM 10
#define DIR_X 0
#define DIR_Y 1
#define DIR_Z 2

extern void err_quit(const char *, ...);
extern void err_sys(const char * , ...);
extern void err_doit(int errnoflag, int level, const char * fmt, va_list ap);
