#ifndef LOGGER_H
#define LOGGER_H

#include "Timer.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

struct Logger {
	int status;
	int rank;
	FILE* file;
};

static struct Logger *logger;

void logger_init(int status);

void logger_set_on(void);

void logger_set_off(void);

void logger_set_rank(int rank);

void logger_log(char *symbol, char *msg);

void logger_info(char *symbol, char *fmt,  ...);
	
void logger_end(void);

#endif
