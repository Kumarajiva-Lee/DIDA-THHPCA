#include "logger.h"


void logger_init(int status) {
	if (logger != NULL)
		free(logger);
	logger = malloc(sizeof(logger));
	logger->status = status;
}

void logger_set_off(void) {
	if (logger == NULL)
		return;
	logger->status = 0;
}

void logger_set_on(void) {
	if (logger == NULL)
		return;
	logger->status = 1;
}

void logger_set_rank(int rank) {
	if (logger == NULL || !logger->status) {
		return;
	}

	logger->rank = rank;
	char filename[20];
	sprintf(filename, "logs/log.%d", rank);
	logger->file = fopen(filename, "w");
}

void logger_log(char *symbol, char *msg) {
	if (logger == NULL || !logger->status) {
		return;
	}

	double now;
	get_time(now);
	if (logger->file == NULL) {
		printf("/logs not found or rank not specified!");
		return;
	}
	fprintf(logger->file, "time: %lf tag:%s,%s \n", now, symbol, msg);
	fflush(logger->file);
}

void logger_info(char *symbol, char *fmt,  ...) {
	if (logger == NULL || !logger->status) {
		return;
	}

	if (logger->file == NULL) {
		printf("rank not specified!");
		return;
	}
	double now;
	get_time(now);
	fprintf(logger->file, "time: %lf tag:%s ", now, symbol);
	
	va_list _va_list;
	va_start(_va_list, fmt);
	vfprintf(logger->file, fmt, _va_list);
	va_end(_va_list);
	fprintf(logger->file, "\n");
	fflush(logger->file);
}

void logger_end() {
	if (!logger->status) {
		return;
	}

	if (logger->file != NULL) {
		fclose(logger->file);
	}
}

