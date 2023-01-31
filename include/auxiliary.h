/**
 *
 * @file auxiliary.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon auxiliary header
 *
 * @version 1.2.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @author Mathieu Faverge
 * @date 2022-11-09
 *
 */
#ifndef _chameleon_auxiliary_h_
#define _chameleon_auxiliary_h_

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "chameleon/struct.h"
#include "chameleon/tasks.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Get environment variable
 */
#if defined(CHAMELEON_OS_WINDOWS)

static inline int
chameleon_setenv( const char *var, const char *value, int overwrite ) {
    return !(SetEnvironmentVariable( var, value ));
}

static inline char *
chameleon_getenv( const char *var ) {
    char *str;
    int len = 512;
    int rc;
    str = (char*)malloc(len * sizeof(char));
    rc = GetEnvironmentVariable(var, str, len);
    if (rc == 0) {
        free(str);
        str = NULL;
    }
    return str;
}

static inline void
chameleon_cleanenv( char *str ) {
    if (str != NULL) free(str);
}

#else /* Other OS systems */

static inline int
chameleon_setenv(const char *var, const char *value, int overwrite) {
    return setenv(var, value, overwrite);
}

static inline char *
chameleon_getenv(const char *var) {
    return getenv(var);
}

static inline void
chameleon_cleanenv(char *str) {
    (void) str;
}

#endif


static inline int
chameleon_env_is_set_to(char *str, char *value) {
    char *val;
    if ((val = chameleon_getenv(str)) &&
        !strcmp(val, value))
        return 1;
    return 0;
}

static inline int
chameleon_env_is_on(char *str) {
    return chameleon_env_is_set_to(str, "1");
}

static inline int
chameleon_env_is_off(char *str) {
    return chameleon_env_is_set_to(str, "0");
}

static inline int
chameleon_getenv_get_value_int(char *string, int default_value) {
    long int ret;
    char *str = chameleon_getenv(string);
    if (str == NULL) return default_value;

    if (sscanf(str, "%ld", &ret) != 1) {
        perror("sscanf");
        return default_value;
    }

    return (int) ret;
}

/**
 *  Internal routines
 */
void chameleon_warning(const char *func_name, const char *msg_text);

void chameleon_error(const char *func_name, const char *msg_text);

void chameleon_fatal_error(const char *func_name, const char *msg_text);

int chameleon_rank(CHAM_context_t *chamctxt);

int chameleon_tune(cham_tasktype_t func, int M, int N, int NRHS);

#ifdef __cplusplus
}
#endif

#endif /* _chameleon_auxiliary_h_ */
