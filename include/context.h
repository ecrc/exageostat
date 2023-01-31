/**
 *
 * @file context.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon context header
 *
 * @version 1.2.0
 * @author Jakub Kurzak
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2022-11-09
 *
 */
#ifndef _chameleon_context_h_
#define _chameleon_context_h_

#include "chameleon/struct.h"

/**
 *  Routines to handle threads context
 */
#ifdef __cplusplus
extern "C" {
#endif

CHAM_context_t *chameleon_context_create();

CHAM_context_t *chameleon_context_self();

int chameleon_context_destroy();

#ifdef __cplusplus
}
#endif

#endif /* _chameleon_context_h_ */
