/**
 * General definitions for the NSCS code
 */

#ifndef H_NSCS_COMMON
#define H_NSCS_COMMON
#include <stdbool.h>
#include <stdlib.h>
//The type of the integer used for indexing.
//For small problems, up to c.a. 2G use int 
//alternatively use ptrdiff_t.
#define csi int 

//
#define CONE_TYPES 1

//Validation results 
//
#define OK 0
#define VALIDATION_OK 0
#define VALIDATION_ERROR 1
#define BACKTRACK_FAIL 2

//Other errors
#define CENTER_ITERS_EXCEEDED -2
#define OUT_OF_MEMORY -99
#define INTERNAL_ERROR -98
#endif
