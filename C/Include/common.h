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
#define INVALID_PARAMETER 3
#define MISSING_ARGUMENT 4

//Other errors
#define OUT_OF_MEMORY -99
#define INTERNAL_ERROR -98

//Stopping reasons
#define END_CONTINUE  20
#define END_OPTIMAL   21
#define END_ILL_POSED 22
#define END_P_INFEAS  23
#define END_D_INFEAS  24
#define END_MAX_CENTER_ITER 25
#define END_MAX_ITER 26

#endif
