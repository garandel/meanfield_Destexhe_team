#ifndef _MATHSBOX_H_
#define _MATHSBOX_H_


#include <math.h>
#include <stdfix-exp.h>
#include <stdfix-full-iso.h>
//#include "../../common/maths-util.h"




struct mathsbox_t;

typedef struct mathsbox_t {

    REAL error_func_sample;
    
    REAL err_func;
    
    
    
}mathsbox_t;


void error_function(REAL argument, mathsbox_t *restrict mathsbox);

// TO DO: replace midpoint by Simpson methods in a second time.

/****************************************************************************
     expk take  : ~ 570 bytes
      sqrtk take : ~1250 bytes
 *****************************************************************************/



#endif