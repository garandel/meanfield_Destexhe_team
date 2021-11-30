#ifndef _MATHSBOX_H_
#define _MATHSBOX_H_

//#include <sqrt.h>
#include <math.h>
#include <stdfix-exp.h>
#include <stdfix-full-iso.h>
//#include "../../common/maths-util.h"




struct mathsbox_t;

typedef struct mathsbox_t {

    REAL error_func_sample;
    
    REAL err_func;
    
    REAL var_sqrt;
}mathsbox_t;

//typedef struct mathsbox_params_t* mathsbox_pointer_t;

void error_function(REAL argument, mathsbox_t *restrict mathsbox);
#endif