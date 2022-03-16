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
//REAL cycles_numbre;
//REAL var_sqrt;
//typedef struct mathsbox_params_t* mathsbox_pointer_t;

void error_function(REAL argument, mathsbox_t *restrict mathsbox);

//static inline REAL square_root_of(REAL number);

/****************************************************************************
 *   Error function with integral computing by midpoint method OK
 *   Will do the Simpson if ITCM is ok
 *   
 *   Sampling of error function is maybe connected to the time_step need to
 *   investigate.
 *  
 *    expk take  : ~ 570 bytes
      sqrtk take : ~1250 bytes
 *****************************************************************************/
/*
void error_function(REAL argument, mathsbox_t *restrict mathsbox){

    mathsbox->err_func = 0.;
    REAL step = argument/mathsbox->error_func_sample;
    REAL x;
    REAL t;
    //REAL Pi = REAL_CONST(3.1415927);// here was a k
    REAL two_over_sqrt_Pi = REAL_CONST(1.128379167); //APPROXIMATION
    REAL Erf = ZERO;
    REAL Erfc = ZERO;
    
    for(x=0; x<=argument; x+=step){
        
        //Erfc +=  factor*(2/sqrtk(Pi))*expk(-(t*t)); // the real one overflowed ITCM because of expk and sqrtk
        t = x + REAL_HALF(step);
        Erf +=  step*two_over_sqrt_Pi*expk(-(t*t)); //working like this one
        //Erf +=  step*two_over_sqrt_Pi*(-(t*t));//TEST
        //Erf +=  step*(REAL_CONST(2.)/sqrtk(Pi))*expk(-(t*t)); // TEST sqrtk ONE
    }
    Erfc = ONE-Erf;

    mathsbox->err_func = Erfc;

}

static inline s1615 square_root_of(REAL number)
{
     //!! square root method take from Quake III Arena 
     //! source code, attribute to John Carmack.
     //! Under GNU GPL licence.
     //!
     //! Implement here just to see if something was lighter than sqrt.c
     //!
     //!
  
    REAL x, y, f;
    REAL i;
    f = REAL_CONST(1.5);
    
    x = REAL_HALF(number);
    y = number;
    
    i = *(REAL *) &y;
    i = 0x5f3759df - (i >> 1);
    y = *(REAL *) &i;
    y = y * (f - (x * y * y));
    y = y * (f - (x * y * y));
    
    return number*y;
    
}
*/
#endif