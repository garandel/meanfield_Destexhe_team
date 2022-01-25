/*
 * Copyright (c) 2017-2019 The University of Manchester
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

//! \file
//! \brief Izhekevich neuron implementation
#include "../../meanfield/models/meanfield_model_impl.h"

#include <debug.h>
#include "../../meanfield/models/params_from_network.h"
#include "../../meanfield/models/mathsbox.h"
#include "../../meanfield/models/P_fit_polynomial.h"
//#include "../../common/maths-util.h" // i.o to use SQRT(x) and SQR(a)

//! The global parameters of the Izhekevich neuron model
static const global_neuron_params_t *global_params;

/*! \brief For linear membrane voltages, 1.5 is the correct value. However
 * with actual membrane voltage behaviour and tested over an wide range of
 * use cases 1.85 gives slightly better spike timings.
 */
//static const REAL SIMPLE_TQ_OFFSET = REAL_CONST(1.85);


//! \brief The original model uses 0.04, but this (1 ULP larger?) gives better
//! numeric stability.
//!
//! Thanks to Mantas Mikaitis for this!
//static const REAL MAGIC_MULTIPLIER = REAL_CONST(0.040008544921875);

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

/*
//! \brief Multiplies an accum by an accum giving an integer answer.
//! \param[in] y: An s16.15 accum
//! \param[in] x: An s16.15 accum
//! \return The integer part of y*x.

static inline int mulkk(
    s1615 y,
    s1615 x)
{
    return __I((__I64(bitsk(y)) * __I64(bitsk(x))) >> 15);    
}
*/


/*
//! \brief Multiplies a frac by an accum giving an integer answer.
//! \param[in] x: An s0.31 frac
//! \param[in] y: An s16.15 accum
//! \return The integer part of y*x.
static inline int mulrk(
    s031 x,
    s1615 y)
{
    return __I((__I32(bitslr(x)) * __I32(bitsk(y))) >> 8);    
}
*/

/*
//! \brief convert an long frac in integer for rounding giving an accum answer.
//! \param[in] y: An s16.15 accum
//! \param[in] x: An s16.15 accum
//! \return The integer part of y*x.

static inline int_k_t lrkbits( const s031 f)
{
    union { s1615 r; int_k_t inter_k; int_lr_t inter_lr; s031 fx; } x;
    
    x.fx = f;
    int32_t round = __stdfix_round_s32(x.inter_lr, 15);
    x.inter_k = round;
      
    return x.inter_k;
}
*/
    




void threshold_func(ParamsFromNetwork_t *restrict pNetwork, pFitPolynomial_t *restrict Pfit)
{
    
    //!    threshold function coming from :
    //!    Neural Computation 31, 653–680 (2019) doi:10.1162/neco_a_01173
    //!    where P's are polynomials parameters involve in a 
    //!    voltage-effective threshold.
    
    /*
    setting by default to True the square
    because when use by external modules, coeff[5:]=np.zeros(3)
    in the case of a linear threshold
    */
    /*
        If add other inv_constante, this go on DTCM so a priori no problem.
        
    
    */
   
    REAL muV0 = pNetwork->muV0;
    REAL iDmuV0 = pNetwork->iDmuV0;

    REAL sV0 = pNetwork->sV0;
    REAL iDsV0 = pNetwork->iDsV0;

    REAL TvN0 = pNetwork->TvN0;
    REAL iDTvN0 = pNetwork->iDTvN0;

    REAL muV = pNetwork->muV;
    REAL sV = pNetwork->sV;
 
    REAL TvN = pNetwork->TvN;
    REAL Vthre = pNetwork->Vthre;
    

    
    REAL P0 = Pfit->P0;
    REAL P1 = Pfit->P1;
    REAL P2 = Pfit->P2;
    REAL P3 = Pfit->P3;
    REAL P5 = Pfit->P5;
    REAL P6 = Pfit->P6;
    REAL P7 = Pfit->P7;
    REAL P8 = Pfit->P8;
    REAL P9 = Pfit->P9;
    REAL P10 = Pfit->P10;

    Vthre = P0\
        + P1*(muV-muV0)*iDmuV0\
        + P2*(sV-sV0)*iDsV0\
        + P3*(TvN-TvN0)*iDTvN0\
        + P5*((muV-muV0)*iDmuV0)*((muV-muV0)*iDmuV0)\
        + P6*((sV-sV0)*iDsV0)*((sV-sV0)*iDsV0)\
        + P7*((TvN-TvN0)*iDTvN0)*((TvN-TvN0)*iDTvN0)\
        + P8*((muV-muV0)*iDmuV0)*((sV-sV0)*iDsV0)\
        + P9*((muV-muV0)*iDmuV0)*((TvN-TvN0)*iDTvN0)\
        + P10*((sV-sV0)*iDsV0)*((TvN-TvN0)*iDTvN0);
    /*
    int_k_t test;
    test = P1*(muV-muV0)*iDmuV0;
    log_info("Vthre = %3.4k", test);
    */

  
    pNetwork->Vthre = Vthre;

    }
    

void get_fluct_regime_varsup(REAL Ve, REAL Vi, REAL W, ParamsFromNetwork_t *restrict pNetwork)
{
    // Need some comments
    
    REAL gei = pNetwork->gei;
    REAL pconnec = pNetwork->pconnec;
    REAL Ntot = pNetwork->Ntot;
    REAL Qe = pNetwork->Qe;
    REAL Qi = pNetwork->Qi;
    REAL Te = pNetwork->Te;
    REAL Ti = pNetwork->Ti;
    REAL Gl = pNetwork->Gl;
    REAL El = pNetwork->El;
    REAL Ei = pNetwork->Ei;
    REAL Ee = pNetwork->Ee;
    REAL Cm = pNetwork->Cm;
    
    //log_info("gei_k = %1.2k", gei);
    //log_info("Gl_k = %2.4k", Gl);
    //log_info("El_k = %2.4k", El);
        
    //log_info("Ve_k=%3.6k", Ve_k);
    //log_info("gei=%1.2k", gei);
    //log_info("pconnec=%2.6k", pconnec);
    //log_info("Ntot=%6.1k", Ntot);
    
    
    REAL Fe;
    Fe = Ve * (1-gei)*pconnec*Ntot; // default is 1 !!
    REAL Fi;
    Fi = Vi * gei*pconnec*Ntot;
    //log_info("Fe=%3.6k", Fe);
    
    /* normaly = Ve*Te*Qe*Ke with Ke = p*N_exc what it is?
        -> here N_exc = (1-gei)*Ntot*pconnec
        So give the same
        
        Avec muGi=Qi*Ti*Fi+1 ça marche!!????!!  No comprendo!! bcs multiplication give zero
    */
    REAL muGe;
    muGe = Qe*Te*Fe;
    REAL muGi;
    muGi = Qi*Ti*Fi;

    REAL muG;
    muG = Gl + muGe + muGi;
    
    //SOME ERRORS

    //log_info("muG = %6.6k", muG);
    /*
    if (muG<1 && muG>=0){
        muG = 1;
    }  
    else if (muG>-1 && muG<0)
    {
        muG = -1;
    }
    */
    
    /*
    if (muG==0)
    {
        log_error("muG env 0");
        muG=1;
    }
    */
    
    /*
    REAL test = Gl*El;
    log_info("%11.4k",test);
    */
    /*if (test==0)
    {
        log_error("test=0");
        log_info("%11.4k",test); 
    }
    */
    
    
    //int_k_t var_test = __stdfix_smul_k(Gl,El); //(muGe*Ee + muGi*Ei + Gl*El - W_k)/muG;
    //int_k_t muV_k;
    //muV_k = __stdfix_smul_k(Gl,gei);// Gl*gei;//bitsk(pNetwork->gei*pNetwork->Gl);//(muGe*Ee + muGi*Ei + Gl*El)/muG;
    //int_k_t muV_kk = var_test + muV_k;
    //REAL var_test = pNetwork->Gl * pNetwork->gei;
    /*
    if (muV_k>-1 && muV_k <1)
    {
        log_error("muV_k env 0");
    }*/
    //log_info("muV_k = %3.4k", muV_k);
    //log_info("var_test = %3.4k", var_test);
    
    REAL muV  = (muGe*Ee + muGi*Ei + Gl*El)/muG;
    pNetwork->muV = muV ;
    /*
    if (muV_k==0)
    {
        log_error("muV_k Gl = 0");
        //muV_k=1;
    }
    */
    log_info("muV = %11.4k", pNetwork->muV);


    
    //double problem muV_k et  Qe*(Ee+1)/muG
    
    pNetwork->muGn = muG/Gl;
    REAL Tm = Cm/muG;
    REAL Ue = Qe*(Ee-muV)/muG;
    REAL Ui = Qi*(Ei-muV)/muG;
    /*
    if (Fe<1 && Fe>0)//just to insure a non zero division,
    {
        Fe = 1;
    }
    else if (Fe>-1 && Fe<0)
    {
        Fe = -1 ;
    }
    
    if (Fi<1 && Fi>0) 
    {
        Fi = 1;    
    }
    else if (Fi>-1 && Fi<0)
    {
        Fi = -1;
    }
    */
    
    if (Fe==0)
    {
        Fe=1;
    }
    else if (Fi==0)
    {
        Fi=1;
    }
    
    //problem is multiplication '*' that give 0 bcs not saturated arithm
    REAL Tv_num = Fe*(Ue*Te)*(Ue*Te) + Fi*(Ti*Ui)*(Ti*Ui);//too long
    REAL Tv_denom = Fe*(Ue*Te)*(Ue*Te)/(Te+Tm) + Fi*(Ti*Ui)*(Ti*Ui)/(Ti+Tm) ;// too long
    /*
    if (Tv_denom<1){
        Tv_denom += 1;
    }
    
    if (Tv_num<1 && Tv_num>0)
    {
        Tv_num=1;
    }
    if (Tv_num>-1 && Tv_num<0)
    {
        Tv_num=-1;
    }*/
    
    //Tv_num give Error so maybe egal to zero
    //Tv_denom aswell

    REAL Tv = Tv_num / Tv_denom;
    //int_k_t TvN_k = Tv*Gl/Cm;
    
    //log_info("Tv=%11.4k", Tv);
    
    //ERROR : With this method TvN is egal to zero Tv*Gl/Cm
        
    pNetwork->TvN = Tv*Gl/Cm; // TvN is adimensional so usefull var
    
    REAL sV_sqr = REAL_CONST(0.50000)*(Tv_denom);
    
    pNetwork->sV = square_root_of(sV_sqr); // //  sqrtk(sV_sqr);//  kbits(sV_sqr);

}


void TF(REAL Ve, REAL Vi, REAL W,
        ParamsFromNetwork_t *restrict pNetwork,
        pFitPolynomial_t *restrict Pfit,
        mathsbox_t *restrict mathsbox){

/*
    State-variables are directly connected to the struct
    parameters are put in local in order to make the code clear.

*/
    
    
    if (pNetwork->Fout_th != ZERO){
        pNetwork->Fout_th = ACS_DBL_TINY;
    }

    if (Ve < ACS_DBL_TINY){
        Ve += ACS_DBL_TINY;
    }
    if (Vi < ACS_DBL_TINY){
        Vi += ACS_DBL_TINY;
    }

    get_fluct_regime_varsup(Ve, Vi, W, pNetwork);
    threshold_func(pNetwork, Pfit);

    
    /*
    normalement sqrt:
        argument = (pNetwork->Vthre - pNetwork->muV)/sqrtk(REAL_CONST(2.))/pNetwork->sV;

    */
    if (pNetwork->sV<ACS_DBL_TINY){
        pNetwork->sV += ACS_DBL_TINY;
    }
    //factor = REAL_HALF(Gl/(pNetwork->TvN * Cm));
    //REAL argument = (pNetwork->Vthre - pNetwork->muV)/(REAL_CONST(1.4142137)*pNetwork->sV);
    REAL argument = (pNetwork->Vthre - pNetwork->muV)/(REAL_CONST(1.4142137)+pNetwork->sV);
    

    error_function(argument, mathsbox);

    
    REAL Gl = pNetwork->Gl;
    REAL Cm = pNetwork->Cm;
    /*
    pNetwork->Fout_th = (HALF*Gl) * mathsbox->err_func / (Cm*pNetwork->TvN);// In fact = 1/(2.*Tv) * err_func , that's it'!!!
    If remove that's will do less instruction->NOP
    
    Put TvN<-:Tv because Tv not in pNetwork
    REMOVE this correction bcs TvN adimensional so usefull
    
    Some problem with sqrt type give a DIVBY0 error when compil with python
    
    pNetwork->Fout_th = (HALF*Gl) * mathsbox->err_func / (Cm*pNetwork->TvN);
    */
    pNetwork->Fout_th = (HALF*Gl) * mathsbox->err_func / (Cm*pNetwork->TvN);


    if (pNetwork->Fout_th < ACS_DBL_TINY){
        pNetwork->Fout_th += ACS_DBL_TINY;
    }
    
}


void RK2_midpoint_MF(REAL h, meanfield_t *meanfield,
                     ParamsFromNetwork_t *restrict pNetwork,
                     pFitPolynomial_t *restrict Pfit_exc,
                     pFitPolynomial_t *restrict Pfit_inh,
                     mathsbox_t *restrict mathsbox) {
    
    /* Propose for now a=0
    *
    */

    REAL lastVe = meanfield->Ve;
    REAL lastVi = meanfield->Vi;
    REAL lastW = meanfield->w;
    
    REAL tauw = meanfield->tauw;
    REAL T_inv = meanfield->Timescale_inv;
    REAL b = meanfield->b;
               
    
    TF(lastVe, lastVi, lastW, pNetwork, Pfit_exc, mathsbox);    
    REAL lastTF_exc = pNetwork->Fout_th;
    
    
    TF(lastVe, lastVi, lastW, pNetwork, Pfit_inh, mathsbox);
    REAL lastTF_inh = pNetwork->Fout_th;
    
/******************************************************
 *   EULER Explicite method
 *   It's very instable if for now h<0.2 for 20ms
 *   
 *   NEED to give also the error of the method here :
 *   0.5*h^2*u''(t_n) + o(h^2)
 *******************************************************/
    
    /*
    
    REAL k1_exc = (lastTF_exc - lastVe)*T_inv;
    meanfield->Ve += h * k1_exc ;
    
    REAL k1_inh = (lastTF_inh - lastVi)*T_inv;
    meanfield->Vi += h * k1_inh ;
    
    REAL k1_W = -lastW/tauw + meanfield->b * lastVe;
    meanfield->w += h * k1_W;
    
    */
    
/***********************************
 *  RUNGE-KUTTA 2nd order Midpoint *
 ***********************************/
    
    REAL k1_exc = (lastTF_exc - lastVe)*T_inv;
    REAL alpha_exc = lastVe + h*k1_exc;
    REAL k2_exc = (lastTF_exc - alpha_exc )*T_inv;
    
    meanfield->Ve += REAL_HALF(h*(k1_exc + k2_exc));
        
    REAL k1_inh = (lastTF_inh - lastVi)*T_inv;
    REAL alpha_inh = lastVi + h*k1_inh;
    REAL k2_inh = (lastTF_inh - alpha_inh)*T_inv;
    
    meanfield->Vi += REAL_HALF(h*(k1_inh + k2_inh));
    
    REAL k1_W = -lastW/tauw + b * lastVe;
    REAL alpha_w = lastW + h*k1_W;
    REAL k2_W = -alpha_w/tauw + b * lastVe;
 
    meanfield->w += REAL_HALF(h*(k1_W+k2_W));


}

void meanfield_model_set_global_neuron_params(
        const global_neuron_params_t *params) {
    global_params = params;
}

/*perhaps when we will do more than one MF we could uses "num_excitatory_inputs" like the number of ex MF and in MF?
  and maybe is there some contamanation from the neightbourest neighbour MF!
*/
state_t meanfield_model_state_update(
    meanfield_t *restrict meanfield,
    ParamsFromNetwork_t *restrict pNetwork,
    pFitPolynomial_t *restrict Pfit_exc,
    pFitPolynomial_t *restrict Pfit_inh,
    mathsbox_t *restrict mathsbox){
    /*
        uint16_t num_excitatory_inputs, const input_t *exc_input,
		uint16_t num_inhibitory_inputs, const input_t *inh_input,
		input_t external_bias, meanfield_t *restrict meanfield,
        ParamsFromNetwork_t *restrict pNetwork) {
    REAL total_exc = 0;
    REAL total_inh = 0;

    for (int i =0; i<num_excitatory_inputs; i++) {
        total_exc += exc_input[i];
    }
    for (int i =0; i<num_inhibitory_inputs; i++) {
        total_inh += inh_input[i];
    }

    //input_t input_this_timestep = total_exc - total_inh
    //        + external_bias + neuron->I_offset;
    */

    // the best AR update so far
    RK2_midpoint_MF(meanfield->this_h,
                    meanfield,
                    pNetwork,
                    Pfit_exc,
                    Pfit_inh,
                    mathsbox);
    meanfield->this_h = global_params->machine_timestep_ms;

    return meanfield->Ve;
}



void neuron_model_has_spiked(meanfield_t *restrict meanfield) {
    log_debug("in neuron_model_has_spiked, time is ",
              global_params->machine_timestep_ms);
    // reset membrane voltage
    //neuron->V = neuron->C;

    // offset 2nd state variable
    //neuron->U += neuron->D;

    // simple threshold correction - next timestep (only) gets a bump
    //neuron->this_h = global_params->machine_timestep_ms * SIMPLE_TQ_OFFSET;
}

//change name neuron -> meanfield and membrane -> rate
state_t meanfield_model_get_firing_rate_Ve(const meanfield_t *meanfield) {
    return meanfield->Ve;
}

state_t meanfield_model_get_firing_rate_Vi(const meanfield_t *meanfield) {
    return meanfield->Vi;
}

state_t meanfield_model_get_adaptation_W(const meanfield_t *meanfield){
    return meanfield->w;
}


void meanfield_model_print_state_variables(const meanfield_t *meanfield) {
    log_debug("Ve = %11.4k ", meanfield->Ve);
    log_debug("Vi = %11.4k ", meanfield->Vi);
    log_debug("W = %11.4k ", meanfield->w);
}

void meanfield_model_print_parameters(const meanfield_t *meanfield) {
    //log_debug("Ve = %11.4k ", meanfield->Ve);
    //log_debug("Vi = %11.4k ", meanfield->Vi);
    //log_debug("B = %11.4k ", neuron->B);
    //log_debug("C = %11.4k ", neuron->C);
    //log_debug("D = %11.4k ", neuron->D);

    //log_debug("I = %11.4k \n", neuron->I_offset);
}
