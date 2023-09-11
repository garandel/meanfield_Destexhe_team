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
#include "../../meanfield_and_syn/models/meanfield_model_impl.h"

#include <debug.h>
//#include <math.h>

//#include "hls_math.h"
//#inlcude "hls_erf.h"
//#include "ap_fixed.h"

//#include <stdfix-exp.h>
#include "stdfix-erfc.h"
//#include <polynomial.h>

#include "../../../src/common/maths-util.h"
//#include "../../common/maths-util.h" // i.o to use SQRT(x) and SQR(a)

#include "../../meanfield_and_syn/models/params_from_network.h"
#include "../../meanfield_and_syn/models/P_fit_polynomial.h"


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

static inline REAL square_root_of(REAL number)
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

//static int poly[3] = {1,5,3};
/*    
static union {
    s1615 as_s1615;
    input_t as_real;
} number;


static inline REAL erfc_test(REAL x)
{
    
    number.as_real = x;//*exc_syn_values;//
    s1615 x_s1615 = number.as_s1615; 
    
    number.as_s1615 = expk(x_s1615);
    input_t result = number.as_real;
       
    
    //REAL result = erfc(x); //region `ITCM' overflowed by 3448 bytes
        
    log_info("IT'S THE WRONG ONE JUST FOR TEST => NEED ERFC NOT ONLY EXPK ");
    //return __horner_int_b(poly,x,2);
    //return expk(x_s1615);
    return result;
    //return expf(x);
    //return x+1.*20;
        
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
  
    pNetwork->Vthre = Vthre;

    }
    

void get_fluct_regime_varsup(REAL Ve, REAL Vi, REAL W, 
                             ParamsFromNetwork_t *restrict pNetwork)
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
    REAL El_exc = pNetwork->El_exc;
    //REAL El_inh = pNetwork->El_inh;
    REAL Ei = pNetwork->Ei;
    REAL Ee = pNetwork->Ee;
    REAL Cm = pNetwork->Cm;
           
    REAL Fe;
    Fe = Ve * (1-gei)*pconnec*Ntot; // default is 1 !!
    REAL Fi;
    Fi = Vi * gei*pconnec*Ntot;
    
    REAL muGe;
    muGe = Qe*Te*Fe;
    REAL muGi;
    muGi = Qi*Ti*Fi;

    REAL muG;
    muG = Gl + muGe + muGi;
    
    //Average population voltage
    REAL muV  = (muGe*Ee + muGi*Ei + Gl*El_exc - W)/muG;
    
    pNetwork->muV = muV ;

    //REAL muGn = muG/Gl;
    //pNetwork->muGn = muGn; //to uncomment in order to stock it
    REAL Tm = Cm/muG;
    REAL Ue = Qe*(Ee-muV)/muG;
    REAL Ui = Qi*(Ei-muV)/muG;
    
    
    REAL Tv_num_e = Fe*(Ue*Te)*(Ue*Te) ;
    REAL Tv_num_i = Fi*(Ti*Ui)*(Ti*Ui) ;
    REAL Tv_denom_e = Fe*(Ue*Te)*(Ue*Te)/(Te+Tm);
    REAL Tv_denom_i = Fi*(Ti*Ui)*(Ti*Ui)/(Ti+Tm);
    

    REAL Tv = (Tv_num_e + Tv_num_i) / (Tv_denom_e + Tv_denom_i);

    pNetwork->TvN = Tv*Gl/Cm; // TvN is adimensional so usefull var
    
    REAL sV_sqr = REAL_CONST(0.50000)*((Tv_denom_e + Tv_denom_i));
    
    pNetwork->sV =  square_root_of(sV_sqr);  // sV_sqr; //  sqrtk(sV_sqr);//  kbits(sV_sqr);

}


void TF(REAL Ve, REAL Vi, REAL W,
        ParamsFromNetwork_t *restrict pNetwork,
        pFitPolynomial_t *restrict Pfit){

/********************************************************************
 *   State-variables are directly connected to the struct           *
 *   parameters are put in local in order to make the code clear.   *
 *                                                                  *
 ********************************************************************/
    
    
    if (pNetwork->Fout_th != ZERO){
        pNetwork->Fout_th = ACS_DBL_TINY;
    }
    if (pNetwork->muV != ZERO){
        pNetwork->muV = ACS_DBL_TINY;
    }

    if (Ve < ACS_DBL_TINY){
        Ve += ACS_DBL_TINY;
    }
    if (Vi < ACS_DBL_TINY){
        Vi += ACS_DBL_TINY;
    }

    get_fluct_regime_varsup(Ve, Vi, W, pNetwork);
    
    if (pNetwork->Vthre != ZERO){
        pNetwork->Vthre = ACS_DBL_TINY;
    }
    
    threshold_func(pNetwork, Pfit);

    
    if (pNetwork->sV<ACS_DBL_TINY){
        pNetwork->sV += ACS_DBL_TINY;
    }
    
    REAL argument = (pNetwork->Vthre - \
                     pNetwork->muV)/(REAL_CONST(1.4142137)*pNetwork->sV); 
    
    REAL error_func = argument;//erfc(argument); //EXP(argument);// with EXP is compiling 
    
    
    
    log_info("argument = %5.5k and error_func = %5.5k", argument, error_func);
    
    REAL Gl = pNetwork->Gl;
    REAL Cm = pNetwork->Cm;
    
    
    log_info("Cm = %5.5k AND TvN = %5.5k", Cm, pNetwork->TvN);
    REAL one_over_CmxTvN = (0,02001441); //1/(Cm*pNetwork->TvN);//
    
    pNetwork->Fout_th = error_func * (HALF*Gl) * one_over_CmxTvN;// /(Cm*pNetwork->TvN) ;
    log_info("Fout_th = %5.5k \n", pNetwork->Fout_th);
    
    if (pNetwork->Fout_th < ACS_DBL_TINY){
        pNetwork->Fout_th += ACS_DBL_TINY;
    }
    
}


void RK2_midpoint_MF(REAL h, meanfield_t *meanfield,
                     ParamsFromNetwork_t *restrict pNetwork,
                     pFitPolynomial_t *restrict Pfit_exc,
                     pFitPolynomial_t *restrict Pfit_inh,
                     REAL total_exc, REAL total_inh,
                     REAL input_this_timestep){
    
    /* 
     * will add exc_aff=0, inh_aff=0, pure_exc_aff=0
     */
    REAL a_exc = meanfield->a_exc;
    REAL b_exc = meanfield->b_exc;
    REAL El_exc = pNetwork->El_exc;
    REAL tauw_exc = meanfield->tauw_exc;
    REAL ext_drive = pNetwork->ext_drive;
    
    //log_info("input_this_timestep = %11.4k", input_this_timestep);

    REAL lastVe = meanfield->Ve;
    REAL lastVepExtD = lastVe + input_this_timestep + ext_drive;// + total_exc;//
    
    REAL lastVi = meanfield->Vi;// + total_inh;
    REAL lastWe = meanfield->w_exc;
    REAL lastWi = meanfield->w_inh;
    
    //REAL lastW_exc = meanfield->w_exc;
    //REAL lastW_inh = lastW_exc - b*lastVe;
    
    
    REAL T_inv = meanfield->Timescale_inv;
    
/***********************************************************************
 *   EULER Explicite method
 *   It's very instable if for now h<0.2 for 20ms
 *   
 *   NEED to give also the error of the method here :
 *   0.5*h^2*u''(t_n) + o(h^2)
 * 
 * ERROR : overflowed for 936 ITCM when add adaptation W in 
 *         get_fluct_regime_varsup() and TF()
 *         AND reduce h don't give big changes
 *      IF need memory will investigate but not now in 23 feb 2022
 ***********************************************************************/
    /*
    h=h*0.001;
    REAL k1_exc = (lastTF_exc - lastVe)*T_inv;
    meanfield->Ve =  lastVe + h*k1_exc ;

    REAL k1_inh = (lastTF_inh - lastVi)*T_inv;
    meanfield->Vi = lastVi + h*k1_inh ;
    
    REAL k1_W = -lastW/tauw + b * lastVe;
    meanfield->w = lastW + h*k1_W;
    */
    
/************************************************
 *  RUNGE-KUTTA 2nd order Midpoint Try precision*
 ***********************************************/
    TF(lastVepExtD, lastVi, lastWe, pNetwork, Pfit_exc);    
    REAL lastmuV = pNetwork->muV;
    REAL lastTF_exc_1 = pNetwork->Fout_th;
       
    TF(lastVepExtD, lastVi, lastWi, pNetwork, Pfit_inh);
    REAL lastTF_inh_1 = pNetwork->Fout_th;
    
    h=h*0.001;
        
    REAL alpha_exc_1 = T_inv*(lastTF_exc_1 - lastVepExtD);
    REAL lastVe_n2 = lastVepExtD + REAL_HALF(h*alpha_exc_1);
    
    REAL alpha_inh_1 = T_inv*(lastTF_inh_1 - lastVi);
    REAL lastVi_n2 = lastVi + REAL_HALF(h*alpha_inh_1);
    
    TF(lastVe_n2, lastVi_n2, lastWe, pNetwork, Pfit_exc);
    REAL TF_exc_2 = pNetwork->Fout_th;
    
    TF(lastVe_n2, lastVi_n2, lastWi, pNetwork, Pfit_inh);
    REAL TF_inh_2 = pNetwork->Fout_th;
    
    REAL alpha_exc_2 = T_inv*(TF_exc_2 - lastVe_n2);
    REAL alpha_inh_2 = T_inv*(TF_inh_2 - lastVi_n2);
    
    meanfield->Ve += h*alpha_exc_2;
    meanfield->Vi += h*alpha_inh_2;
    
    // Control if output are realistic
    //log_info("%6.1k  %4.9k  %4.9k", h, meanfield->Ve, meanfield->Vi);
    
    REAL k1_We = -lastWe/tauw_exc + b_exc * lastVepExtD + a_exc*(lastmuV-El_exc);
    REAL alpha_we = lastWe + h*k1_We;
    REAL k2_We = -alpha_we/tauw_exc + b_exc * lastVepExtD + a_exc*(lastmuV-El_exc);
 
    meanfield->w_exc += REAL_HALF(h*(k1_We+k2_We));

}

void meanfield_model_set_global_neuron_params(
        const global_neuron_params_t *params) {
    global_params = params;
}

/************************* IDEA ***************************************************
 * perhaps when we will do more than one MF we could uses "num_excitatory_inputs" *
 * like the number of ex MF and in MF?                                            *
 * and maybe is there some contamanation from the neightbourest neighbour MF!     *
 **********************************************************************************/
state_t meanfield_model_state_update(
    meanfield_t *restrict meanfield,
    ParamsFromNetwork_t *restrict pNetwork,
    pFitPolynomial_t *restrict Pfit_exc,
    pFitPolynomial_t *restrict Pfit_inh,
    input_t external_bias,
    uint16_t num_excitatory_inputs, const input_t *exc_input,
    uint16_t num_inhibitory_inputs, const input_t *inh_input) {
    
    REAL total_exc = 0;
    REAL total_inh = 0;

    for (int i =0; i<num_excitatory_inputs; i++) {
        total_exc += exc_input[i];
        //log_info("exc_inputs = %6.6k",exc_input[i]);
    }
    for (int i =0; i<num_inhibitory_inputs; i++) {
        total_inh += inh_input[i];
    }
    
    log_info("total_exc=%8.6k",total_exc);
    log_info("total_inh=%8.6k",total_inh);
    

    input_t input_this_timestep = external_bias;//total_exc + total_inh + external_bias;// + neuron->I_offset;
    //log_info("input_this_timestep = %11.4k", input_this_timestep);

    // the best AR update so far
    RK2_midpoint_MF(meanfield->this_h,
                    meanfield,
                    pNetwork,
                    Pfit_exc,
                    Pfit_inh,
                    total_exc,
                    total_inh,
                    input_this_timestep);
                    
    //meanfield->this_h = global_params->machine_timestep_ms;
    
    //what is the best output for this function? Output of this function is used normally for thershold
    return meanfield->Ve;
}


/*
void neuron_model_has_spiked() { 
    log_debug("in neuron_model_has_spiked, time is ",
              global_params->machine_timestep_ms);
}
*/

state_t meanfield_model_get_firing_rate_Ve(const meanfield_t *meanfield) {
    return meanfield->Ve;
}

state_t meanfield_model_get_firing_rate_Vi(const meanfield_t *meanfield) {
    return meanfield->Vi;
}

state_t meanfield_model_get_adaptation_W(const meanfield_t *meanfield){
    return meanfield->w_exc;
}

/*
void meanfield_model_print_state_variables(const meanfield_t *meanfield) {
    log_debug("Ve = %11.4k ", meanfield->Ve);
    log_debug("Vi = %11.4k ", meanfield->Vi);
    log_debug("W_exc = %11.4k ", meanfield->w_exc);
}
*/
/*
void meanfield_model_print_parameters() { //const meanfield_t *meanfield
    //log_debug("Ve = %11.4k ", meanfield->Ve);
    //log_debug("Vi = %11.4k ", meanfield->Vi);
    //log_debug("B = %11.4k ", neuron->B);
    //log_debug("C = %11.4k ", neuron->C);
    //log_debug("D = %11.4k ", neuron->D);

    //log_debug("I = %11.4k \n", neuron->I_offset);
}
*/