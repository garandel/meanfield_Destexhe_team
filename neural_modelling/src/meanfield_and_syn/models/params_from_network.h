#ifndef _PARAMS_FROM_NETWORK_H_
#define _PARAMS_FROM_NETWORK_H_

#include "../../meanfield_and_syn/models/meanfield_model.h"

typedef struct ParamsFromNetwork_t {
    // nominally 'fixed' parameters
    REAL pconnec;
    REAL Qe;
    REAL Qi;
    REAL Te;
    REAL Ti;
    REAL Ee;
    REAL Ei;
    REAL Ntot;
    REAL gei;
    REAL ext_drive;
    REAL afferent_exc_fraction;
    
    REAL Gl;
    REAL Cm;
    REAL El_exc;
    REAL El_inh;
    
    REAL muV;
    REAL muV0;
    
    REAL iDmuV0;//inverse one
    
    REAL sV;
    REAL sV0;
    REAL iDsV0;//inverse one
        
    REAL muGn;
    
    REAL TvN;
    REAL TvN0;
    REAL iDTvN0; //inverse one
    
    REAL Vthre;
    
    REAL Fout_th;
    uint32_t exc_neighbour_contribution;
    uint32_t inh_neighbour_contribution;

} ParamsFromNetwork_t;


#endif