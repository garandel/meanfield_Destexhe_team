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
//! \brief Inlined neuron implementation following standard component model
#ifndef _MEANFIELD_IMPL_STD_H_
#define _MEANFIELD_IMPL_STD_H_


#include "meanfield_impl.h"

// Includes for model parts used in this implementation
#include <meanfield_and_syn/models/meanfield_model_impl.h>
#include <meanfield_and_syn/models/params_from_network.h>
#include <meanfield_and_syn/models/P_fit_polynomial.h>
#include <meanfield_and_syn/models/mathsbox.h>

//#include <meanfield/input_types/input_type.h>
//#include <meanfield/additional_inputs/additional_input.h>
//#include <meanfield/threshold_types/threshold_type.h>
#include <meanfield_and_syn/synapse_types/synapse_types.h>

// Further includes
#include <debug.h>
#include <bit_field.h>
#include <recording.h>
#include <stdfix-full-iso.h>




//! Indices for recording of words
enum word_recording_indices {
    //! V (somatic potential) recording index
    VE_RECORDING_INDEX = 0,
    VI_RECORDING_INDEX = 1,
    W_RECORDING_INDEX = 2,
    //! Number of recorded word-sized state variables
    N_RECORDED_VARS = 3
};

//! Indices for recording of bitfields
enum bitfield_recording_indices {
    //! Spike event recording index
    SPIKE_RECORDING_BITFIELD = 0,
    //! Number of recorded bitfields
    N_BITFIELD_VARS = 1
};


// This import depends on variables defined above
#include <meanfield/meanfield_recording.h>

//! Array of meanfield states -> will be change in future , the future is now!!
static meanfield_t *meanfield_array;

static ParamsFromNetwork_t *pNetwork_array;

static mathsbox_t *mathsbox_array;

static pFitPolynomial_t *Pfit_exc_array;
static pFitPolynomial_t *Pfit_inh_array;

/*
//! Threshold states array
static threshold_type_t *threshold_type_array;
*/
//! Global parameters for the neurons
static global_neuron_params_t *global_parameters;

//! The synapse shaping parameters
static synapse_param_t *neuron_synapse_shaping_params;

//! The number of steps to run per timestep
static uint n_steps_per_timestep;



#ifndef SOMETIMES_UNUSED
#define SOMETIMES_UNUSED __attribute__((unused))
#endif // !SOMETIMES_UNUSED

SOMETIMES_UNUSED // Marked unused as only used sometimes
//! \brief Initialise the particular implementation of the data
//! \param[in] n_neurons: The number of neurons
//! \return True if successful
static bool meanfield_impl_initialise(uint32_t n_meanfields) {
    // allocate DTCM for the global parameter details

    if (sizeof(global_neuron_params_t)) {
        global_parameters = spin1_malloc(sizeof(global_neuron_params_t));
        if (global_parameters == NULL) {
            log_error("Unable to allocate global neuron parameters"
                    "- Out of DTCM");
            return false;
        }
    }

    // Allocate DTCM for neuron array
    if (sizeof(meanfield_t)) {
        meanfield_array = spin1_malloc(n_meanfields * sizeof(meanfield_t));
        if (meanfield_array == NULL) {
            log_error("Unable to allocate meanfield array - Out of DTCM");
            return false;
        }
    }

    // Allocate DTCM for config array and copy block of data
    if (sizeof(ParamsFromNetwork_t)) {
        pNetwork_array = spin1_malloc(n_meanfields * sizeof(ParamsFromNetwork_t));
        if (pNetwork_array == NULL) {
            log_error("Unable to allocate config array - Out of DTCM");
            return false;
        }
    }
    
    // Allocate DTCM for P fit from polyomial for exc array and copy block of data
    if (sizeof(pFitPolynomial_t)) {
        Pfit_exc_array = spin1_malloc(n_meanfields * sizeof(pFitPolynomial_t));
        if (Pfit_exc_array == NULL) {
            log_error("Unable to allocate Pfit_exc_array - Out of DTCM");
            return false;
        }
    }
    
    // Allocate DTCM for P fit from polyomial for inh array and copy block of data
    if (sizeof(pFitPolynomial_t)) {
        Pfit_inh_array = spin1_malloc(n_meanfields * sizeof(pFitPolynomial_t));
        if (Pfit_inh_array == NULL) {
            log_error("Unable to allocate Pfit_inh_array - Out of DTCM");
            return false;
        }
    }

/*
    // Allocate DTCM for threshold type array and copy block of data
    if (sizeof(threshold_type_t)) {
        threshold_type_array =
                spin1_malloc(n_meanfields * sizeof(threshold_type_t));
        if (threshold_type_array == NULL) {
            log_error("Unable to allocate threshold type array - Out of DTCM");
            return false;
        }
    }
*/
    // Allocate DTCM for synapse shaping parameters
    if (sizeof(synapse_param_t)) {
        neuron_synapse_shaping_params =
                spin1_malloc(n_meanfields * sizeof(synapse_param_t));
        if (neuron_synapse_shaping_params == NULL) {
            log_error("Unable to allocate synapse parameters array"
                    " - Out of DTCM");
            return false;
        }
    }

    return true;
}

SOMETIMES_UNUSED // Marked unused as only used sometimes
// \brief Will be used for communication btw MFs
//! \brief Add inputs to the neuron
//! \param[in] synapse_type_index: the synapse type (e.g. exc. or inh.)
//! \param[in] neuron_index: the index of the neuron
//! \param[in] weights_this_timestep: weight inputs to be added
static void neuron_impl_add_inputs(
        index_t synapse_type_index, index_t neuron_index,
        input_t weights_this_timestep) {
    // simple wrapper to synapse type input function
    synapse_param_t *parameters =
            &neuron_synapse_shaping_params[neuron_index];
    synapse_types_add_neuron_input(synapse_type_index,
            parameters, weights_this_timestep);
}


//! \brief The number of _words_ required to hold an object of given size
//! \param[in] size: The size of object
//! \return Number of words needed to hold the object (not bytes!)
static uint32_t n_words_needed(size_t size) {
    return (size + (sizeof(uint32_t) - 1)) / sizeof(uint32_t);
}

SOMETIMES_UNUSED // Marked unused as only used sometimes
//! \brief Load in the neuron parameters
//! \param[in] address: SDRAM block to read parameters from
//! \param[in] next: Offset of next address in store
//! \param[in] n_neurons: number of neurons
static void neuron_impl_load_neuron_parameters(
        address_t address, uint32_t next, uint32_t n_meanfields) {
    log_info("reading parameters, next is %u, n_meanfields is %u ",
            next, n_meanfields);

    // Read the number of steps per timestep
    n_steps_per_timestep = address[next++];
    if (n_steps_per_timestep > 1) {
        log_debug("Looping over %u steps each timestep", n_steps_per_timestep);
    } else if (n_steps_per_timestep == 0) {
        log_error("bad number of steps per timestep: 0");
        rt_error(RTE_SWERR);
    }

    if (sizeof(global_neuron_params_t)) {
        log_debug("writing neuron global parameters");
        spin1_memcpy(global_parameters, &address[next],
                sizeof(global_neuron_params_t));
        next += n_words_needed(sizeof(global_neuron_params_t));
    }

    if (sizeof(meanfield_t)) {
        log_debug("reading neuron local parameters");
        spin1_memcpy(meanfield_array, &address[next],
                n_meanfields * sizeof(meanfield_t));
        next += n_words_needed(n_meanfields * sizeof(meanfield_t));
    }

    if (sizeof(ParamsFromNetwork_t)) {
        log_debug("reading config parameters");
        spin1_memcpy(pNetwork_array, &address[next],
                n_meanfields * sizeof(ParamsFromNetwork_t));
        next += n_words_needed(n_meanfields * sizeof(ParamsFromNetwork_t));
    }
    
    if (sizeof(pFitPolynomial_t)) {
        log_debug("reading pFitPolynomial exc parameters");
        spin1_memcpy(Pfit_exc_array, &address[next],
                n_meanfields * sizeof(pFitPolynomial_t));
        next += n_words_needed(n_meanfields * sizeof(pFitPolynomial_t));
    }
    
    if (sizeof(pFitPolynomial_t)) {
        log_debug("reading pFitPolynomial inh parameters");
        spin1_memcpy(Pfit_inh_array, &address[next],
                n_meanfields * sizeof(pFitPolynomial_t));
        next += n_words_needed(n_meanfields * sizeof(pFitPolynomial_t));
    }
    
    if (sizeof(mathsbox_t)) {
        log_debug("reading mathsbox parameters");
        spin1_memcpy(mathsbox_array, &address[next],
                n_meanfields * sizeof(mathsbox_t));
        next += n_words_needed(n_meanfields * sizeof(mathsbox_t));
    }
  
/*
    if (sizeof(threshold_type_t)) {
        log_debug("reading threshold type parameters");
        spin1_memcpy(threshold_type_array, &address[next],
                n_meanfields * sizeof(threshold_type_t));
        next += n_words_needed(n_meanfields * sizeof(threshold_type_t));
    }
*/
    if (sizeof(synapse_param_t)) {
        log_debug("reading synapse parameters");
        spin1_memcpy(neuron_synapse_shaping_params, &address[next],
                n_meanfields * sizeof(synapse_param_t));
        next += n_words_needed(n_meanfields * sizeof(synapse_param_t));
    }

    meanfield_model_set_global_neuron_params(global_parameters);

#if LOG_LEVEL >= LOG_DEBUG
    log_debug("-------------------------------------\n");
    for (index_t n = 0; n < n_meanfields; n++) {
        meanfield_model_print_parameters(&meanfield_array[n]);
    }
    log_debug("-------------------------------------\n");
#endif // LOG_LEVEL >= LOG_DEBUG
}


static union {
    int32_t as_int;
    REAL as_real;
} number;

static inline uint32_t abs_int(int32_t r){
    
    if (r < 0) {
	r = (r == INT32_MIN) ? INT32_MAX : -r;
    }
    
    return r;
}

//! Key from meanfield.c
extern uint32_t key;
extern uint32_t total_neighbour_exc;
extern uint32_t total_neighbour_inh;

SOMETIMES_UNUSED // Marked unused as only used sometimes
static void neuron_impl_do_timestep_update(
        uint32_t timer_count, uint32_t time, uint32_t n_neurons) {

    for (uint32_t meanfield_index = 0; meanfield_index < n_neurons; meanfield_index++) {
        // Get the neuron itself
        meanfield_t *this_meanfield = &meanfield_array[meanfield_index];

        // Get the Params from network and polynomial equation for this neuron
        ParamsFromNetwork_t *pNetwork_types = &pNetwork_array[meanfield_index];
        pFitPolynomial_t *Pfit_exc_types = &Pfit_exc_array[meanfield_index];
        pFitPolynomial_t *Pfit_inh_types = &Pfit_inh_array[meanfield_index];
        
        /*
        // Get threshold and additional input parameters for this neuron
        threshold_type_t *the_threshold_type = &threshold_type_array[meanfield_index];
        additional_input_t *additional_inputs =
                &additional_input_array[meanfield_index];
        */
        synapse_param_t *the_synapse_type =
                &neuron_synapse_shaping_params[meanfield_index];

        // Loop however many times requested; do this in reverse for efficiency,
        // and because the index doesn't actually matter
        for (uint32_t i_step = n_steps_per_timestep; i_step > 0; i_step--) {
            // Get the firing rate
            
            state_t firing_rate_Ve = meanfield_model_get_firing_rate_Ve(
                this_meanfield);
            state_t firing_rate_Vi = meanfield_model_get_firing_rate_Vi(
                this_meanfield);            
            // Get adaptation from excitator
            state_t adaptation_W = meanfield_model_get_adaptation_W(this_meanfield);
                        
            //***********************************************************************
            //!!Add to mimic an input from "synapses"!!!*
            //***********************************************************************
            

            
            uint32_t last_total_neighbour_exc = 
                pNetwork_types->exc_neighbour_contribution;
            uint32_t last_total_neighbour_inh = 
                pNetwork_types->inh_neighbour_contribution;

            uint32_t diff_neighbour_exc = 
                (total_neighbour_exc - last_total_neighbour_exc);
            uint32_t diff_neighbour_inh = 
                (total_neighbour_inh - last_total_neighbour_inh);

            pNetwork_types->exc_neighbour_contribution = total_neighbour_exc;
            pNetwork_types->inh_neighbour_contribution = total_neighbour_inh;
  
            
            number.as_int = diff_neighbour_exc; 
            REAL total_neighbour_exc_real = number.as_real;
            
            number.as_int = diff_neighbour_inh; 
            REAL total_neighbour_inh_real = number.as_real;
            
            input_t external_bias = 
                pNetwork_types->afferent_exc_fraction;
            
            
            
            the_synapse_type->exc.synaptic_input_value = total_neighbour_exc_real; 
            //firing_rate_Ve + total_neighbour_exc_real;
            //total_neighbour_exc_real;
            the_synapse_type->inh.synaptic_input_value = total_neighbour_inh_real;
            //firing_rate_Vi + total_neighbour_inh_real;
            // total_neighbour_inh_real;
            
            // Get the exc and inh values from the synapses
            input_t exc_values[NUM_EXCITATORY_RECEPTORS];
            input_t *exc_syn_values =
                    synapse_types_get_excitatory_input(exc_values, the_synapse_type);
            input_t inh_values[NUM_INHIBITORY_RECEPTORS];
            input_t *inh_syn_values =
                    synapse_types_get_inhibitory_input(inh_values, the_synapse_type);
            
            // Operation post synaptic
            
            // Sum g_syn contributions from all receptors for recording
            /*
            REAL total_exc = 0;
            REAL total_inh = 0;
            
            for (int i = 0; i < NUM_EXCITATORY_RECEPTORS; i++) {
                total_exc += exc_syn_values[i];
                log_info("exc_syn_values=%8.6k",exc_syn_values[i]);
            }
            for (int i = 0; i < NUM_INHIBITORY_RECEPTORS; i++) {
                total_inh += inh_syn_values[i];
                log_info("inh_syn_values=%8.6k",inh_syn_values[i]);
            }
            */

            // Do recording if on the first step 
            if (i_step == n_steps_per_timestep) {
                neuron_recording_record_accum(
                        VE_RECORDING_INDEX, meanfield_index, firing_rate_Ve);
                neuron_recording_record_accum(
                        VI_RECORDING_INDEX, meanfield_index, firing_rate_Vi);
                neuron_recording_record_accum(
                        W_RECORDING_INDEX, meanfield_index, adaptation_W);
            }
            
            number.as_real = firing_rate_Ve;// *exc_syn_values;//
            int32_t firing_rate_exc_int = number.as_int; 
            
            number.as_real = firing_rate_Vi;// *exc_syn_values;//
            int32_t firing_rate_inh_int = number.as_int; 
            
            
            //! FUTUR : ad weight mult by firing rate
            
            
            log_info("firing_int_exc = %d",firing_rate_exc_int);
            
            //log_info("total_neighbour_exc = %d", total_neighbour_exc);
            //log_info("last_total_neighbour_exc = %d", last_total_neighbour_exc);
            //log_info("diff_neighbour_exc = %d", diff_neighbour_exc);
            //log_info("total_neighbour_inh = %d", diff_neightbour_inh);
                                    
            #define MASK_EXC ((uint32_t) 0x0)
            #define MASK_INH ((uint32_t) 0x1)
            //bool inh = 1;
            uint32_t concat_exc = abs_int((firing_rate_exc_int<<1) | MASK_EXC);
            uint32_t concat_inh = abs_int((firing_rate_inh_int<<1) | MASK_INH);

            spin1_send_mc_packet(key, concat_exc, WITH_PAYLOAD);
            spin1_send_mc_packet(key, concat_inh, WITH_PAYLOAD);
            
            // Shape the existing input according to the included rule
            //synapse_types_shape_input(the_synapse_type);
            
            // update neuron parameters
            
            meanfield_model_state_update(this_meanfield,
                                          pNetwork_types,
                                          Pfit_exc_types,
                                          Pfit_inh_types,
                                          external_bias,
                                          NUM_EXCITATORY_RECEPTORS,
                                          exc_syn_values,
                                          NUM_INHIBITORY_RECEPTORS,
                                          inh_syn_values);
            
        }
    

#if LOG_LEVEL >= LOG_DEBUG
        meanfield_model_print_state_variables(this_meanfield);
#endif // LOG_LEVEL >= LOG_DEBUG
    }
}

SOMETIMES_UNUSED // Marked unused as only used sometimes
//! \brief Stores neuron parameters back into SDRAM
//! \param[out] address: the address in SDRAM to start the store
//! \param[in] next: Offset of next address in store
//! \param[in] n_neurons: number of neurons
static void neuron_impl_store_neuron_parameters(
        address_t address, uint32_t next, uint32_t n_meanfields) {
    log_debug("writing parameters");

    // Skip over the steps per timestep
    next += 1;

    if (sizeof(global_neuron_params_t)) {
        log_debug("writing neuron global parameters");
        spin1_memcpy(&address[next], global_parameters,
                sizeof(global_neuron_params_t));
        next += n_words_needed(sizeof(global_neuron_params_t));
    }

    if (sizeof(meanfield_t)) {
        log_debug("writing neuron local parameters");
        spin1_memcpy(&address[next], meanfield_array,
                n_meanfields * sizeof(meanfield_t));
        next += n_words_needed(n_meanfields * sizeof(meanfield_t));
    }
    
    if (sizeof(synapse_param_t)) {
        log_debug("writing synapse parameters");
        spin1_memcpy(&address[next], neuron_synapse_shaping_params,
                n_meanfields * sizeof(synapse_param_t));
        next += n_words_needed(n_meanfields * sizeof(synapse_param_t));
    }

    
    if (sizeof(ParamsFromNetwork_t)) {
        log_debug("writing input type parameters");
        spin1_memcpy(&address[next], pNetwork_array,
                n_meanfields * sizeof(ParamsFromNetwork_t));
        next += n_words_needed(n_meanfields * sizeof(ParamsFromNetwork_t));
    }
    

}


//#if LOG_LEVEL >= LOG_DEBUG
//! \brief Print the inputs to the neurons
//! \param[in] n_neurons: The number of neurons

#if LOG_LEVEL >= LOG_DEBUG
void neuron_impl_print_inputs(uint32_t n_meanfields) {
    log_debug("-------------------------------------\n");
    for (index_t i = 0; i < n_meanfields; i++) {
        meanfield_t *meanfield = &meanfield_array[i];
        log_debug("inputs: %k %k", meanfield->a, meanfield->b);
    }
    log_debug("-------------------------------------\n");
}


//! \brief Print the synapse parameters of the neurons
//! \param[in] n_neurons: The number of neurons
void neuron_impl_print_synapse_parameters(uint32_t n_neurons) {
    log_debug("-------------------------------------\n");
    for (index_t n = 0; n < n_neurons; n++) {
        synapse_types_print_parameters(&neuron_synapse_shaping_params[n]);
    }
    log_debug("-------------------------------------\n");
}

//! \brief Get the synapse type character for a synapse type
//! \param[in] synapse_type: The synapse type
//! \return The descriptor character (sometimes two characters)
const char *neuron_impl_get_synapse_type_char(uint32_t synapse_type) {
    return synapse_types_get_type_char(synapse_type);
}
#endif // LOG_LEVEL >= LOG_DEBUG

#endif // _NEURON_IMPL_STANDARD_H_
