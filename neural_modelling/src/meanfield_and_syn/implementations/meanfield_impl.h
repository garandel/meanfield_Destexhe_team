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

//! \dir
//! \brief Neuron implementations
//! \file
//! \brief General API of a neuron implementation
#ifndef _MEANFIELD_IMPL_H_
#define _MEANFIELD_IMPL_H_

#include <common/neuron-typedefs.h>
#include <meanfield_and_syn/send_spike.h>

//! \brief Initialise the particular implementation of the data
//! \param[in] n_neurons: The number of neurons
//! \return True if successful
static bool meanfield_impl_initialise(uint32_t n_meanfields);

//! \brief Add inputs to the neuron
//! \param[in] synapse_type_index: the synapse type (e.g. exc. or inh.)
//! \param[in] neuron_index: the index of the neuron
//! \param[in] weights_this_timestep: weight inputs to be added
static void neuron_impl_add_inputs(
        index_t synapse_type_index, index_t meanfield_index,
        input_t weights_this_timestep);

//! \brief Load in the neuron parameters
//! \param[in] address: SDRAM block to read parameters from
//! \param[in] next: Offset of next address in store
//! \param[in] n_neurons: The number of neurons
static void neuron_impl_load_neuron_parameters(
        address_t address, uint32_t next, uint32_t n_meanfields);

//! \brief Do the timestep update for the particular implementation
//! \param[in] neuron_index: The index of the neuron to update
//! \param[in] external_bias: External input to be applied to the neuron
//! \return True if a spike has occurred
static void neuron_impl_do_timestep_update(
        uint32_t timer_count, uint32_t time, uint32_t n_neurons);

//! \brief Stores neuron parameters back into SDRAM
//! \param[out] address: the address in SDRAM to start the store
//! \param[in] next: Offset of next address in store
//! \param[in] n_neurons: The number of neurons
static void neuron_impl_store_neuron_parameters(
        address_t address, uint32_t next, uint32_t n_meanfields);

#if LOG_LEVEL >= LOG_DEBUG
//! \brief Print the inputs to the neurons
//! \param[in] n_neurons: The number of neurons
void neuron_impl_print_inputs(uint32_t n_meanfields);

//! \brief Print the synapse parameters of the neurons
//! \param[in] n_neurons: The number of neurons
void neuron_impl_print_synapse_parameters(uint32_t n_meanfields);

//! \brief Get the synapse type character for a synapse type
//! \param[in] synapse_type: The synapse type
//! \return The descriptor character (sometimes two characters)
const char *neuron_impl_get_synapse_type_char(uint32_t synapse_type);
#endif // LOG_LEVEL >= LOG_DEBUG

#endif // _NEURON_IMPL_H_
