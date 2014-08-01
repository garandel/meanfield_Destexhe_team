#!/usr/bin/env python
import unittest
import spynnaker.pyNN as pynn
from pprint import pprint as pp
from spynnaker.pyNN.exceptions import ConfigurationException
#Setup
pynn.setup(timestep=1, min_delay=1, max_delay=10.0)

cell_params_lif = {'cm'  : 0.25,
             'i_offset'  : 0.0,
             'tau_m'     : 20.0,
             'tau_refrac': 2.0,
             'tau_syn_E' : 5.0,
             'tau_syn_I' : 5.0,
             'v_reset'   : -70.0,
             'v_rest'    : -65.0,
             'v_thresh'  : -50.0
             }
spike_array = {'spike_times':[0]}
#/Setup


class TestingFixedProbabilityConnector(unittest.TestCase):
    def test_generate_synapse_list(self):
        number_of_neurons = 5
        first_population=pynn.Population(number_of_neurons,pynn.IF_curr_exp,cell_params_lif,label="One pop")
        second_population=pynn.Population(number_of_neurons,pynn.IF_curr_exp,cell_params_lif,label= "Second pop")
        weight = 2
        delay = 1
        synapse_type = first_population._vertex.get_synapse_id('excitatory')
        connection = pynn.FixedProbabilityConnector(0.5,weight,delay)
        synaptic_list = connection.generate_synapse_list(first_population._vertex,first_population._vertex,1,synapse_type)
        pp(synaptic_list.get_rows())

    def test_generate_synapse_list_probability_zero_percent(self):
        number_of_neurons = 5
        first_population=pynn.Population(number_of_neurons,pynn.IF_curr_exp,cell_params_lif,label="One pop")
        second_population=pynn.Population(number_of_neurons,pynn.IF_curr_exp,cell_params_lif,label= "Second pop")
        weight = 2
        delay = 1
        synapse_type = first_population._vertex.get_synapse_id('excitatory')
        connection = pynn.FixedProbabilityConnector(0,weight,delay)
        synaptic_list = connection.generate_synapse_list(first_population._vertex,first_population._vertex,1,synapse_type)
        pp(synaptic_list.get_rows())

    def test_generate_synapse_list_probability_100_percent(self):
        number_of_neurons = 5
        first_population=pynn.Population(number_of_neurons,pynn.IF_curr_exp,cell_params_lif,label="One pop")
        second_population=pynn.Population(number_of_neurons,pynn.IF_curr_exp,cell_params_lif,label= "Second pop")
        weight = 2
        delay = 1
        synapse_type = first_population._vertex.get_synapse_id('excitatory')
        connection = pynn.FixedProbabilityConnector(1,weight,delay)
        #def generate_synapse_list(self, prevertex, postvertex, delay_scale, synapse_type)
        synaptic_list = connection.generate_synapse_list(first_population._vertex,first_population._vertex,1,synapse_type)
        pp(synaptic_list.get_rows())
        self.assertEqual(synaptic_list.get_max_weight(),weight)
        self.assertEqual(synaptic_list.get_min_weight(),weight)
        self.assertEqual(synaptic_list.get_n_rows(),number_of_neurons)
        self.assertEqual(synaptic_list.get_min_max_delay(),(delay,delay))

    def test_generate_synapse_list_probability_200_percent(self):
        with self.assertRaises(ConfigurationException):
            number_of_neurons = 5
            first_population=pynn.Population(number_of_neurons,pynn.IF_curr_exp,cell_params_lif,label="One pop")
            second_population=pynn.Population(number_of_neurons,pynn.IF_curr_exp,cell_params_lif,label= "Second pop")
            weight = 2
            delay = 1
            synapse_type = first_population._vertex.get_synapse_id('excitatory')
            connection = pynn.FixedProbabilityConnector(2,weight,delay)

    def test_synapse_list_generation_for_negative_sized_populations(self):
        with self.assertRaises(ConfigurationException):
            number_of_neurons = 5
            first_population=pynn.Population(number_of_neurons,pynn.IF_curr_exp,cell_params_lif,label="One pop")
            second_population=pynn.Population(number_of_neurons,pynn.IF_curr_exp,cell_params_lif,label= "Second pop")
            weight = 2
            delay = 1
            synapse_type = first_population._vertex.get_synapse_id('excitatory')
            connection = pynn.FixedProbabilityConnector(-0.5,weight,delay)

    def test_synapse_list_generation_for_different_sized_populations(self):
        number_of_neurons = 10
        first_population=pynn.Population(number_of_neurons,pynn.IF_curr_exp,cell_params_lif,label="One pop")
        second_population=pynn.Population(number_of_neurons + 5,pynn.IF_curr_exp,cell_params_lif,label= "Second pop")
        weight = 2
        delay = 1
        connection = pynn.FixedProbabilityConnector(0.1,weight,delay)
        synaptic_list = connection.generate_synapse_list(first_population._vertex,second_population._vertex,1,0)
        pp(synaptic_list.get_rows())


    def test_allow_self_connections(self):
        number_of_neurons = 5
        first_population=pynn.Population(number_of_neurons,pynn.IF_curr_exp,cell_params_lif,label="One pop")
        second_population=pynn.Population(number_of_neurons,pynn.IF_curr_exp,cell_params_lif,label= "Second pop")
        weight = 2
        delay = 1
        synapse_type = first_population._vertex.get_synapse_id('excitatory')
        connection = pynn.FixedProbabilityConnector(1,weight,delay,allow_self_connections = False)
        synaptic_list = connection.generate_synapse_list(first_population._vertex,first_population._vertex,1,synapse_type)
        pp(synaptic_list.get_rows())

if __name__=="__main__":
    unittest.main()