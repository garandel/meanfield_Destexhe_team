# Copyright (c) 2017-2019 The University of Manchester
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from spinn_utilities.overrides import overrides
from data_specification.enums import DataType
from spinn_front_end_common.utilities.constants import (
    MICRO_TO_MILLISECOND_CONVERSION)
from .abstract_neuron_model import AbstractNeuronModel
from spynnaker.pyNN.models.neuron.implementations import (
    AbstractStandardNeuronComponent)

###--Meanfield Params--###
#NBR = "nbr"
A_EXC = "a_exc"
B_EXC = "b_exc"
TAUW_EXC = "tauw_exc"
A_INH = "a_inh"
B_INH = "b_inh"
TAUW_INH = "tauw_inh"
TREFRAC = "Trefrac"
VRESET = "Vreset"
DELTA_V_EXC = "delta_v_exc"
DELTA_V_INH = "delta_v_inh"
AMPNOISE = "ampnoise"
TIMESCALE_INV = "Timescale_inv"
VE = "Ve"
VI = "Vi"
W_EXC = "w_exc"
W_INH = "w_inh"

UNITS = {
    ###--Meanfield--###
    #NBR: "",
    A_EXC: "nS",
    B_EXC: "nS",
    TAUW_EXC: "ms",
    A_INH: "nS",
    B_INH: "nS",
    TAUW_INH: "ms",
    TREFRAC: "ms",
    VRESET: "mV",
    DELTA_V_EXC: "mV",
    DELTA_V_INH: "mV",
    AMPNOISE: "Hz",
    TIMESCALE_INV: "Hz",
    VE: "Hz",
    VI: "Hz",
    W_EXC: "pA",
    W_INH: "pA",
}


class MeanfieldOfAdexNetwork(AbstractNeuronModel):
    """ Model of meanfield due to A.Destehexe et al
    """
    __slots__ = [
        "_a_exc", "_b_exc", "_tauw_exc",
        "_a_inh", "_b_inh", "_tauw_inh",
        "_Trefrac", "_Vreset", "_delta_v_exc", "_delta_v_inh",
        "_ampnoise", "_Timescale_inv",
        "_Ve_init", "_Vi_init",
        "_w_exc_init", "_w_inh_init"
    ]

    def __init__(self,
                 a_exc, b_exc, tauw_exc,
                 a_inh, b_inh, tauw_inh,
                 Trefrac, Vreset, delta_v_exc, delta_v_inh,
                 ampnoise, Timescale_inv,
                 Ve_init, Vi_init,
                 w_exc_init, w_inh_init):
        """
        :param a: :math:`a`
        :type a: float, iterable(float), ~pyNN.random.RandomDistribution or
            (mapping) function
        
        """
        super().__init__(
            [DataType.S1615, #a_exc
            DataType.S1615, #b_exc
            DataType.S1615, #tauw_exc
            DataType.S1615, #a_inh
            DataType.S1615, #b_inh
            DataType.S1615, #tauw_inh
            DataType.S1615, #Trefrac
            DataType.S1615, #Vreset
            DataType.S1615, #delta_v_exc
            DataType.S1615, #delta_v_inh
            DataType.S1615, #ampnoise
            DataType.S1615, #Timescale_inv
            DataType.S1615, #Ve
            DataType.S1615, #Vi
            DataType.S1615, #W_exc
            DataType.S1615, #W_inh
            DataType.S1615],  # this_h (= machine_time_step)
            [DataType.S1615])  # machine_time_step
        self._a_exc = a_exc
        self._b_exc = b_exc
        self._tauw_exc = tauw_exc
        self._a_inh = a_inh
        self._b_inh = b_inh
        self._tauw_inh = tauw_inh
        self._Trefrac = Trefrac
        self._Vreset =Vreset
        self._delta_v_exc = delta_v_exc
        self._delta_v_inh = delta_v_inh
        self._ampnoise = ampnoise
        self._Timescale_inv = Timescale_inv
        self._Ve_init = Ve_init
        self._Vi_init = Vi_init
        self._w_exc_init = w_exc_init
        self._w_inh_init = w_inh_init

    @overrides(AbstractStandardNeuronComponent.get_n_cpu_cycles)
    def get_n_cpu_cycles(self, n_neurons):
        # A bit of a guess
        return 150 * n_neurons

    @overrides(AbstractStandardNeuronComponent.add_parameters)
    def add_parameters(self, parameters):
        ###--neuron--###
        #parameters[NBR] = self._nbr
        parameters[A_EXC] = self._a_exc
        parameters[B_EXC] = self._b_exc
        parameters[TAUW_EXC] = self._tauw_exc
        parameters[A_INH] = self._a_inh
        parameters[B_INH] = self._b_inh
        parameters[TAUW_INH] = self._tauw_inh
        parameters[TREFRAC] = self._Trefrac
        parameters[VRESET] = self._Vreset
        parameters[DELTA_V_EXC] = self._delta_v_exc
        parameters[DELTA_V_INH] = self._delta_v_inh
        parameters[AMPNOISE] = self._ampnoise
        parameters[TIMESCALE_INV] = self._Timescale_inv

    @overrides(AbstractStandardNeuronComponent.add_state_variables)
    def add_state_variables(self, state_variables):
        state_variables[VE] = self._Ve_init
        state_variables[VI] = self._Vi_init
        state_variables[W_EXC] = self._w_exc_init
        state_variables[W_INH] = self._w_inh_init
        
    @overrides(AbstractStandardNeuronComponent.get_units)
    def get_units(self, variable):
        return UNITS[variable]

    @overrides(AbstractStandardNeuronComponent.has_variable)
    def has_variable(self, variable):
        return variable in UNITS

    @overrides(AbstractNeuronModel.get_global_values)
    def get_global_values(self, ts):
        # pylint: disable=arguments-differ
        return [float(ts) / MICRO_TO_MILLISECOND_CONVERSION]

    @overrides(AbstractStandardNeuronComponent.get_values)
    def get_values(self, parameters, state_variables, vertex_slice, ts):
        """
        :param ts: machine time step
        """
        # pylint: disable=arguments-differ

        # Add the rest of the data
        return [
            parameters[A_EXC],parameters[B_EXC],parameters[TAUW_EXC],
            parameters[A_INH],parameters[B_INH],parameters[TAUW_INH],
            parameters[TREFRAC],parameters[VRESET],
            parameters[DELTA_V_EXC],parameters[DELTA_V_INH],
            parameters[AMPNOISE], parameters[TIMESCALE_INV],
            state_variables[VE],
            state_variables[VI],
            state_variables[W_EXC],
            state_variables[W_INH],
            float(ts) / MICRO_TO_MILLISECOND_CONVERSION
        ]

    @overrides(AbstractStandardNeuronComponent.update_values)
    def update_values(self, values, parameters, state_variables):

        # Decode the values
        #(#_nbr,
        (_a_exc, _b_exc, _tauw_exc,
         _a_inh, _b_inh, _tauw_inh,
        _Trefrac, _Vreset, _delta_v_exc, _delta_v_inh,
        _ampnoise, _Timescale_inv,
         Ve, Vi,
         w_exc, w_inh, _this_h) = values

        # Copy the changed data only
        state_variables[VE] = Ve
        state_variables[VI] = Vi
        state_variables[W_EXC] = w_exc
        state_variables[W_INH] = w_inh
        #state_variables[U] = u

################
###--Meanfield--###
################

   # @property
   # def nbr(self):
    #    return self._nbr


    @property
    def a_exc(self):
        return self._a_exc

    @property
    def b_exc(self):
        return self._b_exc

    @property
    def tauw_exc(self):
        return self._tauw_exc
    
    @property
    def a_inh(self):
        return self._a_inh

    @property
    def b_inh(self):
        return self._b_inh

    @property
    def tauw_inh(self):
        return self._tauw_inh

    @property
    def Trefrac(self):
        return self._Trefrac

    @property
    def Vreset(self):
        return self._Vreset

    @property
    def delta_v_exc(self):
        return self._delta_v_exc
    
    @property
    def delta_v_inh(self):
        return self._delta_v_inh

    @property
    def ampnoise(self):
        return self._ampnoise

    @property
    def Timescale_inv(self):
        return self._Timescale_inv

    @property
    def Ve_init(self):
        """ Settable model parameter: :math:`V_{e}`

        :rtype: float
        """
        return self._Ve_init

    @property
    def Vi_init(self):
        """ Settable model parameter: :math:`V_{i}`

        :rtype: float
        """
        return self._Vi_init    

    @property
    def w_exc_init(self):
        """ Settable model parameter: :math:`w`

        :rtype: float
        """
        return self._w_exc_init    
    
    @property
    def w_inh_init(self):
        """ Settable model parameter: :math:`w`

        :rtype: float
        """
        return self._w_inh_init    

