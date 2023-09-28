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
from .abstract_input_type import AbstractInputType
from spynnaker.pyNN.models.neuron.implementations import (
    AbstractStandardNeuronComponent)

###--Syn and connect params--###
PCONNEC = "pconnec"
Q_EXC = "q_exc"
Q_INH = "q_inh"
TSYN_EXC = "Tsyn_exc"
TSYN_INH = "Tsyn_inh"
EREV_EXC = "Erev_exc"
EREV_INH = "Erev_inh"
NTOT = "Ntot"
GEI = "gei"
EXT_DRIVE = "ext_drive"
AFFERENT_EXC_FRACTION = "afferent_exc_fraction"

GL = "Gl"
CM = "Cm"
EL_EXC = "El_exc"
EL_INH = "El_inh"

MUV = "muV"
MUV0 = "muV0"
ONE_OVER_DMUV0 = "one_over_DmuV0"

SV = "sV"
SV0 = "sV0"
ONE_OVER_DSV0 = "one_over_DsV0"

MUGN = "muGn"

TVN = "TvN"
TVN0 = "TvN0"
ONE_OVER_DTVN0 = "one_over_DTvN0"

VTHRE = "Vthre"

FOUT_TH = "Fout_th"

EXC_NEIGHBOUR_CONTRIBUTION = "exc_neighbour_contribution"
INH_NEIGHBOUR_CONTRIBUTION = "inh_neighbour_contribution"

UNITS = {
    ###--Syn connec--###
    PCONNEC: "",
    Q_EXC: "nS",
    Q_INH: "nS",
    TSYN_EXC: "",
    TSYN_INH: "",
    EREV_EXC: "mV",
    EREV_INH: "mV",
    NTOT: "",
    GEI: "",
    EXT_DRIVE: "",
    AFFERENT_EXC_FRACTION: "",
    GL: "Gl",
    CM: "pF",
    EL_EXC: "mV",
    EL_INH: "mV",
    VTHRE: "mV",
    MUV : "",
    MUV0 : "",
    ONE_OVER_DMUV0 : "",
    SV : "",
    SV0 : "",
    ONE_OVER_DSV0 : "",
    MUGN : "",
    TVN : "",
    TVN0 : "",
    ONE_OVER_DTVN0 : "",
    VTHRE : "",
    FOUT_TH : "",
    EXC_NEIGHBOUR_CONTRIBUTION : "",
    INH_NEIGHBOUR_CONTRIBUTION : "",
}


class ParamsFromNetwork(AbstractInputType):
    """ Model of neuron due to Eugene M. Izhikevich et al
    """
    __slots__ = [
        "__pconnec", "_q_exc", "_q_inh", "_Tsyn_exc", "_Tsyn_inh",
        "_Erev_exc", "_Erev_inh", "_Ntot",
        "_gei", "_ext_drive", "_afferent_exc_fraction",
        "_Gl", "_Cm",
        "_El_exc", "_El_inh",
        "_Vthre", "_muV", "_muV0", "_one_over_DmuV0", "_sV",
        "_sV0", "_one_over_DsV0",
        "_muGn", "_TvN", "_TvN0", "_one_over_DTvN0", "_Vthre", "_Fout_th",
        "_exc_neighbour_contribution", "_inh_neighbour_contribution"
    ]

    def __init__(self, pconnec,
                 q_exc, q_inh,
                 Tsyn_exc, Tsyn_inh,
                 Erev_exc, Erev_inh,
                 Ntot, gei, ext_drive,
                 afferent_exc_fraction,
                 Gl, Cm, El_exc, El_inh,
                 muV, muV0,one_over_DmuV0,
                 sV, sV0, one_over_DsV0,
                 muGn,
                 TvN, TvN0, one_over_DTvN0,
                 Vthre, Fout_th, 
                 exc_neighbour_contribution,
                 inh_neighbour_contribution,):
        """
        :param a: :math:`a`
        :type a: float, iterable(float), ~pyNN.random.RandomDistribution or
            (mapping) function
        
        """
        super().__init__(
            
            [###--syn and connect--###
             DataType.S1615, #pconnec
             DataType.S1615, #q_exc
             DataType.S1615, #q_inh
             DataType.S1615, #Tsyn_exc
             DataType.S1615, #Tsyn_inh
             DataType.S1615, #Erev_exc
             DataType.S1615, #Erev_inh
             DataType.S1615, #Ntot
             DataType.S1615, #gei
             DataType.S1615, #ext_drive
             DataType.S1615, #afferent_exc_fraction
             DataType.S1615, #Gm
             DataType.S1615, #Cl
             DataType.S1615, #El_exc
             DataType.S1615, #El_inh
             DataType.S1615,   # muV
             DataType.S1615,   # muV0
             DataType.S1615,   # one_over_DmuV0
             DataType.S1615,   # sV
             DataType.S1615,   # sV0
             DataType.S1615,   # one_over_DsV0
             DataType.S1615,   # muGn
             DataType.S1615,   # TvN
             DataType.S1615,   # TvN0
             DataType.S1615,   # one_over_DTvN0
             DataType.S1615,   # Vthre
             DataType.S1615,   # Fout_th
             DataType.UINT32,  #exc_neighbour_contribution
             DataType.UINT32,])#inh_neighbour_contribution   
        
        self.__pconnec = pconnec
        self._q_exc = q_exc
        self._q_inh = q_inh
        self._Tsyn_exc = Tsyn_exc
        self._Tsyn_inh = Tsyn_inh
        self._Erev_exc = Erev_exc
        self._Erev_inh = Erev_inh
        self._Ntot = Ntot
        self._gei = gei
        self._ext_drive = ext_drive
        self._afferent_exc_fraction = afferent_exc_fraction
        self._Gl = Gl
        self._Cm = Cm
        self._El_exc = El_exc  
        self._El_inh = El_inh
        self._muV = muV
        self._muV0 = muV0
        self._one_over_DmuV0 = one_over_DmuV0
        self._sV = sV
        self._sV0 = sV0
        self._one_over_DsV0 = one_over_DsV0
        self._muGn = muGn
        self._TvN = TvN
        self._TvN0 = TvN0
        self._one_over_DTvN0 = one_over_DTvN0
        self._Vthre = Vthre
        self._Fout_th = Fout_th
        self._exc_neighbour_contribution = exc_neighbour_contribution
        self._inh_neighbour_contribution = inh_neighbour_contribution

    @overrides(AbstractStandardNeuronComponent.get_n_cpu_cycles)
    def get_n_cpu_cycles(self, n_neurons):
        # A bit of a guess
        return 150 * n_neurons

    @overrides(AbstractStandardNeuronComponent.add_parameters)
    def add_parameters(self, parameters):
        ###--syn and connec--###
        parameters[PCONNEC] = self.__pconnec
        parameters[Q_EXC] = self._q_exc
        parameters[Q_INH] = self._q_inh
        parameters[TSYN_EXC] = self._Tsyn_exc
        parameters[TSYN_INH] = self._Tsyn_inh
        parameters[EREV_EXC] = self._Erev_exc
        parameters[EREV_INH] = self._Erev_inh
        parameters[NTOT] = self._Ntot
        parameters[GEI] = self._gei
        parameters[EXT_DRIVE] = self._ext_drive
        parameters[AFFERENT_EXC_FRACTION] = self._afferent_exc_fraction
        parameters[GL] = self._Gl
        parameters[CM] = self._Cm
        parameters[EL_EXC] = self._El_exc
        parameters[EL_INH] = self._El_inh
        #parameters[MUV] = self._muV
        parameters[MUV0] = self._muV0
        parameters[ONE_OVER_DMUV0] = self._one_over_DmuV0
        parameters[SV0] = self._sV0
        parameters[ONE_OVER_DSV0] = self._one_over_DsV0
        parameters[TVN0] = self._TvN0
        parameters[ONE_OVER_DTVN0] = self._one_over_DTvN0

    @overrides(AbstractStandardNeuronComponent.add_state_variables)
    def add_state_variables(self, state_variables):
        state_variables[MUV] = self._muV
        state_variables[SV] = self._sV
        state_variables[MUGN] = self._muGn
        state_variables[TVN] = self._TvN
        state_variables[VTHRE] = self._Vthre
        state_variables[FOUT_TH] = self._Fout_th
        state_variables[EXC_NEIGHBOUR_CONTRIBUTION] = self._exc_neighbour_contribution
        state_variables[INH_NEIGHBOUR_CONTRIBUTION] = self._inh_neighbour_contribution

    @overrides(AbstractStandardNeuronComponent.get_units)
    def get_units(self, variable):
        return UNITS[variable]

    @overrides(AbstractStandardNeuronComponent.has_variable)
    def has_variable(self, variable):
        return variable in UNITS

    @overrides(AbstractNeuronModel.get_global_values)
    def get_global_values(self, ts):
        # pylint: disable=arguments-differ
        pass

    @overrides(AbstractStandardNeuronComponent.get_values)
    def get_values(self, parameters, state_variables, vertex_slice, ts):
        """
        :param ts: machine time step
        """
        # pylint: disable=arguments-differ

        # Add the rest of the data
        return [parameters[PCONNEC],###syn and connect
                parameters[Q_EXC],
                parameters[Q_INH],
                parameters[TSYN_EXC],
                parameters[TSYN_INH],
                parameters[EREV_EXC],
                parameters[EREV_INH],
                parameters[NTOT],
                parameters[GEI],
                parameters[EXT_DRIVE],
                parameters[AFFERENT_EXC_FRACTION],
                parameters[GL],
                parameters[CM],
                parameters[EL_EXC],
                parameters[EL_INH],
                state_variables[MUV],
                parameters[MUV0],
                parameters[ONE_OVER_DMUV0],
                state_variables[SV],
                parameters[SV0],
                parameters[ONE_OVER_DSV0],
                state_variables[MUGN],
                state_variables[TVN],
                parameters[TVN0],
                parameters[ONE_OVER_DTVN0],
                state_variables[VTHRE],
                state_variables[FOUT_TH],
                state_variables[EXC_NEIGHBOUR_CONTRIBUTION],
                state_variables[INH_NEIGHBOUR_CONTRIBUTION]
        ]

    @overrides(AbstractStandardNeuronComponent.update_values)
    def update_values(self, values, parameters, state_variables):

        # Decode the values
        (__pconnec, _q_exc, _q_inh, _Tsyn_exc, _Tsyn_inh,
        _Erev_exc, _Erev_inh, _Ntot, _gei, _ext_drive,
        _afferent_exc_fraction,
        _Gl, _Cm, _El_exc, _El_inh,
        muV, _muV0, _one_over_DmuV0,
        sV, _sV0, _one_over_DsV0,
        muGn,
        TvN, _TvN0, _one_over_DTvN0,
        Vthre, Fout_th,
        exc_neighbour_contribution,
        inh_neighbour_contribution) = values

        # Copy the changed data only
        state_variables[MUV] = muV
        state_variables[SV] = sV
        state_variables[MUGN] = muGn
        state_variables[TVN] = TvN
        state_variables[VTHRE] = Vthre
        state_variables[FOUT_TH] = Fout_th
        state_variables[EXC_NEIGHBOUR_CONTRIBUTION] = exc_neighbour_contribution
        state_variables[INH_NEIGHBOUR_CONTRIBUTION] = inh_neighbour_contribution
        
    @overrides(AbstractInputType.get_global_weight_scale)
    def get_global_weight_scale(self):
        return 1024.0

########################
###--syn and connec--###
########################
    @property
    def pconnec(self):
        return self.__pconnec
    
    @property
    def q_exc(self):
        return self._q_exc

    @property
    def q_inh(self):
        return self._q_inh

    @property
    def Tsyn_exc(self):
        return self._Tsyn_exc

    @property
    def Tsyn_inh(self):
        return self._Tsyn_inh

    @property
    def Erev_exc(self):
        return self._Erev_exc

    @property
    def Erev_inh(self):
        return self._Erev_inh

    @Erev_inh.setter
    def Erev_inh(self, Erev_inh):
        self._Erev_inh = Erev_inh

    @property
    def Ntot(self):
        return self._Ntot

    @property
    def gei(self):
        return self._gei

    @property
    def ext_drive(self):
        return self._ext_drive

    @property
    def afferent_exc_fraction(self):
        return self._afferent_exc_fraction
    
    @property
    def Gl(self):
        return self._Gl

    @property
    def Cm(self):
        return self._Cm

    @property
    def El_exc(self):
        return self._El_exc
    
    @property
    def El_inh(self):
        return self._El_inh
        
    @property
    def muV(self):
        """ Settable model parameter: :math:`a`

        :rtype: float
        """
        return self._muV

    @property
    def muV0(self):
        """ Settable model parameter: :math:`b`

        :rtype: float
        """
        return self._muV0

    @property
    def one_over_DmuV0(self):
        """ Settable model parameter: :math:`c`

        :rtype: float
        """
        return self._one_over_DmuV0

    @property
    def sV(self):
        """ Settable model parameter: :math:`d`

        :rtype: float
        """
        return self._sV
    
    @property
    def sV0(self):
        """ Settable model parameter: :math:`d`

        :rtype: float
        """
        return self._sV0
    
    @property
    def one_over_DsV0(self):
        """ Settable model parameter: :math:`d`

        :rtype: float
        """
        return self._one_over_DsV0

    @property
    def muGn(self):
        """ Settable model parameter: :math:`d`

        :rtype: float
        """
        return self._muGn
    
    @property
    def TvN(self):
        """ Settable model parameter: :math:`d`

        :rtype: float
        """
        return self._TvN
    
    @property
    def TvN0(self):
        """ Settable model parameter: :math:`d`

        :rtype: float
        """
        return self._TvN0
    
    @property
    def one_over_DTvN0(self):
        """ Settable model parameter: :math:`d`

        :rtype: float
        """
        return self._one_over_DTvN0
    
    @property
    def Vthre(self):
        """ Settable model parameter: :math:`d`

        :rtype: float
        """
        return self._Vthre

    @property
    def Fout_th(self):
        """ Settable model parameter: :math:`d`

        :rtype: float
        """
        return self._Fout_th
    
    @property
    def exc_neighbour_contribution(self):
        """ Settable model parameter: :math:`d`

        :rtype: float
        """
        return self._exc_neighbour_contribution
    
    @property
    def inh_neighbour_contribution(self):
        """ Settable model parameter: :math:`d`

        :rtype: float
        """
        return self._inh_neighbour_contribution