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

from spynnaker.pyNN.models.neuron.input_types import InputTypeConductance
from spynnaker.pyNN.models.neuron.neuron_models import MeanfieldOfAdexNetwork
from spynnaker.pyNN.models.neuron.neuron_models import ParamsFromNetwork
from spynnaker.pyNN.models.neuron.neuron_models import pFitPolynomialExc
from spynnaker.pyNN.models.neuron.neuron_models import pFitPolynomialInh
from spynnaker.pyNN.models.neuron.neuron_models import Mathsbox
from spynnaker.pyNN.models.neuron.synapse_types import SynapseTypeExponential
from spynnaker.pyNN.models.neuron.threshold_types import ThresholdTypeStatic
from spynnaker.pyNN.models.neuron import AbstractPyNNMeanfieldModelStandard
from spynnaker.pyNN.models.defaults import default_initial_values

_IZK_THRESHOLD = 30.0


class MeanfieldBase(AbstractPyNNMeanfieldModelStandard):
    """ Izhikevich neuron model with conductance inputs.

    :param a: :math:`a`
    :type a: float, iterable(float), ~pyNN.random.RandomDistribution
        or (mapping) function
    """

    # noinspection PyPep8Naming
    @default_initial_values({"Ve", "Vi", "w_exc", "w_inh", "Fout_th",
                             "muV", "sV", "muGn",
                             "TvN", "Vthre",
                             "err_func", "isyn_exc", "isyn_inh"})
    def __init__(self,
                 a_exc=0.,
                 b_exc=60.,
                 tauw_exc=100.0,
                 a_inh=0.,
                 b_inh=0.,
                 tauw_inh=1.0,
                 Trefrac=5.0,
                 Vreset=-65.,
                 delta_v_exc=2.0,
                 delta_v_inh=0.5,
                 ampnoise=0.0,
                 Timescale_inv=200.,
                 Ve=9.,
                 Vi=23.,
                 w_exc=0.25,
                 w_inh=0.00001, #arbitraire

                 pconnec=0.05,
                 q_exc=1.5,
                 q_inh=5.,
                 Tsyn_exc=5.0,
                 Tsyn_inh=5.0,
                 Erev_exc=0.,
                 Erev_inh=-80.,
                 Ntot=10000.,
                 gei=0.2,
                 ext_drive=2.5,
                 afferent_exc_fraction=1.,

                 Gl=10.,
                 Cm=200.,
                 El_exc=-70.,
                 El_inh=-65.,

                 muV=0.,
                 muV0=-0.06,
                 DmuV0=0.01,

                 sV=0.,
                 sV0=0.004,
                 DsV0=0.006,

                 muGn=0.,

                 TvN=0.,
                 TvN0=0.5,
                 DTvN0=1.,

                 Vthre=-50.,
                 Fout_th=0.,

                 p0_exc=-0.0515518,
                 p1_exc=0.00455197,
                 p2_exc=-0.00760625,
                 p3_exc=0.00094851,
                 p4_exc=0.001,
                 p5_exc=-0.0009863,
                 p6_exc=-0.0026474,
                 p7_exc=-0.0135417,
                 p8_exc=0.0028742,
                 p9_exc=0.0029213,
                 p10_exc=-0.014084,
                 
                 p0_inh=-0.0496832,
                 p1_inh=0.00412289,
                 p2_inh=-0.0054849,
                 p3_inh=-0.0013451,
                 p4_inh=0.001,
                 p5_inh=-0.0010459,
                 p6_inh=0.00306102,
                 p7_inh=-0.0084485,
                 p8_inh=-0.0025717,
                 p9_inh=0.00179862,
                 p10_inh=-0.013830,


                 tau_syn_E=5.0,
                 tau_syn_I=5.0,
                 e_rev_E=0.0,
                 e_rev_I=-70.0,

                 isyn_exc=0.0,
                 isyn_inh=0.0,

                 sample=1000,
                 err_func=0.):
        
        # pylint: disable=too-many-arguments, too-many-locals
        one_over_DmuV0 = 1/DmuV0
        one_over_DsV0 = 1/DsV0
        one_over_DTvN0 = 1/DTvN0
        one_over_Cm = 1/Cm
        one_over_Gl = 1/Gl
        
        muGe_0 = q_exc * Tsyn_exc * (Ve+ext_drive) * (1.-gei) * pconnec*Ntot
        muGi_0 = q_inh * Tsyn_inh * Vi * gei * pconnec * Ntot
        muG_0 = Gl + muGe_0 + muGi_0
        muVV = ((muGe_0*Erev_exc + muGi_0*Erev_inh + Gl*El_exc\
                 - (Ve)*tauw_exc*(b_exc) + a_exc*El_exc) /muG_0)/ (1+a_exc/muG_0)
        w_exc = Ve*tauw_exc*b_exc - a_exc * (El_exc-muVV)
        w_inh = 0.00001
        
       
        
        neuron_model = MeanfieldOfAdexNetwork(a_exc, b_exc, tauw_exc,
                                              a_inh, b_inh, tauw_inh,
                                              Trefrac, Vreset,
                                              delta_v_exc, delta_v_inh,
                                              ampnoise, Timescale_inv,
                                              Ve, Vi,
                                              w_exc, w_inh)
        params_from_network = ParamsFromNetwork(pconnec, q_exc, q_inh,
                                                Tsyn_exc, Tsyn_inh,
                                                Erev_exc, Erev_inh,
                                                Ntot, gei, ext_drive,
                                                afferent_exc_fraction,
                                                Gl, Cm,
                                                El_exc, El_inh,
                                                muV, muV0, one_over_DmuV0,
                                                sV, sV0, one_over_DsV0,
                                                muGn,
                                                TvN, TvN0, one_over_DTvN0,
                                                Vthre, Fout_th)
        p_fit_polynomial_exc = pFitPolynomialExc(p0_exc, p1_exc, p2_exc,
                                              p3_exc, p4_exc, p5_exc, 
                                              p6_exc, p7_exc, p8_exc,
                                              p9_exc, p10_exc)
        p_fit_polynomial_inh = pFitPolynomialInh(p0_inh, p1_inh, p2_inh,
                                              p3_inh, p4_inh, p5_inh,
                                              p6_inh, p7_inh, p8_inh,
                                              p9_inh, p10_inh)
        mathsbox = Mathsbox(sample, err_func)#, cycles_numbre)
        synapse_type = SynapseTypeExponential(
            tau_syn_E, tau_syn_I, isyn_exc, isyn_inh)
        input_type = InputTypeConductance(e_rev_E, e_rev_I)
        threshold_type = ThresholdTypeStatic(_IZK_THRESHOLD)

        super().__init__(
            model_name="meanfield_model_cond",
            binary="meanfield_model_cond.aplx",
            neuron_model=neuron_model,
            params_from_network=params_from_network,
            p_fit_polynomial_exc = p_fit_polynomial_exc, 
            p_fit_polynomial_inh = p_fit_polynomial_inh, 
            mathsbox=mathsbox,
            input_type=input_type,
            synapse_type=synapse_type,
            threshold_type=threshold_type)
