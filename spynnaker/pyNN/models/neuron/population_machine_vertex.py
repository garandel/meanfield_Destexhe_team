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
from enum import Enum

from pacman.executor.injection_decorator import inject_items
from spinn_front_end_common.interface.simulation import simulation_utilities
from spinn_utilities.overrides import overrides
from pacman.model.graphs.machine import MachineVertex
from spinn_front_end_common.utilities.utility_objs import ProvenanceDataItem
from spinn_front_end_common.interface.provenance import (
    ProvidesProvenanceDataFromMachineImpl)
from spinn_front_end_common.interface.buffer_management.buffer_models import (
    AbstractReceiveBuffersToHost)
from spinn_front_end_common.utilities.helpful_functions import (
    locate_memory_region_for_placement)
from spinn_front_end_common.abstract_models import (
    AbstractHasAssociatedBinary, AbstractSupportsBitFieldGeneration,
    AbstractSupportsBitFieldRoutingCompression,
    AbstractGeneratesDataSpecification, AbstractRewritesDataSpecification)
from spinn_front_end_common.interface.profiling import (
    AbstractHasProfileData, profile_utils)
from spinn_front_end_common.interface.profiling.profile_utils import (
    get_profiling_data)
from spinn_front_end_common.utilities.utility_objs import ExecutableType
from spynnaker.pyNN.models.neuron.synapse_dynamics import (
    AbstractSynapseDynamicsStructural)
from spynnaker.pyNN.utilities import constants, bit_field_utilities
from spynnaker.pyNN.models.abstract_models import (
    AbstractSynapseExpandable, AbstractReadParametersBeforeSet)
from spynnaker.pyNN.utilities.constants import POPULATION_BASED_REGIONS
from spinn_front_end_common.utilities import (
    constants as common_constants, helpful_functions)


class PopulationMachineVertex(
        MachineVertex, AbstractReceiveBuffersToHost,
        AbstractHasAssociatedBinary, ProvidesProvenanceDataFromMachineImpl,
        AbstractHasProfileData, AbstractSupportsBitFieldGeneration,
        AbstractSupportsBitFieldRoutingCompression,
        AbstractGeneratesDataSpecification, AbstractSynapseExpandable,
        AbstractRewritesDataSpecification, AbstractReadParametersBeforeSet):

    __slots__ = [
        "__binary_file_name",
        "__recorded_region_ids",
        "__resources",
        "__on_chip_generatable_offset",
        "__on_chip_generatable_size",
        "__drop_late_spikes",
        "__change_requires_neuron_parameters_reload"]

    class EXTRA_PROVENANCE_DATA_ENTRIES(Enum):
        """ Entries for the provenance data generated by standard neuron \
            models.
        """
        #: The number of pre-synaptic events
        PRE_SYNAPTIC_EVENT_COUNT = 0
        #: The number of times the synapse arithmetic saturated
        SATURATION_COUNT = 1
        #: The number of times there was a buffer overflow
        BUFFER_OVERFLOW_COUNT = 2
        #: The current timer tick
        CURRENT_TIMER_TIC = 3
        #: The number of times the plastic synapses saturated during weight
        #: calculation
        PLASTIC_SYNAPTIC_WEIGHT_SATURATION_COUNT = 4
        GHOST_POP_TABLE_SEARCHES = 5
        FAILED_TO_READ_BIT_FIELDS = 6
        DMA_COMPLETES = 7
        SPIKE_PROGRESSING_COUNT = 8
        INVALID_MASTER_POP_HITS = 9
        BIT_FIELD_FILTERED_COUNT = 10
        N_REWIRES = 11
        # the number of packets that were dropped as they arrived too late
        # to be processed
        N_LATE_SPIKES = 12
        # the max filled size of the input buffer
        INPUT_BUFFER_FILLED_SIZE = 13
        # the number of tdma misses
        TDMA_MISSES = 14
        # the maxmimum number of background tasks queued
        MAX_BACKGROUND_QUEUED = 15
        # the number of times the background queue overloaded
        N_BACKGROUND_OVERLOADS = 16

    SATURATION_COUNT_NAME = "Times_synaptic_weights_have_saturated"
    _SATURATION_COUNT_MESSAGE = (
        "The weights from the synapses for {} on {}, {}, {} saturated "
        "{} times. If this causes issues you can increase the "
        "spikes_per_second and / or ring_buffer_sigma "
        "values located within the .spynnaker.cfg file.")

    INPUT_BUFFER_FULL_NAME = "Times_the_input_buffer_lost_packets"
    _INPUT_BUFFER_FULL_MESSAGE = (
        "The input buffer for {} on {}, {}, {} lost packets on {} "
        "occasions. This is often a sign that the system is running "
        "too quickly for the number of neurons per core.  Please "
        "increase the timer_tic or time_scale_factor or decrease the "
        "number of neurons per core.")

    TOTAL_PRE_SYNAPTIC_EVENT_NAME = "Total_pre_synaptic_events"
    LAST_TIMER_TICK_NAME = "Last_timer_tic_the_core_ran_to"
    N_RE_WIRES_NAME = "Number_of_rewires"

    SATURATED_PLASTIC_WEIGHTS_NAME = (
        "Times_plastic_synaptic_weights_have_saturated")
    _SATURATED_PLASTIC_WEIGHTS_MESSAGE = (
        "The weights from the plastic synapses for {} on {}, {}, {} "
        "saturated {} times. If this causes issue increase the "
        "spikes_per_second and / or ring_buffer_sigma values located "
        "within the .spynnaker.cfg file.")

    _N_LATE_SPIKES_NAME = "Number_of_late_spikes"
    _N_LATE_SPIKES_MESSAGE_DROP = (
        "{} packets from {} on {}, {}, {} were dropped from the input buffer, "
        "because they arrived too late to be processed in a given time step. "
        "Try increasing the time_scale_factor located within the "
        ".spynnaker.cfg file or in the pynn.setup() method.")
    _N_LATE_SPIKES_MESSAGE_NO_DROP = (
        "{} packets from {} on {}, {}, {} arrived too late to be processed in"
        " a given time step. "
        "Try increasing the time_scale_factor located within the "
        ".spynnaker.cfg file or in the pynn.setup() method.")

    _MAX_FILLED_SIZE_OF_INPUT_BUFFER_NAME = "Max_filled_size_input_buffer"

    _BACKGROUND_OVERLOADS_NAME = "Times_the_background_queue_overloaded"
    _BACKGROUND_MAX_QUEUED_NAME = "Max_backgrounds_queued"

    _PROFILE_TAG_LABELS = {
        0: "TIMER",
        1: "DMA_READ",
        2: "INCOMING_SPIKE",
        3: "PROCESS_FIXED_SYNAPSES",
        4: "PROCESS_PLASTIC_SYNAPSES"}

    # x words needed for a bitfield covering 256 atoms
    _WORDS_TO_COVER_256_ATOMS = 8

    # provenance data items
    BIT_FIELD_FILTERED_PACKETS = \
        "How many packets were filtered by the bitfield filterer."
    INVALID_MASTER_POP_HITS = "Invalid Master Pop hits"
    SPIKES_PROCESSED = "how many spikes were processed"
    DMA_COMPLETE = "DMA's that were completed"
    BIT_FIELDS_NOT_READ = "N bit fields not able to be read into DTCM"
    GHOST_SEARCHES = "Number of failed pop table searches"
    PLASTIC_WEIGHT_SATURATION = "Times_plastic_synaptic_weights_have_saturated"
    LAST_TIMER_TICK = "Last_timer_tic_the_core_ran_to"
    TOTAL_PRE_SYNAPTIC_EVENTS = "Total_pre_synaptic_events"
    LOST_INPUT_BUFFER_PACKETS = "Times_the_input_buffer_lost_packets"

    N_ADDITIONAL_PROVENANCE_DATA_ITEMS = len(EXTRA_PROVENANCE_DATA_ENTRIES)

    def __init__(
            self, resources_required, recorded_region_ids, label, constraints,
            app_vertex, vertex_slice, drop_late_spikes, binary_file_name):
        """
        :param ~pacman.model.resources.ResourceContainer resources_required:
        :param iterable(int) recorded_region_ids:
        :param str label:
        :param bool drop_late_spikes: control flag for dropping packets.
        :param list(~pacman.model.constraints.AbstractConstraint) constraints:
        :param AbstractPopulationVertex app_vertex:
            The associated application vertex
        :param ~pacman.model.graphs.common.Slice vertex_slice:
            The slice of the population that this implements
        :param str binary_file_name: binary name to be run for this verte
        """
        super().__init__(label, constraints, app_vertex, vertex_slice)
        self.__binary_file_name = binary_file_name
        self.__recorded_region_ids = recorded_region_ids
        self.__resources = resources_required
        self.__drop_late_spikes = drop_late_spikes
        self.__on_chip_generatable_offset = None
        self.__on_chip_generatable_size = None
        self.__change_requires_neuron_parameters_reload = False

    def set_on_chip_generatable_area(self, offset, size):
        self.__on_chip_generatable_offset = offset
        self.__on_chip_generatable_size = size

    @overrides(AbstractSupportsBitFieldGeneration.bit_field_base_address)
    def bit_field_base_address(self, transceiver, placement):
        return locate_memory_region_for_placement(
            placement=placement, transceiver=transceiver,
            region=POPULATION_BASED_REGIONS.BIT_FIELD_FILTER.value)

    @overrides(AbstractSupportsBitFieldRoutingCompression.
               key_to_atom_map_region_base_address)
    def key_to_atom_map_region_base_address(self, transceiver, placement):
        return locate_memory_region_for_placement(
            placement=placement, transceiver=transceiver,
            region=POPULATION_BASED_REGIONS.BIT_FIELD_KEY_MAP.value)

    @overrides(AbstractSupportsBitFieldGeneration.bit_field_builder_region)
    def bit_field_builder_region(self, transceiver, placement):
        return locate_memory_region_for_placement(
            placement=placement, transceiver=transceiver,
            region=POPULATION_BASED_REGIONS.BIT_FIELD_BUILDER.value)

    @overrides(AbstractSupportsBitFieldRoutingCompression.
               regeneratable_sdram_blocks_and_sizes)
    def regeneratable_sdram_blocks_and_sizes(self, transceiver, placement):
        synaptic_matrix_base_address = locate_memory_region_for_placement(
            placement=placement, transceiver=transceiver,
            region=POPULATION_BASED_REGIONS.SYNAPTIC_MATRIX.value)
        return [(
            self.__on_chip_generatable_offset + synaptic_matrix_base_address,
            self.__on_chip_generatable_size)]

    @property
    @overrides(MachineVertex.resources_required)
    def resources_required(self):
        return self.__resources

    @property
    @overrides(ProvidesProvenanceDataFromMachineImpl._provenance_region_id)
    def _provenance_region_id(self):
        return POPULATION_BASED_REGIONS.PROVENANCE_DATA.value

    @property
    @overrides(ProvidesProvenanceDataFromMachineImpl._n_additional_data_items)
    def _n_additional_data_items(self):
        return len(PopulationMachineVertex.EXTRA_PROVENANCE_DATA_ENTRIES)

    @overrides(ProvidesProvenanceDataFromMachineImpl.
               get_provenance_data_from_machine)
    def get_provenance_data_from_machine(self, transceiver, placement):
        provenance_data = self._read_provenance_data(transceiver, placement)
        provenance_items = self._read_basic_provenance_items(
            provenance_data, placement)
        provenance_data = self._get_remaining_provenance_data_items(
            provenance_data)

        times_timer_tic_overran = 0
        for item in provenance_items:
            if item.names[-1] == self._TIMER_TICK_OVERRUN:
                times_timer_tic_overran = item.value

        n_saturations = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.SATURATION_COUNT.value]
        n_buffer_overflows = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.BUFFER_OVERFLOW_COUNT.value]
        n_pre_synaptic_events = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.PRE_SYNAPTIC_EVENT_COUNT.value]
        last_timer_tick = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.CURRENT_TIMER_TIC.value]
        n_plastic_saturations = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.
            PLASTIC_SYNAPTIC_WEIGHT_SATURATION_COUNT.value]
        n_ghost_searches = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.GHOST_POP_TABLE_SEARCHES.value]
        failed_to_read_bit_fields = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.FAILED_TO_READ_BIT_FIELDS.value]
        dma_completes = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.DMA_COMPLETES.value]
        spike_processing_count = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.SPIKE_PROGRESSING_COUNT.value]
        invalid_master_pop_hits = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.INVALID_MASTER_POP_HITS.value]
        n_packets_filtered_by_bit_field_filter = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.BIT_FIELD_FILTERED_COUNT.value]
        n_rewires = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.N_REWIRES.value]
        n_late_packets = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.N_LATE_SPIKES.value]
        input_buffer_max_filled_size = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.INPUT_BUFFER_FILLED_SIZE.value]
        tdma_misses = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.TDMA_MISSES.value]
        max_background_queued = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.MAX_BACKGROUND_QUEUED.value]
        n_background_overloads = provenance_data[
            self.EXTRA_PROVENANCE_DATA_ENTRIES.N_BACKGROUND_OVERLOADS.value]

        label, x, y, p, names = self._get_placement_details(placement)

        # translate into provenance data items
        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self.SATURATION_COUNT_NAME),
            n_saturations, report=n_saturations > 0,
            message=self._SATURATION_COUNT_MESSAGE.format(
                label, x, y, p, n_saturations)))
        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self.INPUT_BUFFER_FULL_NAME),
            n_buffer_overflows, report=n_buffer_overflows > 0,
            message=self._INPUT_BUFFER_FULL_MESSAGE.format(
                label, x, y, p, n_buffer_overflows)))
        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self.TOTAL_PRE_SYNAPTIC_EVENT_NAME),
            n_pre_synaptic_events))
        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self.LAST_TIMER_TICK_NAME),
            last_timer_tick))
        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self.SATURATED_PLASTIC_WEIGHTS_NAME),
            n_plastic_saturations, report=n_plastic_saturations > 0,
            message=self._SATURATED_PLASTIC_WEIGHTS_MESSAGE.format(
                label, x, y, p, n_plastic_saturations)))
        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self.N_RE_WIRES_NAME), n_rewires))
        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self.GHOST_SEARCHES), n_ghost_searches,
            report=n_ghost_searches > 0,
            message=(
                "The number of failed population table searches for {} on {},"
                " {}, {} was {}. If this number is large relative to the "
                "predicted incoming spike rate, try increasing source and "
                "target neurons per core".format(
                    label, x, y, p, n_ghost_searches))))
        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self.BIT_FIELDS_NOT_READ),
            failed_to_read_bit_fields, report=False,
            message=(
                "The filter for stopping redundant DMA's couldn't be fully "
                "filled in, it failed to read {} entries, which means it "
                "required a max of {} extra bytes of DTCM (assuming cores "
                "have at max 255 neurons. Try reducing neurons per core, or "
                "size of buffers, or neuron params per neuron etc.".format(
                    failed_to_read_bit_fields,
                    failed_to_read_bit_fields *
                    self._WORDS_TO_COVER_256_ATOMS))))
        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self.DMA_COMPLETE), dma_completes))
        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self.SPIKES_PROCESSED),
            spike_processing_count))
        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self.INVALID_MASTER_POP_HITS),
            invalid_master_pop_hits, report=invalid_master_pop_hits > 0,
            message=(
                "There were {} keys which were received by core {}:{}:{} which"
                " had no master pop entry for it. This is a error, which most "
                "likely strives from bad routing.".format(
                    invalid_master_pop_hits, x, y, p))))
        provenance_items.append((ProvenanceDataItem(
            self._add_name(names, self.BIT_FIELD_FILTERED_PACKETS),
            n_packets_filtered_by_bit_field_filter,
            report=(
                n_packets_filtered_by_bit_field_filter > 0 and (
                    n_buffer_overflows > 0 or times_timer_tic_overran > 0)),
            message=(
                "There were {} packets received by {}:{}:{} that were "
                "filtered by the Bitfield filterer on the core. These packets "
                "were having to be stored and processed on core, which means "
                "the core may not be running as efficiently as it could. "
                "Please adjust the network or the mapping so that these "
                "packets are filtered in the router to improve "
                "performance.".format(
                    n_packets_filtered_by_bit_field_filter, x, y, p)))))
        late_message = (
            self._N_LATE_SPIKES_MESSAGE_DROP if self.__drop_late_spikes
            else self._N_LATE_SPIKES_MESSAGE_NO_DROP)
        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self._N_LATE_SPIKES_NAME),
            n_late_packets, report=n_late_packets > 0,
            message=late_message.format(n_late_packets, label, x, y, p)))

        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self._MAX_FILLED_SIZE_OF_INPUT_BUFFER_NAME),
            input_buffer_max_filled_size, report=False))

        provenance_items.append(self._app_vertex.get_tdma_provenance_item(
            names, x, y, p, tdma_misses))

        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self._BACKGROUND_MAX_QUEUED_NAME),
            max_background_queued, report=max_background_queued > 1,
            message=(
                "A maximum of {} background tasks were queued on {} on"
                " {}, {}, {}.  Try increasing the time_scale_factor located"
                " within the .spynnaker.cfg file or in the pynn.setup()"
                " method.".format(max_background_queued, label, x, y, p))))
        provenance_items.append(ProvenanceDataItem(
            self._add_name(names, self._BACKGROUND_OVERLOADS_NAME),
            n_background_overloads, report=n_background_overloads > 0,
            message=(
                "On {} on {}, {}, {}, the background queue overloaded {}"
                " times.  Try increasing the time_scale_factor located within"
                " the .spynnaker.cfg file or in the pynn.setup() method."
                .format(label, x, y, p, n_background_overloads))))
        return provenance_items

    @overrides(AbstractReceiveBuffersToHost.get_recorded_region_ids)
    def get_recorded_region_ids(self):
        return self.__recorded_region_ids

    @overrides(AbstractReceiveBuffersToHost.get_recording_region_base_address)
    def get_recording_region_base_address(self, txrx, placement):
        return locate_memory_region_for_placement(
            placement, POPULATION_BASED_REGIONS.NEURON_RECORDING.value, txrx)

    @overrides(AbstractHasProfileData.get_profile_data)
    def get_profile_data(self, transceiver, placement):
        return get_profiling_data(
            POPULATION_BASED_REGIONS.PROFILING.value,
            self._PROFILE_TAG_LABELS, transceiver, placement)

    @overrides(AbstractHasAssociatedBinary.get_binary_file_name)
    def get_binary_file_name(self):
        return self.__binary_file_name

    @overrides(AbstractHasAssociatedBinary.get_binary_start_type)
    def get_binary_start_type(self):
        return ExecutableType.USES_SIMULATION_INTERFACE

    @inject_items({
        "machine_time_step": "MachineTimeStep",
        "time_scale_factor": "TimeScaleFactor",
        "application_graph": "MemoryApplicationGraph",
        "machine_graph": "MemoryMachineGraph",
        "routing_info": "MemoryRoutingInfos",
        "data_n_time_steps": "DataNTimeSteps",
        "n_key_map": "MemoryMachinePartitionNKeysMap"
    })
    @overrides(
        AbstractGeneratesDataSpecification.generate_data_specification,
        additional_arguments={
            "machine_time_step", "time_scale_factor",
            "application_graph", "machine_graph", "routing_info",
            "data_n_time_steps", "n_key_map"
        })
    def generate_data_specification(
            self, spec, placement, machine_time_step, time_scale_factor,
            application_graph, machine_graph, routing_info, data_n_time_steps,
            n_key_map):
        """
        :param machine_time_step: (injected)
        :param time_scale_factor: (injected)
        :param application_graph: (injected)
        :param machine_graph: (injected)
        :param routing_info: (injected)
        :param data_n_time_steps: (injected)
        :param n_key_map: (injected)
        """
        # pylint: disable=too-many-arguments, arguments-differ

        spec.comment("\n*** Spec for block of {} neurons ***\n".format(
            self._app_vertex.neuron_impl.model_name))

        # Reserve memory regions
        self._reserve_memory_regions(spec, machine_graph, n_key_map)

        # Declare random number generators and distributions:
        # TODO add random distribution stuff
        # self.write_random_distribution_declarations(spec)

        # Get the key
        key = routing_info.get_first_key_from_pre_vertex(
            self, constants.SPIKE_PARTITION_ID)

        # Write the setup region
        spec.switch_write_focus(POPULATION_BASED_REGIONS.SYSTEM.value)
        spec.write_array(simulation_utilities.get_simulation_header_array(
            self.__binary_file_name, machine_time_step, time_scale_factor))

        # Write the neuron recording region
        self._app_vertex.neuron_recorder.write_neuron_recording_region(
            spec, POPULATION_BASED_REGIONS.NEURON_RECORDING.value,
            self.vertex_slice, data_n_time_steps)

        # Write the neuron parameters
        self._write_neuron_parameters(
            spec, key, constants.POPULATION_BASED_REGIONS.NEURON_PARAMS.value)

        # write profile data
        profile_utils.write_profile_region_data(
            spec, POPULATION_BASED_REGIONS.PROFILING.value,
            self._app_vertex.n_profile_samples)

        # Get the weight_scale value from the appropriate location
        weight_scale = self._app_vertex.neuron_impl.get_global_weight_scale()

        # allow the synaptic matrix to write its data spec-able data
        self._app_vertex.synapse_manager.write_data_spec(
            spec, self._app_vertex, self.vertex_slice, self, machine_graph,
            application_graph, routing_info, weight_scale, machine_time_step)
        self.set_on_chip_generatable_area(
            self._app_vertex.synapse_manager.host_written_matrix_size(
                self.vertex_slice),
            self._app_vertex.synapse_manager.on_chip_written_matrix_size(
                self.vertex_slice))

        # write up the bitfield builder data
        bit_field_utilities.write_bitfield_init_data(
            spec, self, machine_graph, routing_info,
            n_key_map, POPULATION_BASED_REGIONS.BIT_FIELD_BUILDER.value,
            POPULATION_BASED_REGIONS.POPULATION_TABLE.value,
            POPULATION_BASED_REGIONS.SYNAPTIC_MATRIX.value,
            POPULATION_BASED_REGIONS.DIRECT_MATRIX.value,
            POPULATION_BASED_REGIONS.BIT_FIELD_FILTER.value,
            POPULATION_BASED_REGIONS.BIT_FIELD_KEY_MAP.value,
            POPULATION_BASED_REGIONS.STRUCTURAL_DYNAMICS.value,
            isinstance(
                self._app_vertex.synapse_manager.synapse_dynamics,
                AbstractSynapseDynamicsStructural))

        # End the writing of this specification:
        spec.end_specification()

    @inject_items({"routing_info": "MemoryRoutingInfos"})
    @overrides(
        AbstractRewritesDataSpecification.regenerate_data_specification,
        additional_arguments={"routing_info"})
    def regenerate_data_specification(self, spec, placement, routing_info):
        # pylint: disable=too-many-arguments, arguments-differ

        # reserve the neuron parameters data region
        self._reserve_neuron_params_data_region(spec)

        # write the neuron params into the new DSG region
        self._write_neuron_parameters(
            key=routing_info.get_first_key_from_pre_vertex(
                self, constants.SPIKE_PARTITION_ID),
            spec=spec,
            region_id=constants.POPULATION_BASED_REGIONS.NEURON_PARAMS.value)

        # close spec
        spec.end_specification()

    @overrides(AbstractRewritesDataSpecification.reload_required)
    def reload_required(self):
        return self.__change_requires_neuron_parameters_reload

    @overrides(AbstractRewritesDataSpecification.set_reload_required)
    def set_reload_required(self, new_value):
        self.__change_requires_neuron_parameters_reload = new_value

    def _reserve_memory_regions(self, spec, machine_graph, n_key_map):
        """ Reserve the DSG data regions.

        :param ~.DataSpecificationGenerator spec:
            the spec to write the DSG region to
        :param ~.MachineGraph machine_graph: machine graph
        :param n_key_map: n key map
        :return: None
        """
        spec.comment("\nReserving memory space for data regions:\n\n")

        # Reserve memory:
        spec.reserve_memory_region(
            region=POPULATION_BASED_REGIONS.SYSTEM.value,
            size=common_constants.SIMULATION_N_BYTES,
            label='System')

        self._reserve_neuron_params_data_region(spec)

        spec.reserve_memory_region(
            region=POPULATION_BASED_REGIONS.NEURON_RECORDING.value,
            size=self._app_vertex.neuron_recorder.get_exact_static_sdram_usage(
                self.vertex_slice),
            label="neuron recording")

        profile_utils.reserve_profile_region(
            spec, POPULATION_BASED_REGIONS.PROFILING.value,
            self._app_vertex.n_profile_samples)

        # reserve bit field region
        bit_field_utilities.reserve_bit_field_regions(
            spec, machine_graph, n_key_map, self,
            POPULATION_BASED_REGIONS.BIT_FIELD_BUILDER.value,
            POPULATION_BASED_REGIONS.BIT_FIELD_FILTER.value,
            POPULATION_BASED_REGIONS.BIT_FIELD_KEY_MAP.value)

        self.reserve_provenance_data_region(spec)

    @staticmethod
    def neuron_region_sdram_address(placement, transceiver):
        return helpful_functions.locate_memory_region_for_placement(
                placement, POPULATION_BASED_REGIONS.NEURON_PARAMS.value,
                transceiver)

    def _reserve_neuron_params_data_region(self, spec):
        """ Reserve the neuron parameter data region.

        :param ~data_specification.DataSpecificationGenerator spec:
            the spec to write the DSG region to
        :return: None
        """
        params_size = self._app_vertex.get_sdram_usage_for_neuron_params(
            self.vertex_slice)
        spec.reserve_memory_region(
            region=POPULATION_BASED_REGIONS.NEURON_PARAMS.value,
            size=params_size, label='NeuronParams')

    def _write_neuron_parameters(self, spec, key, region_id):

        self._app_vertex.set_has_run()

        # pylint: disable=too-many-arguments
        n_atoms = self.vertex_slice.n_atoms
        spec.comment("\nWriting Neuron Parameters for {} Neurons:\n".format(
            n_atoms))

        # Set the focus to the memory region:
        spec.switch_write_focus(region_id)

        # store the tdma data here for this slice.
        data = self._app_vertex.generate_tdma_data_specification_data(
            self._app_vertex.vertex_slices.index(self.vertex_slice))
        spec.write_array(data)

        # Write whether the key is to be used, and then the key, or 0 if it
        # isn't to be used
        if key is None:
            spec.write_value(data=0)
            spec.write_value(data=0)
        else:
            spec.write_value(data=1)
            spec.write_value(data=key)

        # Write the number of neurons in the block:
        spec.write_value(data=n_atoms)

        # Write the number of synapse types
        spec.write_value(
            data=self._app_vertex.neuron_impl.get_n_synapse_types())

        # Write the size of the incoming spike buffer
        spec.write_value(data=self._app_vertex.incoming_spike_buffer_size)

        # Write the neuron parameters
        neuron_data = self._app_vertex.neuron_impl.get_data(
            self._app_vertex.parameters, self._app_vertex.state_variables,
            self.vertex_slice)
        spec.write_array(neuron_data)

    @overrides(AbstractSynapseExpandable.gen_on_machine)
    def gen_on_machine(self):
        return self.app_vertex.synapse_manager.gen_on_machine(
            self.vertex_slice)

    @overrides(AbstractSynapseExpandable.read_generated_connection_holders)
    def read_generated_connection_holders(self, transceiver, placement):
        self._app_vertex.synapse_manager.read_generated_connection_holders(
            transceiver, placement)

    @overrides(AbstractReadParametersBeforeSet.read_parameters_from_machine)
    def read_parameters_from_machine(
            self, transceiver, placement, vertex_slice):

        # locate SDRAM address to where the neuron parameters are stored
        neuron_region_sdram_address = self.neuron_region_sdram_address(
            placement, transceiver)

        # shift past the extra stuff before neuron parameters that we don't
        # need to read
        neuron_parameters_sdram_address = (
            neuron_region_sdram_address +
            self._app_vertex.tdma_sdram_size_in_bytes +
            self._app_vertex.BYTES_TILL_START_OF_GLOBAL_PARAMETERS)

        # get size of neuron params
        size_of_region = self._app_vertex.get_sdram_usage_for_neuron_params(
            vertex_slice)
        size_of_region -= (
            self._app_vertex.BYTES_TILL_START_OF_GLOBAL_PARAMETERS +
            self._app_vertex.tdma_sdram_size_in_bytes)

        # get data from the machine
        byte_array = transceiver.read_memory(
            placement.x, placement.y, neuron_parameters_sdram_address,
            size_of_region)

        # update python neuron parameters with the data
        self._app_vertex.neuron_impl.read_data(
            byte_array, 0, vertex_slice, self._app_vertex.parameters,
            self._app_vertex.state_variables)
