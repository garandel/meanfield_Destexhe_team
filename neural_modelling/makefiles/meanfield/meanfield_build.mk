# See Notes in sPyNNaker/neural_modelling/CHANGES_April_2018

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

# If SPINN_DIRS is not defined, this is an error!
ifndef SPINN_DIRS
    $(error SPINN_DIRS is not set.  Please define SPINN_DIRS (possibly by running "source setup" in the spinnaker package folder))
endif

# If NEURAL_MODELLING_DIRS is not defined, this is an error!
ifndef NEURAL_MODELLING_DIRS
    $(error NEURAL_MODELLING_DIRS is not set.  Please define NEURAL_MODELLING_DIRS (possibly by running "source setup" in the sPyNNaker folder))
endif
#Check NEURAL_MODELLING_DIRS
MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
CHECK_PATH := $(NEURAL_MODELLING_DIRS)/makefiles/meanfield/meanfield_build.mk
ifneq ($(CHECK_PATH), $(MAKEFILE_PATH))
    $(error Please check NEURAL_MODELLING_DIRS as based on that this file is at $(CHECK_PATH) when it is actually at $(MAKEFILE_PATH))
endif

# Set logging levels
ifeq ($(SPYNNAKER_DEBUG), DEBUG)
    NEURON_DEBUG = LOG_DEBUG
endif

ifndef NEURON_DEBUG
    NEURON_DEBUG = LOG_INFO
endif

# Add source directory

# Define the directories
# Path flag to replace with the modified dir  (abspath drops the final /)
MEANFIELD_DIR := $(abspath $(NEURAL_MODELLING_DIRS)/src)
MODIFIED_DIR :=$(dir $(abspath $(MEANFIELD_DIR)))modified_src/
SOURCE_DIRS += $(MEANFIELD_DIR)

# Define a rule to find the source directory of the given file.
# This attempts to find each of SOURCE_DIRS within the given file name; the
# first one that matches is then returned.  If none match, an empty string
# will be returned.
define get_source_dir#(file)
$(firstword $(strip $(foreach d, $(sort $(SOURCE_DIRS)), $(findstring $(d), $(1)))))
endef

# Define rule to strip any SOURCE_DIRS from source_file to allow use via local.mk.
# If no match is found, the value is returned untouched 
# (though this will probably fail later).
define strip_source_dirs#(source_file)
$(or $(patsubst $(call get_source_dir, $(1))/%,%,$(1)), $(1))
endef

# Define a rule to replace any SOURCE_DIRS from header_file with the modified_src folder.
define replace_source_dirs#(header_file)
$(patsubst $(call get_source_dir, $(1))%, $(dir $(call get_source_dir, $(1)))modified_src%, $(1))
endef

# Need to build each neuron seperately or complier gets confused
# BUILD_DIR and APP_OUTPUT_DIR end with a / for historictical/ shared reasons
ifndef BUILD_DIR
    BUILD_DIR := $(NEURAL_MODELLING_DIRS)/builds/$(APP)/
endif
ifndef APP_OUTPUT_DIR
    APP_OUTPUT_DIR :=  $(NEURAL_MODELLING_DIRS)/../spynnaker/pyNN/model_binaries
endif

# Check if the neuron implementation is the default one
ifndef MEANFIELD_IMPL_H
    $(error MEANFIELD_IMPL_H is not set.  Please select a neuron implementation)
else
    MEANFIELD_IMPL := $(call strip_source_dirs,$(MEANFIELD_IMPL_H))
    MEANFIELD_IMPL_H := $(call replace_source_dirs,$(MEANFIELD_IMPL_H))
    MEANFIELD_IMPL_STANDARD := meanfield/implementations/meanfield_impl_standard.h
    MEANFIELD_INCLUDES := -include $(MEANFIELD_IMPL_H)
    ifeq ($(MEANFIELD_IMPL), $(MEANFIELD_IMPL_STANDARD))
        
        # Check required inputs and point them to modified sources
		ifndef ADDITIONAL_INPUT_H
		    ADDITIONAL_INPUT_H = $(MODIFIED_DIR)meanfield/additional_inputs/additional_input_none_impl.h
		else
		    ADDITIONAL_INPUT_H := $(call replace_source_dirs,$(ADDITIONAL_INPUT_H))
		endif
		
		ifndef MEANFIELD_MODEL
		    $(error MEANFIELD_MODEL is not set.  Please choose a neuron model to compile)
		else
		    MEANFIELD_MODEL := $(call strip_source_dirs,$(MEANFIELD_MODEL))
		endif
		
		ifndef MEANFIELD_MODEL_H
		    $(error MEANFIELD_MODEL is not set.  Please choose a neuron model to compile)
		else
		    MEANFIELD_MODEL_H := $(call replace_source_dirs,$(MEANFIELD_MODEL_H))
		endif
		
		
		ifndef PARAMS_FROM_NETWORK_H
		    $(error PARAMS_FROM_NETWORK_H is not set.  Please select an input type header file)
		else
		    PARAMS_FROM_NETWORK_H := $(call replace_source_dirs,$(PARAMS_FROM_NETWORK_H))
		endif
		
		ifndef P_FIT_POLYNOMIAL_H
		    $(error P_FIT_POLYNOMIAL_H is not set.  Please select an input type header file)
		else
		    P_FIT_POLYNOMIAL_H := $(call replace_source_dirs,$(P_FIT_POLYNOMIAL_H))
		endif
		
		ifndef INPUT_TYPE_H
		    $(error INPUT_TYPE_H is not set.  Please select an input type header file)
		else
		    INPUT_TYPE_H := $(call replace_source_dirs,$(INPUT_TYPE_H))
		endif
		
#		ifndef MATHSBOX_H
#		    $(error MATHSBOX_H is not set.  Please select an input type header file)
#		else
#		    MATHSBOX_H := $(call replace_source_dirs,$(MATHSBOX_H))
#		endif
		        
		ifndef THRESHOLD_TYPE_H
		    $(error THRESHOLD_TYPE_H is not set.  Please select a threshold type header file)
		else
		    THRESHOLD_TYPE_H := $(call replace_source_dirs,$(THRESHOLD_TYPE_H))
		endif
		
		ifndef SYNAPSE_TYPE_H
		    $(error SYNAPSE_TYPE_H is not set.  Please select a synapse type header file)
		else
		    SYNAPSE_TYPE_H := $(call replace_source_dirs,$(SYNAPSE_TYPE_H))
		endif
		
		MEANFIELD_INCLUDES := \
	      -include $(MEANFIELD_MODEL_H) \
	      -include $(SYNAPSE_TYPE_H) \
	      -include $(PARAMS_FROM_NETWORK_H)\
          -include $(P_FIT_POLYNOMIAL_H)\
          -include $(MATHSBOX_H) \
	      -include $(INPUT_TYPE_H) \
	      -include $(THRESHOLD_TYPE_H) \
	      -include $(ADDITIONAL_INPUT_H) \
	      -include $(MEANFIELD_IMPL_H)
    endif
endif

OTHER_SOURCES_CONVERTED := $(call strip_source_dirs,$(OTHER_SOURCES))

# List all the sources relative to one of SOURCE_DIRS
SOURCES = meanfield/c_main_neurons.c \
          meanfield/meanfield.c \
          meanfield/meanfield_recording.c \
          meanfield/plasticity/synapse_dynamics_remote.c\
          $(MEANFIELD_MODEL) $(OTHER_SOURCES_CONVERTED)

include $(SPINN_DIRS)/make/local.mk

FEC_OPT = $(OTIME)

TEST_COMPILE = $(CC) -DLOG_LEVEL=$(NEURON_DEBUG) $(CFLAGS)

$(BUILD_DIR)meanfield/c_main_neurons.o: $(MODIFIED_DIR)meanfield/c_main_neurons.c
	#c_main_neurons.c
	-@mkdir -p $(dir $@)
	$(TEST_COMPILE) -o $@ $<
    
#$(BUILD_DIR)meanfield/meanfield.o: $(MODIFIED_DIR)meanfield/meanfield.c
#	#meanfield.c
#	-@mkdir -p $(dir $@)
#	$(TEST_COMPILE) -o $@ $<

$(BUILD_DIR)meanfield/meanfield.o: $(MODIFIED_DIR)meanfield/meanfield.c
	# meanfield.o
	-@mkdir -p $(dir $@)
	$(CC) -DLOG_LEVEL=$(NEURON_DEBUG) $(CFLAGS) $(MEANFIELD_INCLUDES) -o $@ $<

$(BUILD_DIR)meanfield/meanfield_recording.o: $(MODIFIED_DIR)meanfield/meanfield_recording.c
	# meanfield_recording.o
	-@mkdir -p $(dir $@)
	$(CC) -DLOG_LEVEL=$(NEURON_DEBUG) $(CFLAGS) $(MEANFIELD_INCLUDES) -o $@ $<

.PRECIOUS: $(MODIFIED_DIR)%.c $(MODIFIED_DIR)%.h $(LOG_DICT_FILE) $(EXTRA_PRECIOUS)
