#include "connection_generator.h"
#include <spin1_api.h>
#include <debug.h>

#define N_CONNECTION_GENERATORS 0

struct connection_generator {
    uint32_t index;
    void *data;
};

struct connection_generator_info {
    uint32_t hash;
    void* (*initialize)(address_t *region);
    uint32_t (*generate)(
        void *data,  uint32_t pre_slice_start, uint32_t pre_slice_end,
        uint32_t pre_neuron_index, uint32_t post_slice_start,
        uint32_t post_slice_count, uint32_t max_row_length, rng_t rng,
        uint16_t *indices);
    void (*free)(void *data);
};

struct connection_generator_info connection_generators[N_CONNECTION_GENERATORS];

void register_connection_generators() {

    // TODO: Fill in the generators
}

connection_generator_t connection_generator_init(
        uint32_t hash, address_t *in_region) {
    for (uint32_t i = 0; i < N_CONNECTION_GENERATORS; i++) {
        if (hash == connection_generators[i].hash) {

            address_t region = *in_region;
            connection_generator_t generator = spin1_malloc(
                sizeof(connection_generator_t));
            if (generator == NULL) {
                log_error("Could not create generator");
                return NULL;
            }
            generator->index = i;
            generator->data = connection_generators[i].initialize(&region);
            *in_region = region;
            return generator;
        }
    }
    log_error("Connection generator with hash %u not found", hash);
    return NULL;
}

uint32_t connection_generator_generate(
        connection_generator_t generator, uint32_t pre_slice_start,
        uint32_t pre_slice_count, uint32_t pre_neuron_index,
        uint32_t post_slice_start, uint32_t post_slice_count,
        uint32_t max_row_length, rng_t rng, uint16_t *indices) {
    return connection_generators[generator->index].generate(
        generator->data, pre_slice_start, pre_slice_count,
        pre_neuron_index, post_slice_start, post_slice_count,
        max_row_length, rng, indices);
}

void connection_generator_free(connection_generator_t generator) {
    connection_generators[generator->index].free(generator->data);
    sark_free(generator);
}
