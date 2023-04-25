//
// Created by sachetto on 13/10/17.
//

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "../alg/grid/grid.h"
#include "../config/stim_config.h"
#include "../config_helpers/config_helpers.h"

SET_SPATIAL_STIM(stim_base_mouse) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    real_cpu stim_size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, stim_size, config, "stim_size");

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config, "current");

    real stim_value;

    uint32_t i;

    ALLOCATE_STIMS();

    OMP(parallel for private(stim_value))
    for(i = 0; i < n_active; i++) {

        bool stim;
        stim = (ac[i]->center.x >= 3000.0 - stim_size) && (ac[i]->center.x <= 3000.0 + stim_size);
        stim &= (ac[i]->center.y >= 2400.0 - stim_size) && (ac[i]->center.y <= 2400.0 + stim_size);
        stim &= (ac[i]->center.z >= 300 - stim_size) && (ac[i]->center.z <= 300 + stim_size);

        if(stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }
}

SET_SPATIAL_STIM(stim_mouse_spiral) {

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

    real stim_current = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, stim_current, config, "current");

    real stim_value;

    int i;

    ALLOCATE_STIMS();

    OMP(parallel for private(stim_value))
    for(i = 0; i < n_active; i++) {

        bool stim;

        stim = (ac[i]->center.x >= 3000.0) && (ac[i]->center.x <= 6000.0);
        stim &= (ac[i]->center.y >= 1940.0) && (ac[i]->center.y <= 6100.0);
        stim &= (ac[i]->center.z >= 2230.0) && (ac[i]->center.z <= 5800.0);

        if(stim) {
            stim_value = stim_current;
        } else {
            stim_value = 0.0;
        }

        SET_STIM_VALUE(i, stim_value);
    }
}
