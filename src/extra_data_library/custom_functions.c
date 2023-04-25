//
// Created by sachetto on 01/10/17.
//

#include <unistd.h>

#include "../config/extra_data_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../utils/file_utils.h"
#include "../domains_library/mesh_info_data.h"
#include "helper_functions.h"
#include "../domains_library/custom_mesh_info_data.h"

SET_EXTRA_DATA(set_extra_data_for_human_full_mesh) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    struct extra_data_for_fibrosis *extra_data = NULL;
    extra_data = set_common_schemia_data(config, num_active_cells);

    struct cell_node ** ac = the_grid->active_cells;

    real_cpu small_scar_center_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_x, config, "small_scar_center_x");

    real_cpu small_scar_center_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_y, config, "small_scar_center_y");

    real_cpu small_scar_center_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_z, config, "small_scar_center_z");

    real_cpu big_scar_center_x = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_x, config, "big_scar_center_x");

    real_cpu big_scar_center_y = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_y, config, "big_scar_center_y");

    real_cpu big_scar_center_z = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_z, config, "big_scar_center_z");

    real_cpu bz_size_big = 0;
    real_cpu bz_size_small = 0;
    real_cpu dist_big = 0;
    real_cpu dist_small = 0;

    uint32_t i;
    bool fibrotic, border_zone;
    int scar_type;

    OMP(parallel for private(dist_big, dist_small))
    for (i = 0; i < num_active_cells; i++) {

        border_zone = DHZB_MESH_TISSUE_TYPE(ac[i]) == BZ;
        scar_type = DHZB_MESH_LOCATION(ac[i]);

        if (ac[i]->active && border_zone) {
            real_cpu center_x = ac[i]->center.x;
            real_cpu center_y = ac[i]->center.y;
            real_cpu center_z = ac[i]->center.z;
            if(scar_type == BIG_SCAR) {
                dist_big = sqrt((center_x - big_scar_center_x) * (center_x - big_scar_center_x) +
                                (center_y - big_scar_center_y) * (center_y - big_scar_center_y) +
                                (center_z - big_scar_center_z) * (center_z - big_scar_center_z));
                OMP(critical(big))
                if (dist_big > bz_size_big) {
                    bz_size_big = dist_big;
                }
            }
            else if(scar_type == SMALL_SCAR) {
                dist_small = sqrt((center_x - small_scar_center_x) * (center_x - small_scar_center_x) +
                                  (center_y - small_scar_center_y) * (center_y - small_scar_center_y) +
                                  (center_z - small_scar_center_z) * (center_z - small_scar_center_z));
                OMP(critical(small))
                if (dist_small > bz_size_small) {
                    bz_size_small = dist_small;
                }
            }
        }
    }

    OMP(parallel for private(dist_big, dist_small))
    for (i = 0; i < num_active_cells; i++) {

        if (ac[i]->active) {
            fibrotic = DHZB_MESH_TISSUE_TYPE(ac[i]) == SCAR;
            border_zone = DHZB_MESH_TISSUE_TYPE(ac[i]) == BZ;
            scar_type = DHZB_MESH_LOCATION(ac[i]);

            if(fibrotic) {
                extra_data->fibrosis[i] = 0.0f;
            }
            else if (border_zone) {
                real_cpu center_x = ac[i]->center.x;
                real_cpu center_y = ac[i]->center.y;
                real_cpu center_z = ac[i]->center.z;
                if(scar_type == BIG_SCAR) {
                    dist_big = sqrt((center_x - big_scar_center_x) * (center_x - big_scar_center_x) +
                                    (center_y - big_scar_center_y) * (center_y - big_scar_center_y) +
                                    (center_z - big_scar_center_z) * (center_z - big_scar_center_z));

                    extra_data->fibrosis[i] = (real)(dist_big / bz_size_big);

                }
                else if(scar_type == SMALL_SCAR) {
                    dist_small = sqrt((center_x - small_scar_center_x) * (center_x - small_scar_center_x) +
                                      (center_y - small_scar_center_y) * (center_y - small_scar_center_y) +
                                      (center_z - small_scar_center_z) * (center_z - small_scar_center_z));

                    extra_data->fibrosis[i] = (real)(dist_small / bz_size_small);
                }
            }
            else {
                extra_data->fibrosis[i] = 1.0f;
            }
        }
    }

    SET_EXTRA_DATA_SIZE(sizeof(struct extra_data_for_fibrosis));

    return (void*)extra_data;
}

SET_EXTRA_DATA(set_extra_data_for_scar_wedge) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    struct extra_data_for_fibrosis *extra_data = NULL;

    extra_data = set_common_schemia_data(config, num_active_cells);

    struct cell_node ** ac = the_grid->active_cells;

    char *scar_size;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR (scar_size, config, "scar_size");

    uint8_t size_code;

    if(strcmp(scar_size, "big") == 0) {
        size_code = 0;
    }
    else if(strcmp(scar_size, "small") == 0) {
        size_code = 1;
    }
    else {
        printf("Function: set_extra_data_for_scar_edge, invalid scar size %s. Valid sizes are big or small. Exiting!\n", scar_size);
        exit(EXIT_FAILURE);
    }

    real_cpu scar_center_x;
    real_cpu scar_center_y;
    real_cpu scar_center_z;

    ////Fibrosis configuration
    //BIG SCAR
    if(size_code == 0) {
        scar_center_x = 95300;
        scar_center_y = 81600;
        scar_center_z = 36800;
    }
    else {
        scar_center_x = 52469;
        scar_center_y = 83225;
        scar_center_z = 24791;
    }
    real_cpu bz_size = 0.0;
    real_cpu dist;

    uint32_t i;
    bool border_zone, fibrotic;

    OMP(parallel for private(dist))
    for (i = 0; i < num_active_cells; i++) {
        if(ac[i]->active) {
            border_zone = DHZB_MESH_TISSUE_TYPE(ac[i]) == BZ;
            if(border_zone) {
                real_cpu center_x = ac[i]->center.x;
                real_cpu center_y = ac[i]->center.y;
                real_cpu center_z = ac[i]->center.z;
                dist =  sqrt((center_x - scar_center_x)*(center_x - scar_center_x) + (center_y - scar_center_y)*(center_y - scar_center_y)  + (center_z - scar_center_z)*(center_z - scar_center_z)  );
                OMP(critical)
                if(dist > bz_size) {
                    bz_size = dist;
                }
            }

        }
    }

    OMP(parallel for private(dist))
    for (i = 0; i < num_active_cells; i++) {

        if(ac[i]->active) {

            fibrotic = DHZB_MESH_TISSUE_TYPE(ac[i]) == SCAR;
            border_zone = DHZB_MESH_TISSUE_TYPE(ac[i]) == BZ;

            if(fibrotic) {
                extra_data->fibrosis[i] = 0.0;
            }
            else if(border_zone) {
                real_cpu center_x = ac[i]->center.x;
                real_cpu center_y = ac[i]->center.y;
                real_cpu center_z = ac[i]->center.z;
                dist =  sqrt((center_x - scar_center_x)*(center_x - scar_center_x) + (center_y - scar_center_y)*(center_y - scar_center_y)  + (center_z - scar_center_z)*(center_z - scar_center_z)  );
                dist = dist/bz_size;

                extra_data->fibrosis[i] = (real)dist;

            }
            else {
                extra_data->fibrosis[i] = 1.0f;
            }

        }
    }

    SET_EXTRA_DATA_SIZE(sizeof(struct extra_data_for_fibrosis));

    return (void*)extra_data;
}

SET_EXTRA_DATA(set_extra_data_for_scv_mesh) {

    uint32_t num_active_cells = the_grid->num_active_cells;

    *extra_data_size = sizeof(uint32_t)*(num_active_cells);

    uint32_t *mapping = (uint32_t*)malloc(*extra_data_size);

    struct cell_node ** ac = the_grid->active_cells;

    int i;

    OMP(parallel for)
    for (i = 0; i < num_active_cells; i++) {
        mapping[i] = TISSUE_TYPE(ac[i]); //endo 0; mid  1; epi  2
    }

    return (void*)mapping;
}

SET_EXTRA_DATA(set_extra_data_for_scv_mesh_with_torord) {

    uint32_t num_eq = 43;   // ToRORd number of equations
    uint32_t num_active_cells = the_grid->num_active_cells;

    // The percentages were taken from the ToRORd paper (Transmural experiment)
    //real side_length_endo = side_length*0.45;
    //real side_length_mid = side_length_endo + side_length*0.25;
    //real side_length_epi = side_length_mid + side_length*0.3;

    // The extra data size is the initial state vector for each celltype plus the mapping of each cell
    *extra_data_size = sizeof(real)*(num_active_cells) + sizeof(real)*num_eq*3;

    real *extra_data = (real*)malloc(*extra_data_size);

    // Set the initial conditions (celltype=ENDO)
    int offset = 0;
    extra_data[0] = -88.7638;
    extra_data[1] = 0.0111;
    extra_data[2] = 7.0305e-5;
    extra_data[3] = 12.1025;
    extra_data[4] = 12.1029;
    extra_data[5] = 142.3002;
    extra_data[6] = 142.3002;
    extra_data[7] = 1.5211;
    extra_data[8] = 1.5214;
    extra_data[9] = 8.1583e-05;
    extra_data[10] = 8.0572e-4;
    extra_data[11] = 0.8286;
    extra_data[12] = 0.8284;
    extra_data[13] = 0.6707;
    extra_data[14] = 0.8281;
    extra_data[15] = 1.629e-4;
    extra_data[16] = 0.5255;
    extra_data[17] = 0.2872;
    extra_data[18] = 9.5098e-4;
    extra_data[19] = 0.9996;
    extra_data[20] = 0.5936;
    extra_data[21] = 4.8454e-4;
    extra_data[22] = 0.9996;
    extra_data[23] = 0.6538;
    extra_data[24] = 8.1084e-9;
    extra_data[25] = 1.0;
    extra_data[26] = 0.939;
    extra_data[27] = 1.0;
    extra_data[28] = 0.9999;
    extra_data[29] = 1.0;
    extra_data[30] = 1.0;
    extra_data[31] = 1.0;
    extra_data[32] = 6.6462e-4;
    extra_data[33] = 0.0012;
    extra_data[34] = 7.0344e-4;
    extra_data[35] = 8.5109e-4;
    extra_data[36] = 0.9981;
    extra_data[37] = 1.3289e-5;
    extra_data[38] = 3.7585e-4;
    extra_data[39] = 0.248;
    extra_data[40] = 1.7707e-4;
    extra_data[41] = 1.6129e-22;
    extra_data[42] = 1.2475e-20;

    // Set the initial conditions (celltype=EPI)
    offset = num_eq;
    extra_data[0+offset] = -89.1400;
    extra_data[1+offset] = 0.0129;
    extra_data[2+offset] = 5.7672e-05;
    extra_data[3+offset] = 12.8363;
    extra_data[4+offset] = 12.8366;
    extra_data[5+offset] = 142.6951;
    extra_data[6+offset] = 142.6951;
    extra_data[7+offset] = 1.8119;
    extra_data[8+offset] = 1.8102;
    extra_data[9+offset] = 6.6309e-05;
    extra_data[10+offset] = 0.00074303;
    extra_data[11+offset] = 0.8360;
    extra_data[12+offset] = 0.8359;
    extra_data[13+offset] = 0.6828;
    extra_data[14+offset] = 0.8357;
    extra_data[15+offset] = 0.00015166;
    extra_data[16+offset] = 0.5401;
    extra_data[17+offset] = 0.3034;
    extra_data[18+offset] = 0.00092716;
    extra_data[19+offset] = 0.9996;
    extra_data[20+offset] = 0.9996;
    extra_data[21+offset] = 0.0004724;
    extra_data[22+offset] = 0.9996;
    extra_data[23+offset] = 0.9996;
    extra_data[24+offset] = 0;
    extra_data[25+offset] = 1.0;
    extra_data[26+offset] = 0.9485;
    extra_data[27+offset] = 1.0;
    extra_data[28+offset] = 0.9999;
    extra_data[29+offset] = 1.0;
    extra_data[30+offset] = 1.0;
    extra_data[31+offset] = 1.0;
    extra_data[32+offset] = 0.00030853;
    extra_data[33+offset] = 0.00053006;
    extra_data[34+offset] = 0.00067941;
    extra_data[35+offset] = 0.00082869;
    extra_data[36+offset] = 0.9982;
    extra_data[37+offset] = 9.5416e-06;
    extra_data[38+offset] = 0.00027561;
    extra_data[39+offset] = 0.2309;
    extra_data[40+offset] = 0.00016975;
    extra_data[41+offset] = 2.8189e-24;
    extra_data[42+offset] = 0;

    // Set the initial conditions (celltype=MID)
    offset = num_eq*2;
    extra_data[0+offset] = -89.1704;
    extra_data[1+offset] = 0.0192;
    extra_data[2+offset] = 6.5781e-05;
    extra_data[3+offset] = 15.0038;
    extra_data[4+offset] = 15.0043;
    extra_data[5+offset] = 143.0403;
    extra_data[6+offset] = 143.0402;
    extra_data[7+offset] = 1.9557;
    extra_data[8+offset] = 1.9593;
    extra_data[9+offset] = 8.166e-05;
    extra_data[10+offset] = 0.00073818;
    extra_data[11+offset] = 0.8365;
    extra_data[12+offset] = 0.8363;
    extra_data[13+offset] = 0.6838;
    extra_data[14+offset] = 0.8358;
    extra_data[15+offset] = 0.00015079;
    extra_data[16+offset] = 0.5327;
    extra_data[17+offset] = 0.2834;
    extra_data[18+offset] = 0.00092527;
    extra_data[19+offset] = 0.9996;
    extra_data[20+offset] = 0.5671;
    extra_data[21+offset] = 0.00047143;
    extra_data[22+offset] = 0.9996;
    extra_data[23+offset] = 0.6261;
    extra_data[24+offset] = 0;
    extra_data[25+offset] = 1.0;
    extra_data[26+offset] = 0.92;
    extra_data[27+offset] = 1.0;
    extra_data[28+offset] = 0.9998;
    extra_data[29+offset] = 1.0;
    extra_data[30+offset] = 1.0;
    extra_data[31+offset] = 1.0;
    extra_data[32+offset] = 0.00051399;
    extra_data[33+offset] = 0.0012;
    extra_data[34+offset] = 0.00069560;
    extra_data[35+offset] = 0.00082672;
    extra_data[36+offset] = 0.9979;
    extra_data[37+offset] = 1.8784e-05;
    extra_data[38+offset] = 0.00054206;
    extra_data[39+offset] = 0.2653;
    extra_data[40+offset] = 0.00016921;
    extra_data[41+offset] = 0;
    extra_data[42+offset] = 0;

    offset = num_eq*3;

    struct cell_node ** ac = the_grid->active_cells;

    int i;

    OMP(parallel for)
    for (i = 0; i < num_active_cells; i++)
    {
        real tag = TISSUE_TYPE(ac[i]); // Mesh tags: endo 0; mid  1; epi  2

        // Change to the ToRORd tag
        if (tag == 1) tag = 2;
        else if (tag == 2) tag = 1;

        extra_data[i+offset] = tag;
    }

    return (void*)extra_data;

}
