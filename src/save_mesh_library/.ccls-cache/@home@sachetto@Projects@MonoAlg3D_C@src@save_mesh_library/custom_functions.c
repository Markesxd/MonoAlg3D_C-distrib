#include <stdbool.h>
#include <stdlib.h>

#include "../3dparty/sds/sds.h"
#include "../alg/grid/grid.h"
#include "../config/save_mesh_config.h"
#include "../utils/utils.h"

#include "../vtk_utils/vtk_unstructured_grid.h"
#include "../vtk_utils/vtk_polydata_grid.h"
#include "../libraries_common/common_data_structures.h"
#include "../domains_library/mesh_info_data.h"

#include "save_mesh_helper.h"

#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif

static char *file_prefix;
static bool binary = false;
static bool clip_with_plain = false;
static bool clip_with_bounds = false;
static bool save_pvd = true;
static bool save_inactive = false;
static bool compress = false;
static bool save_f = false;
static int compression_level = 3;
char *output_dir;
bool save_visible_mask = true;

static bool initialized = false;

void calculate_activation_time_and_apd (struct time_info *time_info, struct config *config, struct grid *the_grid, const real_cpu time_threshold,\
                                            const real_cpu activation_threshold, const real_cpu apd_threshold);
void write_activation_time_maps (struct config *config, struct grid *the_grid, char *output_dir,\
                                    char *file_prefix, bool clip_with_plain, bool clip_with_bounds, bool binary,\
                                    bool save_pvd, bool compress, int compression_level, bool save_f);
void write_tissue_activation_time_for_each_pulse (struct config *config, struct grid *the_grid,\
                                            char *output_dir, float plain_coords[], float bounds[],\
                                            bool clip_with_plain, bool clip_with_bounds, bool binary, bool save_pvd, bool compress, int compression_level, bool save_f);
void set_tissue_vtk_values_with_activation_time_from_current_pulse (void **persistent_data, struct grid *the_grid, const int cur_pulse);
void write_tissue_apd_map (struct config *config, struct grid *the_grid, char *output_dir,\
                                    char *file_prefix, bool clip_with_plain, bool clip_with_bounds, bool binary,\
                                    bool save_pvd, bool compress, int compression_level, bool save_f);
void set_tissue_vtk_values_with_mean_apd (void **persistent_data, struct grid *the_grid);

static void save_visibility_mask(sds output_dir_with_file, ui8_array visible_cells) {
        sds output_dir_with_new_file = sdsnew(output_dir_with_file);
        output_dir_with_new_file = sdscat(output_dir_with_new_file, ".vis");
		FILE *vis = fopen(output_dir_with_new_file, "wb");
		fwrite(visible_cells, sizeof(uint8_t), arrlen(visible_cells), vis);
        sdsfree(output_dir_with_new_file);
		fclose(vis);
}

struct save_with_activation_times_persistent_data {
    struct vtk_unstructured_grid *grid;
    struct point_hash_entry *last_time_v;
    struct point_hash_entry *num_activations;
    struct point_hash_entry *cell_was_active;
    struct point_voidp_hash_entry *activation_times;
    struct point_voidp_hash_entry *apds;
    bool first_save_call;

};

INIT_SAVE_MESH(init_save_with_activation_times_and_apd) {

    config->persistent_data = calloc(1, sizeof(struct save_with_activation_times_persistent_data));
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->cell_was_active, 0.0);
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->last_time_v, -100.0);
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->num_activations, 0);
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->activation_times, NULL);
    hmdefault(((struct save_with_activation_times_persistent_data*)config->persistent_data)->apds, NULL);

    ((struct save_with_activation_times_persistent_data*)config->persistent_data)->grid = NULL;
    ((struct save_with_activation_times_persistent_data*)config->persistent_data)->first_save_call = true;

}

END_SAVE_MESH(end_save_with_activation_times_and_apd) {
    free(config->persistent_data);
    write_activation_time_maps(config,the_grid,output_dir,file_prefix,clip_with_plain,clip_with_bounds,binary,save_pvd,compress,compression_level,save_f);
    write_tissue_apd_map(config,the_grid,output_dir,file_prefix,clip_with_plain,clip_with_bounds,binary,save_pvd,compress,compression_level,save_f);
}

SAVE_MESH(save_with_activation_times_and_apd) {

    int mesh_output_pr = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, mesh_output_pr, config, "mesh_print_rate");

    int iteration_count = time_info->iteration;

    if(mesh_output_pr) {
        if (iteration_count % mesh_output_pr == 0)
            save_as_text_or_binary(time_info, config, the_grid, ode_solver);
    }

    float time_threshold = 10.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, time_threshold, config, "time_threshold");

    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");

    float activation_threshold = -30.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, activation_threshold, config, "activation_threshold");

    float apd_threshold = -83.0f;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(float, apd_threshold, config, "apd_threshold");

    calculate_activation_time_and_apd(time_info,config,the_grid,time_threshold,activation_threshold,apd_threshold);

}

void calculate_activation_time_and_apd (struct time_info *time_info, struct config *config, struct grid *the_grid, const real_cpu time_threshold,\
                                            const real_cpu activation_threshold, const real_cpu apd_threshold) {
	real_cpu current_t = time_info->current_t;
    real_cpu last_t = time_info->final_t;
    real_cpu dt = time_info->dt;

    struct save_with_activation_times_persistent_data *persistent_data = (struct save_with_activation_times_persistent_data*)config->persistent_data;

    struct cell_node *grid_cell = the_grid->first_cell;

    real_cpu center_x, center_y, center_z, dx, dy, dz;
    real_cpu v;

    while(grid_cell != 0) {

        if( grid_cell->active || ( grid_cell->mesh_extra_info && ( FIBROTIC(grid_cell) || BORDER_ZONE(grid_cell) ) ) ) {

            center_x = grid_cell->center.x;
            center_y = grid_cell->center.y;
            center_z = grid_cell->center.z;

            v = grid_cell->v;

            struct point_3d cell_coordinates;
            cell_coordinates.x = center_x;
            cell_coordinates.y = center_y;
            cell_coordinates.z = center_z;

            //dx = grid_cell->discretization.x / 2.0;
            //dy = grid_cell->discretization.y / 2.0;
            //dz = grid_cell->discretization.z / 2.0;

            int n_activations = 0;
            float *apds_array = NULL;
            float *activation_times_array = NULL;

            if(grid_cell->active) {
                float last_v = hmget(persistent_data->last_time_v, cell_coordinates);

                n_activations = (int) hmget(persistent_data->num_activations, cell_coordinates);
                activation_times_array = (float *) hmget(persistent_data->activation_times, cell_coordinates);
                apds_array = (float *) hmget(persistent_data->apds, cell_coordinates);

                int act_times_len = arrlen(activation_times_array);

                if (current_t == 0.0f) {
                    hmput(persistent_data->last_time_v, cell_coordinates, v);
                } else {
                    if ((last_v < activation_threshold) && (v >= activation_threshold)) {

                        if (act_times_len == 0) {
                            n_activations++;
                            hmput(persistent_data->num_activations, cell_coordinates, n_activations);
                                    arrput(activation_times_array, current_t);
                            float tmp = hmget(persistent_data->cell_was_active, cell_coordinates);
                            hmput(persistent_data->cell_was_active, cell_coordinates, tmp + 1);
                            hmput(persistent_data->activation_times, cell_coordinates, activation_times_array);
                        } else { //This is to avoid spikes in the middle of an Action Potential
                            float last_act_time = activation_times_array[act_times_len - 1];
                            if (current_t - last_act_time > time_threshold) {
                                n_activations++;
                                hmput(persistent_data->num_activations, cell_coordinates, n_activations);
                                        arrput(activation_times_array, current_t);
                                float tmp = hmget(persistent_data->cell_was_active, cell_coordinates);
                                hmput(persistent_data->cell_was_active, cell_coordinates, tmp + 1);
                                hmput(persistent_data->activation_times, cell_coordinates, activation_times_array);
                            }
                        }
                    }

                    //CHECK APD
                    bool was_active = (hmget(persistent_data->cell_was_active, cell_coordinates) != 0.0);
                    if (was_active) {
                        if (v <= apd_threshold || (hmget(persistent_data->cell_was_active, cell_coordinates) == 2.0) || (last_t-current_t) <= dt) {

                            int tmp = (int)hmget(persistent_data->cell_was_active, cell_coordinates);
                            int act_time_array_len = arrlen(activation_times_array);
                            //if this in being calculated because we had a new activation before the cell achieved the rest potential,
                            // we need to get the activation before this one
                            real_cpu last_act_time = activation_times_array[act_time_array_len  - tmp];
                            real_cpu apd = current_t - last_act_time;
                                    arrput(apds_array, apd);
                            hmput(persistent_data->apds, cell_coordinates, apds_array);
                            hmput(persistent_data->cell_was_active, cell_coordinates, tmp - 1);
                        }
                    }

                    hmput(persistent_data->last_time_v, cell_coordinates, v);
                }
            }
        }

        grid_cell = grid_cell->next;
    }
}// end of function

void write_activation_time_maps (struct config *config, struct grid *the_grid, char *output_dir,\
                                    char *file_prefix, bool clip_with_plain, bool clip_with_bounds, bool binary,\
                                    bool save_pvd, bool compress, int compression_level, bool save_f) {
	real_cpu current_t = 0.0;
    if(((struct save_with_activation_times_persistent_data *) config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_f, config, "save_f");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_visible_mask, config, "save_visible_mask");

        ((struct save_with_activation_times_persistent_data *) config->persistent_data)->first_save_call = false;
    }

    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[0], config, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[1], config, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[2], config, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[3], config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[4], config, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[5], config, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[0], config, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[1], config, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[2], config, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[3], config, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[4], config, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[5], config, "max_z");
    }

    write_activation_time_for_each_pulse(config,the_grid,output_dir,plain_coords,bounds,clip_with_plain,clip_with_bounds,binary,save_pvd,compress,compression_level,save_f);
}

void write_tissue_activation_time_for_each_pulse (struct config *config, struct grid *the_grid,\
                                            char *output_dir, float plain_coords[], float bounds[],\
                                            bool clip_with_plain, bool clip_with_bounds, bool binary, bool save_pvd, bool compress, int compression_level, bool save_f) {
	struct save_with_activation_times_persistent_data **data = (struct save_with_activation_times_persistent_data **)&config->persistent_data;
	struct cell_node **grid_cell = the_grid->active_cells;

	real_cpu center_x, center_y, center_z;

	center_x = grid_cell[0]->center.x;
    center_y = grid_cell[0]->center.y;
    center_z = grid_cell[0]->center.z;

    struct point_3d cell_coordinates;
    cell_coordinates.x = center_x;
    cell_coordinates.y = center_y;
    cell_coordinates.z = center_z;

	// Get the number of pulses using one cell of the grid
    int n_pulses = (int) hmget((*data)->num_activations, cell_coordinates);

    // Write the activation time map for each pulse
    for (int cur_pulse = 0; cur_pulse < n_pulses; cur_pulse++){
    	sds base_name = create_base_name("Activation_time_map_pulse", cur_pulse, "vtk");
    	ds output_dir_with_file = sdsnew(output_dir);
        output_dir_with_file = sdscat(output_dir_with_file, "/");
        output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, cur_pulse);

        bool read_only_data = ((struct save_with_activation_times_persistent_data *) config->persistent_data)->grid != NULL;
        new_vtk_unstructured_grid_from_alg_grid(&(((struct save_with_activation_times_persistent_data *) config->persistent_data)->grid), the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data, save_f);

		set_tissue_vtk_values_with_activation_time_from_current_pulse(&config->persistent_data,the_grid,cur_pulse);

        save_vtk_unstructured_grid_as_legacy_vtk(((struct save_with_activation_times_persistent_data *) config->persistent_data)->grid, output_dir_with_file, binary, save_f);

	    if(save_visible_mask) {
			save_visibility_mask(output_dir_with_file, (((struct save_with_activation_times_persistent_data *) config->persistent_data)->grid)->cell_visibility);
		}

	    if(the_grid->adaptive) {
	        free_vtk_unstructured_grid(((struct save_with_activation_times_persistent_data *) config->persistent_data)->grid);
	        ((struct save_with_activation_times_persistent_data *) config->persistent_data)->grid = NULL;
	    }

	    sdsfree(output_dir_with_file);
	    sdsfree(base_name);
    }
}

void set_tissue_vtk_values_with_activation_time_from_current_pulse (void **persistent_data, struct grid *the_grid, const int cur_pulse) {
	struct save_with_activation_times_persistent_data **data = (struct save_with_activation_times_persistent_data **)&config->persistent_data;

	struct cell_node **grid_cell = the_grid->active_cells;
    uint32_t num_active_cells =  the_grid->num_active_cells;
    real_cpu center_x, center_y, center_z;

    // We need to pass through the cells in the same order as we did when the VTU grid was created
    for(int i = 0; i < num_active_cells; i++) {
    	center_x = grid_cell[i]->center.x;
        center_y = grid_cell[i]->center.y;
        center_z = grid_cell[i]->center.z;

        struct point_3d cell_coordinates;
        cell_coordinates.x = center_x;
        cell_coordinates.y = center_y;
        cell_coordinates.z = center_z;

        float *activation_times_array = NULL;
        if(grid_cell[i]->active) {

        	activation_times_array = (float *) hmget((*data)->activation_times, cell_coordinates);

            float at;
            if (activation_times_array == NULL)
                at = -1;
            else
                at = activation_times_array[cur_pulse];

            // Update the scalar value from the "vtk_unstructured_grid"
            (*data)->grid->values[i] = at;
        }
    }

}

void write_tissue_apd_map (struct config *config, struct grid *the_grid, char *output_dir,\
                                    char *file_prefix, bool clip_with_plain, bool clip_with_bounds, bool binary,\
                                    bool save_pvd, bool compress, int compression_level, bool save_f){

	real_cpu current_t = 0.0;
    if(((struct save_with_activation_times_persistent_data *) config->persistent_data)->first_save_call) {
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(output_dir, config, "output_dir");
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(file_prefix, config, "file_prefix");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_plain, config, "clip_with_plain");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(clip_with_bounds, config, "clip_with_bounds");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(binary, config, "binary");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_f, config, "save_f");
        GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(save_visible_mask, config, "save_visible_mask");

        ((struct save_with_activation_times_persistent_data *) config->persistent_data)->first_save_call = false;
    }

    float plain_coords[6] = {0, 0, 0, 0, 0, 0};
    float bounds[6] = {0, 0, 0, 0, 0, 0};

    if(clip_with_plain) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[0], config, "origin_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[1], config, "origin_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[2], config, "origin_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[3], config, "normal_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[4], config, "normal_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, plain_coords[5], config, "normal_z");
    }

    if(clip_with_bounds) {
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[0], config, "min_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[1], config, "min_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[2], config, "min_z");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[3], config, "max_x");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[4], config, "max_y");
        GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(float, bounds[5], config, "max_z");
    }

    sds base_name = create_base_name("Apd_map", 0, "vtk");
    ds output_dir_with_file = sdsnew(output_dir);
    output_dir_with_file = sdscat(output_dir_with_file, "/");
    output_dir_with_file = sdscatprintf(output_dir_with_file, base_name, current_t);

	bool read_only_data = ((struct save_with_activation_times_persistent_data *) config->persistent_data)->grid != NULL;
    new_vtk_unstructured_grid_from_alg_grid(&(((struct save_with_activation_times_persistent_data *) config->persistent_data)->grid), the_grid, clip_with_plain, plain_coords, clip_with_bounds, bounds, read_only_data, save_f);

    set_vtk_values_with_mean_apd(&config->persistent_data,the_grid);
    save_vtk_unstructured_grid_as_legacy_vtk(((struct save_with_activation_times_persistent_data *) config->persistent_data)->grid, output_dir_with_file, binary, save_f);

	if(save_visible_mask) {
		save_visibility_mask(output_dir_with_file, (((struct save_with_activation_times_persistent_data *) config->persistent_data)->grid)->cell_visibility);
	}

	if(the_grid->adaptive) {
	    free_vtk_unstructured_grid(((struct save_with_activation_times_persistent_data *) config->persistent_data)->grid);
	    ((struct save_with_activation_times_persistent_data *) config->persistent_data)->grid = NULL;
	}

    sdsfree(output_dir_with_file);
	sdsfree(base_name);
}

void set_tissue_vtk_values_with_mean_apd (void **persistent_data, struct grid *the_grid) {

	struct save_with_activation_times_persistent_data **data = (struct save_with_activation_times_persistent_data **)&config->persistent_data;

	struct cell_node **grid_cell = the_grid->active_cells;
    uint32_t num_active_cells =  the_grid->num_active_cells;
    real_cpu center_x, center_y, center_z;

    // We need to pass through the cells in the same order as we did when the VTU grid was created
    for(int i = 0; i < num_active_cells; i++) {
    	center_x = grid_cell[i]->center.x;
        center_y = grid_cell[i]->center.y;
        center_z = grid_cell[i]->center.z;

        struct point_3d cell_coordinates;
        cell_coordinates.x = center_x;
        cell_coordinates.y = center_y;
        cell_coordinates.z = center_z;

        float *apds_array = NULL;
        if(grid_cell[i]->active) {

        	apds_array = (float *) hmget((*data)->tissue_apds, cell_coordinates);
        	unsigned long apd_len = arrlen(apds_array);

        	// Calculate the mean APD values
        	float mean_value = 0.0;
        	mean_value = calculate_mean(apds_array,apd_len);

            // Update the scalar value from the "vtk_unstructured_grid"
            (*data)->grid->values[i] = mean_value;
        }
    }
}
