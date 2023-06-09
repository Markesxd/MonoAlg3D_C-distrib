//
// Created by sachetto on 13/10/17.
//

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "../alg/grid/grid.h"
#include "../config/restore_state_config.h"

#include "../3dparty/stb_ds.h"


#ifdef COMPILE_CUDA
#include "../gpu_utils/gpu_utils.h"
#endif

#include "../3dparty/sds/sds.h"

struct mesh_data {
    struct point_3d cube_side_length;
    struct point_3d mesh_side_length;
    uint32_t number_of_cells;
    uint32_t num_active_cells;
} __attribute__((packed));

struct cell_data {
    struct point_3d center;
    real_cpu v;
    real_cpu z_front_flux;
    real_cpu z_back_flux;
    real_cpu y_top_flux;
    real_cpu y_down_flux;
    real_cpu x_right_flux;
    real_cpu x_left_flux;
    real_cpu b;
    bool can_change;
    bool active;
    size_t mesh_extra_info_size;
    void *mesh_extra_info;
} __attribute__((packed));

RESTORE_STATE (restore_simulation_state) {

    // Here we restore the saved domain
    if (the_grid) {

        struct mesh_data mesh_data;

        sds tmp = sdsnew (input_dir);
        tmp = sdscat (tmp, "/grid_checkpoint.dat");

        FILE *input_file = fopen (tmp, "rb");

        if (!input_file) {
            fprintf (stderr, "Error opening %s file for restoring the simulation state\n", tmp);
            sdsfree (tmp);
            return false;
        }

        sdsfree (tmp);

        struct point_voidp_hash_entry *mesh_hash = NULL;
        hmdefault(mesh_hash, NULL);

        fread (&mesh_data, sizeof (struct mesh_data), 1, input_file);

        the_grid->cube_side_length = mesh_data.cube_side_length;
        the_grid->mesh_side_length = mesh_data.mesh_side_length;

        // Read the mesh to a point hash
        for (uint32_t i = 0; i < mesh_data.number_of_cells; i++) {
            //struct cell_data *cell_data = (struct cell_data*) malloc(sizeof(struct cell_data));
            struct cell_data *cell_data = MALLOC_ONE_TYPE(struct cell_data);

            // Read center_x, center_y, center_z
            fread(cell_data, sizeof(struct cell_data) - sizeof(void*), 1, input_file); //we need to allocate the extra_mesh_info before reading it

            if(cell_data->mesh_extra_info_size) {
                cell_data->mesh_extra_info = MALLOC_BYTES(void, cell_data->mesh_extra_info_size);
                fread(cell_data->mesh_extra_info, cell_data->mesh_extra_info_size, 1, input_file);
            }

            hmput(mesh_hash, cell_data->center, cell_data);
        }

        printf ("Restoring grid state...\n");

        if(the_grid->adaptive) {
            construct_grid(the_grid);
        }

        struct cell_data * cell_data;

        struct cell_node *grid_cell = the_grid->first_cell;

        while (grid_cell) {

            if (grid_cell->visited) {
                // This cell is already restored
                grid_cell = grid_cell->next;
            } else {

                struct point_3d mesh_point;
                mesh_point.x = grid_cell->center.x;
                mesh_point.y = grid_cell->center.y;
                mesh_point.z = grid_cell->center.z;

                if ((cell_data = (struct cell_data*)hmget (mesh_hash, mesh_point)) != NULL) {

                    grid_cell->active = cell_data->active;

                    // This grid_cell is already in the mesh. We only need to restore the data associated to it...
                    // If the cell is not active we don't need to recover its state
                    if(grid_cell->active) {
                        grid_cell->v = cell_data->v;
                        grid_cell->front_flux = cell_data->z_front_flux;
                        grid_cell->back_flux = cell_data->z_back_flux;
                        grid_cell->top_flux = cell_data->y_top_flux;
                        grid_cell->down_flux = cell_data->y_down_flux;
                        grid_cell->right_flux = cell_data->x_right_flux;
                        grid_cell->left_flux = cell_data->x_left_flux;
                        grid_cell->b = cell_data->b;
                        grid_cell->can_change = cell_data->can_change;
                    }

                    grid_cell->visited = true;
                    grid_cell = grid_cell->next;

                } else {
                    // This grid_cell is not in the  mesh. We need to refine it and check again and again....
                    refine_grid_cell (the_grid, grid_cell);
                    grid_cell = the_grid->first_cell;
                }
            }
        }

        assert (mesh_data.number_of_cells == the_grid->number_of_cells);

        fclose (input_file);
        for(int i = 0; i < hmlen(mesh_hash); i++) {
            free(mesh_hash[i].value);
        }
        hmfree(mesh_hash);
    }

    if (the_monodomain_solver) {

        printf ("Restoring monodomain solver state...\n");

        sds tmp = sdsnew (input_dir);
        tmp = sdscat (tmp, "/monodomain_solver_checkpoint.dat");

        FILE *input_file = fopen (tmp, "rb");

        if (!input_file) {
            fprintf (stderr, "Error opening %s file for restoring the simulation state\n", tmp);
            sdsfree (tmp);
            return false;
        }

        sdsfree (tmp);

        real_cpu config_final_time = the_monodomain_solver->final_time;
        bool config_abort = the_monodomain_solver->abort_on_no_activity;

        fread(the_monodomain_solver, sizeof (struct monodomain_solver), 1, input_file);

        the_monodomain_solver->final_time = config_final_time;
        the_monodomain_solver->abort_on_no_activity = config_abort;

        fread(time_info, sizeof (struct time_info), 1, input_file);

        fclose (input_file);
    }

    if (the_ode_solver) {

        printf ("Restoring ode solver state...\n");

        sds tmp = sdsnew (input_dir);
        tmp = sdscat (tmp, "/ode_solver_checkpoint.dat");

        FILE *input_file = fopen (tmp, "rb");

        if (!input_file) {
            fprintf (stderr, "Error opening %s file for restoring the simulation state\n", tmp);
            sdsfree (tmp);
            return false;
        }

        sdsfree (tmp);

        fread (&(the_ode_solver->adaptive), sizeof(the_ode_solver->adaptive), 1, input_file);
        fread (&(the_ode_solver->max_dt),   sizeof(the_ode_solver->max_dt),   1, input_file);
        fread (&(the_ode_solver->min_dt),   sizeof(the_ode_solver->min_dt),   1, input_file);
        fread (&(the_ode_solver->rel_tol),  sizeof(the_ode_solver->rel_tol),  1, input_file);
        fread (&(the_ode_solver->abs_tol),  sizeof(the_ode_solver->abs_tol),  1, input_file);

        fread (&(the_ode_solver->num_cells_to_solve), sizeof (the_ode_solver->num_cells_to_solve), 1, input_file);

        bool read_cells_to_solve;
        fread (&read_cells_to_solve, sizeof (read_cells_to_solve), 1, input_file);

        size_t num_cells_to_solve = the_ode_solver->num_cells_to_solve;

        the_ode_solver->cells_to_solve = NULL;

        if(read_cells_to_solve) {
            the_ode_solver->cells_to_solve = malloc(sizeof(the_ode_solver->cells_to_solve[0]) * the_ode_solver->num_cells_to_solve);
            fread(the_ode_solver->cells_to_solve, sizeof(the_ode_solver->cells_to_solve[0]), num_cells_to_solve, input_file);
        }

        fread (&(the_ode_solver->gpu), sizeof (the_ode_solver->gpu), 1, input_file);
        fread (&(the_ode_solver->gpu_id), sizeof (the_ode_solver->gpu_id), 1, input_file);

        fread (&(the_ode_solver->pitch), sizeof (the_ode_solver->pitch), 1, input_file);

        fread (&(the_ode_solver->original_num_cells), sizeof (the_ode_solver->original_num_cells), 1, input_file);

        size_t num_sv_entries = the_ode_solver->model_data.number_of_ode_equations;


        if (the_ode_solver->gpu) {

            #ifdef COMPILE_CUDA
            if(the_ode_solver->adaptive) {
                num_sv_entries = num_sv_entries + 3;
            }

            real *sv_cpu;
            sv_cpu = MALLOC_ARRAY_OF_TYPE(real, the_ode_solver->original_num_cells * num_sv_entries);

            fread (sv_cpu, sizeof (real), the_ode_solver->original_num_cells * num_sv_entries, input_file);

            check_cuda_error(cudaMemcpy2D (the_ode_solver->sv, the_ode_solver->pitch, sv_cpu, the_ode_solver->original_num_cells * sizeof (real),
                          the_ode_solver->original_num_cells * sizeof (real), num_sv_entries, cudaMemcpyHostToDevice));
            #endif
        } else {
            fread (the_ode_solver->sv, sizeof (real), the_ode_solver->original_num_cells * num_sv_entries, input_file);
        }

        fread(&(the_ode_solver->extra_data_size), sizeof(the_ode_solver->extra_data_size), 1, input_file);
        fread (the_ode_solver->ode_extra_data, the_ode_solver->extra_data_size, 1, input_file);

        fclose (input_file);
    }

    return true;
}
