//
// Created by sachetto on 30/10/18.
//

#ifndef MONOALG3D_VTK_UNSTRUCTURED_GRID_H
#define MONOALG3D_VTK_UNSTRUCTURED_GRID_H

#include "../alg/grid/grid.h"
#include "../common_types/common_types.h"

struct vtk_unstructured_grid {
    uint32_t num_points;
    uint32_t num_cells;

    //TODO: we handle only meshes with the same number of points per cell. If we need something different we will need to change this
    uint32_t points_per_cell;

    //TODO: we handle only meshes with the same cell_type. If we need something different we will need to change this
    uint8_t cell_type;

    f32_array values;
    real_cpu **fibers;
    int64_array cells;
    //TODO: create a mask here, so we can draw only the faces that are not occluded
    ui8_array cell_visibility;
    point3d_array points;

    float min_v;
    float max_v;
};

struct vtk_unstructured_grid *new_vtk_unstructured_grid();

void new_vtk_unstructured_grid_from_alg_grid(struct vtk_unstructured_grid **vtk_grid, struct grid *grid, bool clip_with_plain,
                                                                     float *plain_coordinates, bool clip_with_bounds,
                                                                     float *bounds, bool read_only_values, bool read_fibers_f);

void save_vtk_unstructured_grid_as_vtu(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary);
void save_vtk_unstructured_grid_as_vtu_compressed(struct vtk_unstructured_grid *vtk_grid, const char *filename, int compression_level);
void save_vtk_unstructured_grid_as_legacy_vtk(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary, bool save_f);
void save_vtk_unstructured_grid_as_alg_file(struct vtk_unstructured_grid *vtk_grid, char *filename, bool binary);
void free_vtk_unstructured_grid(struct vtk_unstructured_grid *vtk_grid);

struct vtk_unstructured_grid * new_vtk_unstructured_grid_from_file(const char *vtu_file_name);
struct vtk_unstructured_grid * new_vtk_unstructured_grid_from_file_with_index(const char *file_name, uint32_t v_index);
void new_vtk_unstructured_grid_from_string(struct vtk_unstructured_grid **vtk_grid, char* source, size_t source_size, bool binary, bool read_only_values, int v_index);
void new_vtk_unstructured_grid_from_string_with_activation_info(struct vtk_unstructured_grid **vtk_grid, char* source, size_t source_size);

#endif // MONOALG3D_VTK_UNSTRUCTURED_GRID_H
