//
// Created by sachetto on 13/10/17.
//

#include "assembly_matrix_config.h"
#include "../3dparty/stb_ds.h"

void print_assembly_matrix_config_values(struct config* s) {

    if(s == NULL) {
        log_info("[assembly_matrix] No Assembly Matrix configuration.\n");
        return;
    }

    log_info("[assembly_matrix] Assembly Matrix function configuration:\n");

    log_info("[assembly_matrix] Assembly Matrix library: %s\n", s->library_file_path);

    if(shlen(s->config_data) == 1) {
        log_info("[assembly_matrix] Assembly Matrix parameter:\n");
    } else if(shlen(s->config_data) > 1) {
        log_info("[assembly_matrix] Assembly Matrix parameters:\n");
    }

    STRING_HASH_PRINT_KEY_VALUE_LOG(s->config_data);

}
