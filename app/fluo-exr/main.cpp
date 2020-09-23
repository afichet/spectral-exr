#include <cassert>
#include <iostream>
#include <vector>

#include <EXRBiSpectralImage.h>

#include "fluo_data.h"


int main(/*int argc, char *argv[]*/) {
    const size_t width = 150;
    const size_t height = 150;

    // Make it square
    const float wl_start = std::max(wi_start, wo_start);
    assert(wi_inc == wo_inc);
    const float wl_inc = wi_inc;

    const size_t wl_i_idx_start = (wl_start - wi_start) / wi_inc;
    const size_t wl_o_idx_start = (wl_start - wo_start) / wo_inc;

    const size_t wl_size = std::min(wi_size, wo_size);

    std::vector<float> wavelengths_nm(wl_size);

    for (size_t i = 0; i < wl_size; i++) {
      wavelengths_nm[i] = wl_start + i * wl_inc;
    }

    EXRBiSpectralImage fluoImage(width, height, wavelengths_nm);

    size_t x_sz = 50;
    size_t y_sz = 50;

    for (size_t y = 0; y < height; y++) {
        size_t c_id = y / y_sz;

        for (size_t x = 0; x < width; x++) {
            size_t r_id = x / x_sz;

            // Checkerboard patern
            const auto ptr = (((c_id % 2) + (r_id % 2)) % 2 == 0)
                                ? fluorescent_pink
                                : fluorescent_yellow;

            for (size_t wl_i_idx = 0; wl_i_idx < wl_size; wl_i_idx++) {
                const size_t db_i = wl_i_idx + wl_i_idx_start;
                assert(db_i < wi_size);

                for (size_t wl_o_idx = wl_i_idx; wl_o_idx < wl_size; wl_o_idx++) {
                    const size_t db_o = wl_o_idx + wl_o_idx_start;
                    assert(db_o < wo_size);

                    fluoImage(x, y, wl_i_idx, wl_o_idx, 0, 0) = std::max(0.F, ptr[db_o][db_i]);
                }
            }
        }
    }

    fluoImage.save("BiSpectral.exr");
}