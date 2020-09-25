#include <iostream>
#include <cmath>

#include <EXRSpectralImage.h>

#include "macbeth_data.h"

using namespace SEXR;

int main(int argc, char* argv[]) {
    (void)argc; (void)argv;

    const size_t width = 600;
    const size_t height = 400;

    std::vector<float> wavelengths(
        std::begin(macbeth_wavelengths), 
        std::end(macbeth_wavelengths));

    EXRSpectralImage spectralImage(
        width, height, 
        wavelengths, 
        REFLECTIVE);

    for (size_t y = 0; y < height; y++) {
        const float v_idx = 4.F * float(y) / float(height);
        const float f_v = v_idx - std::floor(v_idx);

        for (size_t x = 0; x < width; x++) {
            const float u_idx = 6.F * float(x) / float(width);
            const float f_u = u_idx - std::floor(u_idx);

            const float s_v = 0.1;
            const float s_u = 0.1;

            int idx = int(u_idx) + 6*int(v_idx);

            if (
                (
                    (int(v_idx) == 0 && f_v > s_v && f_v < 1.-s_v/2.) ||
                    (int(v_idx) > 0  && int(v_idx) < 3 && f_v > s_v/2. && f_v < 1.-s_v/2.) ||
                    (int(v_idx) == 3 && f_v > s_v/2. && f_v < 1.-s_v)
                ) 
                &&
                (
                    (int(u_idx) == 0 && f_u > s_u && f_u < 1.-s_u/2.) ||
                    (int(u_idx) > 0  && int(u_idx) < 5 && f_u > s_u/2. && f_u < 1.-s_u/2.) ||
                    (int(u_idx) == 5 && f_u > s_u/2. && f_u < 1.-s_u)
                )
            )
            {
                memcpy(&spectralImage(x, y, 0, 0, 0), &macbeth_patches[idx][0], spectralImage.nSpectralBands() * sizeof(float));
            }
        }
    }

    spectralImage.save("Macbeth.exr");

    return 0;
}