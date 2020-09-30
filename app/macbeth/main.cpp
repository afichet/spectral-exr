/**
 * Copyright (c) 2020 
 *  Alban Fichet <alban.fichet@gmx.fr>, 
 *  Romain Pacanowski <romain.pacanowski@inria.fr>, 
 *  Alexander Wilkie <alexander.wilkie@cgg.mff.cuni.cz>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *  * Neither the name of %ORGANIZATION% nor the names of its contributors may be
 * used to endorse or promote products derived from this software without specific
 * prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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
                memcpy(&spectralImage(x, y, 0), &macbeth_patches[idx][0], spectralImage.nSpectralBands() * sizeof(float));
            }
        }
    }

    spectralImage.save("Macbeth.exr");

    return 0;
}