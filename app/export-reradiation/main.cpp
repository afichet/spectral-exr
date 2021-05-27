/**
 * Copyright (c) 2020 - 2021
 * Alban Fichet, Romain Pacanowski, Alexander Wilkie
 * Institut d'Optique Graduate School, CNRS - Universite de Bordeaux,
 * Inria, Charles University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer in the documentation and/or other materials provided
 * with the distribution.
 *
 *  * Neither the name of Institut d'Optique Graduate School, CNRS -
 * Universite de Bordeaux, Inria, Charles University nor the names of
 * its contributors may be used to endorse or promote products derived
 * from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */


#include <iostream>
#include <fstream>

#include <EXRBiSpectralImage.h>

using namespace SEXR;

int main(int argc, char *argv[])
{
    if (argc < 5) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0]
                  << " <bispectral_exr> <x> <y> [<wavelength_i>] <output_file>"
                  << std::endl
                  << std::endl;

        return 0;
    }

    const bool matrixMode = argc == 5;

    const EXRBiSpectralImage image(argv[1]);
    const size_t             x      = std::stoi(argv[2]);
    const size_t             y      = std::stoi(argv[3]);
    const size_t             wl_idx = (matrixMode) ? 0 : std::stoi(argv[4]);
    std::ofstream            tabularOut(argv[matrixMode ? 4 : 5]);

    if (x >= image.width() || y >= image.height()) {
        std::cerr << "Coordinates out of bounds." << std::endl;
        return 0;
    }

    if (matrixMode) {
        tabularOut << "# ";
        for (size_t wl_i_idx = 0; wl_i_idx < image.nSpectralBands();
             wl_i_idx++) {
            tabularOut << image.wavelength_nm(wl_i_idx) << " ";
        }

        tabularOut << "\n";

        for (size_t wl_i_idx = 0; wl_i_idx < image.nSpectralBands();
             wl_i_idx++) {
            for (size_t wl_o_idx = 0; wl_o_idx < image.nSpectralBands();
                 wl_o_idx++) {
                tabularOut << image.getReflectiveValue(x, y, wl_i_idx, wl_o_idx)
                           << " ";
            }

            tabularOut << "\n";
        }
    } else {
        if (wl_idx >= image.nSpectralBands()) {
            std::cerr << "Wavelength index out of bounds." << std::endl;
            return 0;
        }

        tabularOut << "# Remission for wl_i=" << image.wavelength_nm(wl_idx)
                   << "nm\n";

        for (size_t wl_o_idx = 0; wl_o_idx < image.nSpectralBands();
             wl_o_idx++) {
            tabularOut << image.wavelength_nm(wl_o_idx) << " "
                       << image.getReflectiveValue(x, y, wl_idx, wl_o_idx)
                       << "\n";
        }
    }

    return 0;
}