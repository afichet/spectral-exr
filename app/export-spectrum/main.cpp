/**
 * Copyright (c) 2020 Alban Fichet, Romain Pacanowski, Alexander Wilkie
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
#include <fstream>

#include <EXRSpectralImage.h>

using namespace SEXR;

int main(int argc, char *argv[])
{
    if (argc < 5) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <spectral_exr> <x> <y> <output_file>"
                  << std::endl
                  << std::endl;

        return 0;
    }

    const EXRSpectralImage image(argv[1]);
    const size_t           x = std::stoi(argv[2]);
    const size_t           y = std::stoi(argv[3]);
    std::ofstream          tabularOut(argv[4]);

    tabularOut << "# lambda(nm) ";

    for (size_t s = 0; s < image.nStokesComponents(); s++) {
        tabularOut << "S" << s << " ";
    }

    if (image.isReflective()) {
        tabularOut << "T";
    }

    for (size_t wl_idx = 0; wl_idx < image.nSpectralBands(); wl_idx++) {
        tabularOut << image.wavelength_nm(wl_idx);

        for (size_t s = 0; s < image.nStokesComponents(); s++) {
            tabularOut << " " << image.emissive(x, y, wl_idx, s);
        }

        if (image.isReflective()) {
            tabularOut << " " << image.reflective(x, y, wl_idx);
        }

        tabularOut << "\n";
    }

    return 0;
}