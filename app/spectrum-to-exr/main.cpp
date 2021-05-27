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
#include <regex>
#include <fstream>
#include <vector>

#include <EXRSpectralImage.h>

using namespace SEXR;

int main(int argc, char *argv[])
{
    if (argc < 4) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <spectrum> <type> <output_exr>" << std::endl
                  << std::endl
                  << "<spectrum>   A spectrum in comma separated values with "
                     "wavelength_nm, value."
                  << std::endl
                  << "<type>       Can be \"reflective\" or \"emissive\"."
                  << std::endl
                  << "<output_exr> The path to the spectral EXR to create."
                  << std::endl
                  << std::endl;

        return 0;
    }

    std::vector<float> wavelengths_nm;
    std::vector<float> values;

    // Load the CSV file
    const std::string numberRegex = " *(\\d*\\.?\\d*([Ee][+-]?\\d+)?) *";
    const std::regex  e(numberRegex + "," + numberRegex);

    const std::string fileIn  = argv[1];
    const std::string fileOut = argv[3];
    std::cout << "Reading: [" << fileIn << "]" << std::endl;

    std::ifstream inFile(fileIn);
    std::string   line;

    while (std::getline(inFile, line)) {
        std::smatch matches;
        if (std::regex_search(line, matches, e)) {
            wavelengths_nm.push_back(std::stof(matches[1]));
            values.push_back(std::stof(matches[3]));
        }
    }

    std::cout << "Found " << wavelengths_nm.size() << " samples" << std::endl;
    if (wavelengths_nm.size() == 0) {
        std::cerr << "The provided spectrum is empty!" << std::endl;

        return -1;
    }

    SpectrumType type;

    if (strcmp(argv[2], "reflective") == 0) {
        type = REFLECTIVE;
    } else if (strcmp(argv[2], "emissive") == 0) {
        type = EMISSIVE;
    } else {
        std::cerr << "Invalid argument for spectrum type!" << std::endl;
        std::cerr
          << "The spectrum type can either be \"emissive\" or \"reflective\"."
          << std::endl;

        return -1;
    }

    EXRSpectralImage image(1, 1, wavelengths_nm, type);
    if (type == REFLECTIVE) {
        memcpy(
          &image.reflective(0, 0, 0),
          &values[0],
          values.size() * sizeof(float));
    } else {
        memcpy(
          &image.emissive(0, 0, 0, 0),
          &values[0],
          values.size() * sizeof(float));
    }
    image.save(fileOut);

    std::cout << "File saved as: [" << fileOut << "]" << std::endl;

    return 0;
}