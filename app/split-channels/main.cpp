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
#include <fstream>
#include <sstream>

#include <EXRBiSpectralImage.h>
#include <SpectrumAttribute.h>

using namespace SEXR;

void writeAttributeCSVIfExists(const SpectrumAttribute& attr, const std::string& filename) {
    if (attr.size() > 0) {
        std::cout << "Exporting metadata: [" << filename << "]" << std::endl;
        std::ofstream cameraCSV(filename);

        for (size_t i = 0; i < attr.size(); i++) {
            cameraCSV << attr.wavelength_nm(i) << ","
                      << attr.value(i) << std::endl;
        }
    }
}


int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <spectral_exr> <output_folder>" << std::endl
                  << std::endl
                  << "The <output_folder> must have been created prior to the "
                  << "execution and with the correct rights."
                  << std::endl;

        return 0;
    }

    const std::string outputFolder = argv[2];
    const std::string spectralImageFilename = argv[1];

    EXRBiSpectralImage spectralImage(spectralImageFilename);

    // Export individual channels as separate files
    std::cout << "Writing spectral channels in: [" << outputFolder << "]" << std::endl;
    spectralImage.exportChannels(outputFolder);
    
    // If present, write the additional response and transmission
    writeAttributeCSVIfExists(spectralImage.cameraResponse(), outputFolder + "/camera.csv");
    writeAttributeCSVIfExists(spectralImage.lensTransmission(), outputFolder + "/lens.csv");

    for (size_t wl_idx = 0; wl_idx < spectralImage.nSpectralBands(); wl_idx++) {
        const float& wavelength_nm = spectralImage.wavelength_nm(wl_idx);
        std::stringstream filename;
        filename << outputFolder << "/" << wavelength_nm << ".csv";

        writeAttributeCSVIfExists(spectralImage.channelSensitivity(wl_idx), filename.str());
    }

    // Create a simple text file for additional information
    const std::string additionalInfoFilename(outputFolder + "/" + "info.txt");
    std::cout << "Writing image informations: [" << additionalInfoFilename << "]" << std::endl;
    std::ofstream additionalInfo(additionalInfoFilename);

    additionalInfo << "Spectral Image: " << spectralImageFilename  << " " 
                   << spectralImage.width() << "x" 
                   << spectralImage.height() << "px" << std::endl;


    additionalInfo << "\tType: ";

    switch (spectralImage.type())
    {
    case EMISSIVE:
        additionalInfo << "emissive" << std::endl;
        break;

    case REFLECTIVE:
        additionalInfo << "reflective" << std::endl;
        break;

    default:
        additionalInfo << "unknown" << std::endl;
        break;
    }
    

    additionalInfo << "\tPolarised: ";

    if (spectralImage.polarised()) {
        additionalInfo << "YES" << std::endl;
    } else {
        additionalInfo << "NO" << std::endl;
    }


    additionalInfo << "\tSpectral bands: " << spectralImage.nSpectralBands() << std::endl;

    for (size_t wl_idx = 0; wl_idx < spectralImage.nSpectralBands(); wl_idx++) {
        additionalInfo << "\t\t" << spectralImage.wavelength_nm(wl_idx) << "nm" << std::endl;
    }

    additionalInfo << "Metadata:" << std::endl;
    bool haveMetadata = false;

    if (spectralImage.cameraResponse().size() > 0) {
        additionalInfo << "\tHave camera response information" << std::endl;
        haveMetadata = true;
    }

    if (spectralImage.lensTransmission().size() > 0) {
        additionalInfo << "\tHave lens transmission information" << std::endl;
        haveMetadata = true;
    }

    for (size_t wl_idx = 0; wl_idx < spectralImage.nSpectralBands(); wl_idx++) {
        if (spectralImage.channelSensitivity(wl_idx).size() > 0) {
            const float& wavelength_nm = spectralImage.wavelength_nm(wl_idx);
            additionalInfo << "\tFilter response for " << wavelength_nm << "nm" << std::endl;
            haveMetadata = true;
        }
    }

    if (!haveMetadata) {
        std::cout << "\tNone" << std::endl;
    }

    return 0;
}