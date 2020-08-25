#include <iostream>
#include <fstream>
#include <sstream>

#include <EXRSpectralImage.h>
#include <SpectrumAttribute.h>


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

    EXRSpectralImage spectralImage(spectralImageFilename);

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
    case SpectralImage::EMISSIVE:
        additionalInfo << "emissive" << std::endl;
        break;

    case SpectralImage::REFLECTIVE:
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