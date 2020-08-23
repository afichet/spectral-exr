#include <iostream>
#include <fstream>
#include <sstream>

#include <EXRSpectralImage.h>
#include <SpectrumAttribute.h>


void writeAttributeCSVIfExists(const SpectrumAttribute& attr, const std::string& filename) {
    if (attr.size() > 0) {
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

    EXRSpectralImage spectralImage(argv[1]);
    spectralImage.exportChannels(outputFolder);
    
    writeAttributeCSVIfExists(spectralImage.getCameraResponse(), outputFolder + "/camera.csv");
    writeAttributeCSVIfExists(spectralImage.getLensTransmission(), outputFolder + "/lens.csv");

    for (size_t i = 0; i < spectralImage.nSpectralBands(); i++) {
        const float& wavelength_nm = spectralImage.wavelength_nm(i);
        std::stringstream filename;
        filename << outputFolder << "/" << wavelength_nm << ".csv";

        writeAttributeCSVIfExists(spectralImage.getChannelSensitivity(i), filename.str());
    }

    return 0;
}