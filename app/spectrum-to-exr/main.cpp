#include <iostream>
#include <regex>
#include <fstream>
#include <vector>

#include <EXRSpectralImage.h>


int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <spectrum> <type> <output_exr>" << std::endl
                  << std::endl
                  << "<spectrum>   A spectrum in comma separated values with wavelength_nm, value." << std::endl
                  << "<type>       Ca be \"reflective\" or \"emissive\"." << std::endl
                  << "<output_exr> The path of the spectral EXR to create." << std::endl
                  << std::endl;

        return 0;
    }

    std::vector<float> wavelengths_nm;
    std::vector<float> values;

    // Load the CSV file
    const std::string numberRegex = "(\\d*\\.?\\d*([Ee][+-]?\\d+)?)";
    const std::regex e("^" + numberRegex + "," + numberRegex + "$");

    std::ifstream inFile(argv[1]);
    std::string line;

    while(std::getline(inFile, line)) {
        std::smatch matches;
        if (std::regex_search(line, matches, e)) {
            wavelengths_nm.push_back(std::stof(matches[1]));
            values.push_back(std::stof(matches[3]));
        }
    }

    SpectralImage::SpectrumType type;

    if (strcmp(argv[2], "reflective") == 0) {
        type = SpectralImage::REFLECTIVE;
    } else if (strcmp(argv[2], "emissive") == 0) {
        type = SpectralImage::EMISSIVE;
    } else {
        std::cerr << "Invalid argument for spectrum type!" << std::endl;
        std::cerr << "The spectrum type can either be \"emissive\" or \"reflective\"." << std::endl;

        return -1;
    }

    EXRSpectralImage image(1, 1, wavelengths_nm, type);
    memcpy(&image(0, 0, 0), &values[0], values.size() * sizeof(float));
    image.save(argv[3]);

    return 0;
}