#include <iostream>
#include <fstream>

#include <EXRBiSpectralImage.h>

int main(int argc, char* argv[]) 
{
    if (argc < 5) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <bispectral_exr> <x> <y> <output_file>" << std::endl
                  << std::endl;

        return 0;
    }

    const EXRBiSpectralImage image(argv[1]);
    const size_t x = std::stoi(argv[2]);
    const size_t y = std::stoi(argv[3]);
    std::ofstream tabularOut(argv[4]);


    tabularOut << "# ";
    for (size_t wl_i_idx = 0; wl_i_idx < image.nSpectralBands(); wl_i_idx++) {
        tabularOut << image.wavelength_nm(wl_i_idx) << " ";
    }

    tabularOut << "\n";

    for (size_t wl_i_idx = 0; wl_i_idx < image.nSpectralBands(); wl_i_idx++) {
        for (size_t wl_o_idx = 0; wl_o_idx < image.nSpectralBands(); wl_o_idx++) {
            tabularOut << image.getPixelValue(x, y, wl_i_idx, wl_o_idx) << " ";
        }
        
        tabularOut << "\n";
    }

    return 0;
}