#include <iostream>
#include <fstream>

#include <EXRSpectralImage.h>

int main(int argc, char* argv[]) 
{
    if (argc < 5) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <spectral_exr> <x> <y> <output_file>" << std::endl
                  << std::endl;

        return 0;
    }

    const EXRSpectralImage image(argv[1]);
    const size_t x = std::stoi(argv[2]);
    const size_t y = std::stoi(argv[3]);
    std::ofstream tabularOut(argv[4]);

    tabularOut << "# lambda(nm) ";

    for (size_t s = 0; s < image.nStokesComponents(); s++) {
        tabularOut << "S" << s << " ";
    }

    for (size_t m = 0; m < image.nMuellerComponents(); m++) {
        size_t row, col;
        SpectralImage::componentsFromIndex(m, row, col);
        tabularOut << "M" << row << col << " ";
    }

    for (size_t wl_idx = 0; wl_idx < image.nSpectralBands(); wl_idx++) {
        tabularOut << image.wavelength_nm(wl_idx);

        for (size_t p = 0; p < image.nStokesComponents(); p++) {
            tabularOut << " " << image(x, y, wl_idx, p);
        }

        for (size_t p = 0; p < image.nMuellerComponents(); p++) {
            tabularOut << " " << image(x, y, wl_idx, p);
        }

        tabularOut << "\n";
    }

    return 0;
}