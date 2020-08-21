#include <iostream>

#include <EXRSpectralImage.h>

int main(int argc, char* argv[]) {

    if (argc > 2) {
        EXRSpectralImage spectralImage(argv[1]);
        spectralImage.exportChannels(argv[2]);
    } else {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <spectral_exr> <output_folder>" << std::endl
                  << std::endl
                  << "The <output_folder> must have been created prior to the "
                  << "execution and with the correct rights."
                  << std::endl;
    }

    return 0;
}