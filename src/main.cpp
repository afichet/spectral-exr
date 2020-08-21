#include <iostream>
// #include <filesystem>

#include <EXRSpectralImage.h>

int main(int argc, char* argv[]) {

    if (argc > 2) {
        EXRSpectralImage spectralImage(argv[1]);

        spectralImage.exportChannels(argv[2]);
    }

    return 0;
}