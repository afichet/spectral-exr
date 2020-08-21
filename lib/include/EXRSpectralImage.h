#pragma once

#include <vector>
#include <array>
#include <string>

#include "SpectralImage.h"

class EXRSpectralImage : public SpectralImage {
    public:
        EXRSpectralImage(
            size_t width = 0, size_t height = 0,
            const std::vector<float>& wavelengths_nm = std::vector<float>(),
            SpectrumType type = EMISSIVE_IMAGE,
            bool containsPolarisationData = false
        );

        EXRSpectralImage(
            const std::string& filename
        );

        void save(const std::string& filename) const;

        static bool isSpectralChannel(
            const std::string& channelName, 
            int& stokesComponent, 
            float& wavelengths_nm
        );

        static std::string getChannelName(
            int stokesComponent,
            float wavelength_nm
        );
};
