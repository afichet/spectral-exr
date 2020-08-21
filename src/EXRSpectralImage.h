#pragma once

#include <vector>
#include <array>
#include <string>

#include "SpectralImage.h"

class EXRSpectralImage : public SpectralImage {
    public:
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
