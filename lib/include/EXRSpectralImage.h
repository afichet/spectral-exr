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
            SpectrumType type = EMISSIVE,
            bool containsPolarisationData = false
        );

        EXRSpectralImage(
            const std::string& filename
        );

        void save(const std::string& filename) const;

        static SpectralImage::SpectrumType isSpectralChannel(
            const std::string& channelName, 
            int& polarisationComponent, 
            float& wavelengths_nm
        );

        std::string getChannelName(
            int polarisationComponent,
            float wavelength_nm
        ) const;

    static constexpr const char* SPECTRUM_TYPE_ATTR     = "Spectrum type"; 
    static constexpr const char* LENS_TRANSMISSION_ATTR = "Lens transmission"; 
    static constexpr const char* CAMERA_RESPONSE_ATTR   = "Camera response"; 
};
