#pragma once

#include <vector>
#include <array>
#include <string>

#include "SpectralImage.h"

namespace SEXR {

class EXRSpectralImage : public SpectralImage {
    public:
        EXRSpectralImage(
            size_t width = 0, size_t height = 0,
            const std::vector<float>& wavelengths_nm = std::vector<float>(),
            SpectrumType type = BISPECTRAL
        );

        EXRSpectralImage(
            const std::string& filename
        );

        void save(const std::string& filename) const;

        static SpectrumType channelType(
            const std::string& channelName, 
            int& polarisationComponent, 
            float& wavelengths_nm
        );

        static std::string getStokesChannelName(
            int stokesComponent,
            float wavelength_nm
        );

        static std::string getMuellerChannelName(
            int muellerComponent,
            float wavelength_nm
        );

    static constexpr const char* SPECTRUM_TYPE_ATTR         = "Spectrum type"; 
    static constexpr const char* LENS_TRANSMISSION_ATTR     = "Lens transmission"; 
    static constexpr const char* CAMERA_RESPONSE_ATTR       = "Camera response"; 
    static constexpr const char* EXPOSURE_COMPENSATION_ATTR = "EV";
};

} // namespace SEXR
