#pragma once

#include "BiSpectralImage.h"

namespace SEXR {

class EXRBiSpectralImage: public BiSpectralImage {
    public:
        EXRBiSpectralImage(
            size_t width = 0, size_t height = 0,
            const std::vector<float>& wavelengths_nm = std::vector<float>(),
            SpectrumType type = REFLECTIVE
        );

        EXRBiSpectralImage(
            const std::string& filename
        );

        void save(const std::string& filename) const;

        static SpectrumType channelType(
            const std::string& channelName, 
            int& polarisationComponent,
            float& wavelengths_nm,
            float& reradiation_wavelength_nm
        );

        static std::string getEmissiveChannelName(
            int stokesComponent,
            float wavelength_nm
        );

        static std::string getReflectiveChannelName(
            float wavelength_nm
        );

        static std::string getReradiationChannelName(
            float wavelength_nm,
            float reradiation_wavelength_nm
        );

    static constexpr const char* SPECTRUM_TYPE_ATTR         = "Spectrum type"; 
    static constexpr const char* LENS_TRANSMISSION_ATTR     = "Lens transmission"; 
    static constexpr const char* CAMERA_RESPONSE_ATTR       = "Camera response"; 
    static constexpr const char* EXPOSURE_COMPENSATION_ATTR = "EV";
};

} // namespace SEXR
