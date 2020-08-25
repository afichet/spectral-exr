#pragma once

#include "BiSpectralImage.h"

class EXRBiSpectralImage: public BiSpectralImage {
    enum ChannelType {
        DIAGONAL,
        RERADIATION,
        RGBA,
        OTHER
    };

    public:
        EXRBiSpectralImage(
            size_t width = 0, size_t height = 0,
            const std::vector<float>& wavelengths_nm = std::vector<float>(),
            bool containsPolarisationData = false
        );

        EXRBiSpectralImage(
            const std::string& filename
        );

        void save(const std::string& filename) const;

        ChannelType channelType(
            const std::string& channelName, 
            int& muellerComponent,
            float& wavelengths_nm,
            float& reradiation_wavelength_nm
        ) const;

        std::string getDiagonalChannelName(
            int muellerComponent,
            float wavelength_nm
        ) const;

        std::string getReradiationChannelName(
            float wavelength_nm,
            float reradiation_wavelength_nm
        ) const;


};