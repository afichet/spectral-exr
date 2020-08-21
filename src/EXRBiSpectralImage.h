#pragma once

#include "BiSpectralImage.h"

class EXRBiSpectralImage: public BiSpectralImage {
    enum CHANNEL_TYPE {
        SPECTRAL_DIAGONAL,
        SPECTRAL_RERADIATION,
        RGBA,
        OTHER
    };

    public:
        EXRBiSpectralImage(
            const std::string& filename
        );

        //void save(const std::string& filename) const;

        static CHANNEL_TYPE channelType(
            const std::string& channelName, 
            int& stokesComponent,
            float& wavelengths_nm,
            float& reradiation_wavelength_nm
        );

        static float toWavelength_nm(
            const std::string& value,
            const std::string& multiplier,
            const std::string& unit
        );
        
        static std::string getChannelName(
            int stokesComponent,
            float wavelength_nm
        );

        static std::string getChannelName(
            float wavelength_nm,
            float reradiation_wavelength_nm
        );


};