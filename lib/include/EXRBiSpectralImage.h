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
            const std::string& filename
        );

        //void save(const std::string& filename) const;

        ChannelType channelType(
            const std::string& channelName, 
            int& muellerComponent,
            float& wavelengths_nm,
            float& reradiation_wavelength_nm
        ) const;

        static float toWavelength_nm(
            const std::string& value,
            const std::string& multiplier,
            const std::string& unit
        );
        
        std::string getChannelName(
            int muellerComponent,
            float wavelength_nm
        ) const;

        std::string getChannelName(
            float wavelength_nm,
            float reradiation_wavelength_nm
        ) const;


};