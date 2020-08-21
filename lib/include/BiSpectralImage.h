#pragma once

#include "SpectralImage.h"

class BiSpectralImage: public SpectralImage {
    struct Reradiation {
        bool reradiate;
        float outWavelenght;
        std::vector<float> image;
    };

    public:
        BiSpectralImage(
            size_t width = 0, size_t height = 0,
            const std::vector<float>& wavelengths_nm = std::vector<float>(),
            SpectrumType type = EMISSIVE_IMAGE,
            bool containsPolarisationData = false
        );

        //virtual void save(const std::string& filename) const = 0;

        //virtual void exportChannels(const std::string& path) const;

    protected:
        std::vector<std::pair<float, Reradiation>> _reradiation;
};