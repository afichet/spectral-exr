#pragma once

#include "SpectralImage.h"


class BiSpectralImage: public SpectralImage {
    public:
        BiSpectralImage(
            size_t width = 0, size_t height = 0,
            const std::vector<float>& wavelengths_nm = std::vector<float>(),
            bool containsPolarisationData = false
        );

        //virtual void save(const std::string& filename) const = 0;

        //virtual void exportChannels(const std::string& path) const;

        size_t reradiationSize() const {
            return nSpectralBands() * (nSpectralBands() - 1) / 2;
        }

        int idxFromWavelengthIdx(
            size_t wlFrom_idx,
            size_t wlTo_idx
        ) const {
            if (wlFrom_idx < wlTo_idx) {
                return wlTo_idx * (wlTo_idx - 1) / 2 + wlFrom_idx;
            } else {
                return -1;
            }
        }

    protected:
        // Upper right triangular matrices for each pixel
        // pixel stride is reradiationSize()
        std::vector<float> _reradiation; 
};