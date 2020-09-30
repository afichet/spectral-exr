#pragma once

#include "SpectralImage.h"

namespace SEXR {

class BiSpectralImage: public SpectralImage {
    public:
        BiSpectralImage(
            size_t width = 0, size_t height = 0,
            const std::vector<float>& wavelengths_nm = std::vector<float>(),
            SpectrumType type = REFLECTIVE
        );

        virtual void exportChannels(const std::string& path) const;
        virtual void getRGBImage(std::vector<float>& rgbImage) const;

        size_t reradiationSize() const {
            return nSpectralBands() * (nSpectralBands() - 1) / 2;
        }

        static size_t idxFromWavelengthIdx(
            size_t wlFrom_idx,
            size_t wlTo_idx
        );

        static void wavelengthsIdxFromIdx(
            size_t rerad_idx,
            size_t& wlFrom_idx,
            size_t& wlTo_idx
        );

        virtual float getReflectiveValue(
            size_t x, size_t y, 
            size_t wavelengthFrom_idx, size_t wavelengthTo_idx
        ) const;

        virtual float& operator()(
            size_t x, size_t y, 
            size_t wavelengthFrom_idx, size_t wavelengthTo_idx
        );

        virtual const float& operator()(
            size_t x, size_t y, 
            size_t wavelengthFrom_idx, size_t wavelengthTo_idx
        ) const;


    protected:
        // Upper right triangular matrices for each pixel
        // pixel stride is reradiationSize()
        std::vector<float> _reradiation; 
};

} // namespace SEXR
