#pragma once

#include "SpectralImage.h"


class BiSpectralImage: public SpectralImage {
    public:
        BiSpectralImage(
            size_t width = 0, size_t height = 0,
            const std::vector<float>& wavelengths_nm = std::vector<float>(),
            bool containsPolarisationData = false
        );

        virtual void exportChannels(const std::string& path) const;
        virtual void getRGBImage(std::vector<float>& rgbImage) const;

        size_t reradiationSize() const {
            return nSpectralBands() * (nSpectralBands() - 1) / 2;
        }

        static size_t idxFromWavelengthIdx(
            size_t wlFrom_idx,
            size_t wlTo_idx
        ) {
            if (wlFrom_idx < wlTo_idx) {
                return wlTo_idx * (wlTo_idx - 1) / 2 + wlFrom_idx;
            } else {
                return -1;
            }
        }

        void wavelengthsIdxFromIdx(
            size_t rerad_idx,
            size_t& wlFrom_idx,
            size_t& wlTo_idx
        ) const {
            float k = std::floor((std::sqrt(1.F + 8.F * float(rerad_idx)) - 1.F) / 2.F);
            float j = rerad_idx - k * (k + 1) / 2.F;

            wlFrom_idx = k + 1;
            wlTo_idx = j;
        }

        virtual float getPixelValue(
            size_t x, size_t y, 
            size_t wavelengthFrom_idx, size_t wavelengthTo_idx,
            size_t polarsiationComponent = 0) const;

        virtual float& operator()(
            size_t x, size_t y, 
            size_t wavelengthFrom_idx, size_t wavelengthTo_idx,
            size_t polarsiationComponent = 0);

        virtual const float& operator()(
            size_t x, size_t y, 
            size_t wavelengthFrom_idx, size_t wavelengthTo_idx,
            size_t polarsiationComponent = 0) const;


    protected:
        // Upper right triangular matrices for each pixel
        // pixel stride is reradiationSize()
        std::vector<float> _reradiation; 
};