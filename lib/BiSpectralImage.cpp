#include <BiSpectralImage.h>
#include "SpectrumConverter.h"

#include <cassert>


BiSpectralImage::BiSpectralImage(
    size_t width, size_t height,
    const std::vector<float>& wavelengths_nm,
    bool containsPolarisationData
)
    : SpectralImage(width, height, wavelengths_nm, REFLECTIVE, containsPolarisationData)
{
    for (size_t s = 0; s < nPolarsiationComponents(); s++) {
        _pixelBuffers[s].resize(nSpectralBands() * _width * _height);
    }

    _reradiation.resize(reradiationSize() * _width * _height);
}


void BiSpectralImage::getRGBImage(std::vector<float>& rgbImage) 
const {
    rgbImage.resize(3 * width() * height());
    SpectrumConverter sc(type());
    
    std::array<float, 3> rgb;

    for (size_t i = 0; i < width() * height(); i++) {
        sc.spectrumToRGB(
            _wavelengths_nm, 
            &_pixelBuffers[0][nSpectralBands() * i],
            &_reradiation[reradiationSize() * i],
            rgb
            );

        memcpy(&rgbImage[3 * i], &rgb[0], 3 * sizeof(float));
    }
}

float& BiSpectralImage::operator()(
    size_t x, size_t y, 
    size_t wavelengthFrom_idx, size_t wavelengthTo_idx,
    size_t polarsiationComponent)
{
    if (wavelengthFrom_idx == wavelengthTo_idx) {
        return SpectralImage::operator()(x, y, wavelengthFrom_idx, polarsiationComponent);
    }

    assert(polarsiationComponent == 0);
    
    size_t reradIdx = idxFromWavelengthIdx(wavelengthFrom_idx, wavelengthTo_idx);

    return _reradiation[reradiationSize() * (y * width() + x) + reradIdx];
}


const float& BiSpectralImage::operator()(
    size_t x, size_t y, 
    size_t wavelengthFrom_idx, size_t wavelengthTo_idx,
    size_t polarsiationComponent) 
const {
    if (wavelengthFrom_idx == wavelengthTo_idx) {
        return SpectralImage::operator()(x, y, wavelengthFrom_idx, polarsiationComponent);
    }

    assert(polarsiationComponent == 0);
    
    size_t reradIdx = idxFromWavelengthIdx(wavelengthFrom_idx, wavelengthTo_idx);

    return _reradiation[reradiationSize() * (y * width() + x) + reradIdx];
}