#include <BiSpectralImage.h>

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