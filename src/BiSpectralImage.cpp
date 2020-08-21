#include <BiSpectralImage.h>

BiSpectralImage::BiSpectralImage(
    size_t width, size_t height,
    const std::vector<float>& wavelengths_nm,
    bool containsPolarisationData
)
    : SpectralImage(width, height, wavelengths_nm, containsPolarisationData)
{
    for (size_t s = 0; s < nStokesComponents(); s++) {
        _pixelBuffers[s].resize(nSpectralBands() * _width * _height);
    }
}