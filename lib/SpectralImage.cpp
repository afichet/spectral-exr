#include <SpectralImage.h>

#include <sstream>
#include <cassert>

#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>

#include "SpectrumConverter.h"

SpectralImage::SpectralImage(
    size_t width, size_t height,
    const std::vector<float>& wavelengths_nm,
    SpectrumType type,
    bool containsPolarisationData
)
    : _width(width)
    , _height(height)
    , _wavelengths_nm(wavelengths_nm)
    , _isSpectral(!wavelengths_nm.empty())
    , _containsPolarisationData(containsPolarisationData)
    , _spectrumType(type)

{
    for (size_t s = 0; s < nStokesComponents(); s++) {
        _pixelBuffers[s].resize(nSpectralBands() * _width * _height);
    }
}


void SpectralImage::exportChannels(const std::string& path) 
const {
    const size_t xStride = sizeof(float) * nSpectralBands();
    const size_t yStride = xStride * width();

    for (size_t s = 0; s < nSpectralBands(); s++) {
        for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
            const float& wavelength = _wavelengths_nm[wl_idx];
            std::stringstream filepath;
            filepath << path << "/" << "S" << s << " - " << wavelength << "nm.exr";

            Imf::Header exrHeader(width(), height());
            Imf::ChannelList & exrChannels = exrHeader.channels();
            Imf::FrameBuffer exrFrameBuffer;

            exrChannels.insert("Y", Imf::Channel(Imf::FLOAT));
            exrFrameBuffer.insert("Y", Imf::Slice(Imf::FLOAT, (char*)(&_pixelBuffers[s][wl_idx]), xStride, yStride));
            
            Imf::OutputFile exrOut(filepath.str().c_str(), exrHeader);
            exrOut.setFrameBuffer(exrFrameBuffer);
            exrOut.writePixels(height());
        }
    }
}


void SpectralImage::getRGBImage(std::vector<float>& rgbImage) 
const {
    rgbImage.resize(3 * width() * height());
    SpectrumConverter sc((_spectrumType == EMISSIVE_IMAGE) ? SpectrumConverter::EMISSIVE : SpectrumConverter::REFLECTIVE);
    
    std::array<float, 3> rgb;

    for (size_t i = 0; i < width() * height(); i++) {
        sc.spectrumToRGB(
            _wavelengths_nm, 
            &_pixelBuffers[0][nSpectralBands() * i],
            rgb
            );

        memcpy(&rgbImage[3 * i], &rgb[0], 3 * sizeof(float));
    }
}


float& SpectralImage::operator()(
    size_t x, size_t y,
    size_t wavelength_idx, 
    size_t stokes
) {
    assert(x < width());
    assert(y < height());
    assert(wavelength_idx < nSpectralBands());
    assert(stokes < 4);
    return _pixelBuffers[stokes][nSpectralBands() * (y * width() + x) + wavelength_idx];
}


const float& SpectralImage::operator()(
    size_t x, size_t y,
    size_t wavelength_idx, 
    size_t stokes
) const {
    assert(x < width());
    assert(y < height());
    assert(wavelength_idx < nSpectralBands());
    assert(stokes < 4);
    return _pixelBuffers[stokes][nSpectralBands() * (y * width() + x) + wavelength_idx];
}
