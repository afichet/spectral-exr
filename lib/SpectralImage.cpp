#include <SpectralImage.h>

#include <sstream>
#include <cassert>

#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>


SpectralImage::SpectralImage(
    size_t width, size_t height,
    const std::vector<float>& wavelengths_nm,
    bool containsPolarisationData
)
    : _width(width)
    , _height(height)
    , _wavelengths_nm(wavelengths_nm)
    , _isSpectral(!wavelengths_nm.empty())
    , _containsPolarisationData(containsPolarisationData)
{
    for (size_t s = 0; s < nStokesComponents(); s++) {
        _pixelBuffers[s].resize(nSpectralBands() * _width * _height);
    }
}


void SpectralImage::exportChannels(const std::string& path) const {
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
