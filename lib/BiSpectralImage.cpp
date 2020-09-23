#include <BiSpectralImage.h>

#include <sstream>
#include <cassert>

#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>

#include "SpectrumConverter.h"


BiSpectralImage::BiSpectralImage(
    size_t width, size_t height,
    const std::vector<float>& wavelengths_nm
)
    // TODO: spectrum type shall be set explicitly
    : SpectralImage(width, height, wavelengths_nm, REFLECTIVE)
{
    _reradiation.resize(reradiationSize() * _width * _height);
}


void BiSpectralImage::exportChannels(const std::string& path) 
const {
    // Export the diagonal
    SpectralImage::exportChannels(path);

    const size_t xStride = sizeof(float) * reradiationSize();
    const size_t yStride = xStride * width();

    // Export the reradiation
    for (size_t rr = 0; rr < reradiationSize(); rr++) {
        for (size_t wl_i_idx = 0; wl_i_idx < nSpectralBands(); wl_i_idx++) {
            const float& wavelength_i = _wavelengths_nm[wl_i_idx];

            for (size_t wl_o_idx = wl_i_idx + 1; wl_o_idx < nSpectralBands(); wl_o_idx++) {
                const float& wavelength_o = _wavelengths_nm[wl_o_idx];

                std::stringstream filepath;
                filepath << path << "/" << "M00 - " << wavelength_i << "nm - " << wavelength_o << "nm.exr";

                Imf::Header exrHeader(width(), height());
                Imf::ChannelList & exrChannels = exrHeader.channels();
                Imf::FrameBuffer exrFrameBuffer;

                exrChannels.insert("Y", Imf::Channel(Imf::FLOAT));
                exrFrameBuffer.insert("Y", Imf::Slice(Imf::FLOAT, (char*)(&_reradiation[rr]), xStride, yStride));
                
                Imf::OutputFile exrOut(filepath.str().c_str(), exrHeader);
                exrOut.setFrameBuffer(exrFrameBuffer);
                exrOut.writePixels(height());
            }
        }
    }
}


void BiSpectralImage::getRGBImage(std::vector<float>& rgbImage) 
const {
    rgbImage.resize(3 * width() * height());
    SpectrumConverter sc(emissive()); // TODO
    
    std::array<float, 3> rgb;

    for (size_t i = 0; i < width() * height(); i++) {
        sc.spectrumToRGB(
            _wavelengths_nm, 
            &_reflectivePixelBuffers[0][nSpectralBands() * i],
            &_reradiation[reradiationSize() * i],
            rgb
            );

        memcpy(&rgbImage[3 * i], &rgb[0], 3 * sizeof(float));
    }
}

float BiSpectralImage::getPixelValue(
    size_t x, size_t y, 
    size_t wavelengthFrom_idx, size_t wavelengthTo_idx,
    size_t muellerRow,
    size_t muellerColumn) 
const {
    if (wavelengthFrom_idx > wavelengthTo_idx) {
        return 0.F;
    }

    return (*this)(x, y, wavelengthFrom_idx, wavelengthTo_idx, muellerRow, muellerColumn);
}

float& BiSpectralImage::operator()(
    size_t x, size_t y, 
    size_t wavelengthFrom_idx, size_t wavelengthTo_idx,
    size_t muellerRow,
    size_t muellerColumn)
{
    if (wavelengthFrom_idx == wavelengthTo_idx) {
        return SpectralImage::operator()(x, y, wavelengthFrom_idx, muellerRow, muellerColumn);
    }
        
    assert(muellerRow == 0);
    assert(muellerColumn == 0);
    
    size_t reradIdx = idxFromWavelengthIdx(wavelengthFrom_idx, wavelengthTo_idx);

    return _reradiation[reradiationSize() * (y * width() + x) + reradIdx];
}


const float& BiSpectralImage::operator()(
    size_t x, size_t y, 
    size_t wavelengthFrom_idx, size_t wavelengthTo_idx,
    size_t muellerRow,
    size_t muellerColumn) 
const {
    if (wavelengthFrom_idx == wavelengthTo_idx) {
        return SpectralImage::operator()(x, y, wavelengthFrom_idx, muellerRow, muellerColumn);
    }

    assert(muellerRow == 0);
    assert(muellerColumn == 0);

    size_t reradIdx = idxFromWavelengthIdx(wavelengthFrom_idx, wavelengthTo_idx);

    return _reradiation[reradiationSize() * (y * width() + x) + reradIdx];
}