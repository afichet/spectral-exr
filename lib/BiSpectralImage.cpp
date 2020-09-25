#include <BiSpectralImage.h>

#include <sstream>
#include <cassert>
#include <cmath>

#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>

#include "SpectrumConverter.h"

namespace SEXR {

BiSpectralImage::BiSpectralImage(
    size_t width, size_t height,
    const std::vector<float>& wavelengths_nm,
    SpectrumType type
)
    : SpectralImage(width, height, wavelengths_nm, type)
{
    _reradiation.resize(reradiationSize() * _width * _height);
}


void BiSpectralImage::exportChannels(const std::string& path) 
const {
    // TODO: support the emissive part
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
    SpectrumConverter sc(emissive());
    
    std::array<float, 3> rgb;

    if (emissive() && reflective()) {
        for (size_t i = 0; i < width() * height(); i++) {
            sc.spectraToRGB(
                _wavelengths_nm,
                &_reflectivePixelBuffers[0][nSpectralBands() * i],
                &_reradiation[reradiationSize() * i],
                &_emissivePixelBuffers[0][nSpectralBands() * i],
                rgb
                );

            memcpy(&rgbImage[3 * i], &rgb[0], 3 * sizeof(float));
        }
    } else if (reflective()) {
        for (size_t i = 0; i < width() * height(); i++) {
            sc.spectrumToRGB(
                _wavelengths_nm,
                &_reflectivePixelBuffers[0][nSpectralBands() * i],
                &_reradiation[reradiationSize() * i],
                rgb
                );

            memcpy(&rgbImage[3 * i], &rgb[0], 3 * sizeof(float));
        }
    } else {
        SpectralImage::getRGBImage(rgbImage);
    }
}


size_t BiSpectralImage::idxFromWavelengthIdx(
    size_t wlFrom_idx,
    size_t wlTo_idx
) {
    if (wlFrom_idx < wlTo_idx) {
        return wlTo_idx * (wlTo_idx - 1) / 2 + wlFrom_idx;
    } else {
        return -1;
    }
}


void BiSpectralImage::wavelengthsIdxFromIdx(
    size_t rerad_idx,
    size_t& wlFrom_idx,
    size_t& wlTo_idx
) {
    float k = std::floor((std::sqrt(1.F + 8.F * float(rerad_idx)) - 1.F) / 2.F);
    float j = rerad_idx - k * (k + 1) / 2.F;

    wlFrom_idx = j;
    wlTo_idx = k + 1;
}


float BiSpectralImage::getReflectiveValue(
    size_t x, size_t y, 
    size_t wavelengthFrom_idx, size_t wavelengthTo_idx,
    size_t muellerRow,
    size_t muellerColumn) 
const {
    if (!reflective() || wavelengthFrom_idx > wavelengthTo_idx) {
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

} // namespace SEXR
