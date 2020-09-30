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
    if (bispectral()) {
        _reradiation.resize(reradiationSize() * _width * _height);
    }
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
                filepath << path << "/" << "T - " << wavelength_i << "nm - " << wavelength_o << "nm.exr";

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
                &_reflectivePixelBuffer[nSpectralBands() * i],
                &_reradiation[reradiationSize() * i],
                &_emissivePixelBuffers[0][nSpectralBands() * i],
                rgb
                );

            memcpy(&rgbImage[3 * i], &rgb[0], 3 * sizeof(float));
        }

        // Exposure compensation
        for (size_t i = 0; i < width() * height(); i++) {
            for (size_t c = 0; c < 3; c++) {
                rgbImage[3 * i + c] *= std::pow(2.F, _ev);
            }
        }
    } else if (reflective()) {
        for (size_t i = 0; i < width() * height(); i++) {
            sc.spectrumToRGB(
                _wavelengths_nm,
                &_reflectivePixelBuffer[nSpectralBands() * i],
                &_reradiation[reradiationSize() * i],
                rgb
                );

            memcpy(&rgbImage[3 * i], &rgb[0], 3 * sizeof(float));
        }

        // Exposure compensation
        for (size_t i = 0; i < width() * height(); i++) {
            for (size_t c = 0; c < 3; c++) {
                rgbImage[3 * i + c] *= std::pow(2.F, _ev);
            }
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
    size_t wavelengthFrom_idx, size_t wavelengthTo_idx) 
const {
    assert(x < width());
    assert(y < height());
    assert(wavelengthFrom_idx < nSpectralBands());
    assert(wavelengthTo_idx < nSpectralBands());

    if (!reflective() || wavelengthFrom_idx > wavelengthTo_idx) {
        return 0.F;
    }

    if (!bispectral()) {
        return 0.F;
    }

    return this->operator()(x, y, wavelengthFrom_idx, wavelengthTo_idx);
}


float& BiSpectralImage::operator()(
    size_t x, size_t y, 
    size_t wavelengthFrom_idx, size_t wavelengthTo_idx)
{
    assert(reflective());
    assert(x < width());
    assert(y < height());
    assert(wavelengthFrom_idx < nSpectralBands());
    assert(wavelengthTo_idx < nSpectralBands());
    assert(wavelengthFrom_idx <= wavelengthTo_idx);

    if (wavelengthFrom_idx == wavelengthTo_idx) {
        return SpectralImage::operator()(x, y, wavelengthFrom_idx);
    }
    
    assert(bispectral());
    assert(_reradiation.size() == reradiationSize() * width() * height());

    size_t reradIdx = idxFromWavelengthIdx(wavelengthFrom_idx, wavelengthTo_idx);

    return _reradiation[reradiationSize() * (y * width() + x) + reradIdx];
}


const float& BiSpectralImage::operator()(
    size_t x, size_t y, 
    size_t wavelengthFrom_idx, size_t wavelengthTo_idx) 
const {
    assert(reflective());
    assert(x < width());
    assert(y < height());
    assert(wavelengthFrom_idx < nSpectralBands());
    assert(wavelengthTo_idx < nSpectralBands());

    if (wavelengthFrom_idx == wavelengthTo_idx) {
        return SpectralImage::operator()(x, y, wavelengthFrom_idx);
    }

    assert(bispectral());
    assert(_reradiation.size() == reradiationSize() * width() * height());

    size_t reradIdx = idxFromWavelengthIdx(wavelengthFrom_idx, wavelengthTo_idx);

    return _reradiation[reradiationSize() * (y * width() + x) + reradIdx];
}

} // namespace SEXR
