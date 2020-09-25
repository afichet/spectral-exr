#include <SpectralImage.h>

#include <sstream>
#include <cassert>
#include <functional>

#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>

#include "SpectrumConverter.h"


SpectralImage::SpectralImage(
    size_t width, size_t height,
    const std::vector<float>& wavelengths_nm,
    SpectrumType type
)
    : _width(width)
    , _height(height)
    , _wavelengths_nm(wavelengths_nm)
    , _spectrumType(type)
{
    const size_t buffSize = nSpectralBands() * _width * _height;

    for (size_t s = 0; s < nStokesComponents(); s++) {
        _emissivePixelBuffers[s].resize(buffSize);
    }

    for (size_t m = 0; m < nMuellerComponents(); m++) {
        _reflectivePixelBuffers[m].resize(buffSize);
    }

    _channelSensitivities.resize(nSpectralBands());
}


void SpectralImage::exportChannels(const std::string& path) 
const {
    // Utility function to create an EXR from a monochromatic buffer
    std::function<void(const std::string&, const float*)> writeEXR = 
    [this](const std::string& filename, const float* buffer) {
        const size_t xStride = sizeof(float) * nSpectralBands();
        const size_t yStride = xStride * width();

        Imf::Header exrHeader(width(), height());
        Imf::ChannelList & exrChannels = exrHeader.channels();
        Imf::FrameBuffer exrFrameBuffer;

        exrChannels.insert("Y", Imf::Channel(Imf::FLOAT));
        exrFrameBuffer.insert("Y", Imf::Slice(Imf::FLOAT, (char*)(buffer), xStride, yStride));
        
        Imf::OutputFile exrOut(filename.c_str(), exrHeader);
        exrOut.setFrameBuffer(exrFrameBuffer);
        exrOut.writePixels(height());
    };

    // Export the emissive part
    for (size_t s = 0; s < nStokesComponents(); s++) {
        std::stringstream filePrefix;

        filePrefix << "S" << s;

        for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
            const float& wavelength = _wavelengths_nm[wl_idx];
            std::stringstream filepath;
            filepath << path << "/" << filePrefix.str() << " - " << wavelength << "nm.exr";

            writeEXR(filepath.str(), &_emissivePixelBuffers[s][wl_idx]);
        }
    }

    // Export the reflective part
    for (size_t m = 0; m < nMuellerComponents(); m++) {
        std::stringstream filePrefix;
        size_t row, col;

        componentsFromIndex(m, row, col);
        filePrefix << "M" << row << col;

        for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
            const float& wavelength = _wavelengths_nm[wl_idx];
            std::stringstream filepath;
            filepath << path << "/" << filePrefix.str() << " - " << wavelength << "nm.exr";

            writeEXR(filepath.str(), &_reflectivePixelBuffers[m][wl_idx]);
        }
    }
}


void SpectralImage::getRGBImage(std::vector<float>& rgbImage) 
const {
    rgbImage.resize(3 * width() * height());
    SpectrumConverter sc(emissive());
    
    std::array<float, 3> rgb;

    if (emissive() && reflective()) {
        for (size_t i = 0; i < width() * height(); i++) {
            sc.spectraToRGB(
                _wavelengths_nm,
                &_reflectivePixelBuffers[0][nSpectralBands() * i],
                &_emissivePixelBuffers[0][nSpectralBands() * i],
                rgb
                );

            memcpy(&rgbImage[3 * i], &rgb[0], 3 * sizeof(float));
        }
    } else if (emissive()) {
        for (size_t i = 0; i < width() * height(); i++) {
            sc.spectrumToRGB(
                _wavelengths_nm,
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
                rgb
                );

            memcpy(&rgbImage[3 * i], &rgb[0], 3 * sizeof(float));
        }
    }
}


void SpectralImage::setCameraResponse(
    const std::vector<float>& wavelengths_nm,
    const std::vector<float>& values
) {
    assert(wavelengths_nm.size() == values.size());

    _cameraReponse = SpectrumAttribute(wavelengths_nm, values);
}


const SpectrumAttribute& SpectralImage::cameraResponse() 
const { 
    return _cameraReponse; 
}


void SpectralImage::setLensTransmission(
    const std::vector<float>& wavelengths_nm,
    const std::vector<float>& values
) {
    assert(wavelengths_nm.size() == values.size());

    _lensTransmissionSpectra = SpectrumAttribute(wavelengths_nm, values);
}


const SpectrumAttribute& SpectralImage::lensTransmission() 
const { 
    return _lensTransmissionSpectra; 
}


void SpectralImage::setChannelSensitivity(
    size_t wl_idx,
    const std::vector<float>& wavelengths_nm,
    const std::vector<float>& values
) {
    assert(wl_idx < _emissivePixelBuffers[0].size());
    assert(wavelengths_nm.size() == values.size());
    
    _channelSensitivities[wl_idx] = SpectrumAttribute(wavelengths_nm, values);
}


const std::vector<SpectrumAttribute>& SpectralImage::channelSensitivities() 
const { 
    return _channelSensitivities; 
}


const SpectrumAttribute& SpectralImage::channelSensitivity(size_t wl_idx) 
const {
    assert(wl_idx < _channelSensitivities.size());

    return _channelSensitivities[wl_idx]; 
}

// Access the emissive part

float& SpectralImage::operator()(
    size_t x, size_t y,
    size_t wavelength_idx, 
    size_t stokesComponent
) {
    assert(x < width());
    assert(y < height());
    assert(wavelength_idx < nSpectralBands());
    assert(emissive());
    assert(stokesComponent < nStokesComponents());

    return _emissivePixelBuffers[stokesComponent][nSpectralBands() * (y * width() + x) + wavelength_idx];
}


const float& SpectralImage::operator()(
    size_t x, size_t y,
    size_t wavelength_idx, 
    size_t stokesComponent
) const {
    assert(x < width());
    assert(y < height());
    assert(wavelength_idx < nSpectralBands());
    assert(emissive());
    assert(stokesComponent < nStokesComponents());

    return _emissivePixelBuffers[stokesComponent][nSpectralBands() * (y * width() + x) + wavelength_idx];
}

// Access the reflective part

float& SpectralImage::operator()(
    size_t x, size_t y,
    size_t wavelength_idx, 
    size_t muellerRow,
    size_t muellerColumn
) {
    assert(x < width());
    assert(y < height());
    assert(wavelength_idx < nSpectralBands());
    assert(reflective());
    assert(muellerRow <= (isPolarised() ? 3 : 1));
    assert(muellerColumn <= (isPolarised() ? 3 : 1));

    return _reflectivePixelBuffers[indexFromComponents(muellerRow, muellerColumn)][nSpectralBands() * (y * width() + x) + wavelength_idx];
}


const float& SpectralImage::operator()(
    size_t x, size_t y,
    size_t wavelength_idx, 
    size_t muellerRow,
    size_t muellerColumn
) const {
    assert(x < width());
    assert(y < height());
    assert(wavelength_idx < nSpectralBands());
    assert(reflective());
    assert(muellerRow <= (isPolarised() ? 3 : 1));
    assert(muellerColumn <= (isPolarised() ? 3 : 1));

    return _reflectivePixelBuffers[indexFromComponents(muellerRow, muellerColumn)][nSpectralBands() * (y * width() + x) + wavelength_idx];
}


float SpectralImage::getEmissiveValue(
    size_t x, size_t y, 
    size_t wavelength_idx, 
    size_t stokesComponent
) const {
    if (emissive()) {
        return (*this)(x, y, wavelength_idx, stokesComponent);
    }

    return 0.F;
}


float SpectralImage::getReflectiveValue(
    size_t x, size_t y, 
    size_t wavelength_idx, 
    size_t muellerRow, 
    size_t muellerColumn
) const {
    if (reflective()) {
        return (*this)(x, y, wavelength_idx, muellerRow, muellerColumn);
    }

    return 0.F;
}


const float& SpectralImage::wavelength_nm(size_t wl_idx) 
const { 
    assert(wl_idx < _wavelengths_nm.size());

    return _wavelengths_nm[wl_idx]; 
}


size_t SpectralImage::width()  const { return _width; }
size_t SpectralImage::height() const { return _height; }
size_t SpectralImage::nSpectralBands() const { return _wavelengths_nm.size(); }


size_t SpectralImage::nStokesComponents()  const {
    if (emissive()) {
        if (polarised()) {
            return 4;
        } else {
            return 1;
        }
    }

    return 0;
}


size_t SpectralImage::nMuellerComponents() const {
    if (reflective()) {
        if (polarised()) {
            return 16;
        } else {
            return 1;
        }
    }

    return 0;
}


bool SpectralImage::polarised()    const { return isPolarised(_spectrumType); }
bool SpectralImage::emissive()     const { return isEmissive(_spectrumType); }
bool SpectralImage::reflective()   const { return isReflective(_spectrumType); }
SpectrumType SpectralImage::type() const { return _spectrumType; }


void SpectralImage::componentsFromIndex(
    size_t index,
    size_t& row,
    size_t& col
) {
    row = index % 4;
    col = index / 4;
}


size_t SpectralImage::indexFromComponents(
    size_t row,
    size_t col
) {
    assert(row < 4);
    assert(col < 4);

    return 4 * col + row;
}
