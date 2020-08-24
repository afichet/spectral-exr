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
    , _containsPolarisationData(containsPolarisationData)
    , _spectrumType(type)
{
    for (size_t s = 0; s < nPolarsiationComponents(); s++) {
        _pixelBuffers[s].resize(nSpectralBands() * _width * _height);
    }

    _channelSensitivities.resize(nSpectralBands());
}


void SpectralImage::exportChannels(const std::string& path) 
const {
    const size_t xStride = sizeof(float) * nSpectralBands();
    const size_t yStride = xStride * width();

    for (size_t s = 0; s < nPolarsiationComponents(); s++) {
        std::stringstream filePrefix;

        if (emissive()) {
            filePrefix << "S" << s;
        } else {
            size_t row, col;
            componentsFromIndex(s, row, col);
            filePrefix << "M" << row << col;
        }

        for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
            const float& wavelength = _wavelengths_nm[wl_idx];
            std::stringstream filepath;
            filepath << path << "/" << filePrefix.str() << " - " << wavelength << "nm.exr";

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
    SpectrumConverter sc((_spectrumType == EMISSIVE) ? SpectrumConverter::EMISSIVE : SpectrumConverter::REFLECTIVE);
    
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
    assert(wl_idx < _pixelBuffers[0].size());
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

const float& SpectralImage::wavelength_nm(size_t wl_idx) 
const { 
    assert(wl_idx < _wavelengths_nm.size());

    return _wavelengths_nm[wl_idx]; 
}

size_t SpectralImage::width()  const { return _width; }
size_t SpectralImage::height() const { return _height; }
size_t SpectralImage::nSpectralBands() const { return _wavelengths_nm.size(); }

size_t SpectralImage::nPolarsiationComponents() 
const { 
    if (polarised()) {
        switch(type()) {
            case EMISSIVE:
            return 4;

            case REFLECTIVE:
            return 16;

            case UNDEFINED:
            return 0;
        }
    }

    return 1;
}

bool SpectralImage::polarised()    const { return _containsPolarisationData; }
bool SpectralImage::emissive()     const { return _spectrumType == EMISSIVE; }
bool SpectralImage::reflective()   const { return _spectrumType == REFLECTIVE; }
SpectralImage::SpectrumType SpectralImage::type() const { return _spectrumType; }


void SpectralImage::componentsFromIndex(
    size_t index,
    size_t& row,
    size_t& col
) const {
    switch(type()) {
        case EMISSIVE:
        assert(index < 4);

        row = index;
        col = 0;
        break;

        case REFLECTIVE:
        assert(index < 16);

        row = index % 4;
        col = index / 4;
        break;

        case UNDEFINED:
        assert(0);

        row = 0;
        col = 0;
        break;
    }
}



size_t SpectralImage::indexFromComponents(
    size_t row,
    size_t col
) const {

#ifndef NDEBUG
    switch (type()) {
        case EMISSIVE:
        assert(row < 4);
        assert(col == 0);
        break;

        case REFLECTIVE:
        assert(row < 4);
        assert(col < 4);
        break;

        default:
        assert(0);
        break;
    }
#endif

    return 4 * col + row;
}
