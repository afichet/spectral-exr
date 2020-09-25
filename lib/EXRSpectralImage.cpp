#include <EXRSpectralImage.h>
#include "Util.h"

#include <regex>
#include <algorithm>
#include <sstream>

#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfStringAttribute.h>

namespace SEXR {

EXRSpectralImage::EXRSpectralImage(
    size_t width, size_t height,
    const std::vector<float>& wavelengths_nm,
    SpectrumType type
)
    : SpectralImage(width, height, wavelengths_nm, type)
{
}


EXRSpectralImage::EXRSpectralImage(
    const std::string& filename
) 
    : SpectralImage()
{
    Imf::InputFile exrIn(filename.c_str());
    const Imf::Header& exrHeader = exrIn.header();
    const Imath::Box2i& exrDataWindow = exrHeader.dataWindow();
    
    _width  = exrDataWindow.max.x - exrDataWindow.min.x + 1;
    _height = exrDataWindow.max.y - exrDataWindow.min.y + 1;
    _spectrumType = UNDEFINED;

    // -----------------------------------------------------------------------
    // Determine channels' position
    // -----------------------------------------------------------------------

    const Imf::ChannelList& exrChannels = exrHeader.channels();

    std::array<std::vector<std::pair<float, std::string>>, 4> wavelengths_nm_S;
    std::array<std::vector<std::pair<float, std::string>>, 16> wavelengths_nm_M;

    for (Imf::ChannelList::ConstIterator channel = exrChannels.begin(); 
        channel != exrChannels.end(); channel++) {
        // Check if the channel is a spectral one
        int polarisationComponent;
        float wavelength_nm;
        SpectrumType spectralChanel = channelType(channel.name(), polarisationComponent, wavelength_nm);

        if (spectralChanel != UNDEFINED) {
            _spectrumType = _spectrumType | spectralChanel;
            
            if (spectralChanel & EMISSIVE) {
                wavelengths_nm_S[polarisationComponent].push_back(std::make_pair(wavelength_nm, channel.name()));
            } else if (spectralChanel & REFLECTIVE) {
                wavelengths_nm_M[polarisationComponent].push_back(std::make_pair(wavelength_nm, channel.name()));
            }
        }
    }

    // Sort by ascending wavelengths
    for (size_t s = 0; s < nStokesComponents(); s++) {
        std::sort(wavelengths_nm_S[s].begin(), wavelengths_nm_S[s].end());
    }

    for (size_t s = 0; s < nMuellerComponents(); s++) {
        std::sort(wavelengths_nm_M[s].begin(), wavelengths_nm_M[s].end());
    }

    // -------------------------------------------------------------------------
    // Sanity check
    // -------------------------------------------------------------------------

    if (_spectrumType == SpectrumType::UNDEFINED) {
        // Probably an RGB EXR, not our job to handle it
        throw INCORRECT_FORMED_FILE;
    }

    if (emissive())
    {
        // Check we have the same wavelength for each Stokes component
        // Wavelength vectors must be of the same size
        const float base_size_emissive = wavelengths_nm_S[0].size();

        for (size_t s = 1; s < nStokesComponents(); s++) {
            if (wavelengths_nm_S[s].size() != base_size_emissive) {
                throw INCORRECT_FORMED_FILE;
            }

            // Wavelengths must correspond
            for (size_t wl_idx = 0; wl_idx < base_size_emissive; wl_idx++) { 
                if (wavelengths_nm_S[s][wl_idx].first != wavelengths_nm_S[0][wl_idx].first) {
                    throw INCORRECT_FORMED_FILE;
                }
            }
        }
    }

    if (reflective())
    {
        // Check we have the same wavelength for each Mueller component
        // Wavelength vectors must be of the same size
        const float base_size_reflective = wavelengths_nm_M[0].size();

        for (size_t s = 1; s < nMuellerComponents(); s++) {
            if (wavelengths_nm_M[s].size() != base_size_reflective) {
                throw INCORRECT_FORMED_FILE;
            }

            // Wavelengths must correspond
            for (size_t wl_idx = 0; wl_idx < base_size_reflective; wl_idx++) { 
                if (wavelengths_nm_M[s][wl_idx].first != wavelengths_nm_M[0][wl_idx].first) {
                    throw INCORRECT_FORMED_FILE;
                }
            }
        }
    }

    
    // If both reflective and emissive, we need to perform a last sanity check
    if (emissive() && reflective())
    {
        const size_t n_emissive_wavelengths = wavelengths_nm_S[0].size();
        const size_t n_reflective_wavelengths = wavelengths_nm_M[0].size();

        if (n_emissive_wavelengths != n_reflective_wavelengths)
            throw INCORRECT_FORMED_FILE;

        for (size_t wl_idx = 0; wl_idx < n_emissive_wavelengths; wl_idx++) {
            if (wavelengths_nm_S[0][wl_idx] != wavelengths_nm_M[0][wl_idx])
                throw INCORRECT_FORMED_FILE;
        }
    }
   
    // -----------------------------------------------------------------------
    // Allocate memory
    // -----------------------------------------------------------------------

    // Now, we can populate the local wavelength vector
    if (emissive()) {
        _wavelengths_nm.reserve(wavelengths_nm_S[0].size());

        for (const auto& wl_index: wavelengths_nm_S[0]) {
            _wavelengths_nm.push_back(wl_index.first);
        }
    } else {
        _wavelengths_nm.reserve(wavelengths_nm_M[0].size());

        for (const auto& wl_index: wavelengths_nm_M[0]) {
            _wavelengths_nm.push_back(wl_index.first);
        }
    }

    // We allocate pixel buffers memory
    if (emissive()) {
        for (size_t s = 0; s < nStokesComponents(); s++) {
            _emissivePixelBuffers[s].resize(nSpectralBands() * _width * _height);
        }
    }

    if (reflective()) {
        for (size_t s = 0; s < nMuellerComponents(); s++) {
            _reflectivePixelBuffers[s].resize(nSpectralBands() * _width * _height);
        }
    }

    // -----------------------------------------------------------------------
    // Read the pixel data
    // -----------------------------------------------------------------------

    Imf::FrameBuffer exrFrameBuffer;

    const Imf::PixelType compType = Imf::FLOAT;
    const size_t xStride = sizeof(float) * nSpectralBands();
    const size_t yStride = xStride * _width;

    for (size_t s = 0; s < nStokesComponents(); s++) {
        for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
            char* ptrS = (char*)(&_emissivePixelBuffers[s][wl_idx]);
            exrFrameBuffer.insert(
                wavelengths_nm_S[s][wl_idx].second, 
                Imf::Slice(compType, ptrS, xStride, yStride));
        }
    }

    for (size_t s = 0; s < nMuellerComponents(); s++) {
        for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
            char* ptrS = (char*)(&_reflectivePixelBuffers[s][wl_idx]);
            exrFrameBuffer.insert(
                wavelengths_nm_M[s][wl_idx].second, 
                Imf::Slice(compType, ptrS, xStride, yStride));
        }
    }

    exrIn.setFrameBuffer(exrFrameBuffer);
    exrIn.readPixels(exrDataWindow.min.y, exrDataWindow.max.y);

    // -----------------------------------------------------------------------
    // Read metadata
    // -----------------------------------------------------------------------

    // Lens transmission data
    const Imf::StringAttribute * lensTransmissionAttr = exrHeader.findTypedAttribute<Imf::StringAttribute>(LENS_TRANSMISSION_ATTR);
    
    if (lensTransmissionAttr != nullptr) {
        try {
            _lensTransmissionSpectra = SpectrumAttribute(*lensTransmissionAttr);
        } catch (SpectrumAttribute::Error &e) {
            throw INCORRECT_FORMED_FILE;
        }
    }

    // Camera spectral response
    const Imf::StringAttribute * cameraResponseAttr = exrHeader.findTypedAttribute<Imf::StringAttribute>(CAMERA_RESPONSE_ATTR);
    
    if (cameraResponseAttr != nullptr) {
        try {
            _cameraReponse = SpectrumAttribute(*cameraResponseAttr);
        } catch (SpectrumAttribute::Error &e) {
            throw INCORRECT_FORMED_FILE;
        }
    }

    // Each channel sensitivity
    _channelSensitivities.resize(nSpectralBands());

    for (size_t i = 0; i < wavelengths_nm_S[0].size(); i++) {
        const Imf::StringAttribute * filterTransmissionAttr = exrHeader.findTypedAttribute<Imf::StringAttribute>(wavelengths_nm_S[0][i].second);
        
        if (filterTransmissionAttr != nullptr) {
            try {
                _channelSensitivities[i] = SpectrumAttribute(*filterTransmissionAttr);
            } catch (SpectrumAttribute::Error& e) {
                throw INCORRECT_FORMED_FILE;
            }
        }
    }
}


void EXRSpectralImage::save(const std::string& filename) 
const {
    Imf::Header exrHeader(width(), height());
    Imf::ChannelList & exrChannels = exrHeader.channels();

    // -----------------------------------------------------------------------
    // Write the pixel data
    // -----------------------------------------------------------------------

    // Layout framebuffer
    Imf::FrameBuffer exrFrameBuffer;
    const Imf::PixelType compType = Imf::FLOAT;

    // Write RGB version
    std::vector<float> rgbImage;
    getRGBImage(rgbImage);

    const std::array<std::string, 3> rgbChannels = {"R", "G", "B"};
    const size_t xStrideRGB = sizeof(float) * 3;
    const size_t yStrideRGB = xStrideRGB * width();

    for (size_t c = 0; c < 3; c++) {
        char* ptrRGB = (char*)(&rgbImage[c]);
        exrChannels.insert(rgbChannels[c], Imf::Channel(compType));
        exrFrameBuffer.insert(
                rgbChannels[c], 
                Imf::Slice(compType, ptrRGB, xStrideRGB, yStrideRGB));
    }

    // Write spectral version
    const size_t xStride = sizeof(float) * nSpectralBands();
    const size_t yStride = xStride * width();

    for (size_t s = 0; s < nStokesComponents(); s++) {
        for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
            // Populate channel name
            std::string channelName = getStokesChannelName(s, _wavelengths_nm[wl_idx]);
            exrChannels.insert(channelName, Imf::Channel(compType));

            char* ptrS = (char*)(&_emissivePixelBuffers[s][wl_idx]);
            exrFrameBuffer.insert(
                channelName, 
                Imf::Slice(compType, ptrS, xStride, yStride));
        }
    }

    for (size_t m = 0; m < nMuellerComponents(); m++) {
        for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
            // Populate channel name
            std::string channelName = getMuellerChannelName(m, _wavelengths_nm[wl_idx]);
            exrChannels.insert(channelName, Imf::Channel(compType));

            char* ptrS = (char*)(&_reflectivePixelBuffers[m][wl_idx]);
            exrFrameBuffer.insert(
                channelName, 
                Imf::Slice(compType, ptrS, xStride, yStride));
        }
    }

    // -----------------------------------------------------------------------
    // Write metadata
    // -----------------------------------------------------------------------

    if (lensTransmission().size() > 0) {
        exrHeader.insert(LENS_TRANSMISSION_ATTR, lensTransmission().getAttribute());
    }

    if (cameraResponse().size() > 0) {
        exrHeader.insert(CAMERA_RESPONSE_ATTR, cameraResponse().getAttribute());
    }

    if (channelSensitivities().size() > 0) {
        for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
            if (channelSensitivity(wl_idx).size() > 0) {
                std::string channelName = getStokesChannelName(0, _wavelengths_nm[wl_idx]);

                exrHeader.insert(channelName, channelSensitivity(wl_idx).getAttribute());
            }
        }
    }

    // -----------------------------------------------------------------------
    // Write file
    // -----------------------------------------------------------------------
    
    Imf::OutputFile exrOut(filename.c_str(), exrHeader);
    exrOut.setFrameBuffer(exrFrameBuffer);
    exrOut.writePixels(height());
}


SpectrumType EXRSpectralImage::channelType(
    const std::string& channelName,
    int& polarisationComponent,
    float& wavelength_nm
) {
    const std::regex expr
        ("^((S([0-3]))|(M([0-3])([0-3])))\\.(\\d*,?\\d*([Ee][+-]?\\d+)?)(Y|Z|E|P|T|G|M|k|h|da|d|c|m|u|n|p)?(m|Hz)$");
    std::smatch matches;

    const bool matched = std::regex_search(channelName, matches, expr);

    SpectrumType channelType = SpectrumType::UNDEFINED;

    if (matched) {
        if (matches.size() != 11) {
            // Something went wrong with the parsing. This shall not occur.
            throw INTERNAL_ERROR;
        }

        if (matches[1].str()[0] == 'S') {
            channelType = SpectrumType::EMISSIVE;
            polarisationComponent = std::stoi(matches[3].str());
        } else if (matches[1].str()[0] == 'M') {
            channelType = SpectrumType::REFLECTIVE;
            const size_t row = std::stoi(matches[5].str());
            const size_t col = std::stoi(matches[6].str());
            polarisationComponent = indexFromComponents(row, col);
        } else {
            throw INTERNAL_ERROR;
        }
        
        if (polarisationComponent > 0) {
            channelType = channelType | SpectrumType::POLARISED;
        }

        // Get value
        std::string centralValueStr(matches[7].str());
        std::replace(centralValueStr.begin(), centralValueStr.end(), ',', '.');
        const float value = std::stof(centralValueStr);
        
        // Convert to nanometers
        const std::string prefix = matches[9].str();
        const std::string units = matches[10].str();

        try {
            wavelength_nm = Util::strToNanometers(value, prefix, units);
        } catch (std::out_of_range& exception) {
            // Unknown unit or multiplier
            // Something went wrong with the parsing. This shall not occur.
            throw INTERNAL_ERROR;
        }
    }

    return channelType;
}


std::string EXRSpectralImage::getStokesChannelName(
    int stokesComponent,
    float wavelength_nm
) {
    std::stringstream b;
    std::string wavelengthStr = std::to_string(wavelength_nm);
    std::replace(wavelengthStr.begin(), wavelengthStr.end(), '.', ',');

    b  << "S" << stokesComponent << "." << wavelengthStr << "nm";

    const std::string channelName = b.str();

    // "Pedantic" check
#ifndef NDEBUG
    int polarisationComponentChecked;
    float wavelength_nmChecked;

    if (!(channelType(channelName, polarisationComponentChecked, wavelength_nmChecked) & EMISSIVE)) {
        throw INTERNAL_ERROR;
    }

    if (polarisationComponentChecked != stokesComponent 
     || wavelength_nmChecked != wavelength_nmChecked) {
        throw INTERNAL_ERROR;
    }
#endif

    return channelName;
}


std::string EXRSpectralImage::getMuellerChannelName(
    int polarisationComponent,
    float wavelength_nm
) {
    std::stringstream b;
    std::string wavelengthStr = std::to_string(wavelength_nm);
    std::replace(wavelengthStr.begin(), wavelengthStr.end(), '.', ',');
 
    size_t row, col;
    componentsFromIndex(polarisationComponent, row, col);
    b << "M" << row << col << "." << wavelengthStr << "nm";

    const std::string channelName = b.str();

    // "Pedantic" check
#ifndef NDEBUG
    int polarisationComponentChecked;
    float wavelength_nmChecked;

    if (!(channelType(channelName, polarisationComponentChecked, wavelength_nmChecked) & REFLECTIVE)) {
        throw INTERNAL_ERROR;
    }

    if (polarisationComponentChecked != polarisationComponent 
     || wavelength_nmChecked != wavelength_nmChecked) {
        throw INTERNAL_ERROR;
    }
#endif

    return channelName;
}

} // namespace SEXR
