#include <EXRBiSpectralImage.h>
#include "Util.h"

#include <regex>
#include <algorithm>
#include <sstream>
#include <cassert>

#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfStringAttribute.h>

namespace SEXR {

EXRBiSpectralImage::EXRBiSpectralImage(
    size_t width, size_t height,
    const std::vector<float>& wavelengths_nm,
    SpectrumType type
) 
    : BiSpectralImage(width, height, wavelengths_nm, type)
{}


EXRBiSpectralImage::EXRBiSpectralImage(
    const std::string& filename
)
    : BiSpectralImage()
{
    Imf::InputFile exrIn(filename.c_str());
    const Imf::Header& exrHeader = exrIn.header();
    const Imath::Box2i& exrDataWindow = exrHeader.dataWindow();

    _width  = exrDataWindow.max.x - exrDataWindow.min.x + 1;
    _height = exrDataWindow.max.y - exrDataWindow.min.y + 1;
    _spectrumType = SpectrumType::UNDEFINED;

    // -----------------------------------------------------------------------
    // Determine channels' position
    // -----------------------------------------------------------------------

    const Imf::ChannelList& exrChannels = exrHeader.channels();

    std::array<std::vector<std::pair<float, std::string>>, 4> wavelengths_nm_S;
    std::array<std::vector<std::pair<float, std::string>>, 16> wavelengths_nm_M;
    std::vector<std::pair<std::pair<float, float>, std::string>> reradiation_wavelengths_nm;

    for (Imf::ChannelList::ConstIterator channel = exrChannels.begin(); 
        channel != exrChannels.end(); channel++) {
        // Check if the channel is a spectral or a bispectral one
        int polarisationComponent;
        float in_wavelength_nm, out_wavelength_nm;

        SpectrumType currChannelType = channelType(
            channel.name(), 
            polarisationComponent, 
            in_wavelength_nm, out_wavelength_nm);
        
        if (currChannelType != SpectrumType::UNDEFINED) {
            _spectrumType = _spectrumType | currChannelType;
        
            if (isReflective(currChannelType)) {
                if (!isBispectral(currChannelType)) {
                    assert(polarisationComponent < 4);

                    wavelengths_nm_M[polarisationComponent].push_back(
                        std::make_pair(
                            in_wavelength_nm, 
                            channel.name()));
                } else {
                    assert(polarisationComponent == 0);

                    reradiation_wavelengths_nm.push_back(
                        std::make_pair(
                            std::make_pair(in_wavelength_nm, out_wavelength_nm), 
                            channel.name()));
                }
            } else if (isEmissive(currChannelType)) {
                assert(polarisationComponent < 16);

                wavelengths_nm_S[polarisationComponent].push_back(
                    std::make_pair(
                        in_wavelength_nm, 
                        channel.name()));
            }
        }
    }

    // Sort by ascending wavelengths
    for (size_t s = 0; s < nStokesComponents(); s++) {
        std::sort(wavelengths_nm_S[s].begin(), wavelengths_nm_S[s].end());
    }

    for (size_t m = 0; m < nMuellerComponents(); m++) {
        std::sort(wavelengths_nm_M[m].begin(), wavelengths_nm_M[m].end());
    }

    struct {
        bool operator()(
            std::pair<std::pair<float, float>, std::string> a, 
            std::pair<std::pair<float, float>, std::string> b) const
        {   
            if (a.first.first == b.first.first) {
                return a.first.second < b.first.second;
            } else {
                return a.first.first < b.first.first;
            }
        }
    } sortBispectral;

    std::sort(reradiation_wavelengths_nm.begin(), reradiation_wavelengths_nm.end(), sortBispectral);

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

    // Check every single reradiation have all upper wavelength values
    // Note: this shall not be mandatory in the format but, we only support that for now
    if (bispectral()) {
        if (reradiation_wavelengths_nm.size() != reradiationSize()) {
            std::cerr << "Reradiation is incomplete" << std::endl;
        }

        float currDiagWl = wavelength_nm(0);
        size_t diagonalIdx = 0;
        size_t reradIdx = 1;

        for (size_t rr = 0; rr < reradiation_wavelengths_nm.size(); rr++) {
            const auto& rerad = reradiation_wavelengths_nm[rr];

            if (rerad.first.first > currDiagWl) {
                if (diagonalIdx + reradIdx != nSpectralBands()) {
                    std::cerr << "We need the full spectral reradiation specification" << std::endl;
                }
                diagonalIdx++;
                reradIdx = 1;

                currDiagWl = wavelength_nm(diagonalIdx);

            }

            if (rerad.first.first != wavelength_nm(diagonalIdx)) {
                std::cerr << "We need the full spectral reradiation specification" << std::endl;
            }

            if (rerad.first.second != wavelength_nm(diagonalIdx + reradIdx)) {
                std::cerr << "We need the full spectral reradiation specification" << std::endl;
            }

            reradIdx++;
        }
    }

    // -----------------------------------------------------------------------
    // Allocate memory
    // -----------------------------------------------------------------------

    for (size_t s = 0; s < nStokesComponents(); s++) {
        _emissivePixelBuffers[s].resize(nSpectralBands() * width() * height());
    }

    for (size_t m = 0; m < nMuellerComponents(); m++) {
        _reflectivePixelBuffers[m].resize(nSpectralBands() * width() * height());
    }

    if (bispectral()) {
        _reradiation.resize(reradiationSize() * width() * height());
    }

    // -----------------------------------------------------------------------
    // Read the pixel data
    // -----------------------------------------------------------------------

    Imf::FrameBuffer exrFrameBuffer;

    const Imf::PixelType compType = Imf::FLOAT;
    
    // Set the diagonal for reading
    const size_t xStride = sizeof(float) * nSpectralBands();
    const size_t yStride = xStride * width();

    for (size_t s = 0; s < nStokesComponents(); s++) {
        for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
            char* ptrS = (char*)(&_emissivePixelBuffers[s][wl_idx]);
            exrFrameBuffer.insert(
                wavelengths_nm_S[s][wl_idx].second, 
                Imf::Slice(compType, ptrS, xStride, yStride));
        }
    }

    for (size_t m = 0; m < nMuellerComponents(); m++) {
        for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
            char* ptrS = (char*)(&_reflectivePixelBuffers[m][wl_idx]);
            exrFrameBuffer.insert(
                wavelengths_nm_M[m][wl_idx].second, 
                Imf::Slice(compType, ptrS, xStride, yStride));
        }
    }

    if (bispectral()) {
        // Set the reradiation part fo reading
        const size_t xStrideReradiation = sizeof(float) * reradiationSize();
        const size_t yStrideReradiation = xStrideReradiation * width();

        for (size_t rr = 0; rr < reradiationSize(); rr++) {
            char* framebuffer = (char*)(&_reradiation[rr]);
            exrFrameBuffer.insert(
                reradiation_wavelengths_nm[reradiationSize() - rr - 1].second, 
                Imf::Slice(compType, framebuffer, xStrideReradiation, yStrideReradiation));
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

    // Exposure compensation value
    const Imf::StringAttribute * exposureCompensationAttr = exrHeader.findTypedAttribute<Imf::StringAttribute>(EXPOSURE_COMPENSATION_ATTR);

    if (exposureCompensationAttr != nullptr) {
        try {
            _ev = std::stof(exposureCompensationAttr->value());
        } catch (std::invalid_argument& e) {
            throw INCORRECT_FORMED_FILE;
        }
    }
}


void EXRBiSpectralImage::save(const std::string& filename) 
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
            const std::string channelName = getStokesChannelName(s, _wavelengths_nm[wl_idx]);
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
            const std::string channelName = getMuellerChannelName(m, _wavelengths_nm[wl_idx]);
            exrChannels.insert(channelName, Imf::Channel(compType));

            char* ptrS = (char*)(&_reflectivePixelBuffers[m][wl_idx]);
            exrFrameBuffer.insert(
                channelName, 
                Imf::Slice(compType, ptrS, xStride, yStride));
        }
    }

    if (bispectral()) 
    {
        // Write the reradiation
        const size_t xStrideReradiation = sizeof(float) * reradiationSize();
        const size_t yStrideReradiation = xStrideReradiation * _width;

        for (size_t rr = 0; rr < reradiationSize(); rr++) {
                size_t wlFromIdx, wlToIdx;
                wavelengthsIdxFromIdx(rr, wlFromIdx, wlToIdx);

                const std::string channelName = getReradiationChannelName(
                    _wavelengths_nm[wlFromIdx],
                    _wavelengths_nm[wlToIdx]);

                exrChannels.insert(channelName, Imf::Channel(compType));

                char* framebuffer = (char*)(&_reradiation[rr]);
                exrFrameBuffer.insert(
                    channelName, 
                    Imf::Slice(compType, framebuffer, xStrideReradiation, yStrideReradiation));
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

    exrHeader.insert(EXPOSURE_COMPENSATION_ATTR, Imf::StringAttribute(std::to_string(_ev)));


    Imf::OutputFile exrOut(filename.c_str(), exrHeader);
    exrOut.setFrameBuffer(exrFrameBuffer);
    exrOut.writePixels(height());
}


SpectrumType EXRBiSpectralImage::channelType(
    const std::string& channelName,
    int& polarisationComponent,
    float& wavelength_nm,
    float& reradiation_wavelength_nm
) {
    const std::string exprMueller = "M([0-3])([0-3])";
    const std::string exprStokes  = "S([0-3])";
    const std::string exprPola    = "((" + exprStokes + ")|(" + exprMueller + "))";
    const std::string exprValue   = "(\\d*,?\\d*([Ee][+-]?\\d+)?)";
    const std::string exprUnits   = "(Y|Z|E|P|T|G|M|k|h|da|d|c|m|u|n|p)?(m|Hz)";

    const std::regex exprDiagonal("^" + exprPola + "\\." + exprValue + exprUnits + "$");
    const std::regex exprRerad   ("^" + exprMueller + "\\." + exprValue + exprUnits + "\\." + exprValue + exprUnits + "$");

    std::smatch matches;

    const bool matchedDiagonal = std::regex_search(channelName, matches, exprDiagonal);

    if (matchedDiagonal) {
        if (matches.size() != 11) {
            // Something went wrong with the parsing. This shall not occur.
            throw INTERNAL_ERROR;
        }
        
        SpectrumType type;
        int row, col;

        switch (matches[1].str()[0])
        {
            case 'S':
                type = SpectrumType::EMISSIVE;
                polarisationComponent = std::stoi(matches[3].str());
                break;

            case 'M':
                type = SpectrumType::REFLECTIVE;
                row                   = std::stoi(matches[5].str());
                col                   = std::stoi(matches[6].str());
                polarisationComponent = indexFromComponents(row, col);
                break;
            default:
                return SpectrumType::UNDEFINED;
        }

        if (polarisationComponent > 0) {
            type = type | SpectrumType::POLARISED;
        }

        // Get value illumination
        std::string centralValueStr(matches[7].str());
        std::replace(centralValueStr.begin(), centralValueStr.end(), ',', '.');
        const float value = std::stof(centralValueStr);

        wavelength_nm = Util::strToNanometers(
            value, // Comma separated floating point value
            matches[9].str(), // Unit multiplier
            matches[10].str()  // Units
            );

        return type;
    }

    const bool matchedRerad = std::regex_search(channelName, matches, exprRerad);

    if (matchedRerad) {
        if (matches.size() != 11) {
            // Something went wrong with the parsing. This shall not occur.
            throw INTERNAL_ERROR;
        }

        size_t row = std::stoi(matches[1].str());
        size_t col = std::stoi(matches[2].str());
        
        polarisationComponent = indexFromComponents(row, col);
        
        // Reradiation is depolarising...
        if (polarisationComponent > 0) {
            return SpectrumType::UNDEFINED;
        }

        // Get value illumination
        std::string centralValueStrI(matches[3].str());
        std::replace(centralValueStrI.begin(), centralValueStrI.end(), ',', '.');
        const float value_i = std::stof(centralValueStrI);

        wavelength_nm = Util::strToNanometers(
            value_i,             
            matches[5].str(), // Unit multiplier
            matches[6].str()  // Units
            );

        // Get value reradiation
        std::string centralValueStrO(matches[7].str());
        std::replace(centralValueStrO.begin(), centralValueStrO.end(), ',', '.');
        const float value_o = std::stof(centralValueStrO);

        reradiation_wavelength_nm = Util::strToNanometers(
            value_o,
            matches[9].str(), // Unit multiplier
            matches[10].str() // Units
            );

        return SpectrumType::BISPECTRAL;
    }

    return SpectrumType::UNDEFINED;
}


std::string EXRBiSpectralImage::getStokesChannelName(
    int stokesComponent,
    float wavelength_nm
) {
    assert(stokesComponent < 4);

    std::stringstream b;
    std::string wavelengthStr = std::to_string(wavelength_nm);
    std::replace(wavelengthStr.begin(), wavelengthStr.end(), '.', ',');

    b  << "S" << stokesComponent << "." << wavelengthStr << "nm";

    const std::string channelName = b.str();

#ifndef NDEBUG
    int polarisationComponentChecked;
    float wavelength_nmChecked;
    float wavelength_nmCheckedOut;
    SpectrumType t = channelType(channelName, polarisationComponentChecked, wavelength_nmChecked, wavelength_nmCheckedOut);
    
    assert(isEmissive(t));
    assert(polarisationComponentChecked == stokesComponent);
    assert(wavelength_nmChecked == wavelength_nm);
#endif

    return channelName;
}


std::string EXRBiSpectralImage::getMuellerChannelName(
    int polarisationComponent,
    float wavelength_nm
) {
    assert(polarisationComponent < 16);

    std::stringstream b;
    std::string wavelengthStr = std::to_string(wavelength_nm);
    std::replace(wavelengthStr.begin(), wavelengthStr.end(), '.', ',');
 
    size_t row, col;
    componentsFromIndex(polarisationComponent, row, col);
    b << "M" << row << col << "." << wavelengthStr << "nm";

    const std::string channelName = b.str();

#ifndef NDEBUG
    int polarisationComponentChecked;
    float wavelength_nmChecked;
    float wavelength_nmCheckedOut;
    SpectrumType t = channelType(channelName, polarisationComponentChecked, wavelength_nmChecked, wavelength_nmCheckedOut);
    
    assert(isReflective(t));
    assert(polarisationComponentChecked == polarisationComponent);
    assert(wavelength_nmChecked == wavelength_nm);
#endif

    return channelName;
}


std::string EXRBiSpectralImage::getReradiationChannelName(
    float wavelength_nm,
    float reradiation_wavelength_nm
) {
    std::string reradWavelengthStr = std::to_string(reradiation_wavelength_nm);
    std::replace(reradWavelengthStr.begin(), reradWavelengthStr.end(), '.', ',');

    std::stringstream b;
    b << getMuellerChannelName(0, wavelength_nm) << '.' << reradWavelengthStr << "nm";

    const std::string channelName = b.str();

#ifndef NDEBUG
    int muellerComponent;
    float wavelength_nmChecked;
    float reradiation_wavelength_nmChecked;

    SpectrumType t = channelType(channelName, muellerComponent, wavelength_nmChecked, reradiation_wavelength_nmChecked);
    
    assert(isBispectral(t));
    assert(muellerComponent == 0);
    assert(wavelength_nmChecked == wavelength_nm);
    assert(reradiation_wavelength_nmChecked == reradiation_wavelength_nm);
#endif

    return channelName;
}

} // namespace SEXR
