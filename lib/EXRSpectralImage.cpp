#include <EXRSpectralImage.h>

#include <regex>
#include <algorithm>
#include <sstream>

#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>


EXRSpectralImage::EXRSpectralImage(
    size_t width, size_t height,
    const std::vector<float>& wavelengths_nm,
    bool containsPolarisationData
)
    : SpectralImage(width, height, wavelengths_nm, containsPolarisationData)
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

    // Determine position of the channels
    const Imf::ChannelList& exrChannels = exrHeader.channels();

    std::array<std::vector<std::pair<float, std::string>>, 4> wavelengths_nm_S;

    for (Imf::ChannelList::ConstIterator channel = exrChannels.begin(); 
        channel != exrChannels.end(); channel++) {
        // Check if the channel is a spectral one
        int stokes;
        float wavelength_nm;
        bool spectralChanel = isSpectralChannel(channel.name(), stokes, wavelength_nm);

        if (spectralChanel) {
            _isSpectral = true;

            if (stokes > 0) {
                _containsPolarisationData = true;
            }
            
            wavelengths_nm_S[stokes].push_back(std::make_pair(wavelength_nm, channel.name()));
        }
    }

    // Sort by ascending wavelengths
    for (int s = 0; s < nStokesComponents(); s++) {
        std::sort(wavelengths_nm_S[s].begin(), wavelengths_nm_S[s].end());
    }

    // Check we have the same wavelength for each stoke component
    if (_containsPolarisationData) {
        // Wavelength vectors must be of the same size
        const float base_size = wavelengths_nm_S[0].size();
        for (size_t s = 1; s < 4; s++) {
            if (wavelengths_nm_S[s].size() != base_size) {
                throw INCORECTED_FORMED_FILE;
            }
        }

        // Wavelengths must correspond
        for (size_t wl_idx = 0; wl_idx < wavelengths_nm_S[0].size(); wl_idx++) {
            const float base_wl = wavelengths_nm_S[0][wl_idx].first;
            for (size_t s = 1; s < 4; s++) {
                if (wavelengths_nm_S[s][wl_idx].first != base_wl) {
                    throw INCORECTED_FORMED_FILE;
                }
            }
        }
    }

    // Now, we can populate the local wavelength vector
    _wavelengths_nm.reserve(wavelengths_nm_S[0].size());

    for (const auto& wl_index: wavelengths_nm_S[0]) {
        _wavelengths_nm.push_back(wl_index.first);
    }


    // We allocate _pixelBuffer memory
    for (size_t s = 0; s < nStokesComponents(); s++) {
        _pixelBuffers[s].resize(nSpectralBands() * _width * _height);
    }

    // Then, we read the pixel data
    Imf::FrameBuffer exrFrameBuffer;

    if (_isSpectral) {
        const Imf::PixelType compType = Imf::FLOAT;
        const size_t xStride = sizeof(float) * nSpectralBands();
        const size_t yStride = xStride * _width;

        for (size_t s = 0; s < nStokesComponents(); s++) {
            for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
                char* ptrS = (char*)(&_pixelBuffers[s][wl_idx]);
                exrFrameBuffer.insert(
                    wavelengths_nm_S[s][wl_idx].second, 
                    Imf::Slice(compType, ptrS, xStride, yStride));
            }
        }
    } else {
        // Probably an RGB EXR, not our job to handle it
        throw INCORECTED_FORMED_FILE;
    }

    exrIn.setFrameBuffer(exrFrameBuffer);
    exrIn.readPixels(exrDataWindow.min.y, exrDataWindow.max.y);
}


void EXRSpectralImage::save(const std::string& filename) const {
    Imf::Header exrHeader(width(), height());
    Imf::ChannelList & exrChannels = exrHeader.channels();

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
            std::string channelName = getChannelName(s, _wavelengths_nm[wl_idx]);
            exrChannels.insert(channelName, Imf::Channel(compType));

            char* ptrS = (char*)(&_pixelBuffers[s][wl_idx]);
            exrFrameBuffer.insert(
                channelName, 
                Imf::Slice(compType, ptrS, xStride, yStride));
        }
    }

    Imf::OutputFile exrOut(filename.c_str(), exrHeader);
    exrOut.setFrameBuffer(exrFrameBuffer);
    exrOut.writePixels(height());
}


bool EXRSpectralImage::isSpectralChannel(
    const std::string& channelName,
    int& stokesComponent,
    float& wavelength_nm
) {
    const std::map<std::string, float> unit_prefix = {
        {"Y", 1e24}, {"Z", 1e21}, {"E", 1e18}, {"P", 1e15}, {"T", 1e12},
        {"G", 1e9} , {"M", 1e6} , {"k", 1e3} , {"h", 1e2} , {"da", 1e1},
        {"d", 1e-1}, {"c", 1e-2}, {"m", 1e-3}, {"u", 1e-6}, {"n", 1e-9},
        {"p", 1e-12}
    };
    
    const std::regex expr
        ("^S([0-3])\\.(\\d*,?\\d*([Ee][+-]?\\d+)?)(Y|Z|E|P|T|G|M|k|h|da|d|c|m|u|n|p)?(m|Hz)$");
    std::smatch matches;

    const bool matched = std::regex_search(channelName, matches, expr);

    if (matched) {
        if (matches.size() != 6) {
            // Something went wrong with the parsing. This shall not occur.
            throw INTERNAL_ERROR;
        }
        
        stokesComponent = std::stoi(matches[1].str());

        // Get value
        std::string centralValueStr(matches[2].str());
        std::replace(centralValueStr.begin(), centralValueStr.end(), ',', '.');
        float value = std::stof(centralValueStr);
        
        // Apply multiplier
        const std::string prefix = matches[4].str();
        
        if (prefix.size() > 0) {
            try {
                value *= unit_prefix.at(prefix);
            } catch (std::out_of_range& exception) {
                // Unknown unit multiplier
                // Something went wrong with the parsing. This shall not occur.
                throw INTERNAL_ERROR;
            }
        }

        // Apply units
        const std::string units = matches[5].str();
        
        if (units == "Hz") {
            wavelength_nm = 299792458.F/value * 1e9;
        } else if (units == "m") {
            wavelength_nm = value * 1e9;
        } else {
            // Unknown unit
            // Something went wrong with the parsing. This shall not occur.
            throw INTERNAL_ERROR;
        }		
    }

    return matched;
}


std::string EXRSpectralImage::getChannelName(
    int stokesComponent,
    float wavelength_nm
) {
    std::stringstream b;
    std::string wavelengthStr = std::to_string(wavelength_nm);
    std::replace(wavelengthStr.begin(), wavelengthStr.end(), '.', ',');

    b  << "S" << stokesComponent << "." << wavelengthStr << "nm";

    const std::string channelName = b.str();

    // "Pedantic" check
    int stokesComponentChecked;
    float wavelength_nmChecked;

    if (!isSpectralChannel(channelName, stokesComponentChecked, wavelength_nmChecked)) {
        throw INTERNAL_ERROR;
    }

    if (stokesComponentChecked != stokesComponent 
     || wavelength_nmChecked != wavelength_nmChecked) {
        throw INTERNAL_ERROR;
    }

    return channelName;
}