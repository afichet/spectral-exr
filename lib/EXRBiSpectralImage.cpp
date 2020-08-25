#include <EXRBiSpectralImage.h>

#include <regex>
#include <algorithm>
#include <sstream>

#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>


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

    // -----------------------------------------------------------------------
    // Determine channels' position
    // -----------------------------------------------------------------------

    const Imf::ChannelList& exrChannels = exrHeader.channels();

    std::array<std::vector<std::pair<float, std::string>>, 4> diagonal_wavelengths_nm;
    std::vector<std::pair<std::pair<float, float>, std::string>> reradiation_wavelengths_nm;

    for (Imf::ChannelList::ConstIterator channel = exrChannels.begin(); 
        channel != exrChannels.end(); channel++) {
        // Check if the channel is a spectral or a bispectral one
        int muellerComponent;
        float in_wavelength_nm, out_wavelength_nm;

        ChannelType currChannelType = channelType(
            channel.name(), 
            muellerComponent, 
            in_wavelength_nm, out_wavelength_nm);
        
        if (   currChannelType == DIAGONAL
            || currChannelType == RERADIATION) {
            
            if (muellerComponent > 0) {
                _containsPolarisationData = true;
            }
        }

        switch(currChannelType) {
            case DIAGONAL:
                diagonal_wavelengths_nm[muellerComponent].push_back(
                    std::make_pair(
                        in_wavelength_nm, 
                        channel.name()));
            break;

            case RERADIATION:
                if (muellerComponent != 0) throw INTERNAL_ERROR;

                reradiation_wavelengths_nm.push_back(
                    std::make_pair(
                        std::make_pair(in_wavelength_nm, out_wavelength_nm), 
                        channel.name())
                    );
            break;

            default:
            break;
        }
    }

    // Sort by ascending wavelengths
    for (size_t s = 0; s < nPolarsiationComponents(); s++) {
        std::sort(diagonal_wavelengths_nm[s].begin(), diagonal_wavelengths_nm[s].end());
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

    // Check we have the same wavelength for each stoke component on the main diagonal
    if (_containsPolarisationData) {
        // Wavelength vectors must be of the same size
        const float base_size = diagonal_wavelengths_nm[0].size();
        for (size_t s = 1; s < 4; s++) {
            if (diagonal_wavelengths_nm[s].size() != base_size) {
                throw INCORRECT_FORMED_FILE;
            }
        }

        // Wavelengths must correspond
        for (size_t wl_idx = 0; wl_idx < diagonal_wavelengths_nm[0].size(); wl_idx++) {
            const float base_wl = diagonal_wavelengths_nm[0][wl_idx].first;
            for (size_t s = 1; s < 4; s++) {
                if (diagonal_wavelengths_nm[s][wl_idx].first != base_wl) {
                    throw INCORRECT_FORMED_FILE;
                }
            }
        }
    }

    // Now, we can populate the local wavelength vector
    _wavelengths_nm.reserve(diagonal_wavelengths_nm[0].size());

    for (const auto& wl_index: diagonal_wavelengths_nm[0]) {
        _wavelengths_nm.push_back(wl_index.first);
    }

    // Check every single reradiation have all upper wavelength values
    // Note: this shall not be mandatory in the format but, we only support that for now
    if (reradiation_wavelengths_nm.size() > 0) {
        float currWl = reradiation_wavelengths_nm[0].first.first;
        size_t wlReradCount = 0;
        size_t wlIdx = 0;
        for (const auto& rr: reradiation_wavelengths_nm) {
            if (rr.first.first == currWl) {
                if (rr.first.first != wavelength_nm(wlIdx)) {
                    std::cerr << "We need full spectral reradiation specification" << std::endl;
                }
                wlReradCount++;
            } else {
                if (wlReradCount != nSpectralBands() - wlIdx - 1) {
                    std::cerr << "We need full spectral reradiation specification" << std::endl;
                }
                wlIdx++;
                wlReradCount = 0;
            }
        }
    }

    // Allocate memory

    for (size_t s = 0; s < nPolarsiationComponents(); s++) {
        _pixelBuffers[s].resize(nSpectralBands() * _width * _height);
    }

    _reradiation.resize(reradiationSize() * _width * _height);

    // -----------------------------------------------------------------------
    // Read the pixel data
    // -----------------------------------------------------------------------

    Imf::FrameBuffer exrFrameBuffer;

    const Imf::PixelType compType = Imf::FLOAT;
    
    // Set the diagonal for reading
    const size_t xStrideDiagonal = sizeof(float) * nSpectralBands();
    const size_t yStrideDiagnoal = xStrideDiagonal * _width;

    for (size_t s = 0; s < nPolarsiationComponents(); s++) {
        for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
            char* framebuffer = (char*)(&_pixelBuffers[s][wl_idx]);
            exrFrameBuffer.insert(
                diagonal_wavelengths_nm[s][wl_idx].second, 
                Imf::Slice(compType, framebuffer, xStrideDiagonal, yStrideDiagnoal));
        }
    }

    // Set the reradiation part fo reading
    const size_t xStrideReradiation = sizeof(float) * reradiationSize();
    const size_t yStrideReradiation = xStrideReradiation * _width;

    for (size_t rr = 0; rr < reradiationSize(); rr++) {
            char* framebuffer = (char*)(&_reradiation[rr]);
            exrFrameBuffer.insert(
                reradiation_wavelengths_nm[rr].second, 
                Imf::Slice(compType, framebuffer, xStrideReradiation, yStrideReradiation));
    }

    exrIn.setFrameBuffer(exrFrameBuffer);
    exrIn.readPixels(exrDataWindow.min.y, exrDataWindow.max.y);
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

#warning This shall be adapted for bi-spectral
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
    const size_t xStrideDiagonal = sizeof(float) * nSpectralBands();
    const size_t yStrideDiagonal = xStrideDiagonal * width();

    for (size_t s = 0; s < nPolarsiationComponents(); s++) {
        for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
            // Populate channel name
            std::string channelName = getDiagonalChannelName(s, _wavelengths_nm[wl_idx]);
            exrChannels.insert(channelName, Imf::Channel(compType));

            char* ptrS = (char*)(&_pixelBuffers[s][wl_idx]);
            exrFrameBuffer.insert(
                channelName, 
                Imf::Slice(compType, ptrS, xStrideDiagonal, yStrideDiagonal));
        }
    }

    // Write the reradiation

}


EXRBiSpectralImage::ChannelType EXRBiSpectralImage::channelType(
    const std::string& channelName,
    int& muellerComponent,
    float& wavelength_nm,
    float& reradiation_wavelength_nm
) const {
    const std::string exprStokes = "M([0-3])([0-3])";
    const std::string exprWave   = "(\\d*,?\\d*([Ee][+-]?\\d+)?)";
    const std::string exprUnits  = "(Y|Z|E|P|T|G|M|k|h|da|d|c|m|u|n|p)?(m|Hz)";

    const std::regex exprDiagonal("^" + exprStokes + "\\." + exprWave + exprUnits + "$");
    const std::regex exprRerad   ("^" + exprStokes + "\\." + exprWave + exprUnits + "\\." + exprWave + exprUnits + "$");

    std::smatch matches;

    const bool matchedDiagonal = std::regex_search(channelName, matches, exprDiagonal);

    if (matchedDiagonal) {
        if (matches.size() != 7) {
            // Something went wrong with the parsing. This shall not occur.
            throw INTERNAL_ERROR;
        }
        
        size_t row = std::stoi(matches[1].str());
        size_t col = std::stoi(matches[2].str());
        
        muellerComponent = indexFromComponents(row, col);

        wavelength_nm = toWavelength_nm(
            matches[3].str(), // Comma separated floating point value
            matches[5].str(), // Unit multiplier
            matches[6].str()  // Units
            );

        return DIAGONAL;	
    }

    const bool matchedRerad = std::regex_search(channelName, matches, exprRerad);

    if (matchedRerad) {
        if (matches.size() != 11) {
            // Something went wrong with the parsing. This shall not occur.
            throw INTERNAL_ERROR;
        }

        size_t row = std::stoi(matches[1].str());
        size_t col = std::stoi(matches[2].str());
        
        muellerComponent = indexFromComponents(row, col);
        
        // TODO: double check that. But reradiation shall not be polarising I think...
        if (muellerComponent != 0) return OTHER;

        wavelength_nm = toWavelength_nm(
            matches[3].str(), // Comma separated floating point value
            matches[5].str(), // Unit multiplier
            matches[6].str()  // Units
            );

        reradiation_wavelength_nm = toWavelength_nm(
            matches[7].str(), // Comma separated floating point value
            matches[9].str(), // Unit multiplier
            matches[10].str() // Units
            );

        return RERADIATION;
    }

    return OTHER;
}


float EXRBiSpectralImage::toWavelength_nm(
    const std::string& valueStr,
    const std::string& multiplierStr,
    const std::string& unitStr
) {
    float wavelength_nm;

    // Get value
    std::string localValueStr = valueStr;
    std::replace(localValueStr.begin(), localValueStr.end(), ',', '.');
    float value = std::stof(localValueStr);
    
    // Apply multiplier
    const std::map<std::string, float> unit_prefix = {
        {"Y", 1e24}, {"Z", 1e21}, {"E", 1e18}, {"P", 1e15}, {"T", 1e12},
        {"G", 1e9} , {"M", 1e6} , {"k", 1e3} , {"h", 1e2} , {"da", 1e1},
        {"d", 1e-1}, {"c", 1e-2}, {"m", 1e-3}, {"u", 1e-6}, {"n", 1e-9},
        {"p", 1e-12}
    };

    if (multiplierStr.size() > 0) {
        try {
            value *= unit_prefix.at(multiplierStr);
        } catch (std::out_of_range& exception) {
            // Unknown unit multiplier
            // Something went wrong with the parsing. This shall not occur.
            throw INTERNAL_ERROR;
        }
    }

    // Apply units
    if (unitStr == "Hz") {
        wavelength_nm = 299792458.F/value * 1e9;
    } else if (unitStr == "m") {
        wavelength_nm = value * 1e9;
    } else {
        // Unknown unit
        // Something went wrong with the parsing. This shall not occur.
        throw INTERNAL_ERROR;
    }

    return wavelength_nm;
}

        
std::string EXRBiSpectralImage::getDiagonalChannelName(
    int muellerComponent,
    float wavelength_nm
) const {
    std::stringstream b;
    std::string wavelengthStr = std::to_string(wavelength_nm);
    std::replace(wavelengthStr.begin(), wavelengthStr.end(), '.', ',');

    b  << "S" << muellerComponent << "." << wavelengthStr << "nm";

    const std::string channelName = b.str();

    // "Pedantic" check
    int muellerComponentChecked;
    float wavelength_nmChecked;
    float reradiation_nm;

    if (channelType(
        channelName, 
        muellerComponentChecked, 
        wavelength_nmChecked, 
        reradiation_nm) != DIAGONAL) {
        throw INTERNAL_ERROR;
    }

    if (muellerComponentChecked != muellerComponent 
     || wavelength_nmChecked != wavelength_nmChecked) {
        throw INTERNAL_ERROR;
    }

    return channelName;
}


std::string EXRBiSpectralImage::getReradiationChannelName(
    float wavelength_nm,
    float reradiation_wavelength_nm
) const {
    std::string reradWavelengthStr = std::to_string(reradiation_wavelength_nm);
    std::replace(reradWavelengthStr.begin(), reradWavelengthStr.end(), '.', ',');

    std::stringstream b;
    b << getChannelName(0, wavelength_nm) << '.' << reradWavelengthStr << "nm";

    const std::string channelName = b.str();

    // "Pedantic" check
    int stokesComponentChecked;
    float wavelength_nmChecked;
    float reradiation_wavelength_nmChecked;

    if (channelType(
        channelName, 
        stokesComponentChecked, 
        wavelength_nmChecked, 
        reradiation_wavelength_nmChecked) != RERADIATION) {
        throw INTERNAL_ERROR;
    }

    if (stokesComponentChecked != 0 
     || wavelength_nmChecked != wavelength_nmChecked
     || reradiation_wavelength_nmChecked != reradiation_wavelength_nm) {
        throw INTERNAL_ERROR;
    }

    return channelName;
}