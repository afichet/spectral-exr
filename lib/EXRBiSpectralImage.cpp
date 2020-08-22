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

    // Determine position of the channels
    const Imf::ChannelList& exrChannels = exrHeader.channels();

    std::array<std::vector<std::pair<float, std::string>>, 4> diagonal_wavelengths_nm_S;
    std::vector<std::pair<std::pair<float, float>, std::string>> reradiation_wavelengths_nm;

    for (Imf::ChannelList::ConstIterator channel = exrChannels.begin(); 
        channel != exrChannels.end(); channel++) {
        // Check if the channel is a spectral or a bispectral one
        int stokes;
        float in_wavelength_nm, out_wavelength_nm;

        CHANNEL_TYPE currChannelType = channelType(
            channel.name(), 
            stokes, in_wavelength_nm, out_wavelength_nm);
        
        if (   currChannelType == SPECTRAL_DIAGONAL
            || currChannelType == SPECTRAL_RERADIATION) {
            _isSpectral = true;
            
            if (stokes > 0) {
                _containsPolarisationData = true;
            }
        }

        switch(currChannelType) {
            case SPECTRAL_DIAGONAL:
                diagonal_wavelengths_nm_S[stokes].push_back(
                    std::make_pair(
                        in_wavelength_nm, 
                        channel.name()));
            break;

            case SPECTRAL_RERADIATION:
                if (stokes != 0) throw INTERNAL_ERROR;

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
    for (size_t s = 0; s < nStokesComponents(); s++) {
        std::sort(diagonal_wavelengths_nm_S[s].begin(), diagonal_wavelengths_nm_S[s].end());
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
        const float base_size = diagonal_wavelengths_nm_S[0].size();
        for (size_t s = 1; s < 4; s++) {
            if (diagonal_wavelengths_nm_S[s].size() != base_size) {
                throw INCORRECT_FORMED_FILE;
            }
        }

        // Wavelengths must correspond
        for (size_t wl_idx = 0; wl_idx < diagonal_wavelengths_nm_S[0].size(); wl_idx++) {
            const float base_wl = diagonal_wavelengths_nm_S[0][wl_idx].first;
            for (size_t s = 1; s < 4; s++) {
                if (diagonal_wavelengths_nm_S[s][wl_idx].first != base_wl) {
                    throw INCORRECT_FORMED_FILE;
                }
            }
        }
    }

    // Now, we can populate the local wavelength vector
    _wavelengths_nm.reserve(diagonal_wavelengths_nm_S[0].size());
    _reradiation.reserve(diagonal_wavelengths_nm_S[0].size());

    for (const auto& wl_index: diagonal_wavelengths_nm_S[0]) {
        _wavelengths_nm.push_back(wl_index.first);
    }
}


EXRBiSpectralImage::CHANNEL_TYPE EXRBiSpectralImage::channelType(
    const std::string& channelName,
    int& stokesComponent,
    float& wavelength_nm,
    float& reradiation_wavelength_nm
) {
    const std::string exprStokes = "S([0-3])";
    const std::string exprWave   = "(\\d*,?\\d*([Ee][+-]?\\d+)?)";
    const std::string exprUnits  = "(Y|Z|E|P|T|G|M|k|h|da|d|c|m|u|n|p)?(m|Hz)";

    const std::regex exprDiagonal("^" + exprStokes + "\\." + exprWave + exprUnits + "$");
    const std::regex exprRerad   ("^" + exprStokes + "\\." + exprWave + exprUnits + "\\." + exprWave + exprUnits + "$");

    std::smatch matches;

    const bool matchedDiagonal = std::regex_search(channelName, matches, exprDiagonal);

    if (matchedDiagonal) {
        if (matches.size() != 6) {
            // Something went wrong with the parsing. This shall not occur.
            throw INTERNAL_ERROR;
        }
        
        stokesComponent = std::stoi(matches[1].str());
        wavelength_nm = toWavelength_nm(
            matches[2].str(), // Comma separated floating point value
            matches[4].str(), // Unit multiplier
            matches[5].str()  // Units
            );

        return SPECTRAL_DIAGONAL;	
    }

    const bool matchedRerad = std::regex_search(channelName, matches, exprRerad);

    if (matchedRerad) {
        if (matches.size() != 10) {
            // Something went wrong with the parsing. This shall not occur.
            throw INTERNAL_ERROR;
        }

        stokesComponent = std::stoi(matches[1].str());
        if (stokesComponent != 0) return OTHER;

        wavelength_nm = toWavelength_nm(
            matches[2].str(), // Comma separated floating point value
            matches[4].str(), // Unit multiplier
            matches[5].str()  // Units
            );

        reradiation_wavelength_nm = toWavelength_nm(
            matches[6].str(), // Comma separated floating point value
            matches[8].str(), // Unit multiplier
            matches[9].str()  // Units
            );

        return SPECTRAL_RERADIATION;
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

        
std::string EXRBiSpectralImage::getChannelName(
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
    float reradiation_nm;

    if (channelType(
        channelName, 
        stokesComponentChecked, 
        wavelength_nmChecked, 
        reradiation_nm) != SPECTRAL_DIAGONAL) {
        throw INTERNAL_ERROR;
    }

    if (stokesComponentChecked != stokesComponent 
     || wavelength_nmChecked != wavelength_nmChecked) {
        throw INTERNAL_ERROR;
    }

    return channelName;
}


std::string EXRBiSpectralImage::getChannelName(
    float wavelength_nm,
    float reradiation_wavelength_nm
) {
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
        reradiation_wavelength_nmChecked) != SPECTRAL_RERADIATION) {
        throw INTERNAL_ERROR;
    }

    if (stokesComponentChecked != 0 
     || wavelength_nmChecked != wavelength_nmChecked
     || reradiation_wavelength_nmChecked != reradiation_wavelength_nm) {
        throw INTERNAL_ERROR;
    }

    return channelName;
}