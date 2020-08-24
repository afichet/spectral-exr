#include <SpectrumAttribute.h>

#include "Util.h"

#include <string>
#include <regex>


SpectrumAttribute::SpectrumAttribute() {}


SpectrumAttribute::SpectrumAttribute(
    std::vector<float> wavelengths_nm,
    std::vector<float> values
)
    : _wavelengths_nm(wavelengths_nm)
    , _values(values)
{
    if (_wavelengths_nm.size() != _values.size()) {
        throw NOT_SAME_VECTOR_SIZE;
    }
}


SpectrumAttribute::SpectrumAttribute(
    const Imf::StringAttribute& attributeValue
) {
    std::string attributeValueStr = attributeValue.value();
    
    const std::string floatRegex  = "(\\d*\\.?\\d*([Ee][+-]?\\d+)?)";
    const std::string units       = "(Y|Z|E|P|T|G|M|k|h|da|d|c|m|u|n|p)?(m|Hz)";

    const std::regex e(floatRegex + units + ":" + floatRegex + ";");


    std::vector<std::pair<float, float>> wavelengthValues;

    for (std::smatch matches; 
         std::regex_search(attributeValueStr, matches, e); 
         attributeValueStr = matches.suffix())
    {
        if (matches.size() != 7) {
            throw PARSING_ERROR;
        }

        const std::string& waveValue = matches[1];
        const std::string& prefix    = matches[3];
        const std::string& units     = matches[4];
        const std::string& value     = matches[5];

        wavelengthValues.push_back(
            std::make_pair(
                Util::strToNanometers(std::stof(waveValue), prefix, units),
                std::stof(value)
            )
        );
    }

    // Sort ascending values
    std::sort(wavelengthValues.begin(), wavelengthValues.end());

    // Populate class data
    _wavelengths_nm.reserve(wavelengthValues.size());
    _values.reserve(wavelengthValues.size());

    for (const auto& wavelengthValue: wavelengthValues) {
        _wavelengths_nm.push_back(wavelengthValue.first);
        _values.push_back(wavelengthValue.second);
    }
}


Imf::StringAttribute SpectrumAttribute::getAttribute() 
const {
    std::stringstream attrValue;

    for (size_t i = 0; i < _wavelengths_nm.size(); i++) {
        attrValue << _wavelengths_nm[i] << "nm:" << _values[i] << ";";
    }

    return Imf::StringAttribute(attrValue.str());
}