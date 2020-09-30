/**
 * Copyright (c) 2020 
 *  Alban Fichet <alban.fichet@gmx.fr>, 
 *  Romain Pacanowski <romain.pacanowski@inria.fr>, 
 *  Alexander Wilkie <alexander.wilkie@cgg.mff.cuni.cz>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *  * Neither the name of %ORGANIZATION% nor the names of its contributors may be
 * used to endorse or promote products derived from this software without specific
 * prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <SpectrumAttribute.h>

#include "Util.h"

#include <string>
#include <regex>

namespace SEXR {

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

} // namespace SEXR
