#pragma once

#include <vector>
#include <OpenEXR/ImfStandardAttributes.h>

namespace SEXR {

class SpectrumAttribute {
    public:
        enum Error {
            NOT_SAME_VECTOR_SIZE,
            PARSING_ERROR
        };

        SpectrumAttribute();

        SpectrumAttribute(
            std::vector<float> wavelengths_nm,
            std::vector<float> values
        );

        SpectrumAttribute(
            const Imf::StringAttribute& attributeValue
        );

        Imf::StringAttribute getAttribute() const;

              std::vector<float>& wavelengths_nm()       { return _wavelengths_nm; }
        const std::vector<float>& wavelengths_nm() const { return _wavelengths_nm; }

              std::vector<float>& values()       { return _values; }
        const std::vector<float>& values() const { return _values; }

              float& wavelength_nm(size_t i)       { return _wavelengths_nm[i]; }
        const float& wavelength_nm(size_t i) const { return _wavelengths_nm[i]; }

              float& value(size_t i)       { return _values[i]; }
        const float& value(size_t i) const { return _values[i]; }

        size_t size() const { return _wavelengths_nm.size(); }

    protected:
        std::vector<float> _wavelengths_nm;
        std::vector<float> _values;
};

} // namespace SEXR
