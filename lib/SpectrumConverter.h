#pragma once

#include <array>
#include <vector>

class SpectrumConverter {
    public:
        SpectrumConverter();

        SpectrumConverter(
            const float& cmfFirstWavelength_nm,
            const std::array<std::vector<float>, 3> &xyzCmfs,
            const std::array<float, 9> xyzToRgb
        );

        void spectrumToXyz(
            const std::vector<float>& wavelengths_nm,
            const float* spectrum,
            std::array<float, 3>& XYZ
        ) const;

        void spectrumToRgb(
            const std::vector<float>& wavelengths_nm,
            const float* spectrum,
            std::array<float, 3>& RGB
        ) const;

        float firstWavelength() const { return _cmfFirstWavelength_nm; }
        float lastWavelength() const { return _cmfFirstWavelength_nm + _xyzCmfs[0].size() - 1; }

        size_t cmfWavelengthIndex(float wavelength_nm) const;
        size_t cmfWavelengthValue(size_t index) const;

    protected:
        float _cmfFirstWavelength_nm;
        std::array<std::vector<float>, 3> _xyzCmfs;
        std::array<float, 9> _xyzToRgb;

};
