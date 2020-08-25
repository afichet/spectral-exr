#pragma once

#include <array>
#include <vector>

class SpectrumConverter {
    public:
        enum SPECTRUM_TYPE {
            REFLECTIVE,
            EMISSIVE
        };

        SpectrumConverter(SPECTRUM_TYPE type = EMISSIVE);

        SpectrumConverter(
            const float& cmfFirstWavelength_nm,
            const std::array<std::vector<float>, 3> &xyzCmfs,
            const std::array<float, 9> xyzToRgb
        );

        float firstWavelength() const { return _cmfFirstWavelength_nm; }
        float lastWavelength() const { return _cmfFirstWavelength_nm + _xyzCmfs[0].size() - 1; }

        size_t cmfWavelengthIndex(float wavelength_nm) const;
        size_t cmfWavelengthValue(size_t index) const;

        void spectrumToXYZ(
            const std::vector<float>& wavelengths_nm,
            const float* spectrum,
            std::array<float, 3>& XYZ
        ) const;

        void spectrumToRGB(
            const std::vector<float>& wavelengths_nm,
            const float* spectrum,
            std::array<float, 3>& RGB
        ) const;

        // Bi-spectral
        void spectrumToXYZ(
            const std::vector<float>& wavelengths_nm,
            const float* diagonal,
            const float* reradiation,
            std::array<float, 3>& XYZ
        ) const;

        void spectrumToRGB(
            const std::vector<float>& wavelengths_nm,
            const float* diagonal,
            const float* reradiation,
            std::array<float, 3>& RGB
        ) const;

    protected:
        void emissiveSpectrumToXYZ(
            const std::vector<float>& wavelengths_nm,
            const float* spectrum,
            std::array<float, 3>& XYZ
        ) const;
        
        void reflectiveSpectrumToXYZ(
            const std::vector<float>& wavelengths_nm,
            const float* spectrum,
            std::array<float, 3>& XYZ
        ) const;


        bool _emissiveSpectrum;

        float _illuminantFirstWavelenght_nm;
        std::vector<float> _illuminantSPD;

        float _cmfFirstWavelength_nm;
        std::array<std::vector<float>, 3> _xyzCmfs;
        std::array<float, 9> _xyzToRgb;
};