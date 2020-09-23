#pragma once

#include <array>
#include <vector>

class SpectrumConverter {
    public:
        SpectrumConverter(bool emissiveSpectrum = true);

        SpectrumConverter(
            const float& cmfFirstWavelength_nm,
            const std::array<std::vector<float>, 3> &xyzCmfs,
            const std::array<float, 9> xyzToRgb
        );

        float firstWavelength() const { return _cmfFirstWavelength_nm; }
        float lastWavelength() const { return _cmfFirstWavelength_nm + _xyzCmfs[0].size() - 1; }

        size_t cmfWavelengthIndex(float wavelength_nm) const;
        size_t cmfWavelengthValue(size_t index) const;

        // The spectrum provided must be either emissive or reflective.
        // This depends on the constructor used.
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

        // Here, we provide two spectra, one for the reflective part,
        // the other for the emissive part.
        // The result will be the lighting multiplied by the reflective
        // added to the emissive spectrum.
        void spectraToXYZ(
            const std::vector<float>& wavelengths_nm,
            const float* reflectiveSpectrum,
            const float* emissiveSpectrum,
            std::array<float, 3>& XYZ
        ) const;

        void spectraToRGB(
            const std::vector<float>& wavelengths_nm,
            const float* reflectiveSpectrum,
            const float* emissiveSpectrum,
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