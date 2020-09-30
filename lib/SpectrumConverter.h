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

#pragma once

#include <array>
#include <vector>

namespace SEXR {

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


        // Bi-spectral
        void spectraToXYZ(
            const std::vector<float>& wavelengths_nm,
            const float* diagonal,
            const float* reradiation,
            const float* emissiveSpectrum,
            std::array<float, 3>& XYZ
        ) const;

        void spectraToRGB(
            const std::vector<float>& wavelengths_nm,
            const float* diagonal,
            const float* reradiation,
            const float* emissiveSpectrum,
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

} // namespace SEXR
