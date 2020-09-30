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

#include "BiSpectralImage.h"

namespace SEXR {

class EXRBiSpectralImage: public BiSpectralImage {
    public:
        EXRBiSpectralImage(
            size_t width = 0, size_t height = 0,
            const std::vector<float>& wavelengths_nm = std::vector<float>(),
            SpectrumType type = REFLECTIVE
        );

        EXRBiSpectralImage(
            const std::string& filename
        );

        void save(const std::string& filename) const;

        static SpectrumType channelType(
            const std::string& channelName, 
            int& polarisationComponent,
            float& wavelengths_nm,
            float& reradiation_wavelength_nm
        );

        static std::string getEmissiveChannelName(
            int stokesComponent,
            float wavelength_nm
        );

        static std::string getReflectiveChannelName(
            float wavelength_nm
        );

        static std::string getReradiationChannelName(
            float wavelength_nm,
            float reradiation_wavelength_nm
        );

    static constexpr const char* SPECTRUM_TYPE_ATTR         = "Spectrum type"; 
    static constexpr const char* LENS_TRANSMISSION_ATTR     = "Lens transmission"; 
    static constexpr const char* CAMERA_RESPONSE_ATTR       = "Camera response"; 
    static constexpr const char* EXPOSURE_COMPENSATION_ATTR = "EV";
};

} // namespace SEXR
