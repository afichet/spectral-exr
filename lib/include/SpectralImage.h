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

#include <vector>
#include <array>
#include <string>

#include <SpectrumAttribute.h>
#include <SpectrumType.h>

namespace SEXR {

class SpectralImage {
    public:
        enum Errors {
            UNSUPORTED_FILE,
            INTERNAL_ERROR,
            READ_ERROR,
            WRITE_ERROR,
            INCORRECT_FORMED_FILE
        };

        SpectralImage(
            size_t width = 0, size_t height = 0,
            const std::vector<float>& wavelengths_nm = std::vector<float>(),
            SpectrumType type = EMISSIVE
        );

        virtual void save(const std::string& filename) const = 0;

        virtual void exportChannels(const std::string& path) const;
        virtual void getRGBImage(std::vector<float>& rgbImage) const;

        void setCameraResponse(
            const std::vector<float>& wavelengths_nm,
            const std::vector<float>& values
        );

        const SpectrumAttribute& cameraResponse()   const;


        void setLensTransmission(
            const std::vector<float>& wavelengths_nm,
            const std::vector<float>& values
        );


        const SpectrumAttribute& lensTransmission() const;

        void setChannelSensitivity(
            size_t wl_idx,
            const std::vector<float>& wavelengths_nm,
            const std::vector<float>& values
        );

        
        const std::vector<SpectrumAttribute>& channelSensitivities() const;

        const SpectrumAttribute& channelSensitivity(size_t wl_idx) const;

        
        void setExposureCompensationValue(float ev);

        const float& exposureCompensationValue() const;

        // Access the emissive part
        virtual float& operator()(
            size_t x, size_t y, 
            size_t wavelength_idx, 
            size_t stokesComponent);

        virtual const float& operator()(
            size_t x, size_t y, 
            size_t wavelength_idx, 
            size_t stokesComponent) const;

        // Access the reflective part
        virtual float& operator()(
            size_t x, size_t y, 
            size_t wavelength_idx);

        virtual const float& operator()(
            size_t x, size_t y, 
            size_t wavelength_idx) const;

        // Those are not direct memory access
        // They can be called whatever the image type is

        virtual float getEmissiveValue(
            size_t x, size_t y, 
            size_t wavelength_idx, 
            size_t stokesComponent = 0) const;

        virtual float getReflectiveValue(
            size_t x, size_t y, 
            size_t wavelength_idx) const;

        const float& wavelength_nm(size_t wl_idx) const;

        size_t width()  const;
        size_t height() const;

        size_t nSpectralBands()     const;
        size_t nStokesComponents()  const;
        
        bool polarised()    const;
        bool emissive()     const;
        bool reflective()   const;
        bool bispectral()   const;
        SpectrumType type() const;

    protected:
        size_t _width, _height;
        float _ev;
                
        // We can have up to 6 pixel buffers:
        // - 1 for emissive unpolarised images (S0)
        // - 1 for reflective unpolarised images (RE)
        // - 4 for emissive polarised images (S0, S1, S2, S3)
        
        std::vector<float> _reflectivePixelBuffer;
        std::array<std::vector<float>, 4> _emissivePixelBuffers;

        std::vector<float> _wavelengths_nm;
        SpectrumType _spectrumType;

        SpectrumAttribute _lensTransmissionSpectra;
        SpectrumAttribute _cameraReponse;
        std::vector<SpectrumAttribute> _channelSensitivities;
};

} // namespace SEXR
