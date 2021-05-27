/**
 * Copyright (c) 2020 - 2021
 * Alban Fichet, Romain Pacanowski, Alexander Wilkie
 * Institut d'Optique Graduate School, CNRS - Universite de Bordeaux,
 * Inria, Charles University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer in the documentation and/or other materials provided
 * with the distribution.
 *
 *  * Neither the name of Institut d'Optique Graduate School, CNRS -
 * Universite de Bordeaux, Inria, Charles University nor the names of
 * its contributors may be used to endorse or promote products derived
 * from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */


#pragma once

#include <vector>
#include <array>
#include <string>

#include <SpectrumAttribute.h>
#include <SpectrumType.h>

namespace SEXR
{
    class SpectralImage
    {
      public:
        enum Errors
        {
            UNSUPORTED_FILE,
            INTERNAL_ERROR,
            READ_ERROR,
            WRITE_ERROR,
            INCORRECT_FORMED_FILE
        };

        enum PolarisationHandedness
        {
            LEFT_HANDED,
            RIGHT_HANDED
        };

        /**
         * Creates a new spectral image.
         *
         * @param width width of the image.
         * @param height height of the image.
         * @param wavelengths_nm wavlengths in nanometers of the image.
         * @param type spectrum type represented in the image.
         * @param handedness polarisation handedness convention.
         */
        SpectralImage(
          size_t                    width          = 0,
          size_t                    height         = 0,
          const std::vector<float> &wavelengths_nm = std::vector<float>(),
          SpectrumType              type           = EMISSIVE,
          PolarisationHandedness    handedness     = RIGHT_HANDED);

        /**
         * Saves the image to an EXR file.
         *
         * @param filename path where the image shall be saved.
         */
        virtual void save(const std::string &filename) const = 0;

        /**
         * Export each channel value in an individual EXR image.  For
         * each spectral and bispectral channel, an individual image
         * will be created in the folder pointed by path. The folder
         * must exists prior to the call of this method.
         *
         *  - Reflective / Transmissive channels are prefixed with "T - "
         *  - Emissive channels are prefixed with "S - "
         *
         * Then, the wavelength in nanometer is appended with "nm"
         * unit specified.
         *
         * @param path folder path where to export the images.
         */
        virtual void exportChannels(const std::string &path) const;

        /**
         * Get the RGB version of the image. The rgb image is stored
         * as rgbImage[3 * (row * width + col) + color].
         *
         * @param rgbImage a reference to the pointer where to store
         * the RGB version of this bispectral image.
         */
        virtual void getRGBImage(std::vector<float> &rgbImage) const;

        void setCameraResponse(
          const std::vector<float> &wavelengths_nm,
          const std::vector<float> &values);

        const SpectrumAttribute &cameraResponse() const;

        /**
         * Sets the lens transmissivity curve.
         *
         * @param wavelengths_nm wavelengths of the transmissive spectrum.
         * @param values transmissives values of the lens.
         */
        void setLensTransmission(
          const std::vector<float> &wavelengths_nm,
          const std::vector<float> &values);

        /** Gets the lens transmissivitiy function. */
        const SpectrumAttribute &lensTransmission() const;

        /**
         * Sets a "filter" transmissivity curve for a specific
         * wavelength.
         *
         * @param wl_idx wavelenght index which the filter corresponds to.
         * @param wavelengths_nm wavelengths of the transmissive spectrum.
         * @param values transmissives values of the filter.
         */
        void setChannelSensitivity(
          size_t                    wl_idx,
          const std::vector<float> &wavelengths_nm,
          const std::vector<float> &values);

        /** Gets the transmissivity curves for each wavelength */
        const std::vector<SpectrumAttribute> &channelSensitivities() const;

        /**
         * Gets the transmissivity curve for a specific wavelength
         *
         * @param wl_idx index of the wavlength to get the sensitivity curve from.
         */
        const SpectrumAttribute &channelSensitivity(size_t wl_idx) const;

        /**
         * Sets the exposure compensation values to use for the RGB
         * representation.
         */
        void setExposureCompensationValue(float ev);

        /**
         * Gets the exposure compensation values used by the RGB
         * representation.
         */
        const float &exposureCompensationValue() const;


        // Access the emissive part
        virtual float &emissive(
          size_t x, size_t y, size_t wavelength_idx, size_t stokesComponent);

        virtual const float &emissive(
          size_t x,
          size_t y,
          size_t wavelength_idx,
          size_t stokesComponent) const;

        // Access the reflective part

        /**
         * Gives a reference to the reflective element at location x,
         * y for given a wavelength index.
         *
         * @param x column coordinate in the image in pixels (0 on left).
         * @param y row coordinate in the image in pixels (0 on top).
         * @param wavelength_idx index of the radiating wavelength.
         */
        virtual float &reflective(size_t x, size_t y, size_t wavelength_idx);

        virtual const float &
        reflective(size_t x, size_t y, size_t wavelength_idx) const;

        // Those are not direct memory access
        // They can be called whatever the image type is

        /**
         * Return the emissive value for a given pixel for a given
         * wavelength.
         *
         * @param x column coordinate in the image in pixels (0 on left).
         * @param y row coordinate in the image in pixels (0 on top).
         * @param wavelength_idx index of the wavelength.
         * @param stokesComponent index of the Stokes component.
         */
        virtual float getEmissiveValue(
          size_t x,
          size_t y,
          size_t wavelength_idx,
          size_t stokesComponent = 0) const;

        /**
         * Return the reflective value for a given pixel for a given
         * wavelength.
         *
         * @param x column coordinate in the image in pixels (0 on left).
         * @param y row coordinate in the image in pixels (0 on top).
         * @param wavelength_idx index of the wavelength.
         */
        virtual float
        getReflectiveValue(size_t x, size_t y, size_t wavelength_idx) const;

        /**
         * Gets the wavelength value in nanometer from a specific
         * wavelength index.
         */
        const float &wavelength_nm(size_t wl_idx) const;

        /** Gets the width of the image in pixels. */
        size_t width() const;

        /** Gets the height of the image in pixels. */
        size_t height() const;

        /** Gets the number of spectral bands. */
        size_t nSpectralBands() const;

        /**
         * Gets the number of Stokes components. 1 for a non polarised
         * image, 4 for a polarised image.
         */
        size_t nStokesComponents() const;

        /** True if the image is polarised, False otherwise */
        bool isPolarised() const;

        /** True if the image contains emissive data, False otherwise */
        bool isEmissive() const;

        /** True if the image contains reflective data, False otherwise */
        bool isReflective() const;

        /** True if the image contains bispectral data, False otherwise */
        bool isBispectral() const;

        /** Spectrum type contains at each pixel location in the image */
        SpectrumType type() const;

      protected:
        size_t _width, _height;
        float  _ev;

        // We can have up to 6 pixel buffers:
        // - 1 for emissive unpolarised images (S0)
        // - 1 for reflective unpolarised images (RE)
        // - 4 for emissive polarised images (S0, S1, S2, S3)

        std::vector<float>                _reflectivePixelBuffer;
        std::array<std::vector<float>, 4> _emissivePixelBuffers;

        std::vector<float> _wavelengths_nm;
        SpectrumType       _spectrumType;

        PolarisationHandedness _polarisationHandedness;

        SpectrumAttribute              _lensTransmissionSpectra;
        SpectrumAttribute              _cameraReponse;
        std::vector<SpectrumAttribute> _channelSensitivities;
    };

}   // namespace SEXR
