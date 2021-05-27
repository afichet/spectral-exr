/**
 * Copyright (c) 2020 Alban Fichet, Romain Pacanowski, Alexander Wilkie
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

namespace SEXR
{
    class EXRBiSpectralImage: public BiSpectralImage
    {
      public:
        /**
         * Creates a new spectral or bispectral image.
         *
         * @param width width of the image.
         * @param height height of the image.
         * @param wavelengths_nm wavlengths in nanometers of the image.
         * @param type spectrum type represented in the image.
         */
        EXRBiSpectralImage(
          size_t                    width          = 0,
          size_t                    height         = 0,
          const std::vector<float> &wavelengths_nm = std::vector<float>(),
          SpectrumType              type           = REFLECTIVE,
          PolarisationHandedness    handedness     = RIGHT_HANDED);

        /**
         * Loads a spectral or bispectral image from an EXR file.
         *
         * @param filename path to the image to load.
         */
        EXRBiSpectralImage(const std::string &filename);

        /**
         * Saves the bispectral image to an EXR file.
         *
         * @param filename path where the image shall be saved.
         */
        void save(const std::string &filename) const;

        static SpectrumType channelType(
          const std::string &channelName,
          int &              polarisationComponent,
          double &           wavelengths_nm,
          double &           reradiation_wavelength_nm);

        /**
         * Gets the channel name used in the EXR file for a specific
         * emissive component.
         *
         * @param stokesComponent index of the Stokes component (0-3).
         * @param wavelength_nm wavelength in nanometers.
         *
         * @returns std::string containing the emissive channel name
         * for the given Stokes component at a specific wavelength.
         */
        static std::string
        getEmissiveChannelName(int stokesComponent, double wavelength_nm);

        /**
         * Gets the channel name used in the EXR file for a specific
         * non reradiating component.
         *
         * @param wavelength_nm wavelength in nanometers.
         *
         * @returns std::string containing the reflective channel name
         * at a specific wavelength.
         */
        static std::string getReflectiveChannelName(double wavelength_nm);

        /**
         * Gets the channel name used in the EXR file for a specific
         * reradiating component.
         *
         * @param wavelength_nm radiating wavelength in nanometers.
         * @param reradiation_wavelength_nm reemissive wavelength in
         * nanometers.
         *
         * @returns std::string containing the reflective channel name
         * at a specific radiating wavelength to a specific reemission
         * band.
         */
        static std::string getReradiationChannelName(
          double wavelength_nm, double reradiation_wavelength_nm);

        static constexpr const char *VERSION_ATTR = "spectralLayoutVersion";
        static constexpr const char *SPECTRUM_TYPE_ATTR  = "spectrumType";
        static constexpr const char *EMISSIVE_UNITS_ATTR = "emissiveUnits";
        static constexpr const char *LENS_TRANSMISSION_ATTR
          = "lensTransmission";
        static constexpr const char *CAMERA_RESPONSE_ATTR = "cameraResponse";
        static constexpr const char *EXPOSURE_COMPENSATION_ATTR = "EV";
        static constexpr const char *POLARISATION_HANDEDNESS_ATTR
          = "polarisationHandedness";
    };

}   // namespace SEXR
