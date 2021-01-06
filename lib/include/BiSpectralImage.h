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

#include "SpectralImage.h"

namespace SEXR
{
  class BiSpectralImage: public SpectralImage
  {
  public:
    /**
     * Creates a new spectral or bispectral image.
     *
     * @param width width of the image.
     * @param height height of the image.
     * @param wavelengths_nm wavlengths in nanometers of the image.
     * @param type spectrum type represented in the image.
     * @param handedness polarisation handedness convention.
     */
    BiSpectralImage(
      size_t                    width          = 0,
      size_t                    height         = 0,
      const std::vector<float> &wavelengths_nm = std::vector<float>(),
      SpectrumType              type           = REFLECTIVE,
      PolarisationHandedness    handedness     = RIGHT_HANDED);

    /**
     * Export each channel value in an individual EXR image.
     * For each spectral and bispectral channel, an individual image
     * will be created in the folder pointed by path. The folder must
     * exists prior to the call of this method.
     *  - Reflective / Transmissive channels are prefixed with "T - "
     *  - Emissive channels are prefixed with "S - "
     * Then, the wavelength in nanometer is appended with "nm" unit
     * specified.
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

    /**
     * Number of elements needed to store the reradiation part of
     * the image.
     */
    size_t reradiationSize() const
    {
      return nSpectralBands() * (nSpectralBands() - 1) / 2;
    }

    /**
     * Gives the index where the reradiation is stored from indices
     * of radiating wavelength and reemission wavelength.
     *
     * @param wlFrom_idx index of the radiating wavelength.
     * @param wlTo_idx index of the reemissive wavelength.
     *
     * @returns index where the reradiation is stored.
     */
    static size_t idxFromWavelengthIdx(size_t wlFrom_idx, size_t wlTo_idx);

    /**
     * Gives the radiating and reemissive indices from the index where
     * the reradiation is stored.
     *
     * @param rerad_idx index where the reradiation is stored.
     * @param wlFrom_idx index of the radiating wavelength.
     * @param wlTo_idx index of the reemissive wavelength.
     */
    static void wavelengthsIdxFromIdx(
      size_t rerad_idx, size_t &wlFrom_idx, size_t &wlTo_idx);

    /**
     * Return the reflective value for a given pixel for a given
     * radiating wavelength to a specific reemitting wavelength.
     * If wavelengthFrom_idx == wavelengthTo_idx, this returns the
     * diagonal (non fluorescent part).
     *
     * @param x column coordinate in the image in pixels (0 on left).
     * @param y row coordinate in the image in pixels (0 on top).
     * @param wavelengthFrom_idx index of the radiating wavelength.
     * @param wavelengthTo_idx index of the reemitting wavelength.
     */
    virtual float getReflectiveValue(
      size_t x,
      size_t y,
      size_t wavelengthFrom_idx,
      size_t wavelengthTo_idx) const;

    /**
     * Gives a reference to the reflective element at location x, y
     * for given radiating and reemissive wavelengths indices.
     *
     * @param x column coordinate in the image in pixels (0 on left).
     * @param y row coordinate in the image in pixels (0 on top).
     * @param wavelengthFrom_idx index of the radiating wavelength.
     * @param wavelengthTo_idx index of the reemitting wavelength.
     */
    virtual float &reflective(
      size_t x, size_t y, size_t wavelengthFrom_idx, size_t wavelengthTo_idx);

    virtual const float &reflective(
      size_t x,
      size_t y,
      size_t wavelengthFrom_idx,
      size_t wavelengthTo_idx) const;


  protected:
    // Upper right triangular matrices for each pixel
    // pixel stride is reradiationSize()
    std::vector<float> _reradiation;
  };

}   // namespace SEXR
