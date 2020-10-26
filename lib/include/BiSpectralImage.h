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
    BiSpectralImage(
      size_t                    width          = 0,
      size_t                    height         = 0,
      const std::vector<float> &wavelengths_nm = std::vector<float>(),
      SpectrumType              type           = REFLECTIVE);

    virtual void exportChannels(const std::string &path) const;
    virtual void getRGBImage(std::vector<float> &rgbImage) const;

    size_t reradiationSize() const
    {
      return nSpectralBands() * (nSpectralBands() - 1) / 2;
    }

    static size_t idxFromWavelengthIdx(size_t wlFrom_idx, size_t wlTo_idx);

    static void wavelengthsIdxFromIdx(
      size_t rerad_idx, size_t &wlFrom_idx, size_t &wlTo_idx);

    virtual float getReflectiveValue(
      size_t x,
      size_t y,
      size_t wavelengthFrom_idx,
      size_t wavelengthTo_idx) const;

    virtual float &operator()(
      size_t x, size_t y, size_t wavelengthFrom_idx, size_t wavelengthTo_idx);

    virtual const float &operator()(
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
