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
 *
 *  * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer in the documentation and/or other materials provided
 * with the distribution.
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

#include <BiSpectralImage.h>

#include <sstream>
#include <cassert>
#include <cmath>

#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfFrameBuffer.h>

#include "SpectrumConverter.h"

namespace SEXR
{
    BiSpectralImage::BiSpectralImage(
      size_t                    width,
      size_t                    height,
      const std::vector<float> &wavelengths_nm,
      SpectrumType              type,
      PolarisationHandedness    handedness)
      : SpectralImage(width, height, wavelengths_nm, type, handedness)
    {
        if (isBispectral()) {
            _reradiation.resize(reradiationSize() * _width * _height);
        }
    }


    void BiSpectralImage::exportChannels(const std::string &path) const
    {
        // Export the diagonal
        SpectralImage::exportChannels(path);

        if (isBispectral()) {
            const size_t xStride = sizeof(float) * reradiationSize();
            const size_t yStride = xStride * width();

            // Export the reradiation
            for (size_t rr = 0; rr < reradiationSize(); rr++) {
                for (size_t wl_i_idx = 0; wl_i_idx < nSpectralBands();
                     wl_i_idx++) {
                    const float &wavelength_i = _wavelengths_nm[wl_i_idx];

                    for (size_t wl_o_idx = wl_i_idx + 1;
                         wl_o_idx < nSpectralBands();
                         wl_o_idx++) {
                        const float &wavelength_o = _wavelengths_nm[wl_o_idx];

                        std::stringstream filepath;
                        filepath << path << "/"
                                 << "T - " << wavelength_i << "nm - "
                                 << wavelength_o << "nm.exr";

                        Imf::Header       exrHeader(width(), height());
                        Imf::ChannelList &exrChannels = exrHeader.channels();
                        Imf::FrameBuffer  exrFrameBuffer;

                        exrChannels.insert("Y", Imf::Channel(Imf::FLOAT));
                        exrFrameBuffer.insert(
                          "Y",
                          Imf::Slice(
                            Imf::FLOAT,
                            (char *)(&_reradiation[rr]),
                            xStride,
                            yStride));

                        Imf::OutputFile exrOut(
                          filepath.str().c_str(),
                          exrHeader);
                        exrOut.setFrameBuffer(exrFrameBuffer);
                        exrOut.writePixels(height());
                    }
                }
            }
        }
    }


    void BiSpectralImage::getRGBImage(std::vector<float> &rgbImage) const
    {
        if (!isBispectral()) {
            SpectralImage::getRGBImage(rgbImage);
        } else {
            rgbImage.resize(3 * width() * height());
            SpectrumConverter sc(isEmissive());

            std::array<float, 3> rgb;

            if (isEmissive() && isReflective()) {
                for (size_t i = 0; i < width() * height(); i++) {
                    sc.spectraToRGB(
                      _wavelengths_nm,
                      &_reflectivePixelBuffer[nSpectralBands() * i],
                      &_reradiation[reradiationSize() * i],
                      &_emissivePixelBuffers[0][nSpectralBands() * i],
                      rgb);

                    memcpy(&rgbImage[3 * i], &rgb[0], 3 * sizeof(float));
                }

                // Exposure compensation
                for (size_t i = 0; i < width() * height(); i++) {
                    for (size_t c = 0; c < 3; c++) {
                        rgbImage[3 * i + c] *= std::pow(2.F, _ev);
                    }
                }
            } else if (isReflective()) {
                for (size_t i = 0; i < width() * height(); i++) {
                    sc.spectrumToRGB(
                      _wavelengths_nm,
                      &_reflectivePixelBuffer[nSpectralBands() * i],
                      &_reradiation[reradiationSize() * i],
                      rgb);

                    memcpy(&rgbImage[3 * i], &rgb[0], 3 * sizeof(float));
                }

                // Exposure compensation
                for (size_t i = 0; i < width() * height(); i++) {
                    for (size_t c = 0; c < 3; c++) {
                        rgbImage[3 * i + c] *= std::pow(2.F, _ev);
                    }
                }
            }
        }
    }


    size_t
    BiSpectralImage::idxFromWavelengthIdx(size_t wlFrom_idx, size_t wlTo_idx)
    {
        if (wlFrom_idx < wlTo_idx) {
            return wlTo_idx * (wlTo_idx - 1) / 2 + wlFrom_idx;
        } else {
            return -1;
        }
    }


    void BiSpectralImage::wavelengthsIdxFromIdx(
      size_t rerad_idx, size_t &wlFrom_idx, size_t &wlTo_idx)
    {
        float k
          = std::floor((std::sqrt(1.F + 8.F * float(rerad_idx)) - 1.F) / 2.F);
        float j = rerad_idx - k * (k + 1) / 2.F;

        wlFrom_idx = j;
        wlTo_idx   = k + 1;
    }


    float BiSpectralImage::getReflectiveValue(
      size_t x,
      size_t y,
      size_t wavelengthFrom_idx,
      size_t wavelengthTo_idx) const
    {
        assert(x < width());
        assert(y < height());
        assert(wavelengthFrom_idx < nSpectralBands());
        assert(wavelengthTo_idx < nSpectralBands());

        if (!isReflective() || wavelengthFrom_idx > wavelengthTo_idx) {
            return 0.F;
        }

        if (!isBispectral()) {
            return 0.F;
        }

        return reflective(x, y, wavelengthFrom_idx, wavelengthTo_idx);
    }


    float &BiSpectralImage::reflective(
      size_t x, size_t y, size_t wavelengthFrom_idx, size_t wavelengthTo_idx)
    {
        assert(isReflective());
        assert(x < width());
        assert(y < height());
        assert(wavelengthFrom_idx < nSpectralBands());
        assert(wavelengthTo_idx < nSpectralBands());
        assert(wavelengthFrom_idx <= wavelengthTo_idx);

        if (wavelengthFrom_idx == wavelengthTo_idx) {
            return SpectralImage::reflective(x, y, wavelengthFrom_idx);
        }

        assert(isBispectral());
        assert(_reradiation.size() == reradiationSize() * width() * height());

        size_t reradIdx
          = idxFromWavelengthIdx(wavelengthFrom_idx, wavelengthTo_idx);

        return _reradiation[reradiationSize() * (y * width() + x) + reradIdx];
    }


    const float &BiSpectralImage::reflective(
      size_t x,
      size_t y,
      size_t wavelengthFrom_idx,
      size_t wavelengthTo_idx) const
    {
        assert(isReflective());
        assert(x < width());
        assert(y < height());
        assert(wavelengthFrom_idx < nSpectralBands());
        assert(wavelengthTo_idx < nSpectralBands());

        if (wavelengthFrom_idx == wavelengthTo_idx) {
            return SpectralImage::reflective(x, y, wavelengthFrom_idx);
        }

        assert(isBispectral());
        assert(_reradiation.size() == reradiationSize() * width() * height());

        size_t reradIdx
          = idxFromWavelengthIdx(wavelengthFrom_idx, wavelengthTo_idx);

        return _reradiation[reradiationSize() * (y * width() + x) + reradIdx];
    }

}   // namespace SEXR
