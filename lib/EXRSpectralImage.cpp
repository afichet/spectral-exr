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

#include <EXRSpectralImage.h>
#include "Util.h"

#include <regex>
#include <algorithm>
#include <sstream>
#include <cassert>

#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfFrameBuffer.h>

namespace SEXR
{
  EXRSpectralImage::EXRSpectralImage(
    size_t                    width,
    size_t                    height,
    const std::vector<float> &wavelengths_nm,
    SpectrumType              type,
    PolarisationHandedness    handedness)
    : SpectralImage(width, height, wavelengths_nm, type, handedness)
  {}


  EXRSpectralImage::EXRSpectralImage(const std::string &filename)
    : SpectralImage()
  {
    Imf::InputFile      exrIn(filename.c_str());
    const Imf::Header & exrHeader     = exrIn.header();
    const Imath::Box2i &exrDataWindow = exrHeader.dataWindow();

    _width        = exrDataWindow.max.x - exrDataWindow.min.x + 1;
    _height       = exrDataWindow.max.y - exrDataWindow.min.y + 1;
    _spectrumType = UNDEFINED;

    // -----------------------------------------------------------------------
    // Determine channels' position
    // -----------------------------------------------------------------------

    const Imf::ChannelList &exrChannels = exrHeader.channels();

    std::array<std::vector<std::pair<float, std::string>>, 4> wavelengths_nm_S;
    std::vector<std::pair<float, std::string>> wavelengths_nm_reflective;

    for (Imf::ChannelList::ConstIterator channel = exrChannels.begin();
         channel != exrChannels.end();
         channel++) {
      // Check if the channel is a spectral one
      int          polarisationComponent;
      double       wavelength_nm;
      SpectrumType spectralChannel
        = channelType(channel.name(), polarisationComponent, wavelength_nm);

      if (spectralChannel != SpectrumType::UNDEFINED) {
        _spectrumType = _spectrumType | spectralChannel;

        if (isEmissiveSpectrum(spectralChannel)) {
          wavelengths_nm_S[polarisationComponent].push_back(
            std::make_pair(wavelength_nm, channel.name()));
        } else if (isReflectiveSpectrum(spectralChannel)) {
          wavelengths_nm_reflective.push_back(
            std::make_pair(wavelength_nm, channel.name()));
        }
      }
    }

    // Sort by ascending wavelengths
    for (size_t s = 0; s < nStokesComponents(); s++) {
      std::sort(wavelengths_nm_S[s].begin(), wavelengths_nm_S[s].end());
    }

    if (isReflective()) {
      std::sort(
        wavelengths_nm_reflective.begin(),
        wavelengths_nm_reflective.end());
    }

    // -------------------------------------------------------------------------
    // Sanity check
    // -------------------------------------------------------------------------

    if (_spectrumType == SpectrumType::UNDEFINED) {
      // Probably an RGB EXR, not our job to handle it
      throw INCORRECT_FORMED_FILE;
    }

    if (isEmissive()) {
      // Check we have the same wavelength for each Stokes component
      // Wavelength vectors must be of the same size
      const float base_size_emissive = wavelengths_nm_S[0].size();

      for (size_t s = 1; s < nStokesComponents(); s++) {
        if (wavelengths_nm_S[s].size() != base_size_emissive) {
          throw INCORRECT_FORMED_FILE;
        }

        // Wavelengths must correspond
        for (size_t wl_idx = 0; wl_idx < base_size_emissive; wl_idx++) {
          if (
            wavelengths_nm_S[s][wl_idx].first
            != wavelengths_nm_S[0][wl_idx].first) {
            throw INCORRECT_FORMED_FILE;
          }
        }
      }
    }

    // If both reflective and emissive, we need to perform a last sanity check
    if (isEmissive() && isReflective()) {
      const size_t n_emissive_wavelengths   = wavelengths_nm_S[0].size();
      const size_t n_reflective_wavelengths = wavelengths_nm_reflective.size();

      if (n_emissive_wavelengths != n_reflective_wavelengths)
        throw INCORRECT_FORMED_FILE;

      for (size_t wl_idx = 0; wl_idx < n_emissive_wavelengths; wl_idx++) {
        if (wavelengths_nm_S[0][wl_idx] != wavelengths_nm_reflective[wl_idx])
          throw INCORRECT_FORMED_FILE;
      }
    }

    // -----------------------------------------------------------------------
    // Allocate memory
    // -----------------------------------------------------------------------

    // Now, we can populate the local wavelength vector
    if (isEmissive()) {
      _wavelengths_nm.reserve(wavelengths_nm_S[0].size());

      for (const auto &wl_index : wavelengths_nm_S[0]) {
        _wavelengths_nm.push_back(wl_index.first);
      }
    } else {
      _wavelengths_nm.reserve(wavelengths_nm_reflective.size());

      for (const auto &wl_index : wavelengths_nm_reflective) {
        _wavelengths_nm.push_back(wl_index.first);
      }
    }

    // We allocate pixel buffers memory
    for (size_t s = 0; s < nStokesComponents(); s++) {
      _emissivePixelBuffers[s].resize(nSpectralBands() * width() * height());
    }

    if (isReflective()) {
      _reflectivePixelBuffer.resize(nSpectralBands() * width() * height());
    }

    // -----------------------------------------------------------------------
    // Read the pixel data
    // -----------------------------------------------------------------------

    Imf::FrameBuffer exrFrameBuffer;

    const Imf::PixelType compType = Imf::FLOAT;
    const size_t         xStride  = sizeof(float) * nSpectralBands();
    const size_t         yStride  = xStride * _width;

    for (size_t s = 0; s < nStokesComponents(); s++) {
      for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
        char *ptrS = (char *)(&_emissivePixelBuffers[s][wl_idx]);
        exrFrameBuffer.insert(
          wavelengths_nm_S[s][wl_idx].second,
          Imf::Slice(compType, ptrS, xStride, yStride));
      }
    }

    if (isReflective()) {
      for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
        char *ptrS = (char *)(&_reflectivePixelBuffer[wl_idx]);
        exrFrameBuffer.insert(
          wavelengths_nm_reflective[wl_idx].second,
          Imf::Slice(compType, ptrS, xStride, yStride));
      }
    }

    exrIn.setFrameBuffer(exrFrameBuffer);
    exrIn.readPixels(exrDataWindow.min.y, exrDataWindow.max.y);

    // -----------------------------------------------------------------------
    // Read metadata
    // -----------------------------------------------------------------------

    // Lens transmission data
    const Imf::StringAttribute *lensTransmissionAttr
      = exrHeader.findTypedAttribute<Imf::StringAttribute>(
        LENS_TRANSMISSION_ATTR);

    if (lensTransmissionAttr != nullptr) {
      try {
        _lensTransmissionSpectra = SpectrumAttribute(*lensTransmissionAttr);
      } catch (SpectrumAttribute::Error &e) {
        throw INCORRECT_FORMED_FILE;
      }
    }

    // Camera spectral response
    const Imf::StringAttribute *cameraResponseAttr
      = exrHeader.findTypedAttribute<Imf::StringAttribute>(
        CAMERA_RESPONSE_ATTR);

    if (cameraResponseAttr != nullptr) {
      try {
        _cameraReponse = SpectrumAttribute(*cameraResponseAttr);
      } catch (SpectrumAttribute::Error &e) {
        throw INCORRECT_FORMED_FILE;
      }
    }

    // Each channel sensitivity
    _channelSensitivities.resize(nSpectralBands());

    for (size_t i = 0; i < wavelengths_nm_S[0].size(); i++) {
      const Imf::StringAttribute *filterTransmissionAttr
        = exrHeader.findTypedAttribute<Imf::StringAttribute>(
          wavelengths_nm_S[0][i].second);

      if (filterTransmissionAttr != nullptr) {
        try {
          _channelSensitivities[i] = SpectrumAttribute(*filterTransmissionAttr);
        } catch (SpectrumAttribute::Error &e) {
          throw INCORRECT_FORMED_FILE;
        }
      }
    }

    // Exposure compensation value
    const Imf::StringAttribute *exposureCompensationAttr
      = exrHeader.findTypedAttribute<Imf::StringAttribute>(
        EXPOSURE_COMPENSATION_ATTR);

    if (exposureCompensationAttr != nullptr) {
      try {
        _ev = std::stof(exposureCompensationAttr->value());
      } catch (std::invalid_argument &e) {
        throw INCORRECT_FORMED_FILE;
      }
    }

    // Polarisation handedness
    const Imf::StringAttribute *polarisationHandednessAttr
      = exrHeader.findTypedAttribute<Imf::StringAttribute>(
        POLARISATION_HANDEDNESS_ATTR);

    if (polarisationHandednessAttr != nullptr) {
      if (polarisationHandednessAttr->value() == "left") {
        _polarisationHandedness = LEFT_HANDED;
      } else if (polarisationHandednessAttr->value() == "right") {
        _polarisationHandedness = RIGHT_HANDED;
      } else {
        throw INCORRECT_FORMED_FILE;
      }
    }
  }


  void EXRSpectralImage::save(const std::string &filename) const
  {
    Imf::Header       exrHeader(width(), height());
    Imf::ChannelList &exrChannels = exrHeader.channels();

    // -----------------------------------------------------------------------
    // Write the pixel data
    // -----------------------------------------------------------------------

    // Layout framebuffer
    Imf::FrameBuffer     exrFrameBuffer;
    const Imf::PixelType compType = Imf::FLOAT;

    // Write RGB version
    std::vector<float> rgbImage;
    getRGBImage(rgbImage);

    const std::array<std::string, 3> rgbChannels = {"R", "G", "B"};
    const size_t                     xStrideRGB  = sizeof(float) * 3;
    const size_t                     yStrideRGB  = xStrideRGB * width();

    for (size_t c = 0; c < 3; c++) {
      char *ptrRGB = (char *)(&rgbImage[c]);
      exrChannels.insert(rgbChannels[c], Imf::Channel(compType));
      exrFrameBuffer.insert(
        rgbChannels[c],
        Imf::Slice(compType, ptrRGB, xStrideRGB, yStrideRGB));
    }

    // Write spectral version
    const size_t xStride = sizeof(float) * nSpectralBands();
    const size_t yStride = xStride * width();

    for (size_t s = 0; s < nStokesComponents(); s++) {
      for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
        // Populate channel name
        const std::string channelName
          = getEmissiveChannelName(s, _wavelengths_nm[wl_idx]);
        exrChannels.insert(channelName, Imf::Channel(compType));

        char *ptrS = (char *)(&_emissivePixelBuffers[s][wl_idx]);
        exrFrameBuffer.insert(
          channelName,
          Imf::Slice(compType, ptrS, xStride, yStride));
      }
    }

    if (isReflective()) {
      for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
        // Populate channel name
        const std::string channelName
          = getReflectiveChannelName(_wavelengths_nm[wl_idx]);
        exrChannels.insert(channelName, Imf::Channel(compType));

        char *ptrS = (char *)(&_reflectivePixelBuffer[wl_idx]);
        exrFrameBuffer.insert(
          channelName,
          Imf::Slice(compType, ptrS, xStride, yStride));
      }
    }

    // -----------------------------------------------------------------------
    // Write metadata
    // -----------------------------------------------------------------------

    if (lensTransmission().size() > 0) {
      exrHeader.insert(
        LENS_TRANSMISSION_ATTR,
        lensTransmission().getAttribute());
    }

    if (cameraResponse().size() > 0) {
      exrHeader.insert(CAMERA_RESPONSE_ATTR, cameraResponse().getAttribute());
    }

    if (channelSensitivities().size() > 0) {
      for (size_t wl_idx = 0; wl_idx < nSpectralBands(); wl_idx++) {
        if (channelSensitivity(wl_idx).size() > 0) {
          std::string channelName
            = getEmissiveChannelName(0, _wavelengths_nm[wl_idx]);

          exrHeader.insert(
            channelName,
            channelSensitivity(wl_idx).getAttribute());
        }
      }
    }

    exrHeader.insert(
      EXPOSURE_COMPENSATION_ATTR,
      Imf::StringAttribute(std::to_string(_ev)));

    // Polarisation handedness
    if (isPolarised()) {
      Imf::StringAttribute handednessAtrrValue(
        _polarisationHandedness == LEFT_HANDED ? "left" : "right");

      exrHeader.insert(
        POLARISATION_HANDEDNESS_ATTR,
        Imf::StringAttribute(handednessAtrrValue));
    }

    // -----------------------------------------------------------------------
    // Write file
    // -----------------------------------------------------------------------

    Imf::OutputFile exrOut(filename.c_str(), exrHeader);
    exrOut.setFrameBuffer(exrFrameBuffer);
    exrOut.writePixels(height());
  }


  SpectrumType EXRSpectralImage::channelType(
    const std::string &channelName,
    int &              polarisationComponent,
    double &           wavelength_nm)
  {
    const std::regex expr(
      "^((S([0-3]))|T)\\.(\\d*,?\\d*([Ee][+-]?\\d+)?)(Y|Z|E|P|T|G|M|k|h|"
      "da|d|c|m|u|n|p)?(m|Hz)$");
    std::smatch matches;

    const bool matched = std::regex_search(channelName, matches, expr);

    SpectrumType channelType = SpectrumType::UNDEFINED;

    if (matched) {
      if (matches.size() != 8) {
        // Something went wrong with the parsing. This shall not occur.
        throw INTERNAL_ERROR;
      }

      switch (matches[1].str()[0]) {
        case 'S':
          channelType           = SpectrumType::EMISSIVE;
          polarisationComponent = std::stoi(matches[3].str());
          if (polarisationComponent > 0) {
            channelType = channelType | SpectrumType::POLARISED;
          }
          break;

        case 'T':
          channelType = SpectrumType::REFLECTIVE;
          break;

        default:
          return SpectrumType::UNDEFINED;
      }

      // Get value
      std::string centralValueStr(matches[4].str());
      std::replace(centralValueStr.begin(), centralValueStr.end(), ',', '.');
      const double value = std::stod(centralValueStr);

      // Convert to nanometers
      const std::string prefix = matches[6].str();
      const std::string units  = matches[7].str();

      try {
        wavelength_nm = Util::strToNanometers(value, prefix, units);
      } catch (std::out_of_range &exception) {
        // Unknown unit or multiplier
        // Something went wrong with the parsing. This shall not occur.
        throw INTERNAL_ERROR;
      }
    }

    return channelType;
  }


  std::string EXRSpectralImage::getEmissiveChannelName(
    int stokesComponent, double wavelength_nm)
  {
    assert(stokesComponent < 4);

    std::stringstream b;
    std::string       wavelengthStr = std::to_string(wavelength_nm);
    std::replace(wavelengthStr.begin(), wavelengthStr.end(), '.', ',');

    b << "S" << stokesComponent << "." << wavelengthStr << "nm";

    const std::string channelName = b.str();

#ifndef NDEBUG
    int    polarisationComponentChecked;
    double wavelength_nmChecked;

    SpectrumType t = channelType(
      channelName,
      polarisationComponentChecked,
      wavelength_nmChecked);

    assert(isEmissiveSpectrum(t));
    assert(polarisationComponentChecked == stokesComponent);
    assert(wavelength_nmChecked == wavelength_nm);
#endif

    return channelName;
  }


  std::string EXRSpectralImage::getReflectiveChannelName(double wavelength_nm)
  {
    std::stringstream b;
    std::string       wavelengthStr = std::to_string(wavelength_nm);
    std::replace(wavelengthStr.begin(), wavelengthStr.end(), '.', ',');

    b << "T"
      << "." << wavelengthStr << "nm";

    const std::string channelName = b.str();

#ifndef NDEBUG
    int    polarisationComponentChecked;
    double wavelength_nmChecked;

    SpectrumType t = channelType(
      channelName,
      polarisationComponentChecked,
      wavelength_nmChecked);

    assert(isReflectiveSpectrum(t));
    assert(wavelength_nmChecked == wavelength_nm);
#endif

    return channelName;
  }

}   // namespace SEXR
