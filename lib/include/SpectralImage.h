#pragma once

#include <vector>
#include <array>
#include <string>

#include <SpectrumAttribute.h>


class SpectralImage {
    public:
        enum Errors {
            UNSUPORTED_FILE,
            INTERNAL_ERROR,
            READ_ERROR,
            WRITE_ERROR,
            INCORRECT_FORMED_FILE
        };

        enum SpectrumType {
            REFLECTIVE,
            EMISSIVE
        };

        SpectralImage(
            size_t width = 0, size_t height = 0,
            const std::vector<float>& wavelengths_nm = std::vector<float>(),
            SpectrumType type = EMISSIVE,
            bool containsPolarisationData = false
        );

        virtual void save(const std::string& filename) const = 0;

        virtual void exportChannels(const std::string& path) const;
        virtual void getRGBImage(std::vector<float>& rgbImage) const;

        size_t width()  const { return _width; }
        size_t height() const { return _height; }

        size_t nSpectralBands()    const { return _wavelengths_nm.size(); }
        size_t nStokesComponents() const { return (_containsPolarisationData ? 4 : 1); }
        
        bool polarised()    const { return _containsPolarisationData; }
        bool emissive()     const { return _spectrumType == EMISSIVE; }
        bool reflective()   const { return _spectrumType == REFLECTIVE; }
        SpectrumType type() const { return _spectrumType; }



        const float& wavelength_nm(size_t wl_idx) const { return _wavelengths_nm[wl_idx]; }

        void setCameraResponse(
            const std::vector<float>& wavelengths_nm,
            const std::vector<float>& values
        );

        void setLensTransmission(
            const std::vector<float>& wavelengths_nm,
            const std::vector<float>& values
        );

        void setChannelSensitivity(
            size_t wl_idx,
            const std::vector<float>& wavelengths_nm,
            const std::vector<float>& values
        );

        const SpectrumAttribute& getLensTransmission() const { return _lensTransmissionSpectra; }
        const SpectrumAttribute& getCameraResponse()   const { return _cameraReponse; }
        const SpectrumAttribute& getChannelSensitivity(size_t wl_idx) const { return _channelSensitivity[wl_idx]; }

        float& operator()(
            size_t x, size_t y, 
            size_t wavelength_idx, size_t stokes = 0);

        const float& operator()(
            size_t x, size_t y, 
            size_t wavelength_idx, size_t stokes = 0) const;

    protected:
        size_t _width, _height;
        std::array<std::vector<float>, 4> _pixelBuffers;
        std::vector<float> _wavelengths_nm;
        bool _isSpectral;
        bool _containsPolarisationData;
        SpectrumType _spectrumType;

        SpectrumAttribute _lensTransmissionSpectra;
        SpectrumAttribute _cameraReponse;
        std::vector<SpectrumAttribute> _channelSensitivity;
};
