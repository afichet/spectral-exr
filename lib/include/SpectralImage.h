#pragma once

#include <vector>
#include <array>
#include <string>

#include <SpectrumAttribute.h>
#include <SpectrumType.h>

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
            size_t wavelength_idx, 
            size_t muellerRow, size_t muellerColumn);

        virtual const float& operator()(
            size_t x, size_t y, 
            size_t wavelength_idx, 
            size_t muellerRow, size_t muellerColumn) const;

        // Those are not direct memory access
        // They can be called whatever the image type is

        virtual float getEmissiveValue(
            size_t x, size_t y, 
            size_t wavelength_idx, 
            size_t stokesComponent = 0) const;

        virtual float getReflectiveValue(
            size_t x, size_t y, 
            size_t wavelength_idx, 
            size_t muellerRow = 0, size_t muellerColumn = 0) const;

        const float& wavelength_nm(size_t wl_idx) const;

        size_t width()  const;
        size_t height() const;

        size_t nSpectralBands()     const;
        size_t nStokesComponents()  const;
        size_t nMuellerComponents() const;
        
        bool polarised()    const;
        bool emissive()     const;
        bool reflective()   const;
        SpectrumType type() const;

        static void componentsFromIndex(
            size_t index,
            size_t& row,
            size_t& col
        );

        static size_t indexFromComponents(
            size_t row,
            size_t col
        );

    protected:
        size_t _width, _height;
                
        // We can have up to 20 pixel buffers:
        // - 1 for emissive unpolarised images (S0)
        // - 1 for reflective unpolarised images (M00)
        // - 4 for emissive polarised images (S0, S1, S2, S3)
        // - 16 for reflective polarised images (M00, M01, M02, M03, M10, ... M33)
        
        std::array<std::vector<float>, 4> _reflectivePixelBuffers;
        std::array<std::vector<float>, 16> _emissivePixelBuffers;

        std::vector<float> _wavelengths_nm;
        SpectrumType _spectrumType;

        SpectrumAttribute _lensTransmissionSpectra;
        SpectrumAttribute _cameraReponse;
        std::vector<SpectrumAttribute> _channelSensitivities;
};
