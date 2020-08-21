#pragma once

#include <vector>
#include <array>
#include <string>

class SpectralImage {
    public:
        enum Errors {
            UNSUPORTED_FILE,
            INTERNAL_ERROR,
            READ_ERROR,
            WRITE_ERROR,
            INCORECTED_FORMED_FILE
        };

        SpectralImage(
            size_t width = 0, size_t height = 0,
            const std::vector<float>& wavelengths_nm = std::vector<float>(),
            bool containsPolarisationData = false
        );

        virtual void save(const std::string& filename) const = 0;

        virtual void exportChannels(const std::string& path) const;

        size_t width() const { return _width; }

        size_t height() const { return _height; }

        size_t nSpectralBands() const { return _wavelengths_nm.size(); }

        size_t nStokesComponents() const { return (_containsPolarisationData ? 4 : 1); }

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
};
