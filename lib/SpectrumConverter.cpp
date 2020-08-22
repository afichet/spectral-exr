#include "SpectrumConverter.h"
#include "MathUtil.h"
#include "spectrum_data.h"

#include <cstring>
#include <cassert>
#include <cmath>


SpectrumConverter::SpectrumConverter(SPECTRUM_TYPE type)
    : _emissiveSpectrum(type == EMISSIVE)
    , _cmfFirstWavelength_nm(CIE1931_2DEG_FIRST_WAVELENGTH_NM)
{
    if (!_emissiveSpectrum) {
        _illuminantFirstWavelenght_nm = D_65_FIRST_WAVELENGTH_NM;
        _illuminantSPD = std::vector<float>(std::begin(D_65_SPD), std::end(D_65_SPD));
    }

    _xyzCmfs[0] = std::vector<float>(std::begin(CIE1931_2DEG_X), std::end(CIE1931_2DEG_X));
    _xyzCmfs[1] = std::vector<float>(std::begin(CIE1931_2DEG_Y), std::end(CIE1931_2DEG_Y));
    _xyzCmfs[2] = std::vector<float>(std::begin(CIE1931_2DEG_Z), std::end(CIE1931_2DEG_Z));
       
    memcpy(&_xyzToRgb[0], XYZ_TO_SRGB_D65_MATRIX, 9 * sizeof(float));
}


SpectrumConverter::SpectrumConverter(
    const float& cmfFirstWavelength_nm,
    const std::array<std::vector<float>, 3> &xyzCmfs,
    const std::array<float, 9> xyzToRgb
)
    : _cmfFirstWavelength_nm(cmfFirstWavelength_nm)
    , _xyzCmfs(xyzCmfs)
    , _xyzToRgb(xyzToRgb)
{
}


void SpectrumConverter::spectrumToXYZ(
    const std::vector<float>& wavelengths_nm,
    const float* spectrum,
    std::array<float, 3>& XYZ
) const  {
    if (_emissiveSpectrum) {
        emissiveSpectrumToXYZ(wavelengths_nm, spectrum, XYZ);
    } else {
        reflectiveSpectrumToXYZ(wavelengths_nm, spectrum, XYZ);
    }
}


void SpectrumConverter::spectrumToRGB(
    const std::vector<float>& wavelengths_nm,
    const float* spectrum,
    std::array<float, 3>& RGB
) const {
    std::array<float, 3> XYZ;
    spectrumToXYZ(wavelengths_nm, spectrum, XYZ);

    memset(&RGB[0], 0, 3*sizeof(float));

    // Convert to RGB using the provided matrix
    for (size_t channel = 0; channel < 3; channel++) {
        for (size_t col = 0; col < 3; col++) {
            RGB[channel] += XYZ[col] * _xyzToRgb[3 * channel + col];
        }

        // Ensure RGB values are > 0
        RGB[channel] = std::max(RGB[channel], 0.F);        
    }
}


size_t SpectrumConverter::cmfWavelengthIndex(float wavelength_nm) 
const { 
    assert(wavelength_nm >= firstWavelength());
    assert(wavelength_nm <= lastWavelength());

    float idx_f = float(wavelength_nm - _cmfFirstWavelength_nm);

    assert(idx_f >= 0);
    assert(idx_f < _xyzCmfs[0].size());

    return static_cast<size_t>(std::round(idx_f));
}


size_t SpectrumConverter::cmfWavelengthValue(size_t index) 
const { 
    // assert(index < _xyzCmfs[0].size());
    if (index >= _xyzCmfs[0].size()) return 0;
    return _cmfFirstWavelength_nm + index;
}


void SpectrumConverter::emissiveSpectrumToXYZ(
    const std::vector<float>& wavelengths_nm,
    const float* spectrum,
    std::array<float, 3>& XYZ
) const {
    memset(&XYZ[0], 0, 3*sizeof(float));

    if (wavelengths_nm.size() == 0) {
        return;
    }

    const float start_wavelength = std::max(firstWavelength(), wavelengths_nm.front());
    const float end_wavelength   = std::min(lastWavelength(), wavelengths_nm.back());

    // Early exit, selection out of range
    if (end_wavelength < start_wavelength) {
        return;
    }

    assert(start_wavelength <= end_wavelength);

    for (size_t idx_value = 0; idx_value < wavelengths_nm.size() - 1; idx_value++) {
        float wl_a = wavelengths_nm[idx_value];
        float wl_b = wavelengths_nm[idx_value + 1];

        // We have not reached yet the starting point
        if (start_wavelength > wl_b) {
            continue;
        }

        // We have finished the integration
        if (end_wavelength < wl_a) {
            break;
        }

        if (start_wavelength > wl_a) {
            wl_a = start_wavelength;
        }

        if (end_wavelength < wl_b) {
            wl_b = end_wavelength;
        }

        const size_t idx_cmf_start = cmfWavelengthIndex(wl_a);
        size_t       idx_cmf_end   = cmfWavelengthIndex(wl_b);

        // On last intervall we need to include the last wavelength of the spectrum
        if (idx_value == wavelengths_nm.size() - 2) {
            idx_cmf_end = idx_cmf_end + 1;
        }

        for (size_t idx_cmf = idx_cmf_start; idx_cmf < idx_cmf_end; idx_cmf++) {
            const float curr_wl    = cmfWavelengthValue(idx_cmf);
            const float curr_value = MathUtil::interp(
                curr_wl,
                wavelengths_nm[idx_value],
                wavelengths_nm[idx_value + 1],
                spectrum[idx_value],
                spectrum[idx_value + 1]);

            for (size_t c = 0; c < 3; c++) {
                XYZ[c] += _xyzCmfs[c][idx_cmf] * curr_value;
            }
        }
    }
}


void SpectrumConverter::reflectiveSpectrumToXYZ(
    const std::vector<float>& wavelengths_nm,
    const float* spectrum,
    std::array<float, 3>& XYZ
) const {
    memset(&XYZ[0], 0, 3*sizeof(float));

    if (wavelengths_nm.size() == 0) {
        return;
    }

    const float illuminant_last_wavelength = _illuminantFirstWavelenght_nm + _illuminantSPD.size() - 1;
    const float start_wavelength = std::max(std::max(_illuminantFirstWavelenght_nm, firstWavelength()), wavelengths_nm.front());
    const float end_wavelength   = std::min(std::min(illuminant_last_wavelength, lastWavelength()), wavelengths_nm.back());

    // Early exit, selection out of range
    if (end_wavelength < start_wavelength) {
        return ;
    }

    assert(start_wavelength <= end_wavelength);

    float normalisation_factor(0);

    for (size_t idx_value = 0; idx_value < wavelengths_nm.size() - 1; idx_value++) {
        float wl_a = wavelengths_nm[idx_value];
        float wl_b = wavelengths_nm[idx_value + 1];

        // We have not reached yet the starting point
        if (start_wavelength > wl_b) {
            continue;
        }

        // We have finished the integration
        if (end_wavelength < wl_a) {
            break;
        }

        if (start_wavelength > wl_a) {
            wl_a = start_wavelength;
        }

        if (end_wavelength < wl_b) {
            wl_b = end_wavelength;
        }

        const size_t idx_curve_start = cmfWavelengthIndex(wl_a);
        size_t       idx_curve_end   = cmfWavelengthIndex(wl_b);

        // On last intervall we need to include the last wavelength of the spectrum
        if (idx_value == wavelengths_nm.size() - 2) {
            idx_curve_end = idx_curve_end + 1;
        }

        for (size_t idx_curve = idx_curve_start; idx_curve < idx_curve_end; idx_curve++) {
            const float curr_wl = cmfWavelengthValue(idx_curve);

            const size_t idx_illu_a = curr_wl - _illuminantFirstWavelenght_nm;
            assert(curr_wl >= _illuminantFirstWavelenght_nm);
            assert(idx_illu_a < _illuminantSPD.size());

            const float illu_value = _illuminantSPD[idx_illu_a];

            normalisation_factor += illu_value * _xyzCmfs[1][idx_curve]; // Y

            const float curr_value = 
                illu_value * MathUtil::interp(
                                curr_wl,
                                wavelengths_nm[idx_value],
                                wavelengths_nm[idx_value + 1],
                                spectrum[idx_value],
                                spectrum[idx_value + 1]);

            for (size_t c = 0; c < 3; c++) {
                XYZ[c] += curr_value * _xyzCmfs[c][idx_curve];
            }
        }
    }

    for (size_t c = 0; c < 3; c++) {
        XYZ[c] /= normalisation_factor;
    }
}
