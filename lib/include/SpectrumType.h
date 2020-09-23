#pragma once

enum SpectrumType {
    UNDEFINED  = 0,
    REFLECTIVE = 2,
    EMISSIVE   = 4,
    BISPECTRAL = 8,
    POLARISED  = 16
};


inline SpectrumType operator|(SpectrumType a, SpectrumType b)
{
    return static_cast<SpectrumType>(static_cast<int>(a) | static_cast<int>(b));
}

