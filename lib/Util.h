#pragma once

#include <string>
#include <map>
#include <stdexcept>

class Util {
    public:
        static float interp(float x, float x0, float x1, float y0, float y1)
        {
            return lerp(y0, y1, alpha(x0, x1, x));
        }

        static float alpha(float x0, float x1, float x)
        {
            return (x - x0) / (x1 - x0);
        }

        static float lerp(float a, float b, float t)
        {
            return a + t * (b - a);
        }

        static float strToNanometers(
            const float& value, 
            const std::string& prefix,
            const std::string& units) 
        {
            float wavelength_nm = value;

            const std::map<std::string, float> unit_prefix = {
                {"Y", 1e24}, {"Z", 1e21}, {"E", 1e18}, {"P", 1e15}, {"T", 1e12},
                {"G", 1e9} , {"M", 1e6} , {"k", 1e3} , {"h", 1e2} , {"da", 1e1},
                {"d", 1e-1}, {"c", 1e-2}, {"m", 1e-3}, {"u", 1e-6}, {"n", 1e-9},
                {"p", 1e-12}
            };

            // Apply multiplier
            if (prefix.size() > 0) {
                wavelength_nm *= unit_prefix.at(prefix);
            }

            // Apply units
            if (units == "Hz") {
                wavelength_nm = 299792458.F/wavelength_nm * 1e9;
            } else if (units == "m") {
                wavelength_nm = wavelength_nm * 1e9;
            } else {
                // Unknown unit
                // Something went wrong with the parsing. This shall not occur.
                throw std::out_of_range("Unknown unit");
            }
            
            return wavelength_nm;
        }
};