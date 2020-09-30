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

#include <string>
#include <map>
#include <stdexcept>

namespace SEXR {

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
            if (prefix == "n" && units == "m") return value;
            
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


        static size_t idxFromWavelengthIdx(
            size_t wlFrom_idx,
            size_t wlTo_idx
        ) {
            if (wlFrom_idx < wlTo_idx) {
                return wlTo_idx * (wlTo_idx - 1) / 2 + wlFrom_idx;
            } else {
                return -1;
            }
        }
};

} // namespace SEXR
