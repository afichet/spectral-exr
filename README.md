# OpenEXR Spectral Image
This is an example code for reading and writing OpenEXR spectral images. The main code is placed in `lib` folder.

## Compilation
To compile this code, you need a C++11 compliant compiler, the OpenEXR library installed on your system and CMake.

```bash
mkdir build
cmake ..
make
```

This will compile an example program.

## Sample programs
You can find in `app` several sample programs using spectral OpenEXR.

# Mabeth
Macbeth executable (`macbeth`) generates a Macbeth colour chart using the spectral data from: http://www.babelcolor.com/colorchecker-2.htm. It will generate a file `Macbeth.exr` from the execution place.

```bash
./bin/macbeth
```

# Fluo EXR
Fluo EXR executable (`fluo-exr`) generates a simple checkerboard made of two fluorescent patches: 3M fluorescent yellow Post-It (R) sticker and 3M fluorescent pink Post-It (R) sticker. Data courtesy of Labsphere Inc. It will generate a file `BiSpectral.exr`

```bash
./bin/fluo-exr
```

# Export spectrum
Export spectrum executable (`export-spectrum`) extract from the given pixel location the stored spectrum in an ASCII file. Each columns correspond to a polarisation component if present in the image. $S_0$, $S_1$, $S_2$, $S_3$ for emissive images and $M_{00}$ ... $M_{33}$ for reflective images.

It takes as arguments:
- A spectral EXR
- The x coordinate of the pixel to extract
- The y coordinate of the pixel to extract
- The output file

```bash
./bin/export-spectrum Macbeth.exr 15 15 Macbeth.txt
```

## License
Copyright (c) 2020, <Authors>

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.