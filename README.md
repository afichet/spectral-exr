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

# Sample programs
You can find in `app` several sample programs using spectral OpenEXR.

## Macbeth
Macbeth executable (`macbeth`) generates a Macbeth colour chart using the spectral data from: http://www.babelcolor.com/colorchecker-2.htm. It will generate a file `Macbeth.exr` from the execution place.

```bash
./bin/macbeth
```

## Fluo EXR
Fluo EXR executable (`fluo-exr`) generates a simple checker-board made of two fluorescent patches: 3M fluorescent yellow Post-It (R) sticker and 3M fluorescent pink Post-It (R) sticker. Data courtesy of Labsphere Inc. It will generate a file `BiSpectral.exr`

```bash
./bin/fluo-exr
```

## Export spectrum
Export spectrum executable (`export-spectrum`) extracts from the given pixel location the stored spectrum in an ASCII file. Each columns correspond to a polarisation component if present in the image. $S_0$, $S_1$, $S_2$, $S_3$ for emissive images and $M_{00}$ ... $M_{33}$ for reflective images. First column is the wavelength in manometers. Each column is separated by a space.

It takes as arguments:
- A spectral EXR
- The x coordinate of the pixel to extract
- The y coordinate of the pixel to extract
- The output file

```bash
./bin/export-spectrum Macbeth.exr 15 15 Macbeth.txt
```

You can plot the spectrum using `gnuplot`
```gnuplot
plot "Macbeth.txt" u 1:2 w line
```

## Export reradiation
Export reradiation executable (`export-reradiation`) extracts from the give pixel location the stored reradiation matrix in a ASCII file.

It takes as arguments:
- A bi-spectral EXR
- The x coordinate of the pixel to extract
- The y coordinate of the pixel to extract
- The output file

```bash
./bin/export-spectrum BiSpectral.exr 15 15 reradiation.txt
```

You can plot the reradiation matrix using `gnuplot`
```gnuplot
plot "reradiation.txt" matrix w image
```

## Merge EXR
Merge EXR executable (`merge-exr`) creates an emissive spectral EXR from a folder containing collection of monochromatic EXRs. You can use the Cornell box data as input http://www.graphics.cornell.edu/online/box/data.html.

It takes as arguments:
- Folder path containing the images
- Starting wavelength (in manometers)
- Increment of wavelength between images (in manometers)
- Output file
- Optional additional arguments:
  - Camera response in CSV format (comma separated)
  - Lens transmission in CSV format (comma separated)
  - Each filter corresponding to each channel transmissions in CSV format (comma separated)

To convert the Matlab matrices of the Cornell data in CSV format, we use the following GNU Octave code:
```Matlab
load filters.mat

csvwrite("F1.csv", [wavelen F1])
csvwrite("F2.csv", [wavelen F2])
csvwrite("F3.csv", [wavelen F3])
csvwrite("F4.csv", [wavelen F4])
csvwrite("F5.csv", [wavelen F5])
csvwrite("F6.csv", [wavelen F6])
csvwrite("F7.csv", [wavelen F7])


load lens.mat
csvwrite("lens.csv", [wavelen lens])

load camera.mat
csvwrite("camera.csv", [wavelen_cam.' response.'])
```

Then, we execute the program as follows:
```bash
./bin/merge-exr \
    data/cornell 400 50 output/CornellBox.exr \
    data/cornell/camera.csv \
    data/cornell/lens.csv \
    data/cornell/F1.csv \
    data/cornell/F2.csv \
    data/cornell/F3.csv \
    data/cornell/F4.csv \
    data/cornell/F5.csv \
    data/cornell/F6.csv \
    data/cornell/F7.csv
```

## Spectrum to EXR
Spectrum to EXR executable (`spectrum-to-exr`) creates a 1x1px spectral OpenEXR from a given spectrum.

It takes as arguments:
- A spectrum in CSV format (comma separated)
- Type of spectrum (either `reflective` or `emissive`)
- Output file

For example:
```bash
./bin/spectrum-to-exr data/D65.csv emissive output/D65.exr
```

## Split channels
Split channels executable (`split-channels`) splits a spectral EXR is separate single monochromatic EXR files.

It takes as arguments:
- A spectral EXR file
- An output folder

For example:
```bash
./bin/split-channels Macbeth.exr output
```
