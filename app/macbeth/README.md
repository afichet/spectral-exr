# Macbeth

This sample code exports a reflective spectral OpenEXR image representing a Macbeth colour chart.

It contains both the spectral data for each patch. Data comes from http://www.babelcolor.com/colorchecker-2.htm. It also generates the RGB layers to allow the EXR image to be displayed in a traditionnal graphic software like Gimp.

Since this is a RGB reflective image, the spectral data are multiplied by a D65 illuminant and renormalised. The renormalisation is the divistion of the integral of the the spectral power distribution of CIE D65 with the colour matching function used, CIE 1931 XYZ 2 degress observer in this instance. 