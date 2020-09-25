#include <SpectrumAttribute.h>
#include <EXRSpectralImage.h>

#include <filesystem>
#include <string>
#include <algorithm>
#include <vector>
#include <regex>
#include <fstream>

#include <OpenEXR/ImfRgbaFile.h>
#include <OpenEXR/ImfArray.h>

using namespace SEXR;

bool checkExtension(const std::string& path) 
{
    if (path.length() < 4) { return false; }

    const char* c[2] = {".EXR", ".exr"};
    
    for (int i = 0; i < 5; i++) {
        if (path[path.length() - i] != c[0][4 - i]
         && path[path.length() - i] != c[1][4 - i]) {
            return false;
        }
    }

    return true;
}


void load_csv(
    const std::string& filename,
    std::vector<float>& wavelengths_nm,
    std::vector<float>& values)
{
    wavelengths_nm.clear();
    values.clear();

    const std::string floatRegex  = " *(\\d*\\.?\\d*([Ee][+-]?\\d+)?) *";
    const std::regex e(floatRegex + "," + floatRegex);

    std::ifstream inFile(filename);
    std::string line;

    while(std::getline(inFile, line)) {
        std::smatch matches;
        if (std::regex_search(line, matches, e)) {
            wavelengths_nm.push_back(std::stof(matches[1]));
            values.push_back(std::stof(matches[3]));
        }
    }
}


int main(int argc, char* argv[]) 
{
    if (argc < 5) {
        std::cout << "Usage:" << std::endl
                  << "------" << std::endl
                  << argv[0] << " <folder> <start_wl_nm> <increment_wl_nm> <output_file> <camera_response> <lens_transmission> <channels_sensitivity...>" << std::endl
                  << std::endl;

        return 0;
    }

    // List files in folder
    std::vector<std::string> files;
    const std::string allowed_extension = ".exr";

    for (auto& p: std::filesystem::directory_iterator(argv[1])) {
        std::filesystem::path path = p.path();

        if (std::filesystem::is_regular_file(path) && checkExtension(path)) {
            files.push_back(path);
        }
    }

    // Sort files
    std::sort(files.begin(), files.end());

    // Get wavelength info
    const float start_wl_nm = std::stof(argv[2]);
    const float increment_wl_nm = std::stof(argv[3]);

    size_t width(0), height(0);
    std::vector<Imf::Rgba> pixels;

    std::vector<float> wavelengths;
    wavelengths.reserve(files.size());

    std::vector<float> spectralFramebuffer;

    std::cout << "Using images:" << std::endl;

    for (size_t i = 0; i < files.size(); i++) {
        wavelengths.push_back(start_wl_nm + i * increment_wl_nm);
        Imf::RgbaInputFile file(files[i].c_str());
        std::cout << "\tAt " << wavelengths[i] << "nm: [" << files[i] << "]" << std::endl;

        Imath::Box2i dw = file.dataWindow();
        const size_t c_width  = dw.max.x - dw.min.x + 1;
        const size_t c_height = dw.max.y - dw.min.y + 1;

        if (i == 0) {
            width = c_width;
            height = c_height;
            spectralFramebuffer.resize(width * height * files.size());
        } else if (width != c_width || height != c_height) {
            std::cerr << "Image sizes does not match" << std::endl;
            return -1;
        }

        pixels.resize(height * width);
        file.setFrameBuffer(&pixels[0], 1, width);
        file.readPixels(dw.min.y, dw.max.y);

        // We can't memcpy: half to float conversion and taking every channel
        for (size_t j = 0; j < width * height; j++) {
            spectralFramebuffer[files.size() * j + i] = (pixels[j].r + pixels[j].g + pixels[j].b) / 3.F;
        }
    }

    // Now, create the spectral image
    EXRSpectralImage spectralImage(width, height, wavelengths, EMISSIVE);

    memcpy(&spectralImage(0, 0, 0, 0), &spectralFramebuffer[0], width * height * wavelengths.size() * sizeof(float));

    std::cout << std::endl;

    if (argc >= 6) {
        std::cout << "Using transmission / response spectrum information:" << std::endl;
    }

    if (argc >= 6) {
        std::vector<float> camera_wavelengths_nm;
        std::vector<float> camera_response;

        std::cout << "\tCamera response: [" << argv[5] << "]" << std::endl;

        load_csv(argv[5], camera_wavelengths_nm, camera_response);
        spectralImage.setCameraResponse(camera_wavelengths_nm, camera_response);
    }

    if (argc >= 7) {
        std::vector<float> lens_wavelengths_nm;
        std::vector<float> lens_transmission;

        std::cout << "\tLens transmission: [" << argv[6] << "]" << std::endl;

        load_csv(argv[6], lens_wavelengths_nm, lens_transmission);
        spectralImage.setLensTransmission(lens_wavelengths_nm, lens_transmission);
    }

    if (argc >= 8) {
        for (int c = 7; c < argc; c++) {
            std::vector<float> filter_wavelengths_nm;
            std::vector<float> filter_transmission;

            std::cout << "\tFilter transmission at: " << wavelengths[c - 7] << "nm: [" << argv[c] << "]" << std::endl;

            load_csv(argv[c], filter_wavelengths_nm, filter_transmission);
            spectralImage.setChannelSensitivity(c - 7, filter_wavelengths_nm, filter_transmission);
        }
    }

    spectralImage.save(argv[4]);
    std::cout << std::endl << "File saved as: [" << argv[4] << "]" << std::endl;

    return 0;
}