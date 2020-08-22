#include <SpectrumAttribute.h>


int main(int argc, char* argv[]) {

    SpectrumAttribute sa(Imf::StringAttribute("195Hz,12;512nm,13;800.2E-3kHz,30;200.2E-3kHz,30;"));

    for (size_t i = 0; i < sa.size(); i++) {
        std::cout << sa.wavelength_nm(i) << " ";
        std::cout << sa.value(i) << std::endl;
    }

    return 0;
}