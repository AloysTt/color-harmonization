#include <iostream>
#include "ColorSpaces.h"

int main(int argc, char **argv)
{
	ColorRGB rgb{0.6f, 0.2f, 0.7f};
	ColorYCbCr ycbcr{};
	ColorRGB res{};

	rgb_to_ycbcr(rgb, ycbcr);
	ycbcr_to_rgb(ycbcr, res);

	std::cout << "RGB: " << rgb.r << " " << rgb.g << " " << rgb.b << std::endl;
	std::cout << "YCbCr: " << ycbcr.y << " " << ycbcr.cb << " " << ycbcr.cr << std::endl;
	std::cout << "RGB (converted back): " << res.r << " " << res.g << " " << res.b << std::endl;
	return 0;
}
