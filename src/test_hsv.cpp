#include <iostream>
#include "HSV.h"

int main(int argc, char **argv)
{
	ColorRGB rgb{0.6f, 0.2f, 0.7f};
	ColorHSV hsv{};
	ColorRGB res{};

	rgb_to_hsv(rgb, hsv);
	hsv_to_rgb(hsv, res);

	std::cout << "RGB: " << rgb.r << " " << rgb.g << " " << rgb.b << std::endl;
	std::cout << "HSV: " << hsv.h << " " << hsv.s << " " << hsv.v << std::endl;
	std::cout << "RGB (converted back): " << res.r << " " << res.g << " " << res.b << std::endl;
	return 0;
}
