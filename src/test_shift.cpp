#include <random>
#include <algorithm>
#include "image_ppm.h"
#include "harmonization.h"

typedef unsigned char uchar;

int main(int argc, char **argv)
{

	if (argc != 2)
	{
		printf("Usage: ImageIn.ppm\n");
		return 1;
	}
	std::string imageName{argv[1]};
	std::string baseName = {imageName.substr(0, imageName.size()-4)};

	Distribution distrib = Distribution::create_I_template(60.0f, 60.0f);
	int h, w, size3, size;
	lire_nb_lignes_colonnes_image_ppm(imageName.data(), &h, &w);
	size3 = h * w * 3;
	size = h*w;
	uchar * image = new uchar[size3];
	lire_image_ppm(imageName.data(), image, size);

	uchar *imageOut = new uchar[size3];
	shift(distrib, h, w, image, imageOut);
	ecrire_image_ppm((baseName+"_shift.ppm").data(), imageOut, h, w);


	delete [] image;
	delete [] imageOut;
	return 0;
}

