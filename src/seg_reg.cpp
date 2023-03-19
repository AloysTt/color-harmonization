#include <iostream>
#include <vector>
#include "image_ppm.h"
#include "HSV.h"

using namespace std;
typedef unsigned char uchar;


bool pixelSimilar(uchar* ImgIn, int seed, int pt, int seuil)
{
	ColorRGB rgbPixel{(float)ImgIn[3*pt] / 255.0f, (float)ImgIn[3*pt + 1] / 255.0f, (float)ImgIn[3*pt + 2] / 255.0f};
	ColorHSV hsvPixel{};
	rgb_to_hsv(rgbPixel, hsvPixel);
	ColorRGB rgbSeed{(float)ImgIn[3*seed] / 255.0f, (float)ImgIn[3*seed + 1] / 255.0f, (float)ImgIn[3*seed + 2] / 255.0f};
	ColorHSV hsvSeed{};
	rgb_to_hsv(rgbSeed, hsvSeed);

    int hDiff = min(abs(hsvPixel.h - hsvSeed.h), 180 - abs(hsvPixel.h - hsvSeed.h));
    int sDiff = abs(hsvPixel.s - hsvSeed.s);
    int vDiff = abs(hsvPixel.v - hsvSeed.v);
    int dist = sqrt(hDiff * hDiff + sDiff * sDiff + vDiff * vDiff);

    return (dist < seuil);
}

void regionGrow(uchar* imgIn, uchar* reg, int seed, int seuil, int nH, int nW)
{
    vector<int> points;
    points.push_back(seed);
    int i, j;

    while (!points.empty())
    {
        int pt = points.back();
        points.pop_back();
        i = pt%nW;
        j = pt/nW;

        if (i < 0 || j < 0 || i >= nW || j >= nH || reg[pt] == 255)
            continue;

        if (pixelSimilar(imgIn, seed, pt, seuil))
        {
            reg[pt] = 255;
            points.push_back(j*nW+i+1);
            points.push_back(j*nW+i-1);
            points.push_back((j+1)*nW+i);
            points.push_back((j-1)*nW+i);
        }
    }
}

int main(int argc, char **argv)
{

	if (argc != 2)
	{
		printf("Usage: ImageIn.ppm\n");
		return 1;
	}
	std::string imageName{argv[1]};
	std::string baseName = {imageName.substr(0, imageName.size()-4)};

    int nH, nW, taille, taille3;
    lire_nb_lignes_colonnes_image_ppm(imageName.data(), &nH, &nW);
	taille3 = nH * nW * 3;
	taille = nH*nW;
	uchar * ImgIn = new uchar[taille3];
	uchar * ImgOut = new uchar[taille];
    memset(ImgOut, 0, taille * sizeof(uchar)); // Init de l'image de sortie à 0
	lire_image_ppm(imageName.data(), ImgIn, taille);

    int seed = 50*nW+128; // Point de départ pour la croissance de région
    int seuil = 20; // Seuil de similitude

    regionGrow(ImgIn, ImgOut, seed, seuil, nH, nW); // Croissance de région à partir du seed

    //Erosion/dilatation

	ecrire_image_pgm((baseName+"_seg.ppm").data(), ImgOut, nH, nW);

    delete [] ImgIn; delete [] ImgOut;
    
    return 0;
}