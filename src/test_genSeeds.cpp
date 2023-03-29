#include <iostream>
#include <vector>
#include <cstring>
#include <algorithm>
#include "image_ppm.h"
#include "ColorSpaces.h"
#include "segmentation.h"
#include "util.h"

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned short ushort;
/*using namespace std;
typedef unsigned char uchar;


bool pixelSimilar(uchar* ImgIn, int seed, int pt, int seuil)
{
	ColorRGB rgbPixel{(float)ImgIn[3*pt] / 255.0f, (float)ImgIn[3*pt + 1] / 255.0f, (float)ImgIn[3*pt + 2] / 255.0f};
	ColorHSV hsvPixel{};
	rgb_to_hsv(rgbPixel, hsvPixel);
	ColorRGB rgbSeed{(float)ImgIn[3*seed] / 255.0f, (float)ImgIn[3*seed + 1] / 255.0f, (float)ImgIn[3*seed + 2] / 255.0f};
	ColorHSV hsvSeed{};
	rgb_to_hsv(rgbSeed, hsvSeed);

    float hDiff = min(abs(hsvPixel.h - hsvSeed.h), 180 - abs(hsvPixel.h - hsvSeed.h));
    hDiff /= 360;
    float sDiff = abs(hsvPixel.s - hsvSeed.s);
    float vDiff = abs(hsvPixel.v - hsvSeed.v);
    float dist = sqrt(hDiff * hDiff + sDiff * sDiff + vDiff * vDiff);
    dist*=360;

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

void rgbToYcbcr(uchar* ImgIn, uchar* ImgOut, int nH, int nW){
    int nW3 = nW*3;
    for (int i=0; i<nH; i++){
        for(int j=0; j<nW; j++){
            int place = i*nW3+j*3;
            float R = ImgIn[place];
            float G = ImgIn[place+1];
            float B = ImgIn[place+2];
            ImgOut[place] = 0.299f*R + 0.587*G + 0.114f*B;
            ImgOut[place+1] = -0.1687f*R - 0.3313f*G + 0.5f*B +128.0f;
            ImgOut[place+2] = 0.5f*R - 0.4187f*G - 0.0813f*B + 128.0f;
        }
    }
}

void ycbcrTorgb(uchar* ImgIn, uchar* ImgOut){
    int nW3 = nW*3;
    for (int i=0; i<nH; i++){
        for(int j=0; j<nW; j++){
            int place = i*nW3+j*3;
            float Y = ImgIn[place];
            float Cb = ImgIn[place+1];
            float Cr = ImgIn[place+2];
            ImgOut[place] = clamp(Y + 1.402f*(Cr - 128.0f));
            ImgOut[place+1] = clamp(Y - 0.34414f*(Cb - 128.0f) - 0.71414f*(Cr -128.0f));
            ImgOut[place+2] = clamp(Y + 1.772f*(Cb-128.0f));
        }
    }
}*/


/*
double autoSeuil(uchar* ImgIn, int nH, int nW){
    //int seuil;
    double max=0;
    double otsu;
    for(int i=0; i<256; i++){
        otsu = otsu(ImgIn, i, nH, nW);
        if(otsu>max){
            max = otsu;
            //seuil = i;
        }
    }
    return max;
}

void seedList(uchar* ImgIn, uchar* ImgOut, int nH, int nW){
    double similiSeuil = autoSeuil(ImgIn, nH, nW);
    int taille = nH*nW;
    int nW3 = nW*3;
    int placeY, placeCb, placeCr;
	uchar * ImgCond1 = new uchar[taille];
    double vY, vCb, vCr, moy;
    double tabV[taille];
    double tabVNorm[taille];
    for (int i=1; i<nH-1; i++){
        for(int j=1; j<nW-1; j++){
            placeY = i*nW3+j*3;
            moy=(ImgIn[placeY]
                +ImgIn[placeY+3]
                +ImgIn[placeY-3]
                +ImgIn[placeY-nW3]
                +ImgIn[placeY-nW3-3]
                +ImgIn[placeY-nW3+3]
                +ImgIn[placeY+nW3]
                +ImgIn[placeY+nW3+3]
                +ImgIn[placeY+nW3-3]
                )/9;
            vY=sqrt((pow(ImgIn[placeY]-moy,2)
                    +pow(ImgIn[placeY-3]-moy,2)
                    +pow(ImgIn[placeY+3]-moy,2)
                    +pow(ImgIn[placeY-nW3]-moy,2)
                    +pow(ImgIn[placeY-3-nW3]-moy,2)
                    +pow(ImgIn[placeY+3-nW3]-moy,2)
                    +pow(ImgIn[placeY+nW3]-moy,2)
                    +pow(ImgIn[placeY-3+nW3]-moy,2)
                    +pow(ImgIn[placeY+3+nW3]-moy,2)
                    )/9
                    );
            placeCb = placeY+1;
            moy=(ImgIn[placeCb]
                +ImgIn[placeCb+3]
                +ImgIn[placeCb-3]
                +ImgIn[placeCb-nW3]
                +ImgIn[placeCb-nW3-3]
                +ImgIn[placeCb-nW3+3]
                +ImgIn[placeCb+nW3]
                +ImgIn[placeCb+nW3+3]
                +ImgIn[placeCb+nW3-3]
                )/9;
            vCb=sqrt((pow(ImgIn[placeCb]-moy,2)
                    +pow(ImgIn[placeCb-3]-moy,2)
                    +pow(ImgIn[placeCb+3]-moy,2)
                    +pow(ImgIn[placeCb-nW3]-moy,2)
                    +pow(ImgIn[placeCb-3-nW3]-moy,2)
                    +pow(ImgIn[placeCb+3-nW3]-moy,2)
                    +pow(ImgIn[placeCb+nW3]-moy,2)
                    +pow(ImgIn[placeCb-3+nW3]-moy,2)
                    +pow(ImgIn[placeCb+3+nW3]-moy,2)
                    )/9
                    );
            placeCr = placeY+2;
            moy=(ImgIn[placeCr]
                +ImgIn[placeCr+3]
                +ImgIn[placeCr-3]
                +ImgIn[placeCr-nW3]
                +ImgIn[placeCr-nW3-3]
                +ImgIn[placeCr-nW3+3]
                +ImgIn[placeY+nW3]
                +ImgIn[placeCr+nW3+3]
                +ImgIn[placeCr+nW3-3]
                )/9;
            vCb=sqrt((pow(ImgIn[placeCr]-moy,2)
                    +pow(ImgIn[placeCr-3]-moy,2)
                    +pow(ImgIn[placeCr+3]-moy,2)
                    +pow(ImgIn[placeCr-nW3]-moy,2)
                    +pow(ImgIn[placeCr-3-nW3]-moy,2)
                    +pow(ImgIn[placeCr+3-nW3]-moy,2)
                    +pow(ImgIn[placeCr+nW3]-moy,2)
                    +pow(ImgIn[placeCr-3+nW3]-moy,2)
                    +pow(ImgIn[placeCr+3+nW3]-moy,2)
                    )/9
                    );
            tabV[i*nW+j]= vY + vCb + vCr;
        }
    }
}*/

void getSeedImageFromOtsuSeeds(const uchar * otsu, const float * imageYCbCr, uchar * out, int h, int w, float maxDistance);

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		printf("Usage: ImageIn.ppm\n");
		return 1;
	}
	std::string imageName{argv[1]};
	std::string baseName = {imageName.substr(0, imageName.size()-4)};

    int h, w, size, size3;
    lire_nb_lignes_colonnes_image_ppm(imageName.data(), &h, &w);
	size3 = h * w * 3;
	size = h * w;
	uchar * image = new uchar[size3];
	lire_image_ppm(imageName.data(), image, size);

	// convert to ycbcr float
	float * imageYCbCr = new float[size3];
	rgb_to_ycbcr(image, imageYCbCr, h, w);

	// get std dev for each pixel
	float * imageSimilarity = new float[size];
	computeSimilarityImage(imageYCbCr, imageSimilarity, h, w);

	// DEBUG: write similarity image
	uchar * testImageSimilarity = new uchar[size];
	for (int i=0; i<size; ++i)
		testImageSimilarity[i] = 255.0f * imageSimilarity[i];
	ecrire_image_pgm((baseName+"_similarity.pgm").data(), testImageSimilarity, h, w);
	delete [] testImageSimilarity;

	// find threshold using otsu
	// we first have to split the image into N buckets (here N=100)
	uint * imageOtsu = new uint[size];
	for (int i=0; i<size; ++i)
		imageOtsu[i] = 99.0f*imageSimilarity[i];
	int threshold = 0;
	double thresholdVal = 0;

	// probability distribution
	double * ddp = new double[100];
	create_ddp(imageOtsu, h, w, ddp, 100);
	for (int i=0; i<100; ++i)
	{
		double otsuVal = otsu(i, ddp, 100, OtsuCriterion::BCW_WCV);
		std::cout << "Threshold " << std::to_string(i) << ": " << std::to_string(otsuVal) << std::endl;
		if (otsuVal > thresholdVal)
		{
			threshold = i;
			thresholdVal = otsuVal;
		}
	}
	delete [] ddp;

	// create binary seed image
	uchar * imageSeedOtsu = new uchar[size];
	for (int i=0; i<size; ++i)
		imageSeedOtsu[i] = imageOtsu[i] >= threshold ? 255 : 0;
	std::cout << "Max threshold: " << std::to_string(threshold) << std::endl;

	// 2nd condition : euclidean distance
	uchar * imageSeed = new uchar[size];
	getSeedImageFromOtsuSeeds(imageSeedOtsu, imageYCbCr, imageSeed, h, w, 0.02f);
	ecrire_image_pgm((baseName+"_seeds.pgm").data(), imageSeed, h, w);

	delete [] imageSeed;
	delete [] imageSeedOtsu;
	delete [] imageOtsu;
	delete [] imageSimilarity;
	delete [] imageYCbCr;
	delete [] image;




	/*uchar * ImgOut = new uchar[size];
    memset(ImgOut, 0, size * sizeof(uchar)); // Init de l'image de sortie à 0
    int inLigne, inCol, inSeuil;
    sscanf (argv[2],"%d",&inLigne);
    sscanf (argv[3],"%d",&inCol);
    sscanf (argv[4],"%d",&inSeuil);
    int ligne = (inLigne/100.f)*h;
    int colonne = (inCol/100.f)*w;

    int seed = ligne*w+colonne; // Point de départ pour la croissance de région
    int seuil = inSeuil; // Seuil de similitude

    regionGrow(ImgIn, ImgOut, seed, seuil, h, w); // Croissance de région à partir du seed

    //Erosion/dilatation

	ecrire_image_pgm((baseName+"_seg.ppm").data(), ImgOut, h, w);*/

//    delete [] ImgIn; delete [] ImgOut;
    
    return 0;
}



void getSeedImageFromOtsuSeeds(const uchar * otsu, const float * imageYCbCr, uchar * out, int h, int w, float maxDistance)
{
	for (int row=1; row<h-1; ++row)
	{
		for (int col=1; col<w-1; ++col)
		{
			if (otsu[row*w+col] == 0)
			{
				out[row*w+col] = 0;
				continue;
			}
			float max = 0.0f;
			for (int i=-1; i<=1; ++i)
			{
				for (int j=-1; j<=1; ++j)
				{
					float dist = ycbcr_distance_relative_euclidean(
						imageYCbCr[row*w*3+col*3],
						imageYCbCr[row*w*3+col*3+1],
						imageYCbCr[row*w*3+col*3+2],
						imageYCbCr[(row+i)*w*3+(col+j)*3],
						imageYCbCr[(row+i)*w*3+(col+j)*3+1],
						imageYCbCr[(row+i)*w*3+(col+j)*3+2]
						);
					if (dist > max)
						max = dist;
				}
			}
			if (max > maxDistance)
				out[row*w+col] = 0;
			else
				out[row*w+col] = 255;
		}
	}
}