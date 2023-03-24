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
}

void histo(uchar *ImgIn, int* T, int nH, int nW){
    int taille = nH*nW;
    for (int i=0; i < taille; i++)
	    T[ImgIn[i]]++;
}

double otsu(uchar* ImgIn, int k, int nH, int nW){
    int T[256];
    memset(T, 0, 256);
    histo(ImgIn, T, nH, nW);
    double w0 = 0;
    double u = 0;
    double v0 = 0;
    double v1 = 0;
    double w1, probi, ur, u0, u1;
    taille = nH*nW;

    for(int i=0 i<k; i++){
        probi = T[i]/taille;
        w0+=probi;
        u+=i*probi;
    }
    ur = u;
    for(int i=k; i<taille; i++)
        ur += i*(T[i]/taille);

    w1=1-w0;
    u0 = u/w0;
    u1 = (ur-u)/w1;

    for(int i=0; i<k; i++)
        v0 += (i-u0)*(i-u0)*((T[i]/taille)/w0);
    for(int i=k; i<taille; i++)
        v1 += (i-u1)*(i-u1)*((T[i]/taille)/w1);

    double BCV = w0*w1*(u1-u0)*(u1-u0);
    double WCV = w0*v0+w1*v1;

    return BCV/WCV;
}

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

}

int main(int argc, char **argv)
{

	if (argc != 5)
	{
		printf("Usage: ImageIn.ppm ligne colonne seuil\n");
		return 1;
	}
   //sscanf (argv[3],"%d",&S);
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
    int inLigne, inCol, inSeuil;
    sscanf (argv[2],"%d",&inLigne);
    sscanf (argv[3],"%d",&inCol);
    sscanf (argv[4],"%d",&inSeuil);
    int ligne = (inLigne/100.f)*nH;
    int colonne = (inCol/100.f)*nW;

    int seed = ligne*nW+colonne; // Point de départ pour la croissance de région
    int seuil = inSeuil; // Seuil de similitude

    regionGrow(ImgIn, ImgOut, seed, seuil, nH, nW); // Croissance de région à partir du seed

    //Erosion/dilatation

	ecrire_image_pgm((baseName+"_seg.ppm").data(), ImgOut, nH, nW);

    delete [] ImgIn; delete [] ImgOut;
    
    return 0;
}