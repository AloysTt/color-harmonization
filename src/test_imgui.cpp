// Dear ImGui: standalone example application for GLFW + OpenGL 3, using programmable pipeline
// (GLFW is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan/Metal graphics context creation, etc.)
// If you are new to Dear ImGui, read documentation from the docs/ folder + read the top of imgui.cpp.
// Read online: https://github.com/ocornut/imgui/tree/master/docs

#define IMGUI_IMPL_OPENGL_LOADER_CUSTOM
#include "glad/glad.h"
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <stdio.h>
#include <iostream>
#include <cstring>
#include "image_ppm.h"
#include "ImGuiFileDialog.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#include "harmonization.h"
#include "ColorSpaces.h"
#include "segmentation.h"

typedef unsigned char uchar;

bool startedThread = false;
std::atomic<bool> threadDone = false;
std::jthread * renderThread = nullptr;
std::atomic<bool> regionsDone = false;

ImGuiTextBuffer buf;

static void glfw_error_callback(int error, const char* description)
{
	fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

HarmonyType textToHarmonyType(const std::string & type)
{
	if (type == "Complémentaire")
		return HarmonyType::COMPLEMENTARY;
	if (type == "Triadique")
		return HarmonyType::TRIADIC;
	if (type == "Tétradique carré")
		return HarmonyType::TETRADIC_SQUARE;
	if (type == "Tétradique rectangle")
		return HarmonyType::TETRADIC_RECTANGLE;
	if (type == "Complémentaire-split")
		return HarmonyType::SPLIT_COMPLEMENTARY;
	if (type == "Analogue")
		return HarmonyType::ANALOGOUS;
	return HarmonyType::COMPLEMENTARY;
}

//void stopThread()
//{
//	if (startedThread && renderThread->joinable())
//	{
//		renderThread->request_stop();
//		renderThread->join();
//		delete renderThread;
//		startedThread = false;
//	}
//}


void updatePreviewThreaded(const ImVec4 & clear_color, int w, int h, const unsigned char *imageIn, unsigned char *imageOut,
				   const char *const *items, int selectedHarmony, bool Kmean)
{
	HarmonyType type = textToHarmonyType({items[selectedHarmony]});
	ColorHSV hsv{};
	rgb_to_hsv({clear_color.x, clear_color.y, clear_color.z}, hsv);
	float hue = hsv.h;

	startedThread = true;
	threadDone = false;
	renderThread = new std::jthread([Kmean, type, h, w, imageIn, imageOut, hue]()
		{
			if (Kmean)
			{
				buf.append("Shifting with k-mean.\n");
				KMeanShift(type, h, w, imageIn, imageOut);
			}
			else
			{
				buf.append("Shifting with given color.\n");
				shiftColor(hue, type, h, w, imageIn, imageOut);
			}
			buf.append("Done.\n");
			threadDone = true;
		}
	);
}


void updatePreview(const ImVec4 & clear_color, int w, int h, const unsigned char *imageIn, unsigned char *imageOut,
						   const char *const *items, int selectedHarmony, GLuint textureID, bool Kmean)
{
	if (!startedThread)
	{
		updatePreviewThreaded(clear_color, w, h, imageIn, imageOut, items, selectedHarmony, Kmean);
	}
	else if (!threadDone)
	{
		delete renderThread;
		updatePreviewThreaded(clear_color, w, h, imageIn, imageOut, items, selectedHarmony, Kmean);
	}
}

void computeRegionsThreaded(int h, int w, const unsigned char * image, int * regionIds, Region ** regions, float thresholdAvg, float thresholdSize)
{
	startedThread = true;
	threadDone = false;
	renderThread = new std::jthread([h, w, image, regions, regionIds, thresholdSize, thresholdAvg]()
		{
			buf.clear();
			int size = h*w;
			int size3 = 3*size;
			// convert to ycbcr float
			float * imageYCbCr = new float[size3];
			rgb_to_ycbcr(image, imageYCbCr, h, w);

			// compute seeds
			buf.append("computing seeds\n");
			uchar * imageSeed = new uchar[size];
			std::memset(imageSeed, 0u, sizeof(uchar)*size);
			computeSeeds(imageYCbCr, imageSeed, h, w, 100, 0.02f);

			// create regions based solely on seeds
			buf.append("creating regions\n");
			int regionCount = createRegionIDFromSeeds(h, w, imageSeed, regionIds);

			// compute region size and average
			buf.append("computing region size and average\n");
			*regions = new Region[regionCount];
			computeRegionSizeAndAvg(*regions, regionCount, regionIds, imageYCbCr, h, w);

			// compute frontiers
			buf.append("computing frontier\n");
			std::multiset<PxDist, CmpPxDist> frontier;
			computeFrontier(h, w, imageYCbCr, regionIds, *regions, frontier);

			// grow regions
			buf.append("growing regions\n");
			growRegions(h, w, imageYCbCr, regionIds, *regions, frontier);

			// merge
			buf.append("merging regions\n");
			// compute neighbors of each regions
			std::set<int> * regionNeighbours = new std::set<int>[regionCount];
			computeRegionNeighbors(h, w, regionIds, regionNeighbours);
			int * regionsAssociations = new int[regionCount];
			std::memset(regionsAssociations, -1, sizeof(int)*regionCount);

			// merging by average
			buf.append("\t1. by average color\n");
			constexpr float THRESHOLD_AVG = 0.05f;
			int regionMin;
			mergeByAverage(regionCount, *regions, regionNeighbours, thresholdAvg, regionsAssociations);

			buf.append("\t2. by size\n");
			const int THRESHOLD_SIZE = (int)(size/1000.0f);
			mergeBySize(regionCount, *regions, regionNeighbours, regionsAssociations, regionMin, thresholdSize);

			// update regionsIds
			buf.append("updating region IDs\n");
			updateRegionsID(size, regionIds, regionsAssociations);

			delete [] imageYCbCr;
			delete [] imageSeed;
			delete [] regionNeighbours;
			delete [] regionsAssociations;
			buf.append("\nDone computing regions.\n");
			regionsDone = true;
			threadDone = true;
		}
	);
}


void computeRegions(int h, int w, const unsigned char * image, int * regionIds, Region ** regions, float thresholdAvg, float thresholdSize)
{
	if (!startedThread)
	{
		computeRegionsThreaded(h, w, image, regionIds, regions, thresholdAvg, thresholdSize);
	}
	else if (!threadDone)
	{
		delete renderThread;
		computeRegionsThreaded(h, w, image, regionIds, regions, thresholdAvg, thresholdSize);
	}
}

// Main code
int main(int argc, char** argv)
{
	glfwSetErrorCallback(glfw_error_callback);
	if (!glfwInit())
		return 1;

	// Decide GL+GLSL versions
#if defined(IMGUI_IMPL_OPENGL_ES2)
	// GL ES 2.0 + GLSL 100
    const char* glsl_version = "#version 100";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#elif defined(__APPLE__)
	// GL 3.2 + GLSL 150
    const char* glsl_version = "#version 150";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // Required on Mac
#else
	// GL 3.0 + GLSL 130
	const char* glsl_version = "#version 130";
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
	//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
#endif

	// Create window with graphics context
	GLFWwindow* window = glfwCreateWindow(1280, 720, "Harmonisation", nullptr, nullptr);
	if (window == nullptr)
		return 1;
	glfwMakeContextCurrent(window);
	glfwSwapInterval(1); // Enable vsync

	gladLoadGLLoader((GLADloadproc) glfwGetProcAddress);
	gladLoadGL();

	// Setup Dear ImGui context
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;
	io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
	io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();
	//ImGui::StyleColorsLight();

	// Setup Platform/Renderer backends
	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init(glsl_version);

	// Load Fonts
	// - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
	// - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
	// - If the file cannot be loaded, the function will return a nullptr. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
	// - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling ImFontAtlas::Build()/GetTexDataAsXXXX(), which ImGui_ImplXXXX_NewFrame below will call.
	// - Use '#define IMGUI_ENABLE_FREETYPE' in your imconfig file to use Freetype for higher quality font rendering.
	// - Read 'docs/FONTS.md' for more instructions and details.
	// - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
	// - Our Emscripten build process allows embedding fonts to be accessible at runtime from the "fonts/" folder. See Makefile.emscripten for details.
	//io.Fonts->AddFontDefault();
	//io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\segoeui.ttf", 18.0f);
	//io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
	//io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
	//io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
	//ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f, nullptr, io.Fonts->GetGlyphRangesJapanese());
	//IM_ASSERT(font != nullptr);

	// Our state
	bool show_demo_window = false;
	bool show_another_window = false;
	ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

	//image
	char text[256]= "";
	std::string imageInName;
	std::string imageOutName;
	std::string baseName;
	int w, h, size, size3;
	w = 512;
	h = 512;
	unsigned char* imageIn = nullptr;
	int * regionIDs = nullptr;
	Region * regions = nullptr;
	unsigned char* imageOut = nullptr;
	bool auto_harmo = false;
	const char* items[] = {"Complémentaire", "Triadique", "Tétradique carré", "Tétradique rectangle", "Complémentaire-split", "Analogue"};
	static int selected_item = 0;
	GLuint textureLeft = -1;
	GLuint textureRight = -1;
	int regionsCount = 0;
	float thresholdAvg = 0.1f;
	int thresholdSize = 500;
	ImVec2 imageDisplaySize(500, 500);
	bool computingRegions = false;


	// create texture
	glGenTextures(1, &textureLeft);
	glBindTexture(GL_TEXTURE_2D, textureLeft);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
	glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
	glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
	glGenTextures(1, &textureRight);
	glBindTexture(GL_TEXTURE_2D, textureRight);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
	glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
	glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);


	// Main loop
	while (!glfwWindowShouldClose(window))
	{
		if (startedThread && threadDone)
		{
			delete renderThread;
			startedThread = false;
			if (!computingRegions)
			{
				glBindTexture(GL_TEXTURE_2D, textureRight);
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, imageOut);
			}
			else
			{
				computingRegions = false;

				startedThread = true;
				threadDone = false;

				HarmonyType type = textToHarmonyType({items[selected_item]});
				ColorHSV hsv{};
				rgb_to_hsv({clear_color.x, clear_color.y, clear_color.z}, hsv);
				float hue = hsv.h;

				renderThread = new std::jthread(
					[auto_harmo, type, hue, h, w, imageIn, imageOut, regions, regionIDs]()
					{
						if (auto_harmo)
						{
							buf.append("Shifting with k-mean.\n");
							KMeanShiftRegions(type, h, w, imageIn, imageOut, regions, regionIDs);
						}
						else
						{
							buf.append("Shifting with given color.\n");
							shiftColorRegions(hue, type, h, w, imageIn, imageOut, regions, regionIDs);
						}
						buf.append("\nDone.\n");
						threadDone = true;
					}
				);
			}
		}

		// Poll and handle events (inputs, window resize, etc.)
		// You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
		// - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
		// - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
		// Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
		glfwPollEvents();
		ImVec2 mainWindowSize = ImGui::GetIO().DisplaySize;

		// Start the Dear ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		{
			ImGui::SetNextWindowSize(mainWindowSize);
			ImGui::SetNextWindowPos(ImVec2(0,0));
			ImGui::Begin("Harmonisation d'images", nullptr, ImGuiWindowFlags_NoResize);                          // Create a window called "Hello, world!" and append into it.

			if (ImGui::Button("Ouvrir..."))
			{
				ImGuiFileDialog::Instance()->OpenDialog("ChooseFileDlgKey", "Choisir une image", ".ppm", ".", 1, nullptr, ImGuiFileDialogFlags_Modal);
			}
			if (ImGuiFileDialog::Instance()->Display("ChooseFileDlgKey"))
			{
				// action if OK
				if (ImGuiFileDialog::Instance()->IsOk())
				{
					std::string filePathName = ImGuiFileDialog::Instance()->GetFilePathName();
					std::string filePath = ImGuiFileDialog::Instance()->GetCurrentPath();
					imageInName = filePathName.substr(filePath.size()+1);
//
//					baseName = {imageInName.substr(0, imageInName.size()-4)};
//					lire_nb_lignes_colonnes_image_ppm(filePathName.data(), &h, &w);
//					size = h*w;
//					size3 = size * 3;
//					lire_image_ppm(filePathName.data(), imageIn, size);

					if (imageIn != nullptr)
						free(imageIn);
					imageIn = stbi_load(filePathName.c_str(), &w, &h, nullptr, 3);
					size = h*w;
					size3 = h*w*3;
					glBindTexture(GL_TEXTURE_2D, textureLeft);
					glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, imageIn);

					imageDisplaySize = ImVec2(500, 500.0f*((float)h/(float)w));

					if (imageOut != nullptr)
						delete [] imageOut;
					imageOut = new unsigned char[size3];
					updatePreview(clear_color, w, h, imageIn, imageOut, items, selected_item, textureRight, auto_harmo);

					if (regionIDs != nullptr)
					{
						delete [] regionIDs;
						regions = nullptr;
					}
					if (regions != nullptr)
					{
						delete [] regions;
						regions = nullptr;
					}
					regionsCount = 0;
					regionsDone = false;
				}

				// close
				ImGuiFileDialog::Instance()->Close();
			}
//			if (ImGui::Button("Charger l'image")){
//				imageInName = {text};
//				baseName = {imageInName.substr(0, imageInName.size()-4)};
//				lire_nb_lignes_colonnes_image_ppm(imageInName.data(), &h, &w);
//				size = h*w;
//				size3 = size * 3;
//				imageIn = new unsigned char[size3];
//				lire_image_ppm(imageInName.data(), imageIn, size);
//			}
			ImGui::Image((void*)textureLeft, imageDisplaySize);
			ImGui::SameLine();
			ImGui::SetCursorPosY(300);
			ImGui::Text("=>");
			ImGui::SameLine();
			ImGui::Image((void*)textureRight, imageDisplaySize);
			if (ImGui::ColorEdit3("Couleur", (float*)&clear_color))
			{
				if (!auto_harmo)
				{
					updatePreview(clear_color, w, h, imageIn, imageOut, items, selected_item, textureRight, auto_harmo);
				}
			}
			ImGui::Text("Choix du template");
			ImGui::SetNextItemWidth(200);
			if (ImGui::ListBox("##", &selected_item, items, IM_ARRAYSIZE(items)))
			{
				if (!auto_harmo)
					updatePreview(clear_color, w, h, imageIn, imageOut, items, selected_item, textureRight, auto_harmo);
			}
			ImGui::SameLine();
			if (ImGui::Checkbox("Harmonisation automatique", &auto_harmo))
			{
				updatePreview(clear_color, w, h, imageIn, imageOut, items, selected_item, textureRight, auto_harmo);
			}

			if (ImGui::Button("Calculer l'image finale"))
			{
				if (!regionsDone)
				{
					regionIDs = new int[size];
					computingRegions = true;
					computeRegions(h, w, imageIn, regionIDs, &regions, thresholdAvg, size*1.0f/(float)thresholdSize);
				}
				else
				{
					computingRegions = false;

					startedThread = true;
					threadDone = false;

					HarmonyType type = textToHarmonyType({items[selected_item]});
					ColorHSV hsv{};
					rgb_to_hsv({clear_color.x, clear_color.y, clear_color.z}, hsv);
					float hue = hsv.h;

					renderThread = new std::jthread(
						[auto_harmo, type, hue, h, w, imageIn, imageOut, regions, regionIDs]()
						{
							if (auto_harmo)
							{
								buf.append("Shifting with k-mean.\n");
								KMeanShiftRegions(type, h, w, imageIn, imageOut, regions, regionIDs);
							}
							else
							{
								buf.append("Shifting with given color.\n");
								shiftColorRegions(hue, type, h, w, imageIn, imageOut, regions, regionIDs);
							}
							buf.append("\nDone.\n");
							threadDone = true;
						}
					);
				}
			}
			ImGui::SameLine();
			if (ImGui::Button("Enregistrer sous..."))
			{
				ImGuiFileDialog::Instance()->OpenDialog("SaveFileDlgKey", "Choisir un nom", ".png", ".", 1, nullptr, ImGuiFileDialogFlags_Modal);
//				ecrire_image_ppm(imageOutName.data(), imageOut, h, w);
			}

			if (ImGuiFileDialog::Instance()->Display("SaveFileDlgKey"))
			{
				// action if OK
				if (ImGuiFileDialog::Instance()->IsOk())
				{
					std::string res = ImGuiFileDialog::Instance()->GetFilePathName();

					stbi_write_png(res.c_str(), w, h, 3, (void*)imageOut, 0);
				}

				// close
				ImGuiFileDialog::Instance()->Close();
			}
			bool delRegions = false;
			if (ImGui::SliderFloat("Différence min. entre les régions", &thresholdAvg, 0.01f, 0.2f))
				delRegions = true;
			if (ImGui::SliderInt("Taille min. des régions", &thresholdSize, 2000, 50, "1/%d"))
				delRegions = true;
			if (delRegions)
			{
				if (regionIDs != nullptr)
				{
					delete [] regionIDs;
					regionIDs = nullptr;
				}
				if (regions != nullptr)
				{
					delete [] regions;
					regions = nullptr;
				}
				regionsCount = 0;
				regionsDone = false;
			}


			ImGui::BeginChild("##");
			{
				if (buf.size() > 0)
				{
					ImGui::TextUnformatted(buf.begin());
				}
				ImGui::SetScrollHereY(1.0f);
			}
			ImGui::EndChild();
			//ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
			ImGui::End();
		}

		// Rendering
		ImGui::Render();
		int display_w, display_h;
		glfwGetFramebufferSize(window, &display_w, &display_h);
		glViewport(0, 0, display_w, display_h);
		//glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
		glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
		glClear(GL_COLOR_BUFFER_BIT);
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		glfwSwapBuffers(window);
	}
#ifdef __EMSCRIPTEN__
	EMSCRIPTEN_MAINLOOP_END;
#endif

	// Cleanup
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	free(imageIn);
	delete [] imageOut;
	delete [] regionIDs;
	delete [] regions;

	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}
