add_library(imguifiledialog STATIC
	ImGuiFileDialog.h
	ImGuiFileDialogConfig.h
	ImGuiFileDialog.cpp
	dirent/dirent.h
	stb/stb_image.h
	stb/stb_image_resize.h
	stb/stb_image_write.h
)
target_include_directories(imguifiledialog PUBLIC .)
target_link_libraries(imguifiledialog PRIVATE imgui)
