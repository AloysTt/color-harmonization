add_library(imgui STATIC
	imgui.cpp include/imgui.h include/imconfig.h
	imgui_draw.cpp
	imgui_tables.cpp
	imgui_widgets.cpp
	imgui_demo.cpp
	include/imgui_internal.h
	include/imstb_rectpack.h
	include/imstb_textedit.h
	include/imstb_truetype.h

	imgui_impl_glfw.cpp
	imgui_impl_opengl3.cpp
	include/imgui_impl_glfw.h
	include/imgui_impl_opengl3.h
	include/imgui_impl_opengl3_loader.h
)
target_include_directories(imgui PUBLIC include)
target_link_libraries(imgui PUBLIC glfw)
