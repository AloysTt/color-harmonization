add_library(harmonization STATIC
	harmonization.h harmonization.cpp
	segmentation.h segmentation.cpp
	KMean.h KMean.cpp
	ColorSpaces.h ColorSpaces.cpp
	util.h
	image_ppm.h
)
target_include_directories(harmonization PUBLIC src)

add_executable(main
    main.cpp
)

add_executable(test_hsv
    test_hsv.cpp
	ColorSpaces.h ColorSpaces.cpp
)

add_executable(test_ycbcr
	test_ycbcr.cpp
	ColorSpaces.h ColorSpaces.cpp
)

add_executable(test_sector_borders
	test_sector_borders.cpp
    harmonization.h harmonization.cpp
	KMean.h KMean.cpp
	ColorSpaces.h ColorSpaces.cpp
    image_ppm.h
)
add_custom_command(TARGET test_sector_borders PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/images/ ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(test_select_color
	test_select_color.cpp
	harmonization.h harmonization.cpp
	KMean.h KMean.cpp
	ColorSpaces.h ColorSpaces.cpp
	image_ppm.h
)
add_custom_command(TARGET test_select_color PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
	${CMAKE_SOURCE_DIR}/images/ ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(test_shift
	test_shift.cpp
	harmonization.h harmonization.cpp
	KMean.h KMean.cpp
	ColorSpaces.h ColorSpaces.cpp
	util.h
	image_ppm.h
)
add_custom_command(TARGET test_shift PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
	${CMAKE_SOURCE_DIR}/images/ ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(histogram_h
	histogram_h.cpp
	ColorSpaces.h ColorSpaces.cpp
	image_ppm.h
)
add_custom_command(TARGET histogram_h PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
	${CMAKE_SOURCE_DIR}/images/ ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(test_kmean
	test_kmean.cpp
	harmonization.h harmonization.cpp
	KMean.h KMean.cpp
	ColorSpaces.h ColorSpaces.cpp
	util.h
	image_ppm.h
	)
add_custom_command(TARGET test_kmean PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
	${CMAKE_SOURCE_DIR}/images/ ${CMAKE_CURRENT_BINARY_DIR}
	)


add_executable(test_genSeeds
	test_genSeeds.cpp
	harmonization.h harmonization.cpp
	segmentation.h segmentation.cpp
	KMean.h KMean.cpp
	ColorSpaces.h ColorSpaces.cpp
	util.h
	image_ppm.h
)
add_custom_command(TARGET test_genSeeds PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
	${CMAKE_SOURCE_DIR}/images/ ${CMAKE_CURRENT_BINARY_DIR}
)


add_executable(test_createRegions
	test_createRegions.cpp
	harmonization.h harmonization.cpp
	segmentation.h segmentation.cpp
	KMean.h KMean.cpp
	ColorSpaces.h ColorSpaces.cpp
	util.h
	image_ppm.h
)
add_custom_command(TARGET test_createRegions PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
	${CMAKE_SOURCE_DIR}/images/ ${CMAKE_CURRENT_BINARY_DIR}
)


add_executable(test_mergeRegions
	test_mergeRegions.cpp
	harmonization.h harmonization.cpp
	segmentation.h segmentation.cpp
	KMean.h KMean.cpp
	ColorSpaces.h ColorSpaces.cpp
	util.h
	image_ppm.h
)
add_custom_command(TARGET test_mergeRegions PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
	${CMAKE_SOURCE_DIR}/images/ ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(test_regionsHarmonization
	test_regionsHarmonization.cpp
)
target_link_libraries(test_regionsHarmonization PUBLIC harmonization)
add_custom_command(TARGET test_regionsHarmonization PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
	${CMAKE_SOURCE_DIR}/images/ ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(test_imgui test_imgui.cpp)
target_link_libraries(test_imgui PRIVATE imgui imguifiledialog glfw glad harmonization)
