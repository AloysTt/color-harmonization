add_executable(main
    main.cpp
)

add_executable(test_hsv
    test_hsv.cpp
    HSV.h HSV.cpp
)

add_executable(test_sector_borders
	test_sector_borders.cpp
    harmonization.h harmonization.cpp
    HSV.h HSV.cpp
    image_ppm.h
)
add_custom_command(TARGET test_sector_borders PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/images/ ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(test_select_color
	test_select_color.cpp
	harmonization.h harmonization.cpp
	HSV.h HSV.cpp
	image_ppm.h
)
add_custom_command(TARGET test_select_color PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
	${CMAKE_SOURCE_DIR}/images/ ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(histogram_h
	histogram_h.cpp
	HSV.h HSV.cpp
	image_ppm.h
)
add_custom_command(TARGET histogram_h PRE_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
	${CMAKE_SOURCE_DIR}/images/ ${CMAKE_CURRENT_BINARY_DIR}
)