add_executable(main
    main.cpp
)

add_executable(test_hsv
    test_hsv.cpp
    HSV.h HSV.cpp
)

add_executable(test_harmonization
    test_harmonization.cpp
    harmonization.h harmonization.cpp
    HSV.h HSV.cpp
    image_ppm.h
)
add_custom_command(TARGET test_harmonization PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/images/ ${CMAKE_CURRENT_BINARY_DIR}
)
