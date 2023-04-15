add_library(glad STATIC glad.c include/glad/glad.h include/KHR/khrplatform.h)
target_include_directories(glad PUBLIC include)
