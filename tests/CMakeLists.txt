find_package(PkgConfig REQUIRED)

set(FSLinalg_tests_SRC
	test_aliasing.cpp 
	test_lazy.cpp 
	tests_fslinalg.cpp
	test_chain.cpp
	test_tensor.cpp)

add_executable(tests_fslinalg ${FSLinalg_tests_SRC})

target_include_directories(tests_fslinalg PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(tests_fslinalg PRIVATE FSLinalg gtest)

add_test(NAME AllTests COMMAND tests_fslinalg)

