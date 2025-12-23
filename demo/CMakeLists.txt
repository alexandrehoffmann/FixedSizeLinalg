find_package(PkgConfig REQUIRED)

add_executable(demo_Vector demo_Vector.cpp)
add_executable(demo_Matrix demo_Matrix.cpp)
add_executable(demo_Tensor demo_Tensor.cpp)

target_include_directories(demo_Vector PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(demo_Matrix PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(demo_Tensor PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(demo_Vector PRIVATE FSLinalg)
target_link_libraries(demo_Matrix PRIVATE FSLinalg)
target_link_libraries(demo_Tensor PRIVATE FSLinalg)

