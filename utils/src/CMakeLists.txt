add_executable(leg_model leg_model.cpp)
target_link_libraries(leg_model ${_CORE_LIBRARIES})
message(STATUS "Core libraries: ${_CORE_LIBRARIES}")
