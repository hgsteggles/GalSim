cmake_minimum_required(VERSION 2.6)

set(GALSIM_SRCS
       ${CMAKE_CURRENT_SOURCE_DIR}/Disks2.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/func.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/Galaxy.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/HosTracks.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/KroupaIMF.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/Logger.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/MonteCarloCumulative.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/ProgressBar.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/Random.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/Recipes.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/SpiralArms2.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/Star.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/HII.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/Timer.cpp
       ${CMAKE_CURRENT_SOURCE_DIR}/FileManagement.cpp)

include_directories("${GALSIM_SOURCE_DIR}/lib")
include_directories("${GALSIM_SOURCE_DIR}/lib/lua-5.2.3")
include_directories("${GALSIM_SOURCE_DIR}/src")

add_executable(galsim ${GALSIM_SRCS})
#target_link_libraries(galsim ${LUA_LIBRARIES})
target_link_libraries(galsim ${GALSIM_SOURCE_DIR}/lib/liblua.a dl)
