set(CET_EXPORT EXPORT)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(
  DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic
)
cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

find_package(PandoraSDK 03.04.00 REQUIRED ${CET_EXPORT})

option(PANDORA_MONITORING "Enable Pandora Monitoring" TRUE)
if(PANDORA_MONITORING)
  find_package(PandoraMonitoring 03.05.00 REQUIRED ${CET_EXPORT})
endif()

find_package(Eigen3 3.3 REQUIRED ${CET_EXPORT})

find_package(Torch QUIET ${CET_EXPORT})
if(Torch_FOUND)
  set(PANDORA_LIBTORCH ON)
  message(STATUS "LibTorch found — building DL content (LArPandoraDLContent)")
else()
  set(PANDORA_LIBTORCH OFF)
  message(STATUS "LibTorch not found — skipping DL content build")
endif()

include(${PANDORA_PROJECT_ROOT}/cmake/LArContent_sources.cmake)

if(NOT DEFINED LAR_CONTENT_SRCS OR LAR_CONTENT_SRCS STREQUAL "")
  message(FATAL_ERROR "LAR_CONTENT_SRCS not defined or empty — check LArContent_sources.cmake")
endif()

cet_make_library(
  LIBRARY_NAME LArPandoraContent
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  SOURCE ${LAR_CONTENT_SRCS}
  LIBRARIES
    PUBLIC
      PandoraPFA::PandoraSDK
      $<$<BOOL:${PANDORA_MONITORING}>:PandoraPFA::PandoraMonitoring>
    PRIVATE
      Eigen3::Eigen
)

if(PANDORA_MONITORING)
  target_compile_definitions(LArPandoraContent PUBLIC MONITORING)
endif()

if(PANDORA_LIBTORCH)
  find_package(TBB REQUIRED EXPORT)

  include(${PANDORA_PROJECT_ROOT}/cmake/LArDLContent_sources.cmake)

  if(NOT DEFINED LAR_DL_CONTENT_SRCS OR LAR_DL_CONTENT_SRCS STREQUAL "")
    message(FATAL_ERROR "LAR_DL_CONTENT_SRCS not defined or empty — check LArDLContent_sources.cmake")
  endif()

  cet_make_library(
    LIBRARY_NAME LArPandoraDLContent
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    SOURCE ${LAR_DL_CONTENT_SRCS}
    LIBRARIES
      PUBLIC
        LArPandoraContent
        PandoraPFA::PandoraSDK
        torch
  )

  target_compile_definitions(LArPandoraDLContent PUBLIC PANDORA_LIBTORCH)
endif()

install_source(SUBDIRS larpandoracontent)
install_headers(SUBDIRS larpandoracontent)

cet_cmake_config()
