set(CET_EXPORT EXPORT)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS WERROR NO_UNDEFINED EXTRA_FLAGS -pedantic)
cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

find_package(Torch QUIET ${CET_EXPORT})
if (Torch_FOUND)
  set(PANDORA_LIBTORCH_DEF ON)
else()
  set(PANDORA_LIBTORCH_DEF OFF)
endif()
option(PANDORA_LIBTORCH "Flag for building against LibTorch" ${PANDORA_LIBTORCH_DEF})

find_package(PandoraSDK 03.04.00 REQUIRED ${CET_EXPORT})
option(PANDORA_MONITORING "Enable Pandora Monitoring" TRUE)
if (PANDORA_MONITORING)
  find_package(PandoraMonitoring 03.05.00 REQUIRED ${CET_EXPORT})
endif()

add_subdirectory(larpandoracontent)
if (PANDORA_LIBTORCH)
  add_subdirectory(larpandoradlcontent)
endif()

cet_cmake_config()
