
find_ups_product( libtorch )
find_package(Torch ${CET_PUBLIC} REQUIRED)

include_directories(SYSTEM ${TORCH_INCLUDE_DIRS})

cet_set_compiler_flags(DIAGS CAUTIOUS
      WERROR
      NO_UNDEFINED
      EXTRA_FLAGS
      -pedantic
      ${TORCH_CXX_FLAGS})

set( subdir_list LArControlFlow
                 LArHelpers
                 LArVertex
                 LArSlicing
                 LArSliceTagging
)

cet_make( LIBRARY_NAME LArPandoraDLContent
    SUBDIRS ${subdir_list}
    LIBRARIES
        PANDORASDK
        PANDORAMONITORING
        ${LAR_CONTENT_LIBRARY_NAME}
        torch
        c10
        caffe2
)

install_source( SUBDIRS ${subdir_list} )
install_headers( SUBDIRS ${subdir_list} )
