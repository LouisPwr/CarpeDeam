include_directories(commons)
add_subdirectory(assembler)
add_subdirectory(commons)
add_subdirectory(version)
add_subdirectory(workflow)
add_subdirectory(util)

add_executable(carpedeam
        ${commons_source_files}
        ${carpedeam_assembler_source_files}
        ${carpedeam_workflow_source_files}
        ${util_source_files}
        LocalCommandDeclarations.h
        carpedeam.cpp
        )
mmseqs_setup_derived_target(carpedeam)
target_link_libraries(carpedeam version kerasify flash)

add_dependencies(carpedeam local-generated)

install(TARGETS carpedeam DESTINATION bin)
