 
if(NOT EXISTS "/vol/cloud/louis/apps/CarpeDeam15.3/build/install_manifest.txt")
  message(FATAL_ERROR "Cannot find install manifest: /vol/cloud/louis/apps/CarpeDeam15.3/build/install_manifest.txt")
endif(NOT EXISTS "/vol/cloud/louis/apps/CarpeDeam15.3/build/install_manifest.txt")

file(READ "/vol/cloud/louis/apps/CarpeDeam15.3/build/install_manifest.txt" files)
string(REGEX REPLACE "\n" ";" files "${files}")
foreach(file ${files})
  message(STATUS "Uninstalling $ENV{DESTDIR}${file}")
  if(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
    exec_program(
      "/vol/cloud/louis/miniconda3/envs/all/bin/cmake" ARGS "-E remove \"$ENV{DESTDIR}${file}\""
      OUTPUT_VARIABLE rm_out
      RETURN_VALUE rm_retval
      )
    if(NOT "${rm_retval}" STREQUAL 0)
      message(FATAL_ERROR "Problem when removing $ENV{DESTDIR}${file}")
    endif(NOT "${rm_retval}" STREQUAL 0)
  else(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
    message(STATUS "File $ENV{DESTDIR}${file} does not exist.")
  endif(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
endforeach(file)
