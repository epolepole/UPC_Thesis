file(GLOB src_files *.c)
file(GLOB head_files ${HEAD_DIR} "*.h")
add_executable(LBMSolver ${src_files} ${head_files} )
target_link_libraries(LBMSolver m)