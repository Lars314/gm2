# build libraries
art_make( LIB_LIBRARIES gm2util_blinders_randomlib )

# copy headers
install_headers()

# add subdirectory
add_subdirectory( randomlib )
add_subdirectory( test )
