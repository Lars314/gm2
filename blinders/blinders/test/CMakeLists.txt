# build libs
art_make( 
	MODULE_LIBRARIES
		gm2util_blinders
		gm2util_blinders_randomlib
)

# install all *.fcl files in this directory
file(GLOB fcl_files *.fcl)
install( FILES ${fcl_files}
         DESTINATION ${product}/${version}/fcl )

# Also install in the build area
foreach(aFile ${fcl_files})
  get_filename_component( basename ${aFile} NAME )
  configure_file( ${aFile} ${CMAKE_BINARY_DIR}/${product}/fcl/${basename} COPYONLY )
endforeach(aFile)

### testing routines ###

cet_enable_asserts() # Enable asserts
include(CetTest) # include CET test macro(s)
cet_test_env() # Initialize the testing environment

# add test(s)
cet_test(notblinded HANDBUILT
	TEST_EXEC gm2
	TEST_ARGS --rethrow-all -c notBlinded.fcl --no-output
	DATAFILES notBlinded.fcl
)


