# Captain! If catch2 alpine files do not exist, fetch content
if( NOT EXISTS "usr/include/catch2/" OR 
    NOT EXISTS "usr/lib/libCatch2Main.a" OR
    NOT EXISTS "usr/lib/libCatch2.a" )

  message( "Ahoy! Not found in system locations - installing from source..." )

  include(FetchContent)

  FetchContent_Declare(
    Catch2
      GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.0.1 ) 

  FetchContent_MakeAvailable(Catch2)

else()

find_package( Catch2 3 )

endif()
