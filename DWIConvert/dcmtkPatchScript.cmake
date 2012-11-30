
# file(READ ${dcmtkGenConfig} code)
# string(REPLACE "int i = strerror_r(0, buf, 100);" 
# "buf = strerror_r(0, buf, 100)" code "${code}")
# file(WRITE ${dcmtkGenConfig} "${code}")


file(READ ${dcmtk3rdParty} code)
string(REPLACE "# turn off library if it could not be found
"
"# turn off library if it could not be found
    find_package(JPEG REQUIRED)
    list(APPEND LIBTIFF_LIBS \${JPEG_LIBRARY})
" code "${code}")
file(WRITE ${dcmtk3rdParty} "${code}")
