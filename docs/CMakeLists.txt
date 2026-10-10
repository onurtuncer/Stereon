# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# The Stereon manual: Doxygen XML for the API, rendered with the rest of the
# .rst pages by Sphinx and Breathe. Enabled with -DSTEREON_BUILD_DOCS=ON.
# Python requirements: docs/requirements.txt.
#
#   cmake --build --preset <preset> --target Stereon_Docs

find_package(Doxygen)
find_program(SPHINX_EXECUTABLE NAMES sphinx-build DOC "Sphinx documentation generator")

if(NOT SPHINX_EXECUTABLE)
    message(WARNING "sphinx-build not found; the documentation target is not available. "
                    "Install it with: pip install -r docs/requirements.txt")
    return()
endif()

add_subdirectory(tikz)

set(STEREON_DOCS_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(STEREON_DOCS_HTML ${CMAKE_CURRENT_BINARY_DIR}/html)
set(STEREON_DOCS_STAMP ${CMAKE_CURRENT_BINARY_DIR}/sphinx.stamp)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/xml/index.xml)

file(GLOB_RECURSE STEREON_DOC_PAGES CONFIGURE_DEPENDS
    ${STEREON_DOCS_SOURCE}/*.rst
    ${STEREON_DOCS_SOURCE}/references.bib)

set(_sphinx_depends ${STEREON_DOC_PAGES} ${STEREON_DOCS_SOURCE}/conf.py)

if(DOXYGEN_FOUND)
    file(GLOB_RECURSE STEREON_DOXYGEN_INPUTS CONFIGURE_DEPENDS
        ${PROJECT_SOURCE_DIR}/src/*.hpp
        ${PROJECT_SOURCE_DIR}/src/*.h
        ${PROJECT_SOURCE_DIR}/src/*.cppm
        ${PROJECT_SOURCE_DIR}/src/*.ixx)

    set(DOXYFILE_IN ${STEREON_DOCS_SOURCE}/Doxyfile.in)
    set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)
    file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

    add_custom_command(
        OUTPUT ${DOXYGEN_INDEX_FILE}
        COMMAND Doxygen::doxygen ${DOXYFILE_OUT}
        MAIN_DEPENDENCY ${DOXYFILE_OUT}
        DEPENDS ${STEREON_DOXYGEN_INPUTS} ${DOXYFILE_IN} ${PROJECT_SOURCE_DIR}/README.md
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating Stereon API XML with Doxygen"
        VERBATIM)
    add_custom_target(Stereon_Docs_Doxygen DEPENDS ${DOXYGEN_INDEX_FILE})
    list(APPEND _sphinx_depends ${DOXYGEN_INDEX_FILE})
else()
    message(WARNING "Doxygen not found; the manual will be built without the API reference")
endif()

add_custom_command(
    OUTPUT ${STEREON_DOCS_STAMP}
    COMMAND ${CMAKE_COMMAND} -E env STEREON_DOXYGEN_XML=${DOXYGEN_OUTPUT_DIR}/xml
            ${SPHINX_EXECUTABLE} -b html -q
            -d ${CMAKE_CURRENT_BINARY_DIR}/doctrees
            ${STEREON_DOCS_SOURCE} ${STEREON_DOCS_HTML}
    COMMAND ${CMAKE_COMMAND} -E touch ${STEREON_DOCS_STAMP}
    DEPENDS ${_sphinx_depends}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Building the Stereon manual with Sphinx"
    VERBATIM)

add_custom_target(Stereon_Docs DEPENDS ${STEREON_DOCS_STAMP})
add_dependencies(Stereon_Docs stereon_tikz_figures)

include(GNUInstallDirs)
install(DIRECTORY ${STEREON_DOCS_HTML}/ DESTINATION ${CMAKE_INSTALL_DOCDIR}/html OPTIONAL)
