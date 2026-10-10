# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Placeholder for manual figures drawn in TikZ. Each docs/tikz/<name>.tex is
# a standalone LaTeX document; it is compiled with pdflatex into
# <build>/docs/tikz/<name>.pdf. No figures exist yet; the target is a no-op
# until the first .tex file is added.

file(GLOB STEREON_TIKZ_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/*.tex)
find_program(PDFLATEX_EXECUTABLE pdflatex)

if(NOT STEREON_TIKZ_SOURCES OR NOT PDFLATEX_EXECUTABLE)
    if(STEREON_TIKZ_SOURCES)
        message(WARNING "pdflatex not found; TikZ figures will not be built")
    endif()
    add_custom_target(stereon_tikz_figures)
    return()
endif()

set(_stereon_tikz_pdfs "")
foreach(tex IN LISTS STEREON_TIKZ_SOURCES)
    get_filename_component(fig_name ${tex} NAME_WE)
    set(pdf ${CMAKE_CURRENT_BINARY_DIR}/${fig_name}.pdf)
    add_custom_command(
        OUTPUT ${pdf}
        COMMAND ${PDFLATEX_EXECUTABLE} -interaction=nonstopmode -halt-on-error
                -output-directory ${CMAKE_CURRENT_BINARY_DIR} ${tex}
        DEPENDS ${tex}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Building TikZ figure ${fig_name}.pdf"
        VERBATIM)
    list(APPEND _stereon_tikz_pdfs ${pdf})
endforeach()

add_custom_target(stereon_tikz_figures DEPENDS ${_stereon_tikz_pdfs})
