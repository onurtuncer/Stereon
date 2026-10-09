# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Python bindings (nanobind). Only a version attribute so far: the rule for
# exposing PascalCase C++ names to Python is still open (ADR-0017).
nanobind_add_module(stereon_python StereonModule.cpp)
set_target_properties(stereon_python PROPERTIES OUTPUT_NAME stereon)
target_link_libraries(stereon_python PRIVATE stereon::core stereon_options)
# nanobind_add_module adds these as ordinary include directories; as system
# directories, Stereon's strict warnings do not fire inside them.
target_include_directories(stereon_python SYSTEM PRIVATE
    "${nanobind_SOURCE_DIR}/include" ${Python_INCLUDE_DIRS})

if(STEREON_BUILD_TESTS)
    add_test(NAME "python: import stereon"
        COMMAND Python::Interpreter -c
            "import stereon; assert stereon.__version__ == '${PROJECT_VERSION}', stereon.__version__")
    set_tests_properties("python: import stereon" PROPERTIES
        ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:stereon_python>")
endif()
