# SPDX-FileCopyrightText: 2026 Onur Tuncer and Stereon contributors
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

include(Catch)

# One test executable per library: tests/<library>/*.cpp -> stereon_<library>_tests.
# Test cases are tagged with their library, e.g. [core], plus topic tags.
# Property tests use RapidCheck through rc::prop: include "RapidCheck.hpp"
# from tests/support/.
function(stereon_add_tests library)
    file(GLOB sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${library}/*.cpp")
    add_executable(stereon_${library}_tests ${sources})
    target_link_libraries(stereon_${library}_tests
        PRIVATE stereon::${library} stereon_options Catch2::Catch2WithMain stereon_rapidcheck_catch)
    target_include_directories(stereon_${library}_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/support")
    stereon_anchor_contracts(stereon_${library}_tests)
    target_compile_definitions(stereon_${library}_tests PRIVATE
        STEREON_EXPECTED_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
        STEREON_EXPECTED_VERSION_MINOR=${PROJECT_VERSION_MINOR}
        STEREON_EXPECTED_VERSION_PATCH=${PROJECT_VERSION_PATCH})
    catch_discover_tests(stereon_${library}_tests
        TEST_PREFIX "${library}: "
        DISCOVERY_MODE PRE_TEST)
endfunction()

stereon_add_tests(core)
stereon_exec_sources(core/ExecTests.cpp)
stereon_add_tests(robust)

# Do contracts written in Stereon.Core fire in a unit that imports it? Each
# probe violates one and must print the violation. GCC 16.2 checks the module's
# inline functions only with the contract anchor above, and drops the checks in
# its templates altogether; that test is marked as a known failure and turns
# red when a compiler release fixes it, which is the signal to remove
# stereon_anchor_contracts(). The handle check is reported from a module
# implementation unit and works regardless. Only the reference compiler has
# contracts, and an "ignore" build (the bench preset) checks nothing.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT STEREON_CONTRACT_SEMANTIC STREQUAL "ignore")
    add_executable(stereon_contract_probe probes/ContractProbe.cpp)
    target_link_libraries(stereon_contract_probe PRIVATE stereon::core stereon_options)
    stereon_anchor_contracts(stereon_contract_probe)
    if(STEREON_HOT_CONTRACTS)
        add_test(NAME "core: contracts in module inline functions fire in importers"
            COMMAND stereon_contract_probe inline)
        add_test(NAME "core: contracts in module templates fire in importers (GCC 16.2 drops them)"
            COMMAND stereon_contract_probe template)
        set_tests_properties("core: contracts in module inline functions fire in importers"
            PROPERTIES PASS_REGULAR_EXPRESSION "contract violation")
        set_tests_properties("core: contracts in module templates fire in importers (GCC 16.2 drops them)"
            PROPERTIES PASS_REGULAR_EXPRESSION "contract violation" WILL_FAIL TRUE)
    endif()
    if(STEREON_CHECKED_HANDLES)
        add_test(NAME "core: checked handles report a stale handle"
            COMMAND stereon_contract_probe handle)
        set_tests_properties("core: checked handles report a stale handle"
            PROPERTIES PASS_REGULAR_EXPRESSION "invalid handle Vertex#0@0.*slot released.*contract violation")
    endif()
endif()

# Source checks that need only Python. CI runs the same scripts in its lint job.
find_package(Python3 3.9 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME "tools: library layering"
        COMMAND Python3::Interpreter "${PROJECT_SOURCE_DIR}/tools/check_layering.py"
                --root "${PROJECT_SOURCE_DIR}")
    add_test(NAME "tools: epsilon lint self-test"
        COMMAND Python3::Interpreter "${PROJECT_SOURCE_DIR}/tools/check_epsilon.py" --self-test)
    add_test(NAME "tools: epsilon lint"
        COMMAND Python3::Interpreter "${PROJECT_SOURCE_DIR}/tools/check_epsilon.py"
                --root "${PROJECT_SOURCE_DIR}")
    # The complexity check needs the lizard package; without it the script
    # exits 77, which CTest reports as skipped rather than failed.
    add_test(NAME "tools: complexity self-test"
        COMMAND Python3::Interpreter "${PROJECT_SOURCE_DIR}/tools/check_complexity.py" --self-test)
    add_test(NAME "tools: complexity"
        COMMAND Python3::Interpreter "${PROJECT_SOURCE_DIR}/tools/check_complexity.py"
                --root "${PROJECT_SOURCE_DIR}")
    set_tests_properties("tools: complexity self-test" "tools: complexity"
        PROPERTIES SKIP_RETURN_CODE 77)
endif()
