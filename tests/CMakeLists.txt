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
stereon_add_tests(curve)
stereon_add_tests(surface)
stereon_add_tests(intersect)
stereon_add_tests(topo)
stereon_add_tests(io)

# The reflection serialiser is instantiated in the test unit, which pulls the
# module's <meta> declarations into it. GCC 16.2 then rejects the built-in
# std::contracts::contract_violation that the contract anchor creates
# ("conflicting declaration"; docs/toolchain/gcc16-feature-probe.md), so these
# tests are built without the anchor. They call no inline module function
# that carries a contract, which is what the anchor exists for.
add_executable(stereon_serialize_tests serialize/SerializeTests.cpp)
target_link_libraries(stereon_serialize_tests
    PRIVATE stereon::surface stereon_options Catch2::Catch2WithMain stereon_rapidcheck_catch)
target_include_directories(stereon_serialize_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/support")
catch_discover_tests(stereon_serialize_tests TEST_PREFIX "serialize: " DISCOVERY_MODE PRE_TEST)

# The data behind the manual's fitting figures (docs/figures/): built and
# run with the tests so that it keeps compiling and the figures can always be
# regenerated with docs/figures/make_fitting_figures.py.
add_executable(stereon_figure_data ${PROJECT_SOURCE_DIR}/docs/figures/FigureData.cpp)
target_link_libraries(stereon_figure_data PRIVATE stereon::surface stereon_options)
stereon_anchor_contracts(stereon_figure_data)
add_test(NAME "docs: the figure data generator runs" COMMAND stereon_figure_data)
set_tests_properties("docs: the figure data generator runs"
    PROPERTIES PASS_REGULAR_EXPRESSION "\"surface\":")

# Differential tests against OpenCASCADE (ADR-0015's independent referee).
# OCCT is LGPL-2.1: it is linked into this test executable only, which is
# built on request and never distributed, and no Stereon library depends on
# it (ADR-0008). tools/occt/build_occt.sh builds the pinned version; pass its
# install directory in CMAKE_PREFIX_PATH.
if(STEREON_WITH_OCCT)
    # No version in find_package: OCCT's version file accepts only an exact
    # match, so "8.0" would reject 8.0.1. Any 8.0.x will do.
    find_package(OpenCASCADE CONFIG REQUIRED)
    if(NOT OpenCASCADE_VERSION VERSION_GREATER_EQUAL 8.0 OR OpenCASCADE_VERSION VERSION_GREATER_EQUAL 8.1)
        message(FATAL_ERROR "the OCCT differential tests need OpenCASCADE 8.0.x, found "
                            "${OpenCASCADE_VERSION}; see tools/occt/build_occt.sh")
    endif()
    add_executable(stereon_occt_tests occt/GeometryDifferentialTests.cpp
                                      occt/OcctRobustnessTests.cpp occt/RobustnessHunt.cpp)
    target_link_libraries(stereon_occt_tests
        PRIVATE stereon::surface stereon_options Catch2::Catch2WithMain
                TKernel TKMath TKG3d TKGeomBase TKGeomAlgo)
    target_include_directories(stereon_occt_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/support")
    stereon_anchor_contracts(stereon_occt_tests)
    catch_discover_tests(stereon_occt_tests TEST_PREFIX "occt: " DISCOVERY_MODE PRE_TEST)
endif()

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
    target_link_libraries(stereon_contract_probe PRIVATE stereon::surface stereon_options)
    stereon_anchor_contracts(stereon_contract_probe)
    if(STEREON_HOT_CONTRACTS)
        # The geometry evaluators check their preconditions with plain branches
        # that report from EvaluatorChecks.cpp, so they fire in importers even
        # though the evaluators are templates (hot-path checks: debug presets).
        add_test(NAME "curve: a non-finite parameter is reported"
            COMMAND stereon_contract_probe evaluator-nan)
        add_test(NAME "curve: a parameter outside a bounded domain is reported"
            COMMAND stereon_contract_probe evaluator-domain)
        set_tests_properties("curve: a non-finite parameter is reported"
            PROPERTIES PASS_REGULAR_EXPRESSION "Eval.Curve. called with a non-finite parameter.*contract violation")
        set_tests_properties("curve: a parameter outside a bounded domain is reported"
            PROPERTIES PASS_REGULAR_EXPRESSION "Eval.Curve. called with a parameter outside the domain.*contract violation")
        add_test(NAME "core: contracts in module inline functions fire in importers"
            COMMAND stereon_contract_probe inline)
        add_test(NAME "core: contracts in module templates fire in importers (GCC 16.2 drops them)"
            COMMAND stereon_contract_probe template)
        set_tests_properties("core: contracts in module inline functions fire in importers"
            PROPERTIES PASS_REGULAR_EXPRESSION "contract violation")
        set_tests_properties("core: contracts in module templates fire in importers (GCC 16.2 drops them)"
            PROPERTIES PASS_REGULAR_EXPRESSION "contract violation" WILL_FAIL TRUE)
    endif()
    # Like the arena check, the predicate domain is reported from a module
    # implementation unit, so it fires in every build that evaluates contracts.
    add_test(NAME "robust: coordinates outside the predicate domain are reported"
        COMMAND stereon_contract_probe predicate)
    set_tests_properties("robust: coordinates outside the predicate domain are reported"
        PROPERTIES PASS_REGULAR_EXPRESSION "Orient2D outside the predicate domain.*contract violation")
    # The arena's borrow check is a plain branch that reports from Arena.cpp,
    # so it fires in every build that evaluates contracts.
    add_test(NAME "core: copying an arena inside its edit callback is reported"
        COMMAND stereon_contract_probe arena)
    set_tests_properties("core: copying an arena inside its edit callback is reported"
        PROPERTIES PASS_REGULAR_EXPRESSION "PagedArena copy during an edit callback.*contract violation")
    add_test(NAME "core: taking an arena page inside its edit callback is reported"
        COMMAND stereon_contract_probe arena-page)
    set_tests_properties("core: taking an arena page inside its edit callback is reported"
        PROPERTIES PASS_REGULAR_EXPRESSION "PagedArena GetPage during an edit callback.*contract violation")
    if(STEREON_CHECKED_HANDLES)
        add_test(NAME "core: checked handles report a stale handle"
            COMMAND stereon_contract_probe handle)
        set_tests_properties("core: checked handles report a stale handle"
            PROPERTIES PASS_REGULAR_EXPRESSION "invalid handle Vertex#0@0.*slot released.*contract violation")
    endif()
endif()

# The debugger pretty-printers (tools/debugger/): a target program with one
# value of every printed type, inspected by scripted GDB and LLDB sessions
# that must print the same summaries. Only unoptimised builds keep every
# variable, so the checks run in Debug configurations; macOS is left out until
# its LLDB has been verified on GCC's debug information, and MSVC builds
# because GDB and LLDB cannot read their PDB debug information (the checks
# failed there, LLDB after a 300 s timeout). A missing debugger makes the
# check exit 77, which CTest reports as skipped.
add_executable(stereon_debugger_target debugger/DebuggerTarget.cpp)
target_link_libraries(stereon_debugger_target PRIVATE stereon::surface stereon_options)
stereon_anchor_contracts(stereon_debugger_target)
find_package(Python3 3.9 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND AND CMAKE_BUILD_TYPE STREQUAL "Debug" AND NOT APPLE AND NOT MSVC)
    foreach(debugger gdb lldb)
        add_test(NAME "tools: ${debugger} pretty-printers"
            COMMAND Python3::Interpreter "${CMAKE_CURRENT_SOURCE_DIR}/debugger/check_printers.py"
                    --debugger ${debugger}
                    --target $<TARGET_FILE:stereon_debugger_target>
                    --printers "${PROJECT_SOURCE_DIR}/tools/debugger/stereon_${debugger}.py")
        set_tests_properties("tools: ${debugger} pretty-printers" PROPERTIES SKIP_RETURN_CODE 77)
    endforeach()
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
    add_test(NAME "tools: literal lint self-test"
        COMMAND Python3::Interpreter "${PROJECT_SOURCE_DIR}/tools/check_literals.py" --self-test)
    add_test(NAME "tools: literal lint"
        COMMAND Python3::Interpreter "${PROJECT_SOURCE_DIR}/tools/check_literals.py"
                --root "${PROJECT_SOURCE_DIR}")
    # The complexity check needs the lizard package; without it the script
    # exits 77, which CTest reports as skipped rather than failed.
    add_test(NAME "tools: complexity self-test"
        COMMAND Python3::Interpreter "${PROJECT_SOURCE_DIR}/tools/check_complexity.py" --self-test)
    add_test(NAME "tools: complexity"
        COMMAND Python3::Interpreter "${PROJECT_SOURCE_DIR}/tools/check_complexity.py"
                --root "${PROJECT_SOURCE_DIR}")
    add_test(NAME "tools: benchmark report self-test"
        COMMAND Python3::Interpreter "${PROJECT_SOURCE_DIR}/tools/bench_report.py" --self-test)
    set_tests_properties("tools: complexity self-test" "tools: complexity"
        PROPERTIES SKIP_RETURN_CODE 77)
    # ADR-0008: every dependency's declared licence is on the allow-list,
    # and the licence file each build fetched reads as the declared one.
    add_test(NAME "tools: dependency licences self-test"
        COMMAND Python3::Interpreter "${PROJECT_SOURCE_DIR}/tools/check_dependencies.py"
                --self-test)
    add_test(NAME "tools: dependency licences"
        COMMAND Python3::Interpreter "${PROJECT_SOURCE_DIR}/tools/check_dependencies.py"
                --root "${PROJECT_SOURCE_DIR}" --deps-dir "${FETCHCONTENT_BASE_DIR}")
endif()
