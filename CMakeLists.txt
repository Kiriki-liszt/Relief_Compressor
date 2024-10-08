cmake_minimum_required(VERSION 3.14.0)
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13 CACHE STRING "")

if(SMTG_MAC)
    set(vst3sdk_SOURCE_DIR "/Users/yoon/VST_SDK/vst3sdk")
elseif(SMTG_WIN)
    set(vst3sdk_SOURCE_DIR "C:/VST_SDK/vst3sdk")
endif(SMTG_MAC)

if(GITHUB_ACTIONS)
    set(vst3sdk_SOURCE_DIR ${CMAKE_SOURCE_DIR}/vst3sdk)
endif()

if(NOT vst3sdk_SOURCE_DIR)
    message(FATAL_ERROR "Path to VST3 SDK is empty!")
endif()

project(Relief_Compressor
    # This is your plug-in version number. Change it here only.
    # Version number symbols usable in C++ can be found in
    # source/version.h and ${PROJECT_BINARY_DIR}/projectversion.h.
    VERSION 1.0.0.0 
    DESCRIPTION "Relief_Compressor VST 3 Plug-in"
)

set(SMTG_VSTGUI_ROOT "${vst3sdk_SOURCE_DIR}")

add_subdirectory(${vst3sdk_SOURCE_DIR} ${PROJECT_BINARY_DIR}/vst3sdk)
smtg_enable_vst3_sdk()

smtg_add_vst3plugin(Relief_Compressor
    source/version.h
    source/RLFCMP_cids.h
    source/RLFCMP_shared.h
    source/RLFCMP_shared.cpp
    source/RLFCMP_processor.h
    source/RLFCMP_processor.cpp
    source/RLFCMP_controller.h
    source/RLFCMP_controller.cpp
    source/RLFCMP_entry.cpp
)

#- VSTGUI Wanted ----
if(SMTG_ENABLE_VSTGUI_SUPPORT)
    target_sources(Relief_Compressor
        PRIVATE
            resource/RLFCMP_editor.uidesc
    )
    target_link_libraries(Relief_Compressor
        PRIVATE
            vstgui_support
    )
    smtg_target_add_plugin_resources(Relief_Compressor
        RESOURCES
            "resource/RLFCMP_editor.uidesc"
    )
endif(SMTG_ENABLE_VSTGUI_SUPPORT)
# -------------------

smtg_target_add_plugin_snapshots (Relief_Compressor
    RESOURCES
        resource/C7AA80BF71B65A2A887BCCE96CE488A7_snapshot.png
        resource/C7AA80BF71B65A2A887BCCE96CE488A7_snapshot_2.0x.png
)

target_link_libraries(Relief_Compressor
    PRIVATE
        sdk
)

smtg_target_configure_version_file(Relief_Compressor)

if(SMTG_MAC)
    smtg_target_set_bundle(Relief_Compressor
        BUNDLE_IDENTIFIER io.github.yg331.Relief.Compressor
        COMPANY_NAME "yg331"
    )
    smtg_target_set_debug_executable(Relief_Compressor
        "/Applications/VST3PluginTestHost.app"
        "--pluginfolder;$(BUILT_PRODUCTS_DIR)"
    )
elseif(SMTG_WIN)
    target_sources(Relief_Compressor PRIVATE 
        resource/win32resource.rc
    )
    if(MSVC)
        set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Relief_Compressor)

        smtg_target_set_debug_executable(Relief_Compressor
            "$(ProgramW6432)/Steinberg/VST3PluginTestHost/VST3PluginTestHost.exe"
            "--pluginfolder \"$(OutDir)/\""
        )
    endif()
endif(SMTG_MAC)

# Add an AUv2 target
if (SMTG_MAC AND XCODE AND SMTG_COREAUDIO_SDK_PATH)
	smtg_target_add_auv2(Relief_Compressor-au
		BUNDLE_NAME Relief_Compressor
		BUNDLE_IDENTIFIER io.github.yg331.Relief.Compressor.audiounit
		INFO_PLIST_TEMPLATE resource/au-info.plist
		VST3_PLUGIN_TARGET Relief_Compressor)
endif(SMTG_MAC AND XCODE AND SMTG_COREAUDIO_SDK_PATH)
