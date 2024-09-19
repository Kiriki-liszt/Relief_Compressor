#!/bin/bash
echo "This will codesign .vst and .compenent plugins"
echo "This will execute these lines;"
echo "xattr -r -d com.apple.quarantine /Library/Audio/Plug-Ins/VST3/VLC_EQ.vst3"
echo "xattr -r -d com.apple.quarantine /Library/Audio/Plug-Ins/Components/VLC_EQ.component"
echo "codesign --force --sign - /Library/Audio/Plug-Ins/VST3/VLC_EQ.vst3"
echo "codesign --force --sign - /Library/Audio/Plug-Ins/Components/VLC_EQ.component"
sudo xattr -r -d com.apple.quarantine /Library/Audio/Plug-Ins/VST3/VLC_EQ.vst3
sudo xattr -r -d com.apple.quarantine /Library/Audio/Plug-Ins/Components/VLC_EQ.component
sudo codesign --force --sign - /Library/Audio/Plug-Ins/VST3/VLC_EQ.vst3
sudo codesign --force --sign - /Library/Audio/Plug-Ins/Components/VLC_EQ.component
echo "codesigning ended!"