# Relief Compressor  

Relief Compressor is a compressor for everyday tasks, anywhere from mixing a track to master bus.  

Runs in double precision 64-bit internal processing. Also double precision input / output if supported.  
Selectable Detector path in Peak/RMS and Linear/Logarithmic options.  
Dedicated sidechain-EQ page with 2 band EQ, pass filter and shelf selectable.  
Fixed Look-Ahead at 0.5ms. It uses Kaiser-Bessel FIR as lookahead smoother.  

Windows and Mac, VST3 and AU.  

[![GitHub Release](https://img.shields.io/github/v/release/kiriki-liszt/Relief_Compressor?style=flat-square&label=Get%20latest%20Release)](https://github.com/Kiriki-liszt/Relief_Compressor/releases/latest)
[![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/kiriki-liszt/Relief_Compressor/total?style=flat-square&label=total%20downloads&color=blue)](https://github.com/Kiriki-liszt/Relief_Compressor/releases/latest)  

[![Static Badge](https://img.shields.io/badge/coffee%20maybe%3F%20%3D%5D%20-gray?style=for-the-badge&logo=buy-me-a-coffee)](https://buymeacoffee.com/kirikiaris)  

<img src="https://github.com/Kiriki-liszt/Relief_Compressor/blob/main/screenshot/Relief%20Compressor%20Main%20page.png?raw=true"  width="600"/>
<img src="https://github.com/Kiriki-liszt/Relief_Compressor/blob/main/screenshot/Relief%20Compressor%20SC%20page.png?raw=true"  width="600"/>  

## Windows  

- x64  

## macOS  

- 10.13(High Sierra) to 15.0(Sequoia)  

## How to use  

1. Windows

Unzip Win.zip from latest release and copy to "C:\Program Files\Common Files\VST3".  

2. MacOS(Intel tested, Apple Silicon not tested).  

Unzip MacOS.zip from latest release and copy vst3 to "/Library/Audio/Plug-Ins/VST3" and component to "/Library/Audio/Plug-Ins/Components".  

> If it doesn't go well, configure security options in console as  
>
> ``` console  
> sudo xattr -r -d com.apple.quarantine /Library/Audio/Plug-Ins/VST3/Relief_Compressor.vst3  
> sudo xattr -r -d com.apple.quarantine /Library/Audio/Plug-Ins/Components/Relief_Compressor.component
>
> sudo codesign --force --sign - /Library/Audio/Plug-Ins/VST3/Relief_Compressor.vst3  
> sudo codesign --force --sign - /Library/Audio/Plug-Ins/Components/Relief_Compressor.component
> ```  
>
> tested by @jonasborneland [here](https://github.com/Kiriki-liszt/JS_Inflator_to_VST2_VST3/issues/12#issuecomment-1616671177)

## Licensing  

Relief Compressor is using GPL v3 license.  

### VST  

> Q: I would like to share the source code of my VST 3 plug-in/host on GitHub or other such platform.  
>
> - You can choose the GPLv3 license and feel free to share your plug-ins/host's source code including or referencing the VST 3 SDK's sources on GitHub.  
> - **You are allowed to provide a binary form of your plug-ins/host too, provided that you provide its source code as GPLv3 too.**  
> - Note that you have to follow the Steinberg VST usage guidelines.  
>  
> <https://steinbergmedia.github.io/vst3_dev_portal/pages/FAQ/Licensing.html>  

![VST Logo](https://github.com/Kiriki-liszt/Sky_Blue_EQ4/assets/107096260/142e3c12-cd5f-415d-9b72-8b4f04419633)

VSTSDK 3.7.9 used  
VSTGUI 4.12 used  

## Version logs

v0.1.0 : intial beta.  

## Detectors  

### Peak  

Typical track compressor, widly used in consoles.  
It applies envelope as detected.  

### RMS  

Bit more smooth detector, used in buses and mastering.  
I used a semi-true-RMS detector, appliing envelope in squared manner using 1-pole filter.  
A true RMS would use square-normalize-sum-sqrt.  

### Linear  

Detects level in linear gain scale, and then calculates gain reduction.  
It lets fast peaks go through, sounding more open and punchy.  

### Logarithmic  

Calculates gain reduction first in log dB, and applies level envelope to gain reduction.  
More accurate attack timing and good for peak comtrolling.  

## Lookahead  

Typically lookahead is described as simply delaying main path by some amount.  
It's true but we need more to make it better.  

After delaying main signal path, we can apply smoothing filter to detector path.  
This way, the fastest peak reduction fades-in and reaches maximum reduction in sync with main signal.  
Smoothing can be appied many ways, simple lowpass might work, I choose FIR filtering.  
std::transform_reduce works well, since it supports parallel computation.  

## Why Attack and Release affects Threshold?  

It happens if detected signal is rectified but not smoothed.  
Also it gets more apparent in single decoupled branching detectors, where release cannot smooth out rectified detector signal.  

It can be solved by appling sub-1ms-attack few-ms-release to rectified signal.  

Need to know is, that it changes harmonic pattern generated from compressor.  
Reduction rate of harmonics are more steep.  

## Hilbert detector  

It uses two path with about 90-degree phase difference.  
If both path's signal are square-add-sqrted, it gets ideal level envelope.  
This creates very clean harmonics, So I tuned it to have same amount of low harmonic, but less high harmonics.  

However, it causes overshoots and has issue with square signal.  
It somehow reacts like bandpassed detector signal, so this is given as an option.  

### It sounds more compressed then others  

When compressing, harmonic distortion happens and fills up a removed loudness.  
With Hilbert detector, that harmonic distortion is surpressed, leading to a sence of feel that it grabs more then it should.  
Our ears are used to hear that volume drop compensation by harmonic distortion, so when it doesn't, it feels somewhat unexpected.  
