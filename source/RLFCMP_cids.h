//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/base/funknown.h"
#include "pluginterfaces/vst/vsttypes.h"

namespace yg331 {
//------------------------------------------------------------------------
static const Steinberg::FUID kRLFCMP_ProcessorUID  (0xC7AA80BF, 0x71B65A2A, 0x887BCCE9, 0x6CE488A7);
static const Steinberg::FUID kRLFCMP_ControllerUID (0xA8D88BE9, 0x26535B50, 0xBB180284, 0xF9981F18);

#define RLFCMP_VST3Category "Fx|Dynamics"

// SideChain adaptive oversampling to 192kHz at 48kHz
// start from Filter

enum {
    kParamBypass = 0,
    kParamSoftBypass, 
    kParamZoom,
    kParamOS,          // change it to HQ mode?
    
    kParamSidechainFilter, // BS1770 on/off
    
    // GML 8900 ~ MDWDRC's method is not so straight-forward
    // As patent shows, attack comes from exponent parameter('Exponent'), and adds corrective release('Timing') after
    // Looks like some Airwindows type of algo...
    // maybe, he was right after all
    
    // GML 8900 focuses on Attack blend - Release does not change while blending - for dynamic range control
    // Weiss DS1 focuses on Release blend - Attack does nor change while blending - for getting the right feel when mastering
    
    // GML has fast attack / fast release and slow attack / slow release
    // So, fast release cannot take over slow release
    // resulting 2-stage attack and fixed slow release.
    // Meaning, It generally smooth out dynamic range, but able to catch fast peak if happened
    // truly a dynamic range contoller
    
    // Weiss has same attack, but has fast release and slow release
    // resulting 1-stage attack and 2-stage release.
    // Meaning, fast GR will start to recover fast, and slow down as returning
    // natural? program dependant? auto release - like.
    
    // I think 2-RMS detectors might work better
    
    // Attack - Release - Bias - Blend ?
    
    // Attack is in Log, but Release is in Linear, like MH-CS and GML
    
    // GML style   : Attack - Release - NULL -Bias
    // Weiss style : Attack - Fast R - Slow R - AVG
    
    // Envelope Detection
    kParamAttack,  // in ms, Fast RMS : *= 1/2 , Slow RMS : *= 5
    kParamRelease, // Same?
    kParamUNUSED,  //
    kParamBias,
    // Blend should keep max dB at constant, 0-50 -> 1 : 0-1 / 50-100 -> 1-0 : 1
    // Crossover feel

    // Non-Linear Transfer Curve
    kParamThreshold,   // GML-stlye : input is turned up to make threshold go down, and turned down same amount. Silimar to Auto-gain of Matric Halo CS.
    kParamRatio,       // HOW does ratio affects linear release curve to be log-like in GML 8900 ?????
    kParamKnee,        // HOW does knee affects linear release curve to be log-like in MH-CS ?????
    kParamMakeup,
    
    kParamMix,
    kParamOutput
};

//------------------------------------------------------------------------
} // namespace yg331
