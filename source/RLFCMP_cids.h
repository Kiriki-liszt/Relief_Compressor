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
    kParamZoom,
    kParamOS,          // change it to HQ mode?
    
    kParamScLfIn,
    kParamScLfType,
    kParamScLfFreq,
    kParamScLfGain,
    
    kParamScHfIn,
    kParamScHfType,
    kParamScHfFreq,
    kParamScHfGain,
    
    kParamScListen,
    
    // Envelope Detection
    kParamType,
    kParamAttack,
    kParamRelease,
    kParamLookaheadEnable,

    // Non-Linear Transfer Curve
    kParamThreshold,   // GML-stlye : input is turned up to make threshold go down, and turned down same amount. Silimar to Auto-gain. 
    kParamRatio,
    kParamKnee,
    kParamMakeup,
    
    kParamMix,
    kParamOutput,
    kParamSoftBypass
};

enum
{
    kInMono = 100,
    kInLRMS,
    kInRRMS,
    kInLPeak,
    kInRPeak,
    kOutMono,
    kOutLRMS,
    kOutRRMS,
    kOutLPeak,
    kOutRPeak,
    kGainReduction
};

//------------------------------------------------------------------------
} // namespace yg331
