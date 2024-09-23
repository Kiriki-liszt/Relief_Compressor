//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "RLFCMP_shared.h"
#include "public.sdk/source/vst/vstaudioeffect.h"

namespace yg331 {

//------------------------------------------------------------------------
//  RLFCMP_Processor
//------------------------------------------------------------------------
class RLFCMP_Processor : public Steinberg::Vst::AudioEffect
{
public:
    RLFCMP_Processor ();
    ~RLFCMP_Processor () SMTG_OVERRIDE;

    // Create function
    static Steinberg::FUnknown* createInstance (void* /*context*/)
    {
        return (Steinberg::Vst::IAudioProcessor*)new RLFCMP_Processor;
    }

    //--- ---------------------------------------------------------------------
    // AudioEffect overrides:
    //--- ---------------------------------------------------------------------
    /** Called at first after constructor */
    Steinberg::tresult PLUGIN_API initialize (Steinberg::FUnknown* context) SMTG_OVERRIDE;
    
    /** Called at the end before destructor */
    Steinberg::tresult PLUGIN_API terminate () SMTG_OVERRIDE;
    
    /** Switch the Plug-in on/off */
    Steinberg::tresult PLUGIN_API setActive (Steinberg::TBool state) SMTG_OVERRIDE;

    /** Will be called before any process call */
    Steinberg::tresult PLUGIN_API setupProcessing (Steinberg::Vst::ProcessSetup& newSetup) SMTG_OVERRIDE;
    
    /** Gets the current Latency in samples. */
    // Steinberg::uint32  PLUGIN_API getLatencySamples() SMTG_OVERRIDE;
    
    /** Asks if a given sample size is supported see SymbolicSampleSizes. */
    Steinberg::tresult PLUGIN_API canProcessSampleSize (Steinberg::int32 symbolicSampleSize) SMTG_OVERRIDE;

    /** Here we go...the process call */
    Steinberg::tresult PLUGIN_API process (Steinberg::Vst::ProcessData& data) SMTG_OVERRIDE;
        
    /** For persistence */
    Steinberg::tresult PLUGIN_API setState (Steinberg::IBStream* state) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API getState (Steinberg::IBStream* state) SMTG_OVERRIDE;

//------------------------------------------------------------------------
protected:
    template <typename SampleType>
    void processAudio(SampleType** inputs, SampleType** outputs, int32 numChannels, SampleRate getSampleRate, int32 sampleFrames);
    
    // Parameters
    bool       pBypass     = dftBypass > 0;
    bool       pSoftBypass = dftSoftBypass > 0;
    // pZoom,
    ParamValue pOS = 2.0 / 6.0;
    
    bool       pSidechainFilter = dftSidechainFilter > 0;

    ParamValue pAttack    = paramAttack.ToNormalized(dftAttack);
    ParamValue pRelease   = paramRelease.ToNormalized(dftRelease);
    ParamValue pBias      = paramBias.ToNormalized(dftBias);

    ParamValue pThreshold = paramThreshold.ToNormalized(dftThreshold);
    ParamValue pRatio     = paramRatio.ToNormalized(dftRatio);
    ParamValue pKnee      = paramKnee.ToNormalized(dftKnee);
    ParamValue pMakeup    = paramMakeup.ToNormalized(dftMakeup);
    
    ParamValue pMix       = paramMix.ToNormalized(dftMix);
    ParamValue pOutput    = paramOutput.ToNormalized(dftOutput);
    
    
    double coeff = exp(-1.0 / (2.0 * 0.001 * 48000.0));
    double icoef = 1.0 - coeff;
    double rms_sin[2] = { 0.0, }, rms_cos[2] = { 0.0, }, rms[2] = { 0.0, };
    double peak[2] = {0.0, }, pp = 0.0;
    
    static constexpr int numCoefs = 16; // Number of coefficients, must be even
    static constexpr int numCoefsHalf = 8; // Number of coefficients, must be even
    double transition = 2*20.0/48000; // Sampling frequency is 44.1 kHz. Approx. 90 deg phase difference band is from 20 Hz to 22050 Hz - 20 Hz. The transition bandwidth is twice 20 Hz.

    double coefs[numCoefs];
    
    double c[2][numCoefsHalf] = {
        {0.16514909355907719801,0.73982901254452670958,0.94794090632917971107,0.99120971270525837227},
        {0.48660436861367767358,0.88077943527246449484,0.97793125561632343601,0.99767386185073303473}
    };
    double p1_sx_1[2] = { 0.0, };
    double p1_sx_2[2] = { 0.0, };
    double p1_sy_1[2] = { 0.0, };
    double p1_sy_2[2] = { 0.0, };
    double p2_sx_1[2] = { 0.0, };
    double p2_sx_2[2] = { 0.0, };
    double p2_sy_1[2] = { 0.0, };
    double p2_sy_2[2] = { 0.0, };
    double p3_sx_1[2] = { 0.0, };
    double p3_sx_2[2] = { 0.0, };
    double p3_sy_1[2] = { 0.0, };
    double p3_sy_2[2] = { 0.0, };
    double p4_sx_1[2] = { 0.0, };
    double p4_sx_2[2] = { 0.0, };
    double p4_sy_1[2] = { 0.0, };
    double p4_sy_2[2] = { 0.0, };
    
    // 2-channel, 2-path, 2-state(x, y), 2-memory(x1, x2), numCoefs-stages,
    double state[2][2][2][2][numCoefs] = {0, };

    double sx_1[2] = { 0.0, };
    double sx_2[2] = { 0.0, };
    double sy_1[2] = { 0.0, };
    double sy_2[2] = { 0.0, };
};

//------------------------------------------------------------------------
} // namespace yg331
