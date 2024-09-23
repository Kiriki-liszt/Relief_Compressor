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
    
    // fast RMS min attack 0.5ms
    // slow RMS min attack 5.0ms
    
    double _coeff = exp(-1.0 / (20.0 * 0.001 * 48000.0));
    double _icoef = 1.0 - _coeff;
    
    double coeff = exp(-1.0 / (0.2 * 0.001 * 48000.0));
    double icoef = 1.0 - coeff;
    
    double HT_coef = exp(-1.0 / (10.0 * 0.001 * 48000.0));
    double HT_icof = 1.0 - HT_coef;
    double HT_envl[2];
    
    double rms_sin[2] = { 0.0, }, rms_cos[2] = { 0.0, }, rms[2] = { 0.0, };
    double peak[2] = {0.0, }, pp = 0.0;
    
    static constexpr double logistics_k = 950.0;
    
    // it relates to how Hilbert transforms into RMS in low frequencies.
    // High number will change Hilbert to RMS quickly, low number will change gradually.
    // Also high number will mangle transients more likely
    // it does take effect in highend, but well beyond 20kHz.
    
    // Live and Eco looks like 4
    // Presice looks 4~6
    // Insane looks like 6~8
    static constexpr int numCoefs = 6; // Number of coefficients, must be even
    static constexpr int numCoefsHalf = numCoefs / 2;
    double transition = 2*20.0/48000; // Sampling frequency is 44.1 kHz. Approx. 90 deg phase difference band is from 20 Hz to 22050 Hz - 20 Hz. The transition bandwidth is twice 20 Hz.

    double coefs[numCoefs];
    
    double c[2][numCoefsHalf];

    enum {
        io_x = 0,
        io_y = 1,
        io_num = 2
    };
    enum {
        path_ref = 0,
        path_sft = 1,
        path_num = 2
    };
    // 2-channel, 2-path, 2-state(x, y), 2-memory(x1, x2), numCoefs-stages,
    double state[2][path_num][io_num][2][numCoefs] = {0, };

    double ap_state[2][4][2] = {0, };
    double ap_coef[4] = {1.0, 0.6681786379192988, 0.41421356237309503, 0.198912367379658};
};

//------------------------------------------------------------------------
} // namespace yg331
