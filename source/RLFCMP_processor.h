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
    Steinberg::uint32  PLUGIN_API getLatencySamples() SMTG_OVERRIDE;
    
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
    void processAudio (SampleType** inputs, SampleType** outputs, int32 numChannels, SampleRate getSampleRate, int32 sampleFrames);
    
    inline void sendFloat (Steinberg::Vst::IAttributeList::AttrID aid, double& value);
    
    static constexpr ParamValue fa = 0.25;
    static constexpr ParamValue sa = 9.0;
    static constexpr ParamValue fr = 0.5;
    static constexpr ParamValue sr = 1.0;
    // alpha = 1 - exp(-1 / (sec * 0.001 * SR)) = tan(1 / (2 * sec * 0.001 * SR))
    inline ParamValue getFastAttack  (ParamValue sec, SampleRate SR) { return exp(-1.0 / (fa * sec * 0.001 * SR)); }
    inline ParamValue getSlowAttack  (ParamValue sec, SampleRate SR) { return exp(-1.0 / (sa * sec * 0.001 * SR)); }
    inline ParamValue getFastRelease (ParamValue sec, SampleRate SR) { return exp(-1.0 / (fr * sec * 0.001 * SR)); }
    inline ParamValue getSlowRelease (ParamValue sec, SampleRate SR) { return exp(-1.0 / (sr * sec * 0.001 * SR)); }
    
    inline void call_after_SR_changed ()
    {
        fast_atk  = getFastAttack (paramAttack. ToPlain(pAttack),  projectSR);
        fast_rls  = getFastRelease(paramRelease.ToPlain(pRelease), projectSR);
        slow_atk  = getSlowAttack (paramAttack. ToPlain(pAttack),  projectSR);
        slow_rls  = getSlowRelease(paramRelease.ToPlain(pRelease), projectSR);
        
        for (int32 channel = 0; channel < 2; channel++)
        {
            BS1770_PF[channel]. setSVF(pSidechainFilter, PF_FREQ, PF_Q, PF_dB, simpleSVF::HighShelf, projectSR);
            BS1770_RLB[channel].setSVF(pSidechainFilter, RLB_FREQ, RLB_Q, 0.0, simpleSVF::HighPass,  projectSR);
        }
        
        ParamValue transition = 2.0 * bw / projectSR; // 90 deg phase difference band is from 20 Hz to Nyquist - 20 Hz. The transition bandwidth is twice 20 Hz.

        ParamValue coefs[HT_order];
        
        hiir::PolyphaseIir2Designer::compute_coefs_spec_order_tbw (coefs, HT_order, transition);
        
        // Phase reference path c coefficients
        for (int i = 1, j = 0; i < HT_order; i += 2)
            HT_coefs[path_ref][j++] = coefs[i];

        // +90 deg path c coefficients
        for (int i = 0, j = 0; i < HT_order; i += 2)
            HT_coefs[path_sft][j++] = coefs[i];
    }
    inline void call_after_parameter_changed ()
    {
        fast_atk  = getFastAttack (paramAttack. ToPlain(pAttack),  projectSR);
        fast_rls  = getFastRelease(paramRelease.ToPlain(pRelease), projectSR);
        slow_atk  = getSlowAttack (paramAttack. ToPlain(pAttack),  projectSR);
        slow_rls  = getSlowRelease(paramRelease.ToPlain(pRelease), projectSR);

        bias      = paramBias.ToPlain(pBias);
        biasGain  = DecibelConverter::ToGain(-abs(bias));
        ratio     = paramRatio.ToPlain(pRatio);
        slope     = 1.0 / ratio - 1.0;
        knee      = paramKnee.ToPlain(pKnee);
        kneeHalf  = knee / 2.0;
        threshold = (ratio == 1.0) ? 0.0 : paramThreshold.ToPlain(pThreshold) * (1.0 + (1.0/(ratio - 1.0))); // in dB
        preGain   = (ratio == 1.0) ? 1.0 : DecibelConverter::ToGain(-paramThreshold.ToPlain(pThreshold));    // in Gain
        makeup    = DecibelConverter::ToGain(paramMakeup.ToPlain(pMakeup));
        
        for (int32 channel = 0; channel < 2; channel++)
        {
            BS1770_PF[channel]. setIn(pSidechainFilter);
            BS1770_RLB[channel].setIn(pSidechainFilter);
        }
    }
    
    // Parameters
    bool       pBypass     = dftBypass > 0;
    bool       pSoftBypass = dftSoftBypass > 0;
    // pZoom,
    overSample pOS         = overSample_1x;
    
    bool       pSidechainFilter = dftSidechainFilter > 0;
    
    ParamValue pAttack     = paramAttack.   ToNormalized(dftAttack);
    ParamValue pRelease    = paramRelease.  ToNormalized(dftRelease);
    ParamValue pBias       = paramBias.     ToNormalized(dftBias);
    
    ParamValue pThreshold  = paramThreshold.ToNormalized(dftThreshold);
    ParamValue pRatio      = paramRatio.    ToNormalized(dftRatio);
    ParamValue pKnee       = paramKnee.     ToNormalized(dftKnee);
    ParamValue pMakeup     = paramMakeup.   ToNormalized(dftMakeup);
    
    ParamValue pMix        = paramMix.      ToNormalized(dftMix);
    ParamValue pOutput     = paramOutput.   ToNormalized(dftOutput);
    
    // Values for GUI
    ParamValue Input_L  = 0.0, Input_R  = 0.0;
    ParamValue Output_L = 0.0, Output_R = 0.0;
    ParamValue Gain_Reduction = 0.0;
    
    // Internal plain values
    SampleRate projectSR  = 48000.0;
    SampleRate internalSR = 192000.0;
    ParamValue detectorIndicator = 0.0;
    
    ParamValue fast_atk  = getFastAttack (dftAttack,  projectSR);
    ParamValue fast_rls  = getFastRelease(dftRelease, projectSR);
    ParamValue slow_atk  = getSlowAttack (dftAttack,  projectSR);
    ParamValue slow_rls  = getSlowRelease(dftRelease, projectSR);

    ParamValue bias      = paramBias.ToPlain(pBias);
    ParamValue biasGain  = DecibelConverter::ToGain(-abs(bias));
    ParamValue ratio     = paramRatio.ToPlain(pRatio);
    ParamValue slope     = 1.0 / ratio - 1.0;
    ParamValue knee      = paramKnee.ToPlain(pKnee);
    ParamValue kneeHalf  = knee / 2.0;
    ParamValue threshold = (ratio == 1.0) ? 0.0 : paramThreshold.ToPlain(pThreshold) * (1.0 + (1.0/(ratio - 1.0))); // in dB
    ParamValue preGain   = (ratio == 1.0) ? 1.0 : DecibelConverter::ToGain(-paramThreshold.ToPlain(pThreshold));    // in Gain
    ParamValue makeup    = DecibelConverter::ToGain(paramMakeup.ToPlain(pMakeup));
    
    /*
     numCoefs
     it relates to how Hilbert transforms into RMS in low frequencies.
     High number will change Hilbert to RMS quickly, low number will change gradually.
     Also high number will mangle transients more likely
     it does take effect in highend, but well beyond 20kHz.

     transition acts like crossover between plain RMS ans Hilbert
     so it's not a bad thing with high transition freq.
     It will take care of Square Wave issue.
     Although Higher freq means more aliasing - so we need more cleanup = oversampling
     
     TDR Insane   3x2 coef, bw = 400, k = 950
     TDR Presice  3x2 coef, bw = 300, k = 950
     TDR Live/Eco 2x2 coef, bw = 100, k = 800
     */
    static constexpr ParamValue k_HT = 300.0;
    static constexpr ParamValue k_rms = 600.0; // if small, very short attack makes step in release
    static constexpr ParamValue k_log = 100.0;
    static constexpr ParamValue bw = 300.0; // Hz
    static constexpr int32 HT_stage = 3;
    static constexpr int32 HT_order = HT_stage * 2; // Number of coefficients, must be even

    enum {
        path_ref = 0,
        path_sft = 1,
        path_num = 2
    };
    enum {
        io_x = 0,
        io_y = 1,
        io_num = 2
    };
    enum {
        x1 = 0, x2 = 1,
        y1 = 0, y2 = 1,
        order = 2
    };
    ParamValue HT_coefs[path_num][HT_stage];
    ParamValue state[2][path_num][io_num][order][HT_stage] = {0, }; // 2-channel, 2-path, 2-state(x, y), 2-order(x1, x2), numCoefs-stages,
    ParamValue fast_rms[2], slow_rms[2];
    
    /*
     BS1770 k-weight filter
     Provided filter coefficients in paper are for a sampling rate of 48 kHz, with 'cramp'.
     
     With un-cramped(reference : FabFilter Pro-Q 3) EQ, correct parameters are:
     
     Pre-Filter
     f = 1500.0 / Q = 1.007 / dB = 4.0
     ERROR = +- 0.005dB
     
     LRB
     f = 38.134 / Q = 0.70758 / level = +0.038dB
     
     Pre-Warp works fine, since it's only shlelf
     */
    static constexpr ParamValue PF_FREQ  = 1500.0;
    static constexpr ParamValue PF_Q     = M_SQRT1_2; // == 1.0 ~= 1.007
    static constexpr ParamValue PF_dB    = 4.0;
    static constexpr ParamValue RLB_FREQ = 38.134;
    static constexpr ParamValue RLB_Q    = 0.70758;
    simpleSVF BS1770_PF[2], BS1770_RLB[2];
    
    
    int32 lookaheadSize = 0;
    delayLine lookAheadDelayLine;
    delayLine latencyDelayLine;
    double LAH_coef[256] = {0.0, };
};

//------------------------------------------------------------------------
} // namespace yg331
/*
 BS1770 k-weight filter
 With cramped EQ:
 
 Pre-Filter
 f = 1681.8540338092508599234833913
 q = 1 / sqrt(2) == M_SQRT1_2
 gain = 3.999849576502195080962565
 
 RLB
 f = 38.11054352868593
 q = 0.5
 */
