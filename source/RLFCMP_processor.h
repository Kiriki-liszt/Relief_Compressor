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
    // Dunno why, but these are happy only if these are defines here
    // Otherwise, the capacity of vector gets ridiculusly large
    std::vector<std::deque<ParamValue>*> lookAheadDelayLine;
    std::vector<std::deque<ParamValue>*> latencyDelayLine;
    
    // Internal functions ===========================================================
    template <typename SampleType>
    void processAudio (SampleType** inputs, SampleType** outputs, int32 numChannels, SampleRate getSampleRate, int32 sampleFrames);
    
    // alpha = 1 - exp(-1 / (sec * 0.001 * SR)) = tan(1 / (2 * sec * 0.001 * SR))
    ParamValue getTau  (ParamValue sec, SampleRate SR) { return exp(-1.0 / (sec * 0.001 * SR)); }
    
    void call_after_SR_changed ();
    void call_after_parameter_changed ();
    void sendFloat (Steinberg::Vst::IAttributeList::AttrID aid, double& value);
    
    // Some definitions ===========================================================
    static SMTG_CONSTEXPR ParamValue k_HT  = 300.0;
    static SMTG_CONSTEXPR ParamValue k_rms = 600.0; // if small, very short attack makes step in release
    static SMTG_CONSTEXPR ParamValue k_log = 100.0;
    static SMTG_CONSTEXPR ParamValue bw = 300.0; // Hz
    
    static SMTG_CONSTEXPR int32 HT_stage = 3;
    static SMTG_CONSTEXPR int32 HT_order = HT_stage * 2; // Number of coefficients, must be even
    
    static SMTG_CONSTEXPR int32 path_ref = 0;
    static SMTG_CONSTEXPR int32 path_sft = 1;
    static SMTG_CONSTEXPR int32 path_num = 2;
    
    static SMTG_CONSTEXPR int32 io_x = 0;
    static SMTG_CONSTEXPR int32 io_y = 1;
    static SMTG_CONSTEXPR int32 io_num = 2;
    
    static SMTG_CONSTEXPR int32 x1 = 0;
    static SMTG_CONSTEXPR int32 x2 = 1;
    static SMTG_CONSTEXPR int32 y1 = 0;
    static SMTG_CONSTEXPR int32 y2 = 1;
    static SMTG_CONSTEXPR int32 order = 2;
    
    static SMTG_CONSTEXPR int32 maxLAH = 256;
    
    static SMTG_CONSTEXPR double peakEnvDecay = 0.4; //sec
    static SMTG_CONSTEXPR double peakRMSDecay = 0.3; //sec
    
    // Internal datastructures ===========================================================
    ParamValue HT_coefs[path_num][HT_stage];
    ParamValue HT_state[2][path_num][io_num][order][HT_stage] = {0, }; // 2-channel, 2-path, 2-state(x, y), 2-order(x1, x2), numCoefs-stages,
    
    ParamValue detector_state[2]; // Leaky Integrator, naive one-pole filter, EWMA, etc
    
    PassShelfFilter SC_LF[2], SC_HF[2]; // SideChain Filters, simplified for pass and shelf, 6dB/oct 1p1z filter
    
    int32 lookaheadSize = 0;
    int32 halfTap = lookaheadSize / 2;
    int32 condition = lookaheadSize % 2;
    
    ParamValue LAH_coef[maxLAH] = {0.0, };
    
    LevelEnvelopeFollower VuInputRMS[2] = {
        LevelEnvelopeFollower(LevelEnvelopeFollower::rmsEnv, peakRMSDecay),
        LevelEnvelopeFollower(LevelEnvelopeFollower::rmsEnv, peakRMSDecay)};
    LevelEnvelopeFollower VuOutputRMS[2] = {
        LevelEnvelopeFollower(LevelEnvelopeFollower::rmsEnv, peakRMSDecay),
        LevelEnvelopeFollower(LevelEnvelopeFollower::rmsEnv, peakRMSDecay)};
    LevelEnvelopeFollower VuInputPeak[2] = {
        LevelEnvelopeFollower(LevelEnvelopeFollower::peakEnv, peakEnvDecay),
        LevelEnvelopeFollower(LevelEnvelopeFollower::peakEnv, peakEnvDecay)};
    LevelEnvelopeFollower VuOutputPeak[2] = {
        LevelEnvelopeFollower(LevelEnvelopeFollower::peakEnv, peakEnvDecay),
        LevelEnvelopeFollower(LevelEnvelopeFollower::peakEnv, peakEnvDecay)};

    // Parameters ===========================================================
    bool       pBypass     = dftBypass > 0;
    bool       pSoftBypass = dftSoftBypass > 0;
    // pZoom,
    overSample pOS         = overSample_1x;
    
    bool       pScLfIn     = dftScLfIn > 0;
    ParamValue pScLfType   = PassShelfFilter::rLow;
    ParamValue pScLfFreq   = paramScLfFreq. ToNormalized(dftScLfFreq);
    ParamValue pScLfGain   = paramScLfGain. ToNormalized(dftScLfGain);
    bool       pScHfIn     = dftScHfIn > 0;
    ParamValue pScHfType   = PassShelfFilter::rHigh;
    ParamValue pScHfFreq   = paramScHfFreq. ToNormalized(dftScHfFreq);
    ParamValue pScHfGain   = paramScHfGain. ToNormalized(dftScHfGain);
    bool       pScListen   = false;
    
    ParamValue pType       = paramType.     ToNormalized(dftType);
    ParamValue pAttack     = paramAttack.   ToNormalized(dftAttack);
    ParamValue pRelease    = paramRelease.  ToNormalized(dftRelease);
    bool       pLookaheadEnable = true;
    
    ParamValue pThreshold  = paramThreshold.ToNormalized(dftThreshold);
    ParamValue pRatio      = paramRatio.    ToNormalized(dftRatio);
    ParamValue pKnee       = paramKnee.     ToNormalized(dftKnee);
    ParamValue pMakeup     = paramMakeup.   ToNormalized(dftMakeup);
    
    ParamValue pMix        = paramMix.      ToNormalized(dftMix);
    ParamValue pInput      = paramInput.    ToNormalized(dftInput);
    ParamValue pOutput     = paramOutput.   ToNormalized(dftOutput);
    
    // Values for GUI ========================================================
    ParamValue fInputVuRMS[2] = {0.0, }, fOutputVuRMS[2] = {0.0, };  // for each channel
    ParamValue fInputVuPeak[2] = {0.0, }, fOutputVuPeak[2] = {0.0, };
    ParamValue fGainReduction = 0.0;
    
    // Internal plain values =================================================
    SampleRate projectSR   = 48000.0;
    SampleRate internalSR  = 192000.0;//actually, not used
    ParamValue detectorIndicator = 0.0;
    
    int32      type        = paramType.ToPlainList(pType);
    ParamValue atkCoef     = getTau(dftAttack,  projectSR);
    ParamValue rlsCoef     = getTau(dftRelease, projectSR);

    ParamValue ratio       = paramRatio.ToPlain(pRatio);
    ParamValue slope       = 1.0 / ratio - 1.0;
    ParamValue knee        = paramKnee.ToPlain(pKnee);
    ParamValue kneeHalf    = knee / 2.0;
    ParamValue threshold   = paramThreshold.ToPlain(pThreshold);
    ParamValue makeup      = DecibelConverter::ToGain(paramMakeup.ToPlain(pMakeup));
    
    ParamValue mix         = pMix;
    ParamValue input       = DecibelConverter::ToGain(paramInput. ToPlain(pInput));
    ParamValue output      = DecibelConverter::ToGain(paramOutput.ToPlain(pOutput));
};

//------------------------------------------------------------------------
} // namespace yg331

/*
 BS1770 k-weight filter
 
 Pre-Warped works fine, since it's only shlelf
     Pre-Filter
     f = 1500.0 / Q = 1.007 / dB = 4.0
     
     LRB
     f = 38.134 / Q = 0.70758 / level = +0.038dB
 
 With cramped, untreated EQ:
     Pre-Filter
     f = 1681.8540338092508599234833913
     q = 1 / sqrt(2) == M_SQRT1_2
     gain = 3.999849576502195080962565
     
     RLB
     f = 38.11054352868593
     q = 0.5
 */

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
