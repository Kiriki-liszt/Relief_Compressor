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
    // Internal functions ===========================================================
    template <typename SampleType>
    void processAudio (SampleType** inputs, SampleType** outputs, int32 numChannels, SampleRate getSampleRate, int32 sampleFrames);
    
    // alpha = 1.0 - exp(-1.0 / (sec * 0.001 * SR)) = tan(1.0 / (2.0 * sec * 0.001 * SR))
    static inline ParamValue getTau  (ParamValue msec, SampleRate SR) { return exp(-1.0 / (msec * 0.001 * SR)); }
    
    void call_after_SR_changed ();
    void call_after_parameter_changed ();
    void sendFloat (Steinberg::Vst::IAttributeList::AttrID aid, double& value);
    
    // Some definitions ===========================================================
    static SMTG_CONSTEXPR int32 maxChannel = 2; // Strictly Stereo
    
    static SMTG_CONSTEXPR double dcHz = 0.05;
    
    static SMTG_CONSTEXPR ParamValue k_HT  = 300.0;
    static SMTG_CONSTEXPR ParamValue k_rms = 600.0; // if small, very short attack makes step in release
    static SMTG_CONSTEXPR ParamValue k_detector = 3000.0;
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
    
    // timing of detectors are important then I thought
    // if attack is set to 0, 'tone' of compressor feels like it's smashing at 0 attack, even with attack knob backed off
    static SMTG_CONSTEXPR double detectorAtk = 0.4; //msec, about size of lookahead(fixed to 0.5ms)
    static SMTG_CONSTEXPR double detectorRls = 15.0; //msec ~= HOLD
    
    // Internal datastructures ===========================================================
    std::vector<double> sidechain_EQed[maxChannel];
    
    SVF_12     SC_LF[maxChannel];
    SVF_12     SC_HF[maxChannel];
    
    ParamValue HT_coefs[path_num][HT_stage];
    ParamValue HT_state[maxChannel][path_num][io_num][order][HT_stage] = {0, }; // 2-channel, 2-path, 2-state(x, y), 2-order(x1, x2), numCoefs-stages,
    
    ParamValue hilbert_state[maxChannel];
    ParamValue squared_state[maxChannel];
    ParamValue rectified_state[maxChannel];
    ParamValue envelope_state[maxChannel]; // Leaky Integrator, naive one-pole filter, EWMA, etc

    int32 lookaheadSize = 0;
    int32 halfTap = lookaheadSize / 2;
    int32 condition = lookaheadSize % 2;
    
    ParamValue LAH_coef[maxLAH] = {0.0, };
    
    std::deque<ParamValue> lookAheadDelayLine[maxChannel];
    std::deque<ParamValue> latencyDelayLine[maxChannel];
    
    LevelEnvelopeFollower VuInputRMS[maxChannel] = {
        LevelEnvelopeFollower(LevelEnvelopeFollower::rmsEnv, peakRMSDecay),
        LevelEnvelopeFollower(LevelEnvelopeFollower::rmsEnv, peakRMSDecay)};
    LevelEnvelopeFollower VuOutputRMS[maxChannel] = {
        LevelEnvelopeFollower(LevelEnvelopeFollower::rmsEnv, peakRMSDecay),
        LevelEnvelopeFollower(LevelEnvelopeFollower::rmsEnv, peakRMSDecay)};
    LevelEnvelopeFollower VuInputPeak[maxChannel] = {
        LevelEnvelopeFollower(LevelEnvelopeFollower::peakEnv, peakEnvDecay),
        LevelEnvelopeFollower(LevelEnvelopeFollower::peakEnv, peakEnvDecay)};
    LevelEnvelopeFollower VuOutputPeak[maxChannel] = {
        LevelEnvelopeFollower(LevelEnvelopeFollower::peakEnv, peakEnvDecay),
        LevelEnvelopeFollower(LevelEnvelopeFollower::peakEnv, peakEnvDecay)};

    // Parameters ===========================================================
    bool       pBypass     = dftBypass > 0;
    bool       pSoftBypass = dftSoftBypass > 0;
    // pZoom,
    int        pOS         = overSample_1x;
    
    bool       pScLfIn     = dftScLfIn > 0;
    ParamValue pScLfType   = SVF_6::rLow;
    ParamValue pScLfFreq   = paramScLfFreq. ToNormalized(dftScLfFreq);
    ParamValue pScLfGain   = paramScLfGain. ToNormalized(dftScLfGain);
    bool       pScHfIn     = dftScHfIn > 0;
    ParamValue pScHfType   = SVF_6::rHigh;
    ParamValue pScHfFreq   = paramScHfFreq. ToNormalized(dftScHfFreq);
    ParamValue pScHfGain   = paramScHfGain. ToNormalized(dftScHfGain);
    bool       pScListen   = false;
    
    ParamValue pDType      = paramDetectorType.ToNormalized(dftDetectorType);
    ParamValue pSCTopology = paramSidechainTopology.ToNormalized(dftSidechainTopology);
    bool       pHilbertEnable   = true;
    bool       pLookaheadEnable = true;
    ParamValue pAttack     = paramAttack.   ToNormalized(dftAttack);
    ParamValue pRelease    = paramRelease.  ToNormalized(dftRelease);
    
    ParamValue pThreshold  = paramThreshold.ToNormalized(dftThreshold);
    ParamValue pRatio      = paramRatio.    ToNormalized(dftRatio);
    ParamValue pKnee       = paramKnee.     ToNormalized(dftKnee);
    ParamValue pMakeup     = paramMakeup.   ToNormalized(dftMakeup);
    
    ParamValue pMix        = paramMix.      ToNormalized(dftMix);
    ParamValue pInput      = paramInput.    ToNormalized(dftInput);
    ParamValue pOutput     = paramOutput.   ToNormalized(dftOutput);
    
    // Values for GUI ========================================================
    ParamValue fInputVuRMS[maxChannel] = {0.0, }, fOutputVuRMS[maxChannel] = {0.0, };
    ParamValue fInputVuPeak[maxChannel] = {0.0, }, fOutputVuPeak[maxChannel] = {0.0, };
    ParamValue fGainReduction = 0.0;
    
    // Internal plain values =================================================
    SampleRate projectSR   = 48000.0;
    SampleRate internalSR  = 192000.0;//actually, not used
    ParamValue detectorIndicator = 0.0;
    
    ParamValue dtrAtkCoef  = getTau(detectorAtk,  projectSR);
    ParamValue dtrRlsCoef  = getTau(detectorRls,  projectSR);
    
    int32      dType       = paramDetectorType.ToPlainList(pDType);
    int32      scTopology  = paramSidechainTopology.ToPlainList(pSCTopology);
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
    
    // DC Blocker
    // double DC_state_x[2];
    // double DC_state_y[2];
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
