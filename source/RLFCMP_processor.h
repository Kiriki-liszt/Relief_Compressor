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
    
    /** We want to receive message. */
    Steinberg::tresult PLUGIN_API notify (Steinberg::Vst::IMessage* message) SMTG_OVERRIDE;
    
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
    
    static SMTG_CONSTEXPR int HT_STATE_RESET_SIZE = maxChannel * path_num * io_num * order * HT_stage;
    
    static SMTG_CONSTEXPR int32 maxLAH = 256;
    
    static SMTG_CONSTEXPR double peakEnvDecay = 0.4; //sec, for VU meter
    static SMTG_CONSTEXPR double peakRMSDecay = 0.3; //sec, for VU meter
    
    // timing of detectors are important then I thought
    // if attack is set to 0, 'tone' of compressor feels like it's smashing at 0 attack, even with attack knob backed off
    static SMTG_CONSTEXPR double detectorAtk = 0.3; //msec
    static SMTG_CONSTEXPR double detectorRls = 2.0; //msec ~= HOLD
    static SMTG_CONSTEXPR double detectorHlbtAtk = 0.8; //msec
    static SMTG_CONSTEXPR double detectorHlbtRls = 2.0; //msec ~= HOLD
    
    // Internal datastructures ===========================================================
    std::vector<double> sidechain_EQed[maxChannel];
    std::vector<double> level_vec[maxChannel];
    std::deque<ParamValue> lookAheadDelayLine[maxChannel];
    std::deque<ParamValue> latencyDelayLine[maxChannel];
    
    SVF_12     SC_LF[maxChannel];
    SVF_12     SC_HF[maxChannel];
    
    ParamValue HT_coefs[path_num][HT_stage];
    ParamValue HT_state[maxChannel][path_num][io_num][order][HT_stage] = {0, }; // 2-channel, 2-path, 2-state(x, y), 2-order(x1, x2), numCoefs-stages,
    
    ParamValue detectorAtkState[maxChannel] = {0.0, 0.0};
    ParamValue detectorRlsState[maxChannel] = {0.0, 0.0};
    ParamValue envelope_state[maxChannel]   = {0.0, 0.0}; // Leaky Integrator, naive one-pole filter, EWMA, etc

    int32 lookaheadSize = 24;
    
    ParamValue LAH_coef[maxLAH] = {0.0, };
    
    static SMTG_CONSTEXPR int32 _rms  = LevelEnvelopeFollower::rmsEnv;
    static SMTG_CONSTEXPR int32 _peak = LevelEnvelopeFollower::peakEnv;
    LevelEnvelopeFollower VuInputRMS[maxChannel]   = { LevelEnvelopeFollower(_rms, peakRMSDecay),  LevelEnvelopeFollower(_rms, peakRMSDecay)};
    LevelEnvelopeFollower VuOutputRMS[maxChannel]  = { LevelEnvelopeFollower(_rms, peakRMSDecay),  LevelEnvelopeFollower(_rms, peakRMSDecay)};
    LevelEnvelopeFollower VuInputPeak[maxChannel]  = { LevelEnvelopeFollower(_peak, peakEnvDecay), LevelEnvelopeFollower(_peak, peakEnvDecay)};
    LevelEnvelopeFollower VuOutputPeak[maxChannel] = { LevelEnvelopeFollower(_peak, peakEnvDecay), LevelEnvelopeFollower(_peak, peakEnvDecay)};
    LevelEnvelopeFollower VuDetector = LevelEnvelopeFollower(_peak, peakEnvDecay);
    
    
    // Values for GUI ========================================================
    ParamValue fInputVuRMS [maxChannel] = {0.0, }, fOutputVuRMS [maxChannel] = {0.0, };
    ParamValue fInputVuPeak[maxChannel] = {0.0, }, fOutputVuPeak[maxChannel] = {0.0, };
    ParamValue fDetectorLevel = 0.0;
    ParamValue fGainReduction = 0.0;
    bool       sendUI = false;
    
    // Internal plain values =================================================
    SampleRate projectSR    = 48000.0;
    // SampleRate internalSR   = 192000.0; //actually, not used
    
    ParamValue hilbertDtrAtkCoef = std::sqrt(getTau(detectorHlbtAtk, projectSR));
    ParamValue hilbertDtrRlsCoef = getTau(detectorHlbtRls, projectSR) * getTau(detectorHlbtRls, projectSR);
    
    bool       bypass       = dftBypass;
    int32      OS           = overSample_1x;
    
    bool       scLfIn       = dftScLfIn;
    int32      scLfType     = ScTypePass;
    ParamValue scLfFreq     = dftScLfFreq;
    ParamValue scLfGain     = dftScLfGain;
    bool       scHfIn       = dftScHfIn;
    int32      scHfType     = ScTypeShelf;
    ParamValue scHfFreq     = dftScHfFreq;
    ParamValue scHfGain     = dftScHfGain;
    bool       scListen     = dftScListen;
    
    int32      dType        = dftDetectorType;
    int32      scTopology   = dftSidechainTopology;
    bool       hilbertEnable    = dftHilbertEnable;
    bool       lookaheadEnable  = dftLookaheadEnable;
    ParamValue attack       = dftAttack;
    ParamValue release      = dftRelease;
    ParamValue dtrAtkCoef   = getTau(dftAttack * 0.01,  projectSR);
    ParamValue dtrRlsCoef   = getTau(detectorRls,  projectSR);
    ParamValue atkCoef      = getTau(dftAttack * 0.99,  projectSR);
    ParamValue rlsCoef      = getTau(dftRelease, projectSR);

    ParamValue threshold    = dftThreshold;
    ParamValue ratio        = dftRatio;
    ParamValue slope        = 1.0 / ratio - 1.0;
    ParamValue knee         = dftKnee;
    ParamValue kneeHalf     = knee / 2.0;
    ParamValue makeup       = DecibelConverter::ToGain(dftMakeup);

    ParamValue mix          = dftMix/maxMix;
    ParamValue inputGain    = DecibelConverter::ToGain(dftInput);
    ParamValue outputGain   = DecibelConverter::ToGain(dftOutput);
    bool       softBypass   = dftSoftBypass;
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
