//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#include "RLFCMP_processor.h"
#include "RLFCMP_cids.h"

#include "base/source/fstreamer.h"
#include "pluginterfaces/vst/ivstparameterchanges.h"
#include "public.sdk/source/vst/vstaudioprocessoralgo.h"
#include "public.sdk/source/vst/vsthelpers.h"

using namespace Steinberg;

namespace yg331 {
//------------------------------------------------------------------------
// RLFCMP_Processor
//------------------------------------------------------------------------
RLFCMP_Processor::RLFCMP_Processor ()
{
    //--- set the wanted controller for our processor
    setControllerClass (kRLFCMP_ControllerUID);
}

//------------------------------------------------------------------------
RLFCMP_Processor::~RLFCMP_Processor ()
{}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Processor::initialize (FUnknown* context)
{
    // Here the Plug-in will be instantiated
    
    //---always initialize the parent-------
    tresult result = AudioEffect::initialize (context);
    // if everything Ok, continue
    if (result != kResultOk)
    {
        return result;
    }

    //--- create Audio IO ------
    addAudioInput  (STR16 ("Stereo In"),  Steinberg::Vst::SpeakerArr::kStereo);
    addAudioOutput (STR16 ("Stereo Out"), Steinberg::Vst::SpeakerArr::kStereo);

    /* If you don't need an event bus, you can remove the next line */
    // addEventInput (STR16 ("Event In"), 1);
    
    // JUST for my sanity
    call_after_parameter_changed ();

    return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Processor::terminate ()
{
    // Here the Plug-in will be de-instantiated, last possibility to remove some memory!
    
    //---do not forget to call parent ------
    return AudioEffect::terminate ();
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Processor::setActive (TBool state)
{
    //--- called when the Plug-in is enable/disable (On/Off) -----
    return AudioEffect::setActive (state);
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Processor::process (Vst::ProcessData& data)
{
    Vst::IParameterChanges* paramChanges = data.inputParameterChanges;

    if (paramChanges)
    {
        int32 numParamsChanged = paramChanges->getParameterCount();

        for (int32 index = 0; index < numParamsChanged; index++)
        {
            Vst::IParamValueQueue* paramQueue = paramChanges->getParameterData(index);

            if (paramQueue)
            {
                Vst::ParamValue value;
                int32 sampleOffset;
                int32 numPoints = paramQueue->getPointCount();

                if (paramQueue->getPoint(numPoints - 1, sampleOffset, value) == kResultTrue) {
                    switch (paramQueue->getParameterId()) {
                        case kParamBypass:            bypass     = (value > 0.5); break;
                        // case kParamZoom:       pZoom       = value; break;
                        // case kParamOS:         pOS         = value; break;
                        case kParamScLfIn:            scLfIn        = (value > 0.5); break;
                        case kParamScLfType:          scLfType      = paramScLfType.ToPlainList(value); break;
                        case kParamScLfFreq:          scLfFreq      = paramScLfFreq.ToPlain(value); break;
                        case kParamScLfGain:          scLfGain      = paramScLfGain.ToPlain(value); break;
                        case kParamScHfIn:            scHfIn        = (value > 0.5); break;
                        case kParamScHfType:          scHfType      = paramScHfType.ToPlainList(value); break;
                        case kParamScHfFreq:          scHfFreq      = paramScHfFreq.ToPlain(value); break;
                        case kParamScHfGain:          scHfGain      = paramScHfGain.ToPlain(value); break;
                        case kParamScListen:          scListen      = (value > 0.5); break;
                        case kParamDetectorType:      dType         = paramDetectorType.ToPlainList(value); break;
                        case kParamSidechainTopology: scTopology    = paramSidechainTopology.ToPlainList(value); break;
                        case kParamLookaheadEnable:   lookaheadEnable = (value > 0.5); break;
                        case kParamHilbertEnable:     hilbertEnable   = (value > 0.5); break;
                        case kParamAttack:            attack        = paramAttack.ToPlain(value); break;
                        case kParamRelease:           release       = paramRelease.ToPlain(value); break;
                        case kParamThreshold:         threshold     = paramThreshold.ToPlain(value); break;
                        case kParamRatio:             ratio         = paramRatio.ToPlain(value); break;
                        case kParamKnee:              knee          = paramKnee.ToPlain(value); break;
                        case kParamMakeup:            makeup        = DecibelConverter::ToGain(paramMakeup.ToPlain(value)); break;
                        case kParamMix:               mix           = value; break;
                        case kParamInput:             inputGain     = DecibelConverter::ToGain(paramInput.ToPlain(value)); break;
                        case kParamOutput:            outputGain    = DecibelConverter::ToGain(paramOutput.ToPlain(value)); break;
                        case kParamSoftBypass:        softBypass    = (value > 0.5); break;
                        default: break;
                    }
                    call_after_parameter_changed ();
                }
            }
        }
    }


    if (data.numInputs == 0 || data.numOutputs == 0)
    {
        // nothing to do
        return kResultOk;
    }

    // (simplification) we suppose in this example that we have the same input channel count than the output
    int32 numChannels = data.inputs[0].numChannels;

    //---get audio buffers----------------
    uint32 sampleFramesSize = getSampleFramesSizeInBytes(processSetup, data.numSamples);
    void** in  = getChannelBuffersPointer(processSetup, data.inputs[0]);
    void** out = getChannelBuffersPointer(processSetup, data.outputs[0]);
    Vst::SampleRate SampleRate = processSetup.sampleRate;

    //---check if silence---------------
    // check if all channel are silent then process silent
    if (data.inputs[0].silenceFlags == Vst::getChannelMask(data.inputs[0].numChannels))
    {
        // mark output silence too (it will help the host to propagate the silence)
        data.outputs[0].silenceFlags = data.inputs[0].silenceFlags;

        // the plug-in has to be sure that if it sets the flags silence that the output buffer are
        // clear
        for (int32 i = 0; i < numChannels; i++)
        {
            // do not need to be cleared if the buffers are the same (in this case input buffer are
            // already cleared by the host)
            if (in[i] != out[i])
            {
                memset(out[i], 0, sampleFramesSize);
            }
        }
    }
    else {

        data.outputs[0].silenceFlags = data.inputs[0].silenceFlags;
        //---in bypass mode outputs should be like inputs-----
        if (bypass)
        {
            for (int32 channel = 0; channel < numChannels; channel++)
            {
                memcpy (out[channel], in[channel], sampleFramesSize);
            }
        }
        else
        {
            if (data.symbolicSampleSize == Vst::kSample32) {
                processAudio<Vst::Sample32>((Vst::Sample32**)in, (Vst::Sample32**)out, numChannels, SampleRate, data.numSamples);
            }
            else {
                processAudio<Vst::Sample64>((Vst::Sample64**)in, (Vst::Sample64**)out, numChannels, SampleRate, data.numSamples);
            }
        }
    }

    if (sendUI) {
        // can we draw attack-release curve in GUI...?
        //---send a message
        sendFloat(msgInputPeakL,    fInputVuPeak[0]);
        sendFloat(msgInputPeakR,    fInputVuPeak[1]);
        sendFloat(msgInputRMSL,     fInputVuRMS[0]);
        sendFloat(msgInputRMSR,     fInputVuRMS[1]);
        sendFloat(msgOutputPeakL,   fOutputVuPeak[0]);
        sendFloat(msgOutputPeakR,   fOutputVuPeak[1]);
        sendFloat(msgOutputRMSL,    fOutputVuRMS[0]);
        sendFloat(msgOutputRMSR,    fOutputVuRMS[1]);
        sendFloat(msgGainReduction, fGainReduction);
    }

    return kResultOk;
}

inline void RLFCMP_Processor::sendFloat (Steinberg::Vst::IAttributeList::AttrID aid, double& value)
{
    if (IPtr<Vst::IMessage> message = owned (allocateMessage ()))
    {
       message->setMessageID ("GUI");
       message->getAttributes()->setFloat (aid, value);
       sendMessage (message);
    }
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Processor::setupProcessing (Vst::ProcessSetup& newSetup)
{
    /*
     when instantiated, getState happens.
     then setupProcessing.
     then setState, but it not happen if first time instantiating.
     */

    // fprintf (stdout, "setupProcessing\n");
    
    projectSR = newSetup.sampleRate;
    
    // if      (projectSR > 96000.0) internalSR = projectSR;
    // else if (projectSR > 48000.0) internalSR = projectSR * 2.0;
    // else                          internalSR = projectSR * 4.0;
    
    for (int32 channel = 0; channel < maxChannel; channel++)
    {
        sidechain_EQed[channel].resize(newSetup.maxSamplesPerBlock);
        std::fill(sidechain_EQed[channel].begin(), sidechain_EQed[channel].end(), 0.0);
    }
    
    call_after_SR_changed ();
    call_after_parameter_changed (); // in case setState does not happen
    
    //--- called before any processing ----
    return AudioEffect::setupProcessing (newSetup);
}

uint32  PLUGIN_API RLFCMP_Processor::getLatencySamples()
{
    return lookaheadSize;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Processor::canProcessSampleSize (int32 symbolicSampleSize)
{
    // by default kSample32 is supported
    if (symbolicSampleSize == Vst::kSample32)
        return kResultTrue;

    // disable the following comment if your processing support kSample64
    if (symbolicSampleSize == Vst::kSample64)
        return kResultTrue;

    return kResultFalse;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Processor::notify (Vst::IMessage* message)
{
    if (!message)
        return kInvalidArgument;

    if (strcmp (message->getMessageID (), "GUI") == 0)
    {
        Steinberg::int64 data;
        if (message->getAttributes ()->getInt ("start", data) == kResultOk)
        {
            if (data == 1)
                sendUI = true;
            else
                sendUI = false;
            return kResultOk;
        }
    }

    return AudioEffect::notify (message);
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Processor::setState (IBStream* state)
{
    // called when we load a preset, the model has to be reloaded
    if (!state)
        return kResultFalse;
    
    // This happens after setupProcessing, so have to call call_after_parameter_changed ()
    
    // fprintf (stdout, "setState\n");
    IBStreamer streamer (state, kLittleEndian);

    bool    savedBypass     = dftBypass;
    int32   savedOS         = overSample_1x;
    
    bool    savedScLfIn     = dftScLfIn;
    int32   savedScLfType   = ScTypePass;
    double  savedScLfFreq   = dftScLfFreq;
    double  savedScLfGain   = dftScLfGain;
    bool    savedScHfIn     = dftScHfIn;
    int32   savedScHfType   = ScTypeShelf;
    double  savedScHfFreq   = dftScHfFreq;
    double  savedScHfGain   = dftScHfGain;
    bool    savedScListen   = dftScListen;
    
    int32   savedDType           = dftDetectorType;
    int32   savedScTopology      = dftSidechainTopology;
    bool    savedHilbertEnable   = dftHilbertEnable;
    bool    savedLookaheadEnable = dftLookaheadEnable;
    double  savedAttack          = dftAttack;
    double  savedRelease         = dftRelease;
    
    double  savedThreshold       = dftThreshold;
    double  savedRatio           = dftRatio;
    double  savedKnee            = dftKnee;
    double  savedMakeup          = DecibelConverter::ToGain(dftMakeup);
    
    double  savedMix             = dftMix/maxMix;
    double  savedInputGain       = DecibelConverter::ToGain(dftInput);
    double  savedOutputGain      = DecibelConverter::ToGain(dftOutput);
    bool    savedSoftBypass      = dftSoftBypass;

    if (streamer.readBool  (savedBypass)        == false) savedBypass     = dftBypass;
    if (streamer.readInt32 (savedOS)            == false) savedOS         = overSample_1x;
    
    if (streamer.readBool  (savedScLfIn)        == false) savedScLfIn     = dftScLfIn;
    if (streamer.readInt32 (savedScLfType)      == false) savedScLfType   = ScTypePass;
    if (streamer.readDouble(savedScLfFreq)      == false) savedScLfFreq   = dftScLfFreq;
    if (streamer.readDouble(savedScLfGain)      == false) savedScLfGain   = dftScLfGain;
    if (streamer.readBool  (savedScHfIn)        == false) savedScHfIn     = dftScHfIn;
    if (streamer.readInt32 (savedScHfType)      == false) savedScHfType   = ScTypeShelf;
    if (streamer.readDouble(savedScHfFreq)      == false) savedScHfFreq   = dftScHfFreq;
    if (streamer.readDouble(savedScHfGain)      == false) savedScHfGain   = dftScHfGain;
    if (streamer.readBool  (savedScListen)      == false) savedScListen   = dftScListen;
    
    if (streamer.readInt32 (savedDType)             == false) savedDType           = dftDetectorType;
    if (streamer.readInt32 (savedScTopology)        == false) savedScTopology      = dftSidechainTopology;
    if (streamer.readBool  (savedHilbertEnable)     == false) savedHilbertEnable   = dftHilbertEnable;
    if (streamer.readBool  (savedLookaheadEnable)   == false) savedLookaheadEnable = dftLookaheadEnable;
    if (streamer.readDouble(savedAttack)            == false) savedAttack          = dftAttack;
    if (streamer.readDouble(savedRelease)           == false) savedRelease         = dftRelease;
    
    if (streamer.readDouble(savedThreshold)     == false) savedThreshold       = dftThreshold;
    if (streamer.readDouble(savedRatio)         == false) savedRatio           = dftRatio;
    if (streamer.readDouble(savedKnee)          == false) savedKnee            = dftKnee;
    if (streamer.readDouble(savedMakeup)        == false) savedMakeup          = DecibelConverter::ToGain(dftMakeup);
    
    if (streamer.readDouble(savedMix)           == false) savedMix             = dftMix/maxMix;
    if (streamer.readDouble(savedInputGain)     == false) savedInputGain       = DecibelConverter::ToGain(dftInput);
    if (streamer.readDouble(savedOutputGain)    == false) savedOutputGain      = DecibelConverter::ToGain(dftOutput);
    if (streamer.readBool  (savedSoftBypass)    == false) savedSoftBypass      = dftSoftBypass;
    
    bypass      = savedBypass;
    OS          = savedOS;
    
    scLfIn      = savedScLfIn;
    scLfType    = savedScLfType;
    scLfFreq    = savedScLfFreq;
    scLfGain    = savedScLfGain;
    scHfIn      = savedScHfIn;
    scHfType    = savedScHfType;
    scHfFreq    = savedScHfFreq;
    scHfGain    = savedScHfGain;
    scListen    = savedScListen;
    
    dType           = savedDType;
    scTopology      = savedScTopology;
    hilbertEnable   = savedHilbertEnable;
    lookaheadEnable = savedLookaheadEnable;
    attack          = savedAttack;
    release         = savedRelease;
    
    threshold   = savedThreshold;
    ratio       = savedRatio;
    knee        = savedKnee;
    makeup      = savedMakeup;
    
    mix         = savedMix;
    inputGain   = savedInputGain;
    outputGain  = savedOutputGain;
    softBypass  = savedSoftBypass;
    
    call_after_parameter_changed ();
    
    return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Processor::getState (IBStream* state)
{
    // here we need to save the model
    if (!state)
        return kResultFalse;
    
    // fprintf (stdout, "getState\n");
    IBStreamer streamer (state, kLittleEndian);
    
    streamer.writeBool(bypass);
    streamer.writeInt32(OS);
    
    streamer.writeBool(scLfIn);
    streamer.writeInt32(scLfType);
    streamer.writeDouble(scLfFreq);
    streamer.writeDouble(scLfGain);
    streamer.writeBool(scHfIn);
    streamer.writeInt32(scHfType);
    streamer.writeDouble(scHfFreq);
    streamer.writeDouble(scHfGain);
    streamer.writeBool(scListen);
    
    streamer.writeInt32(dType);
    streamer.writeInt32(scTopology);
    streamer.writeBool(hilbertEnable);
    streamer.writeBool(lookaheadEnable);
    streamer.writeDouble(attack);
    streamer.writeDouble(release);
    
    streamer.writeDouble(threshold);
    streamer.writeDouble(ratio);
    streamer.writeDouble(knee);
    streamer.writeDouble(makeup);     // saved in lin gain, not dB
    
    streamer.writeDouble(mix);
    streamer.writeDouble(inputGain);  // saved in lin gain, not dB
    streamer.writeDouble(outputGain); // saved in lin gain, not dB
    streamer.writeBool(softBypass);
    
    return kResultOk;
}

template <typename SampleType>
void RLFCMP_Processor::processAudio(
    SampleType** inputs,
    SampleType** outputs,
    int32 numChannels,
    Vst::SampleRate SampleRate,
    int32 sampleFrames
)
{
    double invNumChannels = 1.0 / numChannels;
    ParamValue GR_Max = 0.0;
    
    ParamValue maxInputPeak[2] = {0.0, };
    ParamValue maxInputRMS[2] = {0.0, };
    ParamValue maxOutputPeak[2] = {0.0, };
    ParamValue maxOutputRMS[2] = {0.0, };
    
    double sqrtDtrAtkCoef = sqrt(dtrAtkCoef);
    double powrDtrRlsCoef = dtrRlsCoef * dtrRlsCoef;
    double sqrtAtkCoef = sqrt(atkCoef);
    double qirtAtkCoef = cbrt(sqrtAtkCoef); // Quintic Root
    double squrRlsCoef = rlsCoef * rlsCoef;
    double cubcRlsCoef = rlsCoef * squrRlsCoef;

    double hilbertDtrAtkCoef = std::sqrt(getTau(0.8, SampleRate));
    double hilbertDtrRlsCoef = getTau(15.0, SampleRate) * getTau(15.0, SampleRate);
    
    // because of locallity, this is faster
    for (int32 channel = 0; channel < numChannels; channel++)
    {
        for (int sample = 0; sample < sampleFrames; sample++)
        {
            // SideChain Filtering ===========================================================
            double t = SC_LF[channel].computeSVF(inputs[channel][sample]);
            t = SC_HF[channel].computeSVF(t);
            sidechain_EQed[channel][sample] = t;
        }
    }
    
    int32 sample = 0;
    // Process ===========================================================
    while (sampleFrames > sample)
    {
        double level[2] = {0.0, 0.0};

        for (int32 channel = 0; channel < numChannels; channel++)
        {
            Vst::Sample64 inputSample = sidechain_EQed[channel][sample];

            inputSample *= inputGain;

            if (hilbertEnable) {
                // Hilbert detector ===========================================================
                for (int path = 0; path < path_num; path++)
                {
                    double nextInput = inputSample;
                    for (int stage = 0; stage < HT_stage; stage++)
                    {
                        // with double-presicion, DF1 works fine...
                        double ret = HT_coefs[path][stage] * (nextInput + HT_state[channel][path][io_y][y2][stage]) - HT_state[channel][path][io_x][x2][stage];
                        HT_state[channel][path][io_x][x2][stage] = HT_state[channel][path][io_x][x1][stage];
                        HT_state[channel][path][io_x][x1][stage] = nextInput;
                        HT_state[channel][path][io_y][y2][stage] = HT_state[channel][path][io_y][y1][stage];
                        HT_state[channel][path][io_y][y1][stage] = ret;
                        nextInput = ret;
                    }
                }
                double rms_s   = HT_state[channel][path_ref][io_y][y2][HT_stage - 1];
                double rms_c   = HT_state[channel][path_sft][io_y][y1][HT_stage - 1];
                
                double vin = rms_s * rms_s + rms_c * rms_c;
                double delta = vin - detectorAtkState[channel];
                double pp = (delta > 0) ? 1.0 : 0.0; // smooth::Logistics(delta, k_detector); // (delta > 0) ? 1.0 : 0.0;
                double nn = 1.0 - pp;
                double g = hilbertDtrAtkCoef;
                detectorAtkState[channel] = g * detectorAtkState[channel] + (1.0 - g) * vin;
                level[channel] = detectorAtkState[channel];
                // std::hypot(x, y) = sqrt(x^2 + y^2)
            }
            else {
                // Square wave rectifier ===========================================================
                double vin = inputSample * inputSample;
                // Smooth decoupled,
                // not Smooth Branched, because it causes detected level drop
                detectorRlsState[channel] = std::max(vin, powrDtrRlsCoef * detectorRlsState[channel] + (1.0 - powrDtrRlsCoef) * vin);
                detectorAtkState[channel] = sqrtDtrAtkCoef * detectorAtkState[channel] + (1.0 - sqrtDtrAtkCoef) * detectorRlsState[channel];
                level[channel] = detectorAtkState[channel];
            }
        }
        
        // Link stereo ===========================================================
        double monoSample   = 0.0;
        for (int32 channel = 0; channel < numChannels; channel++)
        {
            monoSample   += invNumChannels * level[channel];
        }
        for (int32 channel = 0; channel < numChannels; channel++)
        {
            level[channel] = monoSample;
        }
        
        for (int32 channel = 0; channel < numChannels; channel++)
        {
            Vst::Sample64 inputSample = inputs[channel][sample];
            
            inputSample *= inputGain;
            
            double inputPeak = VuInputPeak[channel].processSample(inputSample);
            double inputRMS  = VuInputRMS[channel].processSample(inputSample);
            if(maxInputPeak[channel] < inputPeak) maxInputPeak[channel] = inputPeak;
            if(maxInputRMS [channel] < inputRMS ) maxInputRMS [channel] = inputRMS;
            
            double env = 0.0;
            double gain = 1.0;
            
            // Linear scale detection ===========================================================
            if (scTopology == ScTopologyLin)
            {
                // Envelope Detection ===========================================================
                switch (dType) {
                    case detectorPeak: default:  // Metric Halo ChannelStrip MIO Comp, and Weiss DS1-MK3 if atk*3 & rls*2
                    {
                        // double r = getTau(paramRelease.ToPlain(pRelease) - detectorRls, SampleRate);
                        double vin = sqrt(level[channel]);
                        double delta = vin - envelope_state[channel];
                        double pp = (delta > 0) ? 1.0 : 0.0; // smooth::Logistics(delta, k_rms); // (delta > 0) ? 1.0 : 0.0;
                        double nn = 1.0 - pp;
                        double g = (pp * sqrtAtkCoef + nn * squrRlsCoef);
                        envelope_state[channel] = g * envelope_state[channel] + (1.0 - g) * vin;
                        env = (envelope_state[channel]);
                        break;
                    }
                    case detectorRMS: // Semi-true-RMS. a true RMS would use square-sum-normalize-sqrt, not this 1p-filter
                    {
                        double vin = level[channel];
                        double delta = vin - envelope_state[channel];
                        double pp = smooth::Logistics(delta, k_rms); // (delta > 0) ? 1.0 : 0.0;
                        double nn = 1.0 - pp;
                        double g = (pp * qirtAtkCoef + nn * cubcRlsCoef);
                        envelope_state[channel] = g * envelope_state[channel] + (1.0 - g) * vin;
                        env = sqrt(envelope_state[channel]);
                        break;
                    }
                }
                
                // Transfer Curve ===========================================================
                env = DecibelConverter::ToDecibel(env);
                
                double overshoot = env - (threshold);
                
                if (overshoot <= -kneeHalf)
                    gain = 0.0;
                else if (overshoot > -kneeHalf && overshoot <= kneeHalf)
                    gain = 0.5 * slope * ((overshoot + kneeHalf) * (overshoot + kneeHalf)) / knee;
                else
                    gain = slope * overshoot;
                
                // Look-ahead FIR ===========================================================
                if (lookaheadEnable)
                {
                    lookAheadDelayLine[channel].push_back(gain);
                    // effectively parallelized version of the default std::inner_product
                    // lookAheadDelayLine[channel].end() - 1 - (lookaheadSize - 1)
                    gain = std::transform_reduce(LAH_coef, LAH_coef + lookaheadSize, lookAheadDelayLine[channel].end() - lookaheadSize, 0.0);
                    lookAheadDelayLine[channel].pop_front();
                }
                else
                {
                    lookAheadDelayLine[channel].push_back(gain);
                    gain = *(lookAheadDelayLine[channel].end() - 1 - lookaheadSize);
                    lookAheadDelayLine[channel].pop_front();
                }
                
                if (gain < GR_Max) GR_Max = gain;
                gain = DecibelConverter::ToGain(gain);
            }
            // Logarithmic level detection ===========================================================
            else
            {
                // Transfer Curve ===========================================================
                env = DecibelConverter::ToDecibel(sqrt(level[channel]));
                
                double overshoot = env - (threshold);
                
                if (overshoot <= -kneeHalf)
                    gain = 0.0;
                else if (overshoot > -kneeHalf && overshoot <= kneeHalf)
                    gain = 0.5 * slope * ((overshoot + kneeHalf) * (overshoot + kneeHalf)) / knee;
                else
                    gain = slope * overshoot;
                
                switch (dType) {
                    case detectorPeak: default: // Sonnox Oxford Dynamics
                    {
                        // Envelope Detection ===========================================================
                        double vin = gain;
                        double delta = vin - envelope_state[channel]; // 0.1 == -120 < 0.2 == -60
                        delta *= -1.0;
                        double pp = (delta > 0) ? 1.0 : 0.0; // smooth::Logistics(delta, k_rms); // (delta > 0) ? 1.0 : 0.0;
                        double nn = 1.0 - pp;
                        double g = (pp * atkCoef + nn * rlsCoef);
                        envelope_state[channel] = g * envelope_state[channel] + (1.0 - g) * vin;
                        // state[channel] = std::min(gain, rlsCoef * state[channel] + (1.0 - rlsCoef) * gain);
                        // envelope_state[channel] = atkCoef * envelope_state[channel] + (1.0 - atkCoef) * state[channel];
                        gain = envelope_state[channel];
                        break;
                    }
                    case detectorRMS : // Semi-true-RMS. a true RMS would use square-sum-normalize-sqrt, not this 1p-filter
                    {
                        // Envelope Detection ===========================================================
                        /*
                         gain^2 == dB * 2
                         -16dB -> 0.158489
                         -32dB -> 0.02511876312 = 0.158489 * 0.158489
                         */
                        // Why? if envelope_state is stored in linear gain, volume drops if switched to other detector
                        double vin = DecibelConverter::ToGain(gain * 2.0);
                        envelope_state[channel] = DecibelConverter::ToGain(envelope_state[channel]);
                        double delta = vin - envelope_state[channel];
                        delta *= -1.0;
                        double pp = (delta > 0) ? 1.0 : 0.0; // smooth::Logistics(delta, k_rms); // (delta > 0) ? 1.0 : 0.0;
                        double nn = 1.0 - pp;
                        double g = (pp * atkCoef + nn * sqrt(rlsCoef));
                        envelope_state[channel] = g * envelope_state[channel] + (1.0 - g) * vin;
                        envelope_state[channel] = DecibelConverter::ToDecibel(envelope_state[channel]);
                        gain = envelope_state[channel] * 0.5;
                        break;
                    }
                }
                
                // Look-ahead FIR ===========================================================
                if (lookaheadEnable)
                {
                    lookAheadDelayLine[channel].push_back(gain);
                    // effectively parallelized version of the default std::inner_product
                    // lookAheadDelayLine[channel].end() - 1 - (lookaheadSize - 1)
                    gain = std::transform_reduce(LAH_coef, LAH_coef + lookaheadSize, lookAheadDelayLine[channel].end() - lookaheadSize, 0.0);
                    lookAheadDelayLine[channel].pop_front();
                }
                else
                {
                    lookAheadDelayLine[channel].push_back(gain);
                    gain = *(lookAheadDelayLine[channel].end() - 1 - lookaheadSize);
                    lookAheadDelayLine[channel].pop_front();
                }
                
                if (gain < GR_Max) GR_Max = gain;
                gain = DecibelConverter::ToGain(gain);
            }
                
            if (scListen)
            {
                inputSample = sidechain_EQed[channel][sample];
                gain = 1.0;
                makeup = 1.0;
            }

            latencyDelayLine[channel].push_back(inputSample);
            inputSample = *(latencyDelayLine[channel].end() - 1 - lookaheadSize);
            latencyDelayLine[channel].pop_front();

            double dry = inputSample;
            
            inputSample *= gain;  // apply gain reduction
            
            inputSample *= makeup; // apply makeup gain before mix
            
            inputSample = mix * inputSample + (1.0 - mix) * dry;
            
            inputSample *= outputGain; // apply final output gain
            
            double outputPeak = VuOutputPeak[channel].processSample(inputSample);
            double outputRMS  = VuOutputRMS[channel].processSample(inputSample);
            if(maxOutputPeak[channel] < outputPeak) maxOutputPeak[channel] = outputPeak;
            if(maxOutputRMS[channel]  < outputRMS ) maxOutputRMS[channel]  = outputRMS;
            
            if (softBypass) inputSample = inputs[channel][sample];
            
            outputs[channel][sample] = (SampleType)(inputSample);
        }
        sample++;
    }
    if (numChannels == 1)
    {
        fInputVuPeak [0] = DecibelConverter::ToDecibel(maxInputPeak [0]);
        fOutputVuPeak[0] = DecibelConverter::ToDecibel(maxOutputPeak[0]);
        fInputVuRMS  [0] = DecibelConverter::ToDecibel(maxInputRMS  [0]);
        fOutputVuRMS [0] = DecibelConverter::ToDecibel(maxOutputRMS [0]);
        
        fInputVuPeak [1] = fInputVuPeak [0];
        fOutputVuPeak[1] = fOutputVuPeak[0];
        fInputVuRMS  [1] = fInputVuRMS  [0];
        fOutputVuRMS [1] = fOutputVuRMS [0];
    }
    else
    {
        fInputVuPeak [0] = DecibelConverter::ToDecibel(maxInputPeak [0]);
        fOutputVuPeak[0] = DecibelConverter::ToDecibel(maxOutputPeak[0]);
        fInputVuRMS  [0] = DecibelConverter::ToDecibel(maxInputRMS  [0]);
        fOutputVuRMS [0] = DecibelConverter::ToDecibel(maxOutputRMS [0]);
        
        fInputVuPeak [1] = DecibelConverter::ToDecibel(maxInputPeak [1]);
        fOutputVuPeak[1] = DecibelConverter::ToDecibel(maxOutputPeak[1]);
        fInputVuRMS  [1] = DecibelConverter::ToDecibel(maxInputRMS  [1]);
        fOutputVuRMS [1] = DecibelConverter::ToDecibel(maxOutputRMS [1]);
    }
    fGainReduction = GR_Max;
    return;
}
    


void RLFCMP_Processor::call_after_SR_changed ()
{
    call_after_parameter_changed ();
    
    lookaheadSize = std::min((int)(0.5 * 0.001 * projectSR), maxLAH); // fixed lookahead at 0.5ms

    for (auto& iter : lookAheadDelayLine) iter.resize(maxLAH);
    for (auto& iter : latencyDelayLine)   iter.resize(maxLAH);

    Kaiser::calcFilter2(lookaheadSize, 3.0, LAH_coef); // ((alpha * pi)/0.1102) + 8.7, alpha == 3 -> -94.22 dB
    
    ParamValue transition = 2.0 * bw / projectSR;
    ParamValue coefs[HT_order];
    hiir::PolyphaseIir2Designer::compute_coefs_spec_order_tbw (coefs, HT_order, transition);
    for (int i = 1, j = 0; i < HT_order; i += 2) HT_coefs[path_ref][j++] = coefs[i];
    for (int i = 0, j = 0; i < HT_order; i += 2) HT_coefs[path_sft][j++] = coefs[i];
    
    for (int32 channel = 0; channel < maxChannel; channel++)
    {
        VuInputRMS[channel].prepare(projectSR);
        VuOutputRMS[channel].prepare(projectSR);
        VuInputPeak[channel].prepare(projectSR);
        VuOutputPeak[channel].prepare(projectSR);
    }
}

void RLFCMP_Processor::call_after_parameter_changed ()
{
    dtrAtkCoef  = getTau(detectorAtk,  projectSR); 
    dtrRlsCoef  = getTau(detectorRls,  projectSR);
    
    atkCoef     = getTau(attack - detectorAtk,  projectSR);
    rlsCoef     = getTau(release - detectorRls, projectSR);
        
    for (int32 channel = 0; channel < maxChannel; channel++)
    {
        auto LfType = [](int t) -> int {
            switch (t) {
                case 0: default: return SVF_12::tHighPass;
                case 1: return SVF_12::tLowShelf;
            }
        };
        auto HfType = [](int t) -> int {
            switch (t) {
                case 0: default: return SVF_12::tLowPass;
                case 1: return SVF_12::tHighShelf;
            }
        };
        SC_LF[channel].setSVF(scLfIn, scLfFreq, scLfGain, M_SQRT1_2, LfType(scLfType), projectSR);
        SC_HF[channel].setSVF(scHfIn, scHfFreq, scHfGain, M_SQRT1_2, HfType(scHfType), projectSR);
    }

    slope     = 1.0 / ratio - 1.0;
    kneeHalf  = knee / 2.0;
}

//------------------------------------------------------------------------
} // namespace yg331



/*
 Hilbert Detector
 
 It consists of two paths with 90 degree phase difference.
 We can say these two paths now have a sin <-> cos relation in Real-Imaginary plane
 So if both paths are squared and added, it becomes a perfect envelope of original signal.
 
 Idealy, 90 degree phase shift of Square wave is going inf.
 However the digital method of generating two path with 90 degree phase difference is by
 generating two NEW path with 90 degree by phase all pass filters,
 not Original and Phase-shifted path.
 So, both paths are not going crazy at square input.
 
 One thing to know, is that square wave going through allpass filter does change shape of square.
 */
/*
 Envelope Detector
 
 Envelop(Level) Detector consists of Rectifier + Capacitor.
 Capacitor is the Leaky integrator, converting AC to DC by reducing the ripple.
 Attack and Release are then applied with RC circuit after this Detector.
 
 I wanted to use [ ZDF == TPT == Trapezoidal integration ] one pole filter for enevelope follower,
 but it drops to -inf at nyquist, changing harmonic response of compressor.
 So I had to use naive implementation of one pole filter.
 correct way would be oversampling by x4 it with TPT, but I'm aiming for a near-zero-latency for realtime use.
 Maybe, MZTi kernal with very small size might work, too.
 It kinda looks like high shelf, but changes center freq and gain dB continuously.
 */

/*
 What I want to do - ADAA to rectifier and coefficient switching
 
 "sqared = sideChain * sideChain" is a non-linear func and ADAA can be applied.
 Also, standard if-else is massive aliasing algo, Logistics function is little less, but not perfect.
 
 AA-IIR compensated will work good, but I havn't done math so long, I forgot ALL.
 
 Antiderivative Antialiasing with Frequency Compensation for Stateful Systems
 https://www.dafx.de/paper-archive/2022/papers/DAFx20in22_paper_4.pdf
 */
/*
 Look-ahead with FIR smoothing
 
 'proper' way, some say, to implement lookahead.
 Kaiser window for best harmonic rejection at low tap size(0.5ms, 24 taps at 48000)
 
 Also, transform_reduce is literally convolution.
 It's faster then naive or symmetrical implementation.
 Manual SIMD might be better, but not bothered and quite happy with this.
 
 Pro-C 2 uses A-symmetrical shape. does it help?
 */



/*
 
// Allpass peak interpolator - Maybe good for simple Peak detector or RMS?
if (false)
{
    double ap_state[2][4][2] = {0, };
    double ap_coef[4] = {1.0, 0.6681786379192988, 0.41421356237309503, 0.198912367379658};
    for (int path = 0; path < 4; path++)
    {
        double nextInput = inputSample;
        double ret = ap_coef[path] * (nextInput - ap_state[channel][path][1]) + ap_state[channel][path][0];
        ap_state[channel][path][0] = nextInput;
        ap_state[channel][path][1] = ret;
    }
    for (int path = 1; path < 4; path++)
    {
        if (inputSample < ap_state[channel][path][1]) inputSample = ap_state[channel][path][1];
    }
}
*/
/*
if (false)
{
    double v0 = inputSample;
    double v1 = (g * v0 + ic1eq[channel]) / (1.0 + g);
    ic1eq[channel] = 2*v1 - ic1eq[channel];
    *ptrOut = (SampleType)(v1);
}
else
{
    ic1eq[channel] = g * ic1eq[channel] + (1.0 - g) * inputSample;
    *ptrOut = (SampleType)(ic1eq[channel]);
}
if(false)
{
    double G = g/(1.0+g);
    double v = (HT_dtct - ic1eq[channel]) * G;
    double y = v + ic1eq[channel];
    ic1eq[channel] = y + v;
    *ptrOut = (SampleType)(y);
}
else
{
    double vin = HT_dtct;
    double t0 = vin - ic1eq[channel];
    double v0 = t0 / (1.0 + g);
    double t2 = g * v0; // v = (HT_dtct - ic1eq[channel]) * G, G = g/(1.0+g)
    double v2 = ic1eq[channel] + t2; // y = v + ic1eq[channel];
    ic1eq[channel] += 2.0 * t2; // ic1eq[channel] = y + v;
    *ptrOut = (SampleType)(sqrt(v2));
}
 */
/*
 typedef struct _state_variable{
     double vin = 0.0;
     double t0 = 0.0;
     double v0 = 0.0;
     double t2 = 0.0;
     double v2 = 0.0;
     double ic1eq = 0.0;
 } state_variable;
 state_variable fast_rms[2], slow_rms[2];
 // slow detector, Hilbert
 double delta = HT_dtct - slow_rms[channel].v2;
 double pp = smooth::Logistics(delta, k_HT); // (delta > 0) ? 1.0 : 0.0;
 double nn = 1.0 - pp;
 double g = (pp * slow_atk + nn * slow_rls);
 
 slow_rms[channel].vin = HT_dtct;
 slow_rms[channel].t0 = slow_rms[channel].vin - slow_rms[channel].ic1eq;
 slow_rms[channel].v0 = slow_rms[channel].t0 / (1.0 + g);
 slow_rms[channel].t2 = g * slow_rms[channel].v0;
 slow_rms[channel].v2 = slow_rms[channel].ic1eq + slow_rms[channel].t2;
 slow_rms[channel].ic1eq += 2.0 * slow_rms[channel].t2;
 
 double slow_env = sqrt(slow_rms[channel].v2);
 */
// Naive FIR convolution takes 63% of CPU
// for (int i = 0; i < lookaheadSize; i++)
//    gg += LAH_coef[i] * lookAheadDelayLine.getSample(channel, i);

// if lookaheadSize == 24 -> halfTap == 12, 0~11,     12~23, condition == 0
// if lookaheadSize == 25 -> halfTap == 12, 0~11, 12, 13~24, condition == 1


/*
 double gg = 0.0;
for (int i = 0; i < halfTap; i++)
    gg += LAH_coef[i] * (lookAheadDelayLine[channel].getSample(channel, i) + lookAheadDelayLine[channel].getSample(channel, lookaheadSize - i));
if (condition) // if lookaheadSize is odd
    gg += LAH_coef[halfTap] * lookAheadDelayLine[channel][halfTap];
 */
