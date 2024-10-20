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
                        case kParamBypass:          pBypass     = (value > 0.5); break;
                        // case kParamZoom:       pZoom       = value; break;
                        // case kParamOS:         pOS         = value; break;
                        case kParamScLfIn:            pScLfIn     = (value > 0.5); break;
                        case kParamScLfType:          pScLfType   = value; break;
                        case kParamScLfFreq:          pScLfFreq   = value; break;
                        case kParamScLfGain:          pScLfGain   = value; break;
                        case kParamScHfIn:            pScHfIn     = (value > 0.5); break;
                        case kParamScHfType:          pScHfType   = value; break;
                        case kParamScHfFreq:          pScHfFreq   = value; break;
                        case kParamScHfGain:          pScHfGain   = value; break;
                        case kParamScListen:          pScListen   = (value > 0.5); break;
                        case kParamDetectorType:      pDType      = value; break;
                        case kParamSidechainTopology: pSCTopology = value; break;
                        case kParamLookaheadEnable:   pLookaheadEnable = (value > 0.5); break;
                        case kParamHilbertEnable:     pHilbertEnable   = (value > 0.5); break;
                        case kParamAttack:            pAttack     = value; break;
                        case kParamRelease:           pRelease    = value; break;
                        case kParamThreshold:         pThreshold  = value; break;
                        case kParamRatio:             pRatio      = value; break;
                        case kParamKnee:              pKnee       = value; break;
                        case kParamMakeup:            pMakeup     = value; break;
                        case kParamMix:               pMix        = value; break;
                        case kParamInput:             pInput      = value; break;
                        case kParamOutput:            pOutput     = value; break;
                        case kParamSoftBypass:        pSoftBypass = (value > 0.5); break;
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
        if (pBypass)
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
    // This happens BEFORE setState
    // fprintf (stdout, "setupProcessing\n");
    
    projectSR = newSetup.sampleRate;
    
    if      (projectSR > 96000.0) internalSR = projectSR;
    else if (projectSR > 48000.0) internalSR = projectSR * 2.0;
    else                          internalSR = projectSR * 4.0;
    
    call_after_SR_changed ();
    call_after_parameter_changed (); // in case setState does not happen
    
    //--- called before any processing ----
    return AudioEffect::setupProcessing (newSetup);
}

uint32  PLUGIN_API RLFCMP_Processor::getLatencySamples()
{
    return lookaheadSize-1;
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
tresult PLUGIN_API RLFCMP_Processor::setState (IBStream* state)
{
    // called when we load a preset, the model has to be reloaded
    if (!state)
        return kResultFalse;
    
    // fprintf (stdout, "setState\n");
    IBStreamer streamer (state, kLittleEndian);
    
    // SAVE IN PLAIN VALUE

    int32           savedBypass     = 0;
    // Vst::ParamValue savedOS         = 0.0;
    int32           savedSidechainFilter = 0.0;
    Vst::ParamValue savedAttack     = 0.0;
    Vst::ParamValue savedRelease    = 0.0;
    Vst::ParamValue savedBias       = 0.0;
    Vst::ParamValue savedThreshold  = 0.0;
    Vst::ParamValue savedRatio      = 0.0;
    Vst::ParamValue savedKnee       = 0.0;
    Vst::ParamValue savedMakeup     = 0.0;
    Vst::ParamValue savedMix        = 0.0;
    Vst::ParamValue savedOutput     = 0.0;
    int32           savedSoftBypass = 0.0;
    
    if (streamer.readInt32 (savedBypass)     == false) savedBypass     = 0;
    // if (streamer.readDouble(savedOS)         == false) savedOS         = 0.0;
    if (streamer.readDouble(savedAttack)     == false) savedAttack     = paramAttack.ToNormalized(dftAttack);
    if (streamer.readDouble(savedRelease)    == false) savedRelease    = paramRelease.ToNormalized(dftRelease);
    if (streamer.readDouble(savedThreshold)  == false) savedThreshold  = paramThreshold.ToNormalized(dftThreshold);
    if (streamer.readDouble(savedRatio)      == false) savedRatio      = paramRatio.ToNormalized(dftRatio);
    if (streamer.readDouble(savedKnee)       == false) savedKnee       = paramKnee.ToNormalized(dftKnee);
    if (streamer.readDouble(savedMakeup)     == false) savedMakeup     = paramMakeup.ToNormalized(dftMakeup);
    if (streamer.readDouble(savedMix)        == false) savedMix        = paramMix.ToNormalized(dftMix);
    if (streamer.readDouble(savedOutput)     == false) savedOutput     = paramOutput.ToNormalized(dftOutput);
    if (streamer.readInt32 (savedSoftBypass) == false) savedSoftBypass = 0;
    
    pBypass = savedBypass > 0;
    // pOS = savedOS;
    pAttack     = savedAttack;
    pRelease    = savedRelease;
    pThreshold  = savedThreshold;
    pRatio      = savedRatio;
    pKnee       = savedKnee;
    pMakeup     = savedMakeup;
    pMix        = savedMix;
    pOutput     = savedOutput;
    pSoftBypass = savedSoftBypass > 0;
    
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
    
    streamer.writeInt32(pBypass ? 1 : 0);
    // streamer.writeDouble(Steinberg::ToNormalized<ParamValue> (static_cast<ParamValue>(pOS), overSample_num));
    streamer.writeDouble(pAttack);
    streamer.writeDouble(pRelease);
    streamer.writeDouble(pThreshold);
    streamer.writeDouble(pRatio);
    streamer.writeDouble(pKnee);
    streamer.writeDouble(pMakeup);
    streamer.writeDouble(pMix);
    streamer.writeDouble(pOutput);
    streamer.writeInt32(pSoftBypass ? 1 : 0);
    
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
    if (sidechain_EQed[0].size() != sampleFrames) sidechain_EQed[0].resize(sampleFrames);
    if (sidechain_EQed[1].size() != sampleFrames) sidechain_EQed[1].resize(sampleFrames);
    
    double invNumChannels = 1.0 / numChannels;
    ParamValue GR_Max = 0.0;
    
    // Dunno why, but uses 0.1 less cpu[PluginDoctor]
    int lookAhead_local = 0.5 * 0.001 * SampleRate;
    if (lookAhead_local > maxLAH) lookAhead_local = maxLAH;
    
    int32 sample = 0;
    
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

    //double K = tan(M_PI * 0.15 / SampleRate); //lowpass
    double vvv = sqrt(getTau(1.0, SampleRate));
    
    // because of locallity, this is faster
    for (int32 channel = 0; channel < numChannels; channel++)
    {
        for (int sample = 0; sample < sampleFrames; sample++)
        {
            // SideChain Filtering ===========================================================
            //DC_state_y[channel] = (1.0 - K) * DC_state_y[channel] + inputSample - DC_state_x[channel];
            //DC_state_x[channel] = inputSample;
            //sideChain[channel] = DC_state_y[channel];
            double t = SC_LF[channel].computeSVF(inputs[channel][sample]);
            t = SC_HF[channel].computeSVF(t);
            sidechain_EQed[channel][sample] = t;
        }
    }

    // Process ===========================================================
    while (sampleFrames > sample)
    {
        double sideChain[2] = {0.0, 0.0};
        double HT_dtct[2]   = {0.0, 0.0};
        double squared[2]   = {0.0, 0.0};
        double rectified[2] = {0.0, 0.0};
        
        for (int32 channel = 0; channel < numChannels; channel++)
        {
            Vst::Sample64 inputSample = sidechain_EQed[channel][sample];

            inputSample *= input;

            if (pHilbertEnable) {
                // Hilbert detector ===========================================================
                for (int path = 0; path < path_num; path++)
                {
                    double nextInput = inputSample;
                    for (int stage = 0; stage < HT_stage; stage++)
                    {
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
                HT_dtct[channel] = rms_s * rms_s + rms_c * rms_c; // Ideal input level, squared
                
                double vin = HT_dtct[channel];
                double delta = vin - rectified_state[channel];
                double pp = (delta > 0) ? 1.0 : 0.0; // smooth::Logistics(delta, k_detector); // (delta > 0) ? 1.0 : 0.0;
                double nn = 1.0 - pp;
                //double g = (pp * 0.0 + nn * powrDtrRlsCoef);
                //double g = (pp * sqrt(getTau(0.1, SampleRate)) + nn * powrDtrRlsCoef);
                double g = vvv; //sqrt(getTau(1.0, SampleRate));
                hilbert_state[channel] = g * hilbert_state[channel] + (1.0 - g) * vin;
                HT_dtct[channel] = hilbert_state[channel];
                // std::hypot(x, y) = sqrt(x^2 + y^2)
            }
            switch (dType) {
                case detectorPeak: default:
                {
                    // Full wave rectifier ===========================================================
                    double vin = std::abs(inputSample);
                    double delta = vin - rectified_state[channel];
                    double pp = (delta >= 0) ? 1.0 : 0.0; // smooth::Logistics(delta, k_detector); // (delta > 0) ? 1.0 : 0.0;
                    double nn = 1.0 - pp;
                    double g = (pp * 0.0 + nn * dtrRlsCoef);
                    rectified_state[channel] = g * rectified_state[channel] + (1.0 - g) * vin;
                    rectified[channel] = rectified_state[channel];
                    break;
                }
                case detectorRMS :
                {
                    // Square wave rectifier ===========================================================
                    double vin = inputSample * inputSample;
                    double delta = vin - squared_state[channel];
                    double pp = (delta >= 0) ? 1.0 : 0.0; // smooth::Logistics(delta, k_detector*k_detector); // (delta > 0) ? 1.0 : 0.0;
                    double nn = 1.0 - pp;
                    double g = (pp * 0.0 + nn * powrDtrRlsCoef);
                    squared_state[channel] = g * squared_state[channel] + (1.0 - g) * vin;
                    squared[channel] = squared_state[channel];
                    break;
                }
            }
        }
        
        // Link stereo ===========================================================
        double boldMonoSample   = 0.0;
        double smoothMonoSample = 0.0;
        double cleanMonoSample  = 0.0;
        for (int32 channel = 0; channel < numChannels; channel++)
        {
            boldMonoSample   += invNumChannels * rectified[channel];
            smoothMonoSample += invNumChannels * squared[channel];
            cleanMonoSample  += invNumChannels * HT_dtct[channel];
        }
        for (int32 channel = 0; channel < numChannels; channel++)
        {
            rectified[channel] = boldMonoSample;
            squared[channel]   = smoothMonoSample;
            HT_dtct[channel]   = cleanMonoSample;
        }
        
        for (int32 channel = 0; channel < numChannels; channel++)
        {
            Vst::Sample64 inputSample = inputs[channel][sample];
            inputSample *= input;
            
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
                    case detectorPeak: default:  // It matches Metric Halo ChannelStrip MIO Comp, and Weiss DS1-MK3 if atk*3 & rls*2
                    {
                        // double r = getTau(paramRelease.ToPlain(pRelease) - detectorRls, SampleRate);
                        double vin = rectified[channel];
                        double delta = vin - envelope_state[channel];
                        double pp = smooth::Logistics(delta, k_rms); // (delta > 0) ? 1.0 : 0.0;
                        double nn = 1.0 - pp;
                        double g = (pp * sqrtAtkCoef + nn * squrRlsCoef);
                        envelope_state[channel] = g * envelope_state[channel] + (1.0 - g) * vin;
                        env = (envelope_state[channel]);
                        break;
                    }
                    case detectorRMS: // Semi-true-RMS. a true RMS would use square-sum-normalize-sqrt, not this 1p-filter
                    {
                        double vin = squared[channel];
                        if (pHilbertEnable) vin = HT_dtct[channel];
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
                if (pLookaheadEnable)
                {
                    lookAheadDelayLine[channel].push_back(gain);
                    gain = std::transform_reduce(lookAheadDelayLine[channel].end() - lookAhead_local, lookAheadDelayLine[channel].end(), LAH_coef, 0.0);
                    lookAheadDelayLine[channel].pop_front();
                }
                else
                {
                    lookAheadDelayLine[channel].push_back(gain);
                    gain = *(lookAheadDelayLine[channel].end() - lookAhead_local);
                    lookAheadDelayLine[channel].pop_front();
                }
                
                if (gain < GR_Max) GR_Max = gain;
                gain = DecibelConverter::ToGain(gain);
            }
            // Logarithmic level detection ===========================================================
            else
            {
                switch (dType) {
                    case detectorPeak: default: // It matches Sonnox Oxford Dynamics, Cenozoix if atk*2, rls*2
                    {
                        // Transfer Curve ===========================================================
                        env = DecibelConverter::ToDecibel(rectified[channel]);
                        
                        double overshoot = env - (threshold);
                        
                        if (overshoot <= -kneeHalf)
                            gain = 0.0;
                        else if (overshoot > -kneeHalf && overshoot <= kneeHalf)
                            gain = 0.5 * slope * ((overshoot + kneeHalf) * (overshoot + kneeHalf)) / knee;
                        else
                            gain = slope * overshoot;
                        
                        // Envelope Detection ===========================================================
                        double vin = gain;
                        double delta = vin - envelope_state[channel]; // 0.1 == -120 < 0.2 == -60
                        delta *= -1.0;
                        double pp = smooth::Logistics(delta, k_rms); // (delta > 0) ? 1.0 : 0.0;
                        double nn = 1.0 - pp;
                        double g = (pp * atkCoef + nn * rlsCoef);
                        envelope_state[channel] = g * envelope_state[channel] + (1.0 - g) * vin;
                        gain = (envelope_state[channel]);
                        break;
                    }
                    case detectorRMS : // Ozone RMS
                    {
                        // Transfer Curve ===========================================================
                        env = DecibelConverter::ToDecibel(sqrt(squared[channel]));
                        
                        double overshoot = env - (threshold);
                        
                        if (overshoot <= -kneeHalf)
                            gain = 0.0;
                        else if (overshoot > -kneeHalf && overshoot <= kneeHalf)
                            gain = 0.5 * slope * ((overshoot + kneeHalf) * (overshoot + kneeHalf)) / knee;
                        else
                            gain = slope * overshoot;
                        
                        // Envelope Detection ===========================================================
                        gain = DecibelConverter::ToGain(gain);
                        gain *= gain;
                        // gain = DecibelConverter::ToDecibel(gain);
                        double vin = gain;
                        double delta = vin - envelope_state[channel];
                        delta *= -1.0;
                        double pp = smooth::Logistics(delta, k_rms); // (delta > 0) ? 1.0 : 0.0;
                        double nn = 1.0 - pp;
                        double g = (pp * atkCoef + nn * sqrt(rlsCoef));
                        envelope_state[channel] = g * envelope_state[channel] + (1.0 - g) * vin;
                        //gain = DecibelConverter::ToDecibel(sqrt(DecibelConverter::ToGain(envelope_state[channel])));
                        gain = DecibelConverter::ToDecibel(sqrt(envelope_state[channel]));
                        break;
                    }
                }
                
                // Look-ahead FIR ===========================================================
                if (pLookaheadEnable)
                {
                    lookAheadDelayLine[channel].push_back(gain);
                    gain = std::transform_reduce(lookAheadDelayLine[channel].end() - lookAhead_local, lookAheadDelayLine[channel].end(), LAH_coef, 0.0);
                    lookAheadDelayLine[channel].pop_front();
                }
                else
                {
                    lookAheadDelayLine[channel].push_back(gain);
                    gain = *(lookAheadDelayLine[channel].end() - lookAhead_local);
                    lookAheadDelayLine[channel].pop_front();
                }
                
                if (gain < GR_Max) GR_Max = gain;
                gain = DecibelConverter::ToGain(gain);
            }
                
            if (pScListen)
            {
                inputSample = sidechain_EQed[channel][sample];
                gain = 1.0;
                makeup = 1.0;
            }

            latencyDelayLine[channel].push_back(inputSample);
            inputSample = *(latencyDelayLine[channel].end() - lookAhead_local);
            latencyDelayLine[channel].pop_front();

            double dry = inputSample;
            
            inputSample *= gain;
            
            inputSample *= makeup;
            
            inputSample = mix * inputSample + (1.0 - mix) * dry;
            
            inputSample *= output;
            
            double outputPeak = VuOutputPeak[channel].processSample(inputSample);
            double outputRMS  = VuOutputRMS[channel].processSample(inputSample);
            if(maxOutputPeak[channel] < outputPeak) maxOutputPeak[channel] = outputPeak;
            if(maxOutputRMS[channel]  < outputRMS ) maxOutputRMS[channel]  = outputRMS;
            
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
    lookaheadSize = std::min((int)(0.5 * 0.001 * projectSR), maxLAH); // fixed lookahead at 0.5ms
    halfTap       = lookaheadSize / 2;
    condition     = lookaheadSize % 2;
    
    for (auto& iter : lookAheadDelayLine)
        iter.resize(maxLAH);
    
    for (auto& iter : latencyDelayLine)
        iter.resize(maxLAH);

    Kaiser::calcFilter2(lookaheadSize, 3.0, LAH_coef); // ((alpha * pi)/0.1102) + 8.7, alpha == 3 -> -94.22 dB
    
    dtrAtkCoef  = getTau(detectorAtk,  projectSR);
    dtrRlsCoef  = getTau(detectorRls,  projectSR);
    
    atkCoef    = getTau(paramAttack. ToPlain(pAttack),  projectSR);
    rlsCoef    = getTau(paramRelease.ToPlain(pRelease), projectSR);
    
    for (int32 channel = 0; channel < maxChannel; channel++)
    {
        // DC_state_x[channel] = 0.0;
        // DC_state_y[channel] = 0.0;
        
        int LfType = paramScLfType.ToPlainList(pScLfType);
        switch (LfType) {
            case 0: LfType = SVF_12::tHighPass; break;
            case 1: LfType = SVF_12::tLowShelf; break;
            default: break;
        }
        int HfType = paramScHfType.ToPlainList(pScHfType);
        switch (HfType) {
            case 0: HfType = SVF_12::tLowPass; break;
            case 1: HfType = SVF_12::tHighShelf; break;
            default: break;
        }
        SC_LF[channel].setSVF(pScLfIn, paramScLfFreq.ToPlain(pScLfFreq), paramScLfGain.ToPlain(pScLfGain), M_SQRT1_2, LfType, projectSR);
        SC_HF[channel].setSVF(pScHfIn, paramScHfFreq.ToPlain(pScHfFreq), paramScHfGain.ToPlain(pScHfGain), M_SQRT1_2, HfType, projectSR);
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
    dType      = paramDetectorType.ToPlainList(pDType);
    scTopology = paramSidechainTopology.ToPlainList(pSCTopology);
    
    atkCoef    = getTau(paramAttack. ToPlain(pAttack),  projectSR);
    rlsCoef    = getTau(paramRelease.ToPlain(pRelease), projectSR);
    
    for (int32 channel = 0; channel < maxChannel; channel++)
    {
        int LfType = paramScLfType.ToPlainList(pScLfType);
        switch (LfType) {
            case 0: LfType = SVF_12::tHighPass; break;
            case 1: LfType = SVF_12::tLowShelf; break;
            default: break;
        }
        int HfType = paramScHfType.ToPlainList(pScHfType);
        switch (HfType) {
            case 0: HfType = SVF_12::tLowPass; break;
            case 1: HfType = SVF_12::tHighShelf; break;
            default: break;
        }
        SC_LF[channel].setSVF(pScLfIn, paramScLfFreq.ToPlain(pScLfFreq), paramScLfGain.ToPlain(pScLfGain), M_SQRT1_2, LfType, projectSR);
        SC_HF[channel].setSVF(pScHfIn, paramScHfFreq.ToPlain(pScHfFreq), paramScHfGain.ToPlain(pScHfGain), M_SQRT1_2, HfType, projectSR);
    }

    ratio     = paramRatio.ToPlain(pRatio);
    slope     = 1.0 / ratio - 1.0;
    knee      = paramKnee.ToPlain(pKnee);
    kneeHalf  = knee / 2.0;
    threshold = paramThreshold.ToPlain(pThreshold);
    makeup    = DecibelConverter::ToGain(paramMakeup.ToPlain(pMakeup));
    
    mix       = pMix;
    input     = DecibelConverter::ToGain(paramInput. ToPlain(pInput));
    output    = DecibelConverter::ToGain(paramOutput.ToPlain(pOutput));
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
