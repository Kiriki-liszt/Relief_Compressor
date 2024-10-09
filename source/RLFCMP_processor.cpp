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
    addAudioInput  (STR16 ("Stereo In"), Steinberg::Vst::SpeakerArr::kStereo);
    addAudioOutput (STR16 ("Stereo Out"), Steinberg::Vst::SpeakerArr::kStereo);

    /* If you don't need an event bus, you can remove the next line */
    // addEventInput (STR16 ("Event In"), 1);
    
    // JUST for my sanity
    call_after_parameter_changed ();
    SC_LF[0].setRange(0);
    SC_LF[1].setRange(0);
    SC_HF[0].setRange(1);
    SC_HF[1].setRange(1);

    return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Processor::terminate ()
{
    // Here the Plug-in will be de-instantiated, last possibility to remove some memory!
    for (auto& iter : latencyDelayLine)
    {
        delete iter;
    }
    
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
                        case kParamScLfIn:          pScLfIn     = (value > 0.5); break;
                        case kParamScLfType:        pScLfType   = value; break;
                        case kParamScLfFreq:        pScLfFreq   = value; break;
                        case kParamScLfGain:        pScLfGain   = value; break;
                        case kParamScHfIn:          pScHfIn     = (value > 0.5); break;
                        case kParamScHfType:        pScHfType   = value; break;
                        case kParamScHfFreq:        pScHfFreq   = value; break;
                        case kParamScHfGain:        pScHfGain   = value; break;
                        case kParamScListen:        pScListen   = (value > 0.5); break;
                        case kParamType:            pType       = value; break;
                        case kParamAttack:          pAttack     = value; break;
                        case kParamRelease:         pRelease    = value; break;
                        case kParamLookaheadEnable: pLookaheadEnable = (value > 0.5); break;
                        case kParamThreshold:       pThreshold  = value; break;
                        case kParamRatio:           pRatio      = value; break;
                        case kParamKnee:            pKnee       = value; break;
                        case kParamMakeup:          pMakeup     = value; break;
                        case kParamMix:             pMix        = value; break;
                        case kParamInput:           pInput      = value; break;
                        case kParamOutput:          pOutput     = value; break;
                        case kParamSoftBypass:      pSoftBypass = (value > 0.5); break;
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
    sendFloat("Input L", Input_L);
    sendFloat("Input R", Input_R);
    sendFloat("Output L", Output_L);
    sendFloat("Output R", Output_R);
    sendFloat("Gain Reduction", Gain_Reduction);
    sendFloat("detectorIndicator", detectorIndicator);

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
    
    Vst::SpeakerArrangement arr;
    getBusArrangement (Steinberg::Vst::BusDirections::kOutput, 0, arr);
    auto numChannels = Vst::SpeakerArr::getChannelCount (arr);
    lookaheadSize = std::min((int)(0.5 * 0.001 * newSetup.sampleRate), maxLAH); // fixed lookahead at 0.5ms
    halfTap = lookaheadSize / 2;
    condition = lookaheadSize % 2;
    
    lookAheadDelayLine.resize(numChannels);
    for (auto& iter : lookAheadDelayLine)
    {
        iter = new std::deque<double>;
        iter->resize(maxLAH);
    }
    
    latencyDelayLine.resize(numChannels);
    for (auto& iter : latencyDelayLine)
    {
        iter = new double[maxLAH];
    }

    Kaiser::calcFilter2(lookaheadSize, 3.0, LAH_coef); // ((alpha * pi)/0.1102) + 8.7, alpha == 3 -> -94.22 dB
    
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
    ParamValue GR_Max = 0.0;
    
    // Dunno why, but lets Auto-Vectorization of transform_reduce
    int lookAhead_local = 0.5 * 0.001 * SampleRate;
    if (lookAhead_local > maxLAH) lookAhead_local = maxLAH;
    
    int32 sample = 0;
    
    ParamValue In_Max[2] = {0.0, };
    ParamValue Out_Max[2] = {0.0, };
    
    while (sampleFrames > sample)
    {
        double sideChain[2] = {0.0, 0.0};
        double HT_dtct[2] = {0.0, 0.0};
        double sqared[2]  = {0.0, 0.0};
        double rectified[2]  = {0.0, 0.0};
        
        for (int32 channel = 0; channel < numChannels; channel++)
        {
            Vst::Sample64 inputSample = inputs[channel][sample];
            
            /*
             input gain
             optional gain stage for headroom
             */
            
            inputSample *= input;
            
            if(In_Max[channel] < inputSample) In_Max[channel] = inputSample;
            
            /*
             Pre Gain
             Acts like auto gain, by appling abs gain of threshold
             */
            inputSample *= preGain;
            
            /*
             Sidechain Filter - BS1770 filter
             
             It simulates how human ear percives sound.
             It helps compressor to behave more natual (to human ear).
             */
            sideChain[channel] = SC_LF[channel].computeSVF(inputSample);
            sideChain[channel] = SC_HF[channel].computeSVF(sideChain[channel]);
            
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
            for (int path = 0; path < path_num; path++) // two paths for 0 and +90
            {
                double nextInput = sideChain[channel]; // spits weird values if abs or squaring input
                for (int stage = 0; stage < HT_stage; stage++)
                {
                    double ret = HT_coefs[path][stage] * (nextInput + state[channel][path][io_y][y2][stage]) - state[channel][path][io_x][x2][stage];
                    state[channel][path][io_x][x2][stage] = state[channel][path][io_x][x1][stage];
                    state[channel][path][io_x][x1][stage] = nextInput;
                    state[channel][path][io_y][y2][stage] = state[channel][path][io_y][y1][stage];
                    state[channel][path][io_y][y1][stage] = ret;
                    nextInput = ret;
                }
            }
            double rms_s   = state[channel][path_ref][io_y][y2][HT_stage - 1];
            double rms_c   = state[channel][path_sft][io_y][y1][HT_stage - 1];
            HT_dtct[channel] = rms_s * rms_s + rms_c * rms_c; // Ideal input level, squared
            sqared[channel] = sideChain[channel] * sideChain[channel]; // input squared
            rectified[channel] = std::abs(sideChain[channel]);
        }
        
        double monoSample = 0.0;
        for (int32 channel = 0; channel < numChannels; channel++)
            monoSample += HT_dtct[channel];
        monoSample /= numChannels;
        for (int32 channel = 0; channel < numChannels; channel++)
            HT_dtct[channel] = monoSample;
        
        monoSample = 0.0;
        for (int32 channel = 0; channel < numChannels; channel++)
            monoSample += sqared[channel];
        monoSample /= numChannels;
        for (int32 channel = 0; channel < numChannels; channel++)
            sqared[channel] = monoSample;
        
        monoSample = 0.0;
        for (int32 channel = 0; channel < numChannels; channel++)
            monoSample += rectified[channel];
        monoSample /= numChannels;
        for (int32 channel = 0; channel < numChannels; channel++)
            rectified[channel] = monoSample;
        
        for (int32 channel = 0; channel < numChannels; channel++)
        {
            Vst::Sample64 inputSample = inputs[channel][sample];
            inputSample *= input;
            inputSample *= preGain;
            
            /*
             Envelope Detector
             
             Envelop(Level) Detector consists of Rectifier + Capacitor.
             Capacitor is the Leaky integrator, converting AC to DC by reducing the ripple.
             Attack and Release are then applied with RC circuit after this Detector.
             
             I wanted to use [ ZDF == TPT == Trapezoidal integration ] one pole filter for enevelope follower,
             but it drops to -inf at nyquist, changing harmonic response of compressor.
             So I had to use naive implementation of one pole filter.
             correct way would be oversampling by x4 it with TPT, but I'm aiming for a Zero-latency for realtime use.
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
             Question is, typical compressors have effective threshold to be inconsistant by time constants.
             The only compressor that I know works consistantly is Sonnox Oxford Dynamics by Paul Frindle.
             suggestion is that this is done by two-stage detector, like analog compressors.
             first stage acts as zero-attack moderate-release, and second stage sets actual attack/release.
             */
            
            double env = 0.0;
            
            switch (type) {
                case detectorType::Bold :
                {
                    double vin = rectified[channel];
                    // fast detector, plain rms
                    double delta = vin - detector_state[channel];
                    double pp = smooth::Logistics(delta, k_rms); // (delta > 0) ? 1.0 : 0.0;
                    double nn = 1.0 - pp;
                    double g = (pp * atkCoef + nn * rlsCoef);
                    detector_state[channel] = g * detector_state[channel] + (1.0 - g) * vin;
                    env = (detector_state[channel]);
                }
                    break;
                    
                case detectorType::Smooth :
                {
                    double vin = sqared[channel];
                    // fast detector, plain rms
                    double delta = vin - detector_state[channel];
                    double pp = smooth::Logistics(delta, k_rms); // (delta > 0) ? 1.0 : 0.0;
                    double nn = 1.0 - pp;
                    double g = (pp * sqrt(atkCoef) + nn * (rlsCoef * rlsCoef));
                    detector_state[channel] = g * detector_state[channel] + (1.0 - g) * vin;
                    env = sqrt(detector_state[channel]);
                }
                    break;
                    
                case detectorType::Clean :
                {
                    double vin = HT_dtct[channel];
                    // slow detector, Hilbert
                    double delta = vin - detector_state[channel];
                    double pp = smooth::Logistics(delta, k_HT); // (delta > 0) ? 1.0 : 0.0;
                    double nn = 1.0 - pp;
                    double g = (pp * sqrt(atkCoef) + nn * (rlsCoef * rlsCoef));
                    detector_state[channel] = g * detector_state[channel] + (1.0 - g) * vin;
                    env = sqrt(detector_state[channel]);
                }
                    break;
                    
                default:
                {
                    double vin = HT_dtct[channel];
                    // slow detector, Hilbert
                    double delta = vin - detector_state[channel];
                    double pp = smooth::Logistics(delta, k_HT); // (delta > 0) ? 1.0 : 0.0;
                    double nn = 1.0 - pp;
                    double g = (pp * sqrt(atkCoef) + nn * (rlsCoef * rlsCoef));
                    detector_state[channel] = g * detector_state[channel] + (1.0 - g) * vin;
                    env = sqrt(detector_state[channel]);
                }
                    break;
            }
            
            env = DecibelConverter::ToDecibel(env);
            
            double overshoot = env - (threshold);
            
            double gain = 1.0;
            
            if (overshoot <= -kneeHalf)
                gain = 0.0;
            else if (overshoot > -kneeHalf && overshoot <= kneeHalf)
                gain = 0.5 * slope * ((overshoot + kneeHalf) * (overshoot + kneeHalf)) / knee;
            else
                gain = slope * overshoot;
            
            /*
             Look-ahead with FIR smoothing
             
             'proper' way, some say, to implement lookahead.
             Kaiser window for best harmonic rejection at low tap size(0.5ms, 24 taps at 48000)
             
             Also, transform_reduce is literally convolution.
             It's faster then naive or symmetrical implementation.
             Manual SIMD might be better, but not bothered and quite happy with this.
             
             Pro-C 2 uses A-symmetrical shape. does it help?
             */
            if (pLookaheadEnable)
            {
                lookAheadDelayLine[channel]->push_back(gain);
                gain = std::transform_reduce(lookAheadDelayLine[channel]->end() - lookAhead_local, lookAheadDelayLine[channel]->end(), LAH_coef, 0.0);
                lookAheadDelayLine[channel]->pop_front();
            }
            else
            {
                lookAheadDelayLine[channel]->push_back(gain);
                gain = *(lookAheadDelayLine[channel]->end() - lookAhead_local);
                lookAheadDelayLine[channel]->pop_front();
            }
            
            if (gain < GR_Max) GR_Max = gain;
            gain = DecibelConverter::ToGain(gain);
            
            if (pScListen)
            {
                inputSample = sideChain[channel];
                inputSample /= preGain;
                gain = 1.0;
                makeup = 1.0;
            }

            // latencyDelayLine[channel]->pushSample(inputSample);
            latencyDelayLine[channel][writePos] = inputSample;
            writePos = (writePos + maxLAH - 1) % maxLAH;
            // inputSample = latencyDelayLine[channel]->getSample(lookAhead_local);
            int index = (readPos + lookAhead_local) % maxLAH;
            inputSample = latencyDelayLine[channel][index];
            readPos = (readPos + maxLAH - 1) % maxLAH;

            double dry = inputSample;
            
            inputSample *= gain;
            
            inputSample *= makeup;
            
            inputSample = mix * inputSample + (1.0 - mix) * dry;
            
            inputSample *= output;
            
            if (Out_Max[channel] < inputSample) Out_Max[channel] = inputSample;
            
            outputs[channel][sample] = (SampleType)(inputSample);
        }
        sample++;
    }
    if (numChannels == 1)
    {
        Input_L = DecibelConverter::ToDecibel(In_Max[0]);
        Output_L = DecibelConverter::ToDecibel(Out_Max[0]);
        Input_R = Input_L;
        Output_R = Output_L;
    }
    if (numChannels == 2)
    {
        Input_L = DecibelConverter::ToDecibel(In_Max[0]);
        Output_L = DecibelConverter::ToDecibel(Out_Max[0]);
        Input_R =  DecibelConverter::ToDecibel(In_Max[1]);
        Output_R = DecibelConverter::ToDecibel(Out_Max[1]);
    }
    Gain_Reduction = GR_Max;
    return;
}
//------------------------------------------------------------------------
} // namespace yg331

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
