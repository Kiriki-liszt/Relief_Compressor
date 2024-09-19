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
                        case kParamBypass:     pBypass     = (value > 0.5); break;
                        case kParamSoftBypass: pSoftBypass = (value > 0.5); break;
                        // case kParamZoom:       pZoom       = value; break;
                        case kParamOS:         pOS         = value; break;
                        case kParamSidechainFilter: pSidechainFilter = value; break;
                        case kParamAttack:     pAttack     = value; break;
                        case kParamRelease:    pRelease    = value; break;
                        case kParamBias:       pBias       = value; break;
                        case kParamThreshold:  pThreshold  = value; break;
                        case kParamRatio:      pRatio      = value; break;
                        case kParamKnee:       pKnee       = value; break;
                        case kParamMakeup:     pMakeup     = value; break;
                        case kParamMix:        pMix        = value; break;
                        case kParamOutput:     pOutput     = value; break;
                        
                        default: break;
                    }
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

	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Processor::setupProcessing (Vst::ProcessSetup& newSetup)
{
    transition = 2*10.0/newSetup.sampleRate;
    FDebugPrint("newSetup.sampleRate = %f\n", newSetup.sampleRate);
    hiir::PolyphaseIir2Designer::compute_coefs_spec_order_tbw (coefs, numCoefs, transition);
    // Phase reference path c coefficients
    for (int i = 1, j = 0; i < numCoefs; i += 2) {
        c[1][j++] = coefs[i];
        FDebugPrint("Ref | c[0][j++] = %f\n", coefs[i]);
    }
    // +90 deg path c coefficients
    for (int i = 0, j = 0; i < numCoefs; i += 2) {
        c[0][j++] = coefs[i];
    }
    
	//--- called before any processing ----
	return AudioEffect::setupProcessing (newSetup);
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Processor::canProcessSampleSize (int32 symbolicSampleSize)
{
	// by default kSample32 is supported
	if (symbolicSampleSize == Vst::kSample32)
		return kResultTrue;

	// disable the following comment if your processing support kSample64
	/* if (symbolicSampleSize == Vst::kSample64)
		return kResultTrue; */

	return kResultFalse;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Processor::setState (IBStream* state)
{
	// called when we load a preset, the model has to be reloaded
	IBStreamer streamer (state, kLittleEndian);
	
	return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Processor::getState (IBStream* state)
{
	// here we need to save the model
	IBStreamer streamer (state, kLittleEndian);

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
    for (int32 channel = 0; channel < numChannels; channel++)
    {
        int32 samples = sampleFrames;

        SampleType* ptrIn  = (SampleType*)inputs[channel];
        SampleType* ptrOut = (SampleType*)outputs[channel];

        while (--samples >= 0)
        {
            Vst::Sample64 inputSample = *ptrIn;

            // two paths for 0 and +90
            for (int path = 0; path < 2; path++)
            {
                double nextInput = inputSample;
                for (int stage = 0; stage < numCoefs/2; stage++)
                {
                    double ret = c[path][stage] * (nextInput + state[channel][path][1][1][stage]) - state[channel][path][0][1][stage];
                    state[channel][path][0][1][stage] = state[channel][path][0][0][stage];
                    state[channel][path][0][0][stage] = nextInput;
                    state[channel][path][1][1][stage] = state[channel][path][1][0][stage];
                    state[channel][path][1][0][stage] = ret;
                    nextInput = ret;
                }
            }
            
            Steinberg::Vst::Sample64 tmp1 = abs(inputSample) * abs(inputSample);
            
            double env;
            if (channel == 0)
            {
                peak[channel] = (peak[channel] * coeff) + (icoef * tmp1);
                env = sqrt(peak[channel]) * M_SQRT2;
                rms_sin[channel] = (rms_sin[channel] * coeff) + (icoef * state[channel][1][1][1][numCoefs/2 - 1] * state[channel][1][1][1][numCoefs/2 - 1]);
                rms_cos[channel] = (rms_cos[channel] * coeff) + (icoef * state[channel][0][1][0][numCoefs/2 - 1] * state[channel][0][1][0][numCoefs/2 - 1]);
                env = sqrt(rms_sin[channel] + rms_cos[channel]);
            }
            else
            {
                rms_sin[channel] = (rms_sin[channel] * coeff) + (icoef * state[channel][1][1][1][numCoefs/2 - 1] * state[channel][1][1][1][numCoefs/2 - 1]);
                rms_cos[channel] = (rms_cos[channel] * coeff) + (icoef * state[channel][0][1][0][numCoefs/2 - 1] * state[channel][0][1][0][numCoefs/2 - 1]);
                env = sqrt(rms_sin[channel] + rms_cos[channel]);
            }
            
            env = 20 * log10(env);
            double gain = 1.0;
            if (env > -15.0) // env == -9
            {
                gain = env - (-15.0);
                double slope = 1.0 / 4.0 - 1.0;
                gain *= slope;
                gain = pow(10, gain * 0.05);
            }

            *ptrOut = (SampleType)(inputSample * gain);
            
            ptrIn++;
            ptrOut++;
        }
    }
    return;
}
//------------------------------------------------------------------------
} // namespace yg331
