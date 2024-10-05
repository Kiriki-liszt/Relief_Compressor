//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "pluginterfaces/vst/vsttypes.h"
// #include "pluginterfaces/base/futils.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <algorithm>
#include <vector>
#include <queue>
#include <numeric>
#include <complex>

namespace yg331 {
//------------------------------------------------------------------------
using ParamValue = Steinberg::Vst::ParamValue;
using SampleRate = Steinberg::Vst::SampleRate;
using int32      = Steinberg::int32;
using uint32     = Steinberg::uint32;

typedef enum {
    overSample_1x,
    overSample_2x,
    overSample_4x,
    overSample_8x,
    overSample_num = 3
} overSample;

//------------------------------------------------------------------------
//  Class for converter
//------------------------------------------------------------------------
class DecibelConverter
{
public:
    static ParamValue ToGain (ParamValue dB)
    {
        return std::pow(10.0, dB * 0.05);
    }

    static ParamValue ToDecibel (ParamValue gain)
    {
        return (gain > 0.0) ? (20.0 * std::log10(gain)) : (-100.0);
    }
private:
    DecibelConverter() = delete;
};

class ParameterConverter
{
public:
    enum paramType {
        range = 0,
        log,
        list
    };
    
    ParameterConverter (ParamValue minValue, ParamValue maxValue, paramType type = paramType::range, int32 numSteps = -1)
    : minValue (minValue), maxValue (maxValue), type(type), numSteps (numSteps) {};

    ParamValue ToNormalized (ParamValue plain) const
    {
        switch (type) {
            case paramType::range: return ToNormalizedRange(plain); break;
            case paramType::log:   return ToNormalizedLog(plain); break;
            case paramType::list:  return ToNormalizedList(plain); break;
            default: break;
        }
        return ToNormalizedRange(plain);
    }

    ParamValue ToPlain (ParamValue normalized) const
    {
        switch (type) {
            case paramType::range: return ToPlainRange(normalized); break;
            case paramType::log:   return ToPlainLog(normalized); break;
            case paramType::list:  return ToPlainList(normalized); break;
            default: break;
        }
        return ToNormalizedRange(normalized);
    }
    
    ParamValue ToNormalizedRange (ParamValue plain) const
    {
        if (maxValue == minValue) return 0.0;
        return (plain - minValue) / (maxValue - minValue);
    }

    ParamValue ToPlainRange (ParamValue normalized) const
    {
        return normalized * (maxValue - minValue) + minValue;
    }
    
    /* Log Conversion
     
      x - x0    log(y) - log(y0)
     ------- = -----------------
     x1 - x0   log(y1) - log(y0)
     
     x0 = 0, x1 = 1 : normalized parameter
     
     norm  = ( log(plain) - log(minValue) ) / ( log(maxValue) - log(minValue) )
           = log(plain / minValue) / log(maxValue / minValue)
     
     log(plain / minValue) = norm * log(maxValue / minValue)
     plain / minValue = exp(norm * log(maxValue / minValue))
     plain = minValue * exp(norm * log(maxValue / minValue))
     
     */
    
    ParamValue ToNormalizedLog (ParamValue plain) const
    {
        if (minValue == 0.0) return 0.0;
        if (!((plain    / minValue) > 0.0)) return 0.0;
        if (!((maxValue / minValue) > 0.0)) return 0.0;
        return std::log(plain / minValue) / std::log(maxValue / minValue);
    }

    ParamValue ToPlainLog (ParamValue normalized) const
    {
        if (minValue == 0.0) return 0.0;
        if (!((maxValue / minValue) > 0.0)) return 0.0;
        return minValue * std::exp(normalized * std::log(maxValue / minValue));
    }
    
    ParamValue ToNormalizedList (int32 plain) const
    {
        // return (Steinberg::ToNormalized<type> (v, step));
        if (numSteps <= 0) return 0;
        return plain / ParamValue (numSteps);
    }

    int32 ToPlainList (ParamValue normalized) const
    {
        // return (Steinberg::FromNormalized<type> (v, step));
        if (numSteps <= 0) return 0;
        int32 cmp(normalized * (numSteps + 1));
        return numSteps < cmp ? numSteps : cmp;
    }

    void setRange (ParamValue _min, ParamValue _max)
    {
        minValue = _min;
        maxValue = _max;
    }
    
    void setSteps (int32 _steps)
    {
        numSteps = _steps;
    }
    
private:
    ParamValue minValue = 0.0;
    ParamValue maxValue = 1.0;
    paramType type = paramType::range;
    int32 numSteps = -1;
};

// GML 8900 ~ MDWDRC's method is not so straight-forward
// As patent shows, attack comes from exponent parameter('Exponent'), and adds corrective release('Timing') after
// Looks like some Airwindows type of algo...
// maybe, he was right after all

// GML 8900 focuses on Attack blend - Release does not change while blending - for dynamic range control
// Weiss DS1 focuses on Release blend - Attack does nor change while blending - for getting the right feel when mastering

// GML has fast attack / fast release and slow attack / slow release
// So, fast release cannot take over slow release
// resulting 2-stage attack and fixed slow release.
// Meaning, It generally smooth out dynamic range, but able to catch fast peak if happened
// truly a dynamic range contoller

// Weiss has same attack, but has fast release and slow release
// resulting 1-stage attack and 2-stage release.
// Meaning, fast GR will start to recover fast, and slow down as returning
// natural? program dependant? auto release - like.

// I think 2-RMS detectors might work better

// Attack - Release - Bias - Blend ?

// Attack is in Log, but Release is in Linear, like MH-CS and GML

// GML style   : Attack - Release - NULL -Bias
// Weiss style : Attack - Fast R - Slow R - AVG

//------------------------------------------------------------------------
//  Min, Max, Default of Parameters
//------------------------------------------------------------------------
static constexpr ParamValue dftBypass          = 0.0;
static constexpr ParamValue dftSoftBypass      = 0.0;
static constexpr ParamValue dftSidechainFilter = 1.0;

static constexpr ParamValue minAttack    = 1.2;
static constexpr ParamValue maxAttack    = 70.0;
static constexpr ParamValue dftAttack    = 8.0;

static constexpr ParamValue minRelease   = 10.0;
static constexpr ParamValue maxRelease   = 2000.0;
static constexpr ParamValue dftRelease   = 80.0;

static constexpr ParamValue minBias      = -18.0;
static constexpr ParamValue maxBias      = 18.0;
static constexpr ParamValue dftBias      = 0.0;

static constexpr ParamValue minThreshold = -40.0;
static constexpr ParamValue maxThreshold = 0.0;
static constexpr ParamValue dftThreshold = -20.0;

static constexpr ParamValue minRatio     = 1.0;
static constexpr ParamValue maxRatio     = 20.0;
static constexpr ParamValue dftRatio     = 4.0;

static constexpr ParamValue minKnee      = 0.0;
static constexpr ParamValue maxKnee      = 20.0;
static constexpr ParamValue dftKnee      = 5.0;

static constexpr ParamValue minMakeup    = 0.0;
static constexpr ParamValue maxMakeup    = 24.0;
static constexpr ParamValue dftMakeup    = 0.0;

static constexpr ParamValue minMix       = 0.0;
static constexpr ParamValue maxMix       = 100.0;
static constexpr ParamValue dftMix       = 100.0;

static constexpr ParamValue minOutput    = -12.0;
static constexpr ParamValue maxOutput    = 12.0;
static constexpr ParamValue dftOutput    = 0.0;

static const ParameterConverter paramAttack   (minAttack,    maxAttack,    ParameterConverter::paramType::log);
static const ParameterConverter paramRelease  (minRelease,   maxRelease,   ParameterConverter::paramType::log);
static const ParameterConverter paramBias     (minBias,      maxBias,      ParameterConverter::paramType::range);
static const ParameterConverter paramThreshold(minThreshold, maxThreshold, ParameterConverter::paramType::range);
static const ParameterConverter paramRatio    (minRatio,     maxRatio,     ParameterConverter::paramType::range);
static const ParameterConverter paramKnee     (minKnee,      maxKnee,      ParameterConverter::paramType::range);
static const ParameterConverter paramMakeup   (minMakeup,    maxMakeup,    ParameterConverter::paramType::range);
static const ParameterConverter paramMix      (minMix,       maxMix,       ParameterConverter::paramType::range);
static const ParameterConverter paramOutput   (minOutput,    maxOutput,    ParameterConverter::paramType::range);

//------------------------------------------------------------------------
//  Hilbert Transform Designer, modified
//------------------------------------------------------------------------
namespace hiir {
class PolyphaseIir2Designer {
public:
    static void    compute_coefs_spec_order_tbw (double coef_arr [], int nbr_coefs, double transition) noexcept
    {
        double         k;
        double         q;
        compute_transition_param (k, q, transition);
        const int      order = nbr_coefs * 2 + 1;

        for (int index = 0; index < nbr_coefs; ++index)
        {
            coef_arr [index] = compute_coef (index, k, q, order);
        }
    };

private:
    static void    compute_transition_param (double &k, double &q, double transition) noexcept
    {
        k  = tan ((1 - transition * 2) * M_PI / 4);
        k *= k;

        double         kksqrt = pow (1 - k * k, 0.25);
        const double   e = 0.5 * (1 - kksqrt) / (1 + kksqrt);
        const double   e2 = e * e;
        const double   e4 = e2 * e2;
        q = e * (1 + e4 * (2 + e4 * (15 + 150 * e4)));
    };
    static double  compute_coef (int index, double k, double q, int order) noexcept
    {
        const int      c    = index + 1;
        const double   num  = compute_acc_num (q, order, c) * pow (q, 0.25);
        const double   den  = compute_acc_den (q, order, c) + 0.5;
        const double   ww   = num / den;
        const double   wwsq = ww * ww;

        const double   x    = sqrt ((1 - wwsq * k) * (1 - wwsq / k)) / (1 + wwsq);
        const double   coef = (1 - x) / (1 + x);

        return coef;
    };
    static double  compute_acc_num (double q, int order, int c) noexcept
    {
        int            i   = 0;
        int            j   = 1;
        double         acc = 0;
        double         q_ii1;
        do
        {
            q_ii1  = pow (q, i * (i + 1)); // ipowp
            q_ii1 *= sin ((i * 2 + 1) * c * M_PI / order) * j;
            acc   += q_ii1;

            j = -j;
            ++i;
        }
        while (fabs (q_ii1) > 1e-150);

        return acc;
    };
    static double  compute_acc_den (double q, int order, int c) noexcept
    {
        int            i   =  1;
        int            j   = -1;
        double         acc =  0;
        double         q_i2;
        do
        {
            q_i2  = pow (q, i * i); // ipowp
            q_i2 *= cos (i * 2 * c * M_PI / order) * j;
            acc  += q_i2;

            j = -j;
            ++i;
        }
        while (fabs (q_i2) > 1e-150);

        return acc;
    };

/*\\\ FORBIDDEN MEMBER FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/

private:

                   PolyphaseIir2Designer ()  = delete;
                   ~PolyphaseIir2Designer () = delete;
                   PolyphaseIir2Designer (const PolyphaseIir2Designer &other) = delete;
                   PolyphaseIir2Designer (PolyphaseIir2Designer &&other)      = delete;
    PolyphaseIir2Designer &
                   operator = (const PolyphaseIir2Designer &other)  = delete;
    PolyphaseIir2Designer &
                   operator = (PolyphaseIir2Designer &&other)       = delete;
    bool           operator == (const PolyphaseIir2Designer &other) = delete;
    bool           operator != (const PolyphaseIir2Designer &other) = delete;

}; // class PolyphaseIir2Designer
}  // namespace hiir

//------------------------------------------------------------------------
//  Kaiser-Bessel Filter designser for Oversampling
//------------------------------------------------------------------------
class Kaiser {
public:
    static constexpr int32 maxTap = 512;
    
    static inline double Ino(double x)
    {
        double d = 0, ds = 1, s = 1;
        do
        {
            d += 2;
            ds *= x * x / (d * d);
            s += ds;
        } while (ds > s * 1e-50);
        return s;
    };

    static void calcFilter(double Fs, double Fa, double Fb, int M, double Att, double dest [])
    {
        // Kaiser windowed FIR filter "DIGITAL SIGNAL PROCESSING, II" IEEE Press pp 123-126.

        int Np = (M - 1) / 2;
        double A[maxTap] = { 0, };
        double Alpha; // == pi * alpha (wikipedia)
        double Inoalpha;

        A[0] = 2 * (Fb - Fa) / Fs;

        for (int j = 1; j <= Np; j++)
            A[j] = (sin(2.0 * j * M_PI * Fb / Fs) - sin(2.0 * j * M_PI * Fa / Fs)) / (j * M_PI);

        if (Att < 21.0)
            Alpha = 0;
        else if (Att > 50.0)
            Alpha = 0.1102 * (Att - 8.7);
        else
            Alpha = 0.5842 * pow((Att - 21), 0.4) + 0.07886 * (Att - 21);

        Inoalpha = Ino(Alpha);

        for (int j = 0; j <= Np; j++)
        {
            dest[Np + j] = A[j] * Ino(Alpha * std::sqrt(1.0 - ((double)(j * j) / (double)(Np * Np)))) / Inoalpha;
        }
        dest[Np + Np] = A[Np] * Ino(0.0) / Inoalpha; // ARM with optimizer level O3 returns NaN == sqrt(1.0 - n/n), while x64 does not...
        for (int j = 0; j < Np; j++)
        {
            dest[j] = dest[M - 1 - j];
        }
    }
    
    static void calcFilter2(int length, double alpha, double dest [])
    {
        int N = length - 1;
        float Alpha = M_PI * alpha;
        float Inoalpha;

        Inoalpha = Ino(Alpha);

        for (int n = 0; n <= N; n++)
        {
            dest[n] = Ino(Alpha * sqrt(1.0 - (2.0 * (double)n / (double)N - 1.0) * (2.0 * (double)n / (double)N - 1.0))) / Inoalpha;
        }
        dest[0] = Ino(0.0) / Inoalpha; // ARM with optimizer level O3 returns NaN == sqrt(1.0 - n/n), while x64 does not...
        dest[N] = Ino(0.0) / Inoalpha;

        float pwr = 0.0;
        for (int i = 0; i < length; i++) {
            pwr += dest[i];
        }
        pwr = 1.0 / pwr;
        for (int i = 0; i < length; i++) {
            dest[i] *= pwr;
        }
    }
};

class Flt {
public:
    double coef alignas(16)[Kaiser::maxTap] = { 0, };
    double buff alignas(16)[Kaiser::maxTap] = { 0, };
    int now = 0;
    int size = Kaiser::maxTap;
    double* buff_ptr = buff;
    void acc() {
        now++;
        if (now == size) {
            now = 0;
        }
    }
    int get_nth(int n) {
        if (now + n >= size) {
            return now + n - size;
        }
        else {
            return now + n;
        }
    }
};

//------------------------------------------------------------------------
//  Some functions for smoothing corners
//------------------------------------------------------------------------
class smooth {
public:
    static inline double ABS(double in, double smoothness)
    {
        double err = sqrt(1.0 + smoothness * smoothness) - 1.0;
        return sqrt((in * in) + smoothness * smoothness) - err; //sqrt(smoothness);
    }

    static inline double Max(double a, double b, double smoothness)
    {
        return (a + b + ABS(a - b, smoothness)) * 0.5;
    }

    static inline double Logistics(double delta, double k)
    {
        // double init = 1.0 / (1.0 + exp(k * -delta)); // This was way to CPU intensive
        return 0.5 + 0.5 * tanh(k * delta); // Still CPU heavy, but less then exp()
    }
};

/*CrestFactor Class:
* Calculates the average Crest-Factor for a given buffer.
* Crest-Factor is time-variable value calculated from the ratio between peak and rms of the signal
*/
class CrestFactor
{
public:

    CrestFactor() = default;

    // Prepares processor with ProcessSpec-Object and recalculates coefficients for current ballistics
    void prepare(const double& fs) {
        sampleRate = fs;
        //Calculate alpha for release time of 200ms, same release time for peak & rms detector
        a1 = exp(-1.0 / (sampleRate * 0.2));
        b1 = 1 - a1;
    };

    // Calculates Crest-Factor for given buffer
    void process(const double* src, int numSamples) {
        //Init accumulators
        if (!peakState) peakState = src[0];
        if (!rmsState) rmsState = src[0];

        //Reset avg attack/release
        avgAttackTime = 0.0;
        avgReleaseTime = 0.0;

        //Calculate averages of auto - attack/release times for a single buffer
        for (int i = 0; i < numSamples; i++)
        {
            //Square of input signal
            const double s = static_cast<double>(src[i]) * static_cast<double>(src[i]);

            //Update peak state
            peakState =  std::max (s, a1 * peakState + b1 * s);

            //Update rms state
            rmsState = a1 * rmsState + b1 * s;

            //calculate squared crest factor
            const double c = peakState / rmsState;
            cFactor = c > 0.0 ? c : 0.0;

            //calculate ballistics
            if (cFactor > 0.0)
            {
                attackTimeInSeconds = 2 * (maxAttackTime / cFactor);
                releaseTimeInSeconds = 2 * (maxReleaseTime / cFactor) - attackTimeInSeconds;

                //Update avg ballistics
                avgAttackTime += attackTimeInSeconds;
                avgReleaseTime += releaseTimeInSeconds;
            }
        }

        // Calculate average ballistics & crest factor
        avgAttackTime /= numSamples;
        avgReleaseTime /= numSamples;
    };

    // Get average calculated attack time of a buffer, call after proces()
    double getAvgAttack() { return avgAttackTime; };

    // Get average calculated release time of a buffer, call after process()
    double getAvgRelease() { return avgReleaseTime; };

private:
    double attackTimeInSeconds{ 0.0 };
    double releaseTimeInSeconds{ 0.14 };
    double avgAttackTime{ 0.0 };
    double avgReleaseTime{ 0.14 };

    double peakState{ 0.0 };
    double rmsState{ 0.0 };
    double a1{ 0.0 }, b1{ 0.0 };
    double sampleRate{ 0.0 };
    double maxAttackTime{ 0.08 }, maxReleaseTime{ 1.0 }; //respective 8ms and 1sec
    double cFactor{ 0.0 };
};


/*Simple exponential moving average filter, also known as 1-pole iir filter
* This class can be used to smooth values over a certain time frame
*/
class SmoothingFilter
{
public:

    SmoothingFilter() = default;

    // Prepares the SmoothingFilter with a sampleRate
    void prepare(const double& fs) {
        sampleRate = fs;
        a1 = 1;
        b1 = 1 - a1;
    };

    // Processes a given sample
    void process(const double& sample) {
        if (first)
        {
            state = sample;
            first = false;
        }
        state = a1 * sample + b1 * state;
    };

    // Sets coefficient manually
    void setAlpha(double a) {
        a1 = a;
        b1 = 1 - a1;
    };

    // Set time-frame in seconds, recalculates needed coefficients
    void setAlphaWithTime(double timeInSeconds) {
        a1 = exp(-1.0 / (sampleRate * timeInSeconds));
        b1 = 1 - a1;
    };

    // Gets current value
    double getState() {
        return state;
    };

private:
    double a1{ 1.0 }, b1{ 0.0 };
    double state{ 0.0 };
    double sampleRate{ 0.0 };
    bool first{ true };
};


/*LevelDetector Class:
* Used to have a smooth representation of the level Might be used in linear or log.
* Domain In this compressor implementation it's used in log.
* Domain after the gain computer to smooth the calculated attenuations,
* Therefore the detector does not have to work on the whole dynamic range of the input signal
*/
class LevelDetector
{
public:
    LevelDetector() = default;

    // Prepares LevelDetector with a ProcessSpec-Object containing samplerate, blocksize and number of channels
    void prepare(const double& fs) {
        sampleRate = fs;
        crestFactor.prepare(fs);
        attackSmoothingFilter.prepare(fs);
        releaseSmoothingFilter.prepare(fs);

        alphaAttack = exp(-1.0 / (sampleRate * attackTimeInSeconds));
        alphaRelease = exp(-1.0 / (sampleRate * releaseTimeInSeconds));
        state01 = 0.0;
        state02 = 0.0;
    };

    // Sets attack time constant
    void setAttack(const double& attack) {
        if (attack != attackTimeInSeconds)
        {
            attackTimeInSeconds = attack; //Time it takes to reach 1-1/e = 0.63
            alphaAttack = exp(-1.0 / (sampleRate * attackTimeInSeconds)); //aA = e^(-1/TA*fs)

            double w = 1 / (2.0 * attackTimeInSeconds * sampleRate);
            g_a = tan(w);
            gt0_a = 1 / (1 + g_a);
            gt1_a = g_a * gt0_a;
        }
    };

    // Sets release time constant
    void setRelease(const double& release) {
        if (release != releaseTimeInSeconds)
        {
            releaseTimeInSeconds = release; //Time it takes to reach 1 - (1-1/e) = 0.37
            alphaRelease = exp(-1.0 / (sampleRate * releaseTimeInSeconds)); //aR = e^(-1/TR*fs)

            // a1 = d, tau = (sampleRate * timeInSeconds), f = 1 / 2pi*tau
            //double tau = releaseTimeInSeconds;
            //double f = 1 / 2 * M_PI * tau;
            //double w = f * M_PI / sampleRate;
            // w = 2 * pi * f = 1 / tau
            double w = 1 / (2.0 * releaseTimeInSeconds * sampleRate);
            // g = tan (pi * f / sample rate)

            g_r = tan(w);
            gt0_r = 1 / (1 + g_r);
            gt1_r = g_r * gt0_r;
        }
    };

    // Sets auto attack to enabled/disabled
    void setAutoAttack(bool isEnabled) {  autoAttack = isEnabled; };

    // Sets auto release to enabled/disabled
    void setAutoRelease(bool isEnabled) { autoRelease = isEnabled; };

    // Gets current attack time constant
    double getAttack() {  return attackTimeInSeconds; };

    // Gets current release time constant
    double getRelease() {  return releaseTimeInSeconds; };

    // Gets calculated attack coefficient
    double getAlphaAttack() { return alphaAttack; };

    // gets calculated release coefficient
    double getAlphaRelease() { return alphaRelease; };

    // Processes a sample with smooth branched peak detector
    // 1. Update state with coef branch, instead of switching state.
    // 2. Smooth transition on switching branches.
    double processPeakBranched(const double& in) {
        
        double d = state01 - in;
        double pp = smooth::Logistics(d, 100.0); // (delta > 0) ? 1.0 : 0.0;
        pp *= pp;
        double nn = 1.0 - pp;
        state01 = (pp * alphaAttack + nn * alphaRelease) * state01
                + (1 - (pp * alphaAttack + nn * alphaRelease)) * in;
        return static_cast<double>(state01); //y_L
        

        /*
        double delta = v2 - in;
        double p = smooth::Step(delta, 100.0); // (delta > 0) ? 1.0 : 0.0;
        p *= p;
        double n = 1.0 - p;
        v2 = (p * gt1_a + n * gt1_r) * in
           + (p * gt0_a + n * gt0_r) * ic1eq;
        ic1eq += 2.0 * (p * g_a + n * g_r) * (in - v2);
        //v2 = gt1 * in + gt0 * ic1eq;
        //ic1eq += 2.0 * g * (in - v2);
        return static_cast<double>(v2);
        */

        //Smooth branched peak detector
        if (in < state01)
            state01 = alphaAttack * state01 + (1 - alphaAttack) * in;
        else
            state01 = alphaRelease * state01 + (1 - alphaRelease) * in;

        return static_cast<double>(state01); //y_L
    };

    // Processes a sample with smooth decoupled peak detector
    double processPeakDecoupled(const double& in) {

        const double input = static_cast<double>(in);
        /*
        v1 = (std::max)(input, gt1_r * input + gt0_r * ic2eq);
        ic2eq += 2.0 * g_r * (input - v1);
        v2 = gt1_a * v1 + gt0_a * ic1eq;
        ic1eq += 2.0 * g_a * (v1 - v2);
        return static_cast<double>(v2); //y_L
        */
        //float smooth_max(float a, float b) { return (a + b + smooth_ABS(a - b)) * 0.5f; };
        //Smooth decoupled peak detector
        //const double input = static_cast<double>(in);
        state02 = (std::max)(input, alphaRelease * state02 + (1 - alphaRelease) * input);
        state01 = alphaAttack * state01 + (1 - alphaAttack) * state02;
        return static_cast<double>(state01);
    };

    // Applies ballistics to given buffer
    void applyBallistics(double* src, int numSamples) {
        // Apply ballistics to src buffer
        for (int i = 0; i < numSamples; i++) {
            //src[i] = processPeakDecoupled(src[i]);
            src[i] = processPeakBranched(src[i]);
        }
    };

    // Processes crest factor and sets ballistics accordingly
    void processCrestFactor(const double* src, int numSamples) {
        if (autoAttack || autoRelease)
        {
            //Crest factor calculation
            crestFactor.process(src, numSamples);
            attackSmoothingFilter.process(crestFactor.getAvgAttack());
            releaseSmoothingFilter.process(crestFactor.getAvgRelease());
            if (autoAttack) setAttack(attackSmoothingFilter.getState());
            if (autoRelease) setRelease(releaseSmoothingFilter.getState());
        }
    };

private:
    CrestFactor crestFactor;
    SmoothingFilter attackSmoothingFilter;
    SmoothingFilter releaseSmoothingFilter;

    double attackTimeInSeconds{ 0.01 }, alphaAttack{ 0.0 };
    double releaseTimeInSeconds{ 0.14 }, alphaRelease{ 0.0 };
    double state01{ 0.0 }, state02{ 0.0 };
    double sampleRate{ 0.0 };
    bool autoAttack{ false };
    bool autoRelease{ false };

    double g_a = 0.5, g_r = 0.5;
    double gt0_a = 1 / (g_a), gt0_r = 1 / (g_r);
    double gt1_a = g_a * gt0_a;
    double gt1_r = g_r * gt0_r;
    double v1 = 0, ic1eq = 0;
    double v2 = 0, ic2eq = 0;
};

/* GainComputer Class:
 * Calculates the needed attenuation to compress a signal with given characteristics
 */
class GainComputer
{
public:

    GainComputer() {
        threshold = -20.0;
        ratio = 2.0;
        slope = 1.0 / ratio - 1.0;
        knee = 6.0;
        kneeHalf = 3.0;
    };

    // Sets the threshold in dB
    void setThreshold(double newTreshold)
    {
        threshold = newTreshold;
    };

    // Sets the ratio in dB
    void setRatio(double newRatio)
    {
        if (ratio != newRatio)
        {
            ratio = newRatio;
            if (ratio > 23.9) ratio = -std::numeric_limits<double>::infinity();
            slope = 1.0 / newRatio - 1.0;
        }
    };

    // Sets the knee-width in dB (if > 0, 2nd order interpolation for soft knee)
    void setKnee(double newKnee)
    {
        if (newKnee != knee)
        {
            knee = newKnee;
            kneeHalf = newKnee / 2.0;
        }
    };

    // Applies characteristics to a given sample, 2nd order spline interpolation
    // returns attenuation Xl == Xg - Yg
    double applyCompression(double& input)
    {
        
        const double overshoot = input - threshold;

        if (overshoot <= -kneeHalf)
            return 0.0;
        if (overshoot > -kneeHalf && overshoot <= kneeHalf)
            return 0.5 * slope * ((overshoot + kneeHalf) * (overshoot + kneeHalf)) / knee;

        return slope * overshoot;
    };

    void applyCompressionToBuffer(double* src, int numSamples)
    {
        for (int i = 0; i < numSamples; i++)
        {
            //const double level = std::max(abs(src[i]), 1e-6);
            double levelInDecibels = DecibelConverter::ToDecibel(src[i]);
            src[i] = applyCompression(levelInDecibels);
        }
    };

private:
    double threshold{ -20.0f };
    double ratio{ 2.0f };
    double knee{ 6.0f }, kneeHalf{ 3.0f };
    double slope{ -0.5f };
};


/* Compressor-Class:
 * The circruit is modeled after the "ideal" VCA-Compressor
 * based on the paper "Digital Dynamic Range Compressor Design �  Tutorial and Analysis"
 * by Giannoulis, Massberg & Reiss
 */
class Compressor
{
public:
    Compressor() = default;
    ~Compressor() {
        rawSidechainSignal = nullptr;
    };

    // Prepares compressor with a ProcessSpec-Object containing samplerate, blocksize and number of channels
    void prepare(const double& sampleRate, const double& maximumBlockSize) {
        _sampleRate = sampleRate;
        _maximumBlockSize = maximumBlockSize;
        ballistics.prepare(2 * sampleRate);
        originalSignal[0].resize(maximumBlockSize, 0.0);
        originalSignal[1].resize(maximumBlockSize, 0.0);
        sidechainSignal.resize(maximumBlockSize, 0.0);
        rawSidechainSignal = sidechainSignal.data();
        smoothedAutoMakeup.prepare(sampleRate);
        smoothedAutoMakeup.setAlpha(0.03);

        for (int ch = 0; ch < 2; ch++) {
            up[ch].size = fir_size;
            dn[ch].size = fir_size;
            Kaiser::calcFilter(96000.0, 0.0, 24000.0, fir_size, 180.0, up[ch].coef);
            Kaiser::calcFilter(96000.0, 0.0, 24000.0, fir_size, 180.0, dn[ch].coef);
            for (int i = 0; i < fir_size; i++) up[ch].coef[i] *= 2.0;
        }
        upsample_buff[0].resize(2 * maximumBlockSize, 0.0);
        upsample_buff[1].resize(2 * maximumBlockSize, 0.0);
        dnsample_buff[0].resize(2 * maximumBlockSize, 0.0);
        dnsample_buff[1].resize(2 * maximumBlockSize, 0.0);
    };

    // Sets compressor to bypassed/not bypassed
    void setPower(bool newPower) { bypassed = newPower; };

    // Sets input in dB
    void setInput(double newInput) { input = newInput; };

    // Sets threshold in dB
    void setThreshold(double thresholdInDb) { gainComputer.setThreshold(thresholdInDb); };

    // Sets ratio in dB
    void setRatio(double rat) { gainComputer.setRatio(rat); };

    // Sets knee-width in dB (> 0 = soft knee)
    void setKnee(double kneeInDb) { gainComputer.setKnee(kneeInDb); };


    void setLookAhead(double _LookAhead) { lookaheadDelay = _LookAhead; };
    void setLookAheadEnable(bool _LookAheadEnable) { LookAheadEnable = _LookAheadEnable; };

    // Sets make-up gain in dB
    void setMakeup(double makeupGainInDb) { makeup = makeupGainInDb; };

    // Sets mix 0.0f - 1.0f
    void setMix(double newMix) { mix = newMix; };

    // Sets attack time in milliseconds
    void setAttack(double attackTimeInMs) { ballistics.setAttack(attackTimeInMs * 0.001); };

    // Sets release time in milliseconds
    void setRelease(double releaseTimeInMs) { ballistics.setRelease(releaseTimeInMs * 0.001); };

    // Sets auto attack to enabled = true or disabled = false
    void setAutoAttack(bool isEnabled) { ballistics.setAutoAttack(isEnabled); };

    // Sets auto release to enabled = true or disabled = false
    void setAutoRelease(bool isEnabled) { ballistics.setAutoRelease(isEnabled); };

    // Sets auto makeup to enabled = true or disabled = false
    void setAutoMakeup(bool isEnabled) { autoMakeupEnabled = isEnabled; };

    // Gets current make-up gain value
    double getMakeup() { return makeup; };

    // Return current sampleRate
    double getSampleRate() { return _sampleRate; };

    double getMaxGainReduction() { return maxGainReduction; };
    
    
    
    
    // Processes input sample
    template <typename SampleType>
    void process(
        SampleType* inputs,
        SampleType* outputs,
        Steinberg::int32 channel)
    {
        // BS1770 Filter = pre-filter + revised low-frequency B curve
        // ABS, smoothed
        // path a. fast RMS detector, using traditional rms
        // path b. slow RMS detector, using Hilbert transform
        // upsample
        // envelope detect to both path, log attack, linear release
        // compute controller gain
        // downsample
        // apply gain
    };
    

    // Processes input buffer
    template <typename SampleType>
    void process(
        SampleType** inputs,
        SampleType** outputs,
        Steinberg::int32 numChannels,
        Steinberg::Vst::SampleRate getSampleRate,
        Steinberg::int32 sampleFrames)
    {
        if (!bypassed)
        {
            // Clear any old samples
            for (int i = 0; i < sampleFrames; i++) {
                originalSignal[0][0] = 0.0;
                originalSignal[1][0] = 0.0;
                rawSidechainSignal[i] = 0.0f;
            }
            maxGainReduction = 0.0f;

            // Apply input gain
            for (int ch = 0; ch < numChannels; ch++) applyInputGain(inputs[ch], sampleFrames);

            // Get max l/r amplitude values and fill sidechain signal
            for (int ch = 0; ch < numChannels; ch++) {
                for (int i = 0; i < sampleFrames; i++) {
                    // Multi-phase all-pass interpolation
                    // https://www.dsprelated.com/freebooks/pasp/First_Order_Allpass_Interpolation.html
                    // p1_sy[ch] = 0.198912367379658 * inputs[ch][i] + p1_sx[ch] - 0.198912367379658 * p1_sy[ch]; // 15000/48000
                    // p1_sx[ch] = inputs[ch][i];
                    // p2_sy[ch] = 0.41421356237309503 * inputs[ch][i] + p2_sx[ch] - 0.41421356237309503 * p2_sy[ch]; // 18000/48000
                    // p2_sx[ch] = inputs[ch][i];
                    // p3_sy[ch] = 0.6681786379192988 * inputs[ch][i] + p3_sx[ch] - 0.6681786379192988 * p3_sy[ch]; // 21000/48000
                    // p3_sx[ch] = inputs[ch][i];

                    //double a = std::max(abs((double)inputs[ch][i]), abs(p1_sx[ch]));
                    //double b = std::max(abs(p2_sy[ch]), abs(p3_sy[ch]));
                    //double c = std::max(abs(p4_sy[ch]), abs(p5_sy[ch]));
                    //double d = std::max(a, std::max(b, c));
                    //double d = std::max(a, b);
                    
                    //rawSidechainSignal[i] = (rawSidechainSignal[i] > d) ? rawSidechainSignal[i] : d;
                    
                    
                    rawSidechainSignal[i] = (rawSidechainSignal[i] > abs((double)inputs[ch][i])) ? rawSidechainSignal[i] : (double)inputs[ch][i];
                    //rawSidechainSignal[i] = std::max(abs(rawSidechainSignal[i]), 1e-6);
                    //rawSidechainSignal[i] = abs(rawSidechainSignal[i]);
                    //rawSidechainSignal[i] = smooth::ABS(rawSidechainSignal[i], 0.0001);

                    //outputs[ch][i] = (SampleType)rawSidechainSignal[i];
                    //outputs[ch][i] = std::abs(inputs[ch][i]);
                }
            }
            if (false) {
                for (int i = 0; i < sampleFrames; i++) {
                    // two paths for 0 and +90
                    for (int p = 0; p < 2; p++) {
                        double ret1 = c[p][0] * (rawSidechainSignal[i] + p1_sy_1[p]) - p1_sx_1[p];
                        p1_sx_1[p] = p1_sx_2[p];
                        p1_sx_2[p] = rawSidechainSignal[i];
                        p1_sy_1[p] = p1_sy_2[p];
                        p1_sy_2[p] = ret1;
                        double ret2 = c[p][1] * (ret1 + p2_sy_1[p]) - p2_sx_1[p];
                        p2_sx_1[p] = p2_sx_2[p];
                        p2_sx_2[p] = ret1; // inputs[ch][i];
                        p2_sy_1[p] = p2_sy_2[p];
                        p2_sy_2[p] = ret2;
                        double ret3 = c[p][2] * (ret2 + p3_sy_1[p]) - p3_sx_1[p];
                        p3_sx_1[p] = p3_sx_2[p];
                        p3_sx_2[p] = ret2; // inputs[ch][i];
                        p3_sy_1[p] = p3_sy_2[p];
                        p3_sy_2[p] = ret3;
                        double ret4 = c[p][3] * (ret3 + p4_sy_1[p]) - p4_sx_1[p];
                        p4_sx_1[p] = p4_sx_2[p];
                        p4_sx_2[p] = ret3; // inputs[ch][i];
                        p4_sy_1[p] = p4_sy_2[p];
                        p4_sy_2[p] = ret4;
                    }
                    double env = sqrt(p4_sy_1[1] * p4_sy_1[1] + p4_sy_2[0] * p4_sy_2[0]); // *0.637?
                    rawSidechainSignal[i] = env;
                }
            }
            //return;

            // upsampling
            for (int ch = 0; ch < 1; ch++) {
                for (int i = 0; i < sampleFrames; i++) {
                    // OS - upsample
                    up[ch].acc();
                    up[ch].buff[up[ch].now] = (double)rawSidechainSignal[i];
                    //up[ch].buff[up[ch].now] = (double)inputs[ch][i];
                    double in0 = 0.0;
                    for (int i = 0; i < tap_hm; i++) {
                        double a = up[ch].buff[up[ch].get_nth(i)];
                        double b = up[ch].buff[up[ch].get_nth(fir_size - 1 - i)];
                        in0 += up[ch].coef[i] * (a + b);
                    }   in0 += up[ch].coef[tap_hm] * up[ch].buff[up[ch].get_nth(tap_hm)];
                    up[ch].acc();
                    up[ch].buff[up[ch].now] = 0.0;
                    double in1 = 0.0;
                    for (int i = 0; i < tap_hm; i++) {
                        double a = up[ch].buff[up[ch].get_nth(i)];
                        double b = up[ch].buff[up[ch].get_nth(fir_size - 1 - i)];
                        in1 += up[ch].coef[i] * (a + b);
                    }   in1 += up[ch].coef[tap_hm] * up[ch].buff[up[ch].get_nth(tap_hm)];

                    upsample_buff[ch][(2 * i) + 0] = in0;
                    upsample_buff[ch][(2 * i) + 1] = in1;
                }
            }

            //int ch = 0;
            /*
            // Get max l/r amplitude values and fill sidechain signal
            for (int ch = 0; ch < numChannels; ch++) {
                for (int i = 0; i < 2*sampleFrames; i++) {
                    upsample_buff[0][i] = (abs(upsample_buff[0][i]) < abs(upsample_buff[ch][i]))
                        ? upsample_buff[0][i] : upsample_buff[ch][i];
                }
            }
            */
            for (int i = 0; i < 2 * sampleFrames; i++)
                upsample_buff[0][i] = smooth::ABS(upsample_buff[0][i], 0.0001);
            ballistics.processCrestFactor(upsample_buff[0].data(), 2 * sampleFrames);
            gainComputer.applyCompressionToBuffer(upsample_buff[0].data(), 2 * sampleFrames);
            ballistics.applyBallistics(upsample_buff[0].data(), 2 * sampleFrames);

            // Do lookahead if enabled
            if (true)//LookAheadEnable
            {
                // Delay input buffer
                //delay.process(buffer);
                int ch = 0;
                int latency = 128; // (int)lookaheadDelay;
                if (latency != delayline[ch].size()) {
                    int32 diff = latency - (int32)delayline[ch].size();
                    if (diff > 0)
                        for (int i = 0; i < diff; i++) delayline[ch].push_back(0.0);
                    else
                        for (int i = 0; i < -diff; i++) delayline[ch].pop_front();
                }

                for (int i = 0; i < 2 * sampleFrames; i++) {
                    double ahead = upsample_buff[0][i];
                    double delayed = delayline[ch].front();
                    delayline[ch].pop_front();
                    delayline[ch].push_back(upsample_buff[0][i]);

                    double acc = std::accumulate(delayline[ch].begin(), delayline[ch].end(), 0.0);
                    acc /= delayline[ch].size();
                    // upsample_buff[0][i] = acc;

                    if (ahead > acc) upsample_buff[0][i] = acc;
                    else upsample_buff[0][i] = acc;

                    //if (acc < delayed) upsample_buff[0][i] = acc;
                    //else upsample_buff[0][i] = delayed;
                }


                // Gain reduction must be smoothed with side chain window width
                // Maybe, crossfading with LAH / ORIG only when LAH > ORIG

                // Process side-chain (delay + gain reduction fade in)

            }

            for (int ch = 0; ch < 1; ch++) {
                for (int i = 0; i < sampleFrames; i++) {
                    // OS - downsample
                    //memmove(dn[ch].buff + 1, dn[ch].buff, sizeof(double) * (maxTap)-1);
                    dn[ch].acc();
                    dn[ch].buff[dn[ch].now] = upsample_buff[ch][(2 * i) + 0];

                    double out = 0.0;

                    for (int i = 0; i < tap_hm; i++) {
                        double a = dn[ch].buff[dn[ch].get_nth(i)];
                        double b = dn[ch].buff[dn[ch].get_nth(fir_size - 1 - i)];
                        out += dn[ch].coef[i] * (a + b);
                    }   out += dn[ch].coef[tap_hm] * dn[ch].buff[dn[ch].get_nth(tap_hm)];

                    //memmove(dn[ch].buff + 1, dn[ch].buff, sizeof(double) * (maxTap)-1);
                    dn[ch].acc();
                    dn[ch].buff[dn[ch].now] = upsample_buff[ch][(2 * i) + 1];

                    rawSidechainSignal[i] = (SampleType)out;
                }
            }

            /* - original

            // Calculate crest factor on max. amplitude values of input buffer
            ballistics.processCrestFactor(rawSidechainSignal, sampleFrames);

            // Compute attenuation - converts side-chain signal from linear to logarithmic domain
            gainComputer.applyCompressionToBuffer(rawSidechainSignal, sampleFrames);
            // Now rawSidechainSignal has ideal GR dB.
            // Smooth attenuation - still logarithmic
            ballistics.applyBallistics(rawSidechainSignal, sampleFrames);

            */

            // Get minimum = max. gain reduction from side chain buffer
            double min = 0;
            for (int i = 0; i < sampleFrames; i++)
                if (rawSidechainSignal[i] < min)
                    min = rawSidechainSignal[i];
            maxGainReduction = min;



            // Calculate auto makeup
            autoMakeup = calculateAutoMakeup(rawSidechainSignal, sampleFrames);

            // Add makeup gain and convert side-chain to linear domain
            for (int i = 0; i < sampleFrames; i++)
                sidechainSignal[i] = DecibelConverter::ToGain(sidechainSignal[i] + makeup + autoMakeup);



            for (int ch = 0; ch < numChannels; ch++)
            {
                int latency = 63 + 64; // (int)lookaheadDelay;
                if (latency != latency_q[ch].size()) {
                    int32 diff = latency - (int32)latency_q[ch].size();
                    if (diff > 0) for (int i = 0; i < diff; i++) latency_q[ch].push(0.0);
                    else for (int i = 0; i < -diff; i++) latency_q[ch].pop();
                }

                for (int j = 0; j < sampleFrames; j++) {
                    // Copy buffer to original signal
                    // originalSignal[ch][j] = (double)inputs[ch][j];

                    originalSignal[ch][j] = latency_q[ch].front();
                    latency_q[ch].pop();
                    latency_q[ch].push((double)inputs[ch][j]);


                    // Multiply attenuation with buffer - apply compression
                    double t = originalSignal[ch][j] * rawSidechainSignal[j];

                    // Mix dry & wet signal
                    outputs[ch][j] = (SampleType)(t * mix + originalSignal[ch][j] * (1.0 - mix));
                }
            }
        }
    };

private:
    template <typename SampleType>
    inline void applyInputGain(SampleType* buffer, int numSamples)
    {
        for (int i = 0; i < numSamples; i++) buffer[i] *= DecibelConverter::ToGain(input);
        /*
        if (prevInput == input)
            buffer.applyGain(0, numSamples, Decibels::decibelsToGain(prevInput));
        else
        {
            buffer.applyGainRamp(0, numSamples, Decibels::decibelsToGain(prevInput), Decibels::decibelsToGain(input));
            prevInput = input;
        }
        */
    };
    inline double calculateAutoMakeup(const double* src, int numSamples)
    {
        double sum = 0.0;
        for (int i = 0; i < numSamples; i++) {
            sum += src[i];
        }

        smoothedAutoMakeup.process(-sum / static_cast<double>(numSamples));
        return autoMakeupEnabled ? static_cast<double>(smoothedAutoMakeup.getState()) : 0.0f;
    };

    //Directly initialize process spec to avoid debugging problems
    double _sampleRate;
    double _maximumBlockSize;

    std::vector<double> originalSignal[2];
    std::vector<double> sidechainSignal;
    double* rawSidechainSignal{ nullptr };

    LevelDetector ballistics;
    GainComputer  gainComputer;
    SmoothingFilter smoothedAutoMakeup;

    double lookaheadDelay{ 0.005 };
    double input{ 0.0 };
    double prevInput{ 0.0 };
    double makeup{ 0.0 };
    double autoMakeup{ 0.0 };
    bool   bypassed{ false };
    bool   autoMakeupEnabled{ false };
    bool   LookAheadEnable{ true };
    double mix{ 1.0 };
    double maxGainReduction{ 0.0 };

    Flt up[2];
    Flt dn[2];
    std::vector<double> upsample_buff[2];
    std::vector<double> dnsample_buff[2];
    std::deque<double> delayline[2];
    const int fir_size = 125;
    const int tap_hm = (fir_size - 1) / 2;
    std::queue<double> latency_q[2];

    double c[2][4] = {
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
};


/* Basic envelope-follwer, to track peak & rms signal level with configurable decay time*/
class LevelEnvelopeFollower
{
public:
    LevelEnvelopeFollower() = default;

    // Prepares envelope follower with given sample rate and recalculates decayInSamples
    // as well as the peak/rms coefficient
    void prepare(const double& fs)
    {
        sampleRate = fs;

        peakDecayInSamples = static_cast<int>(peakDecayInSeconds * sampleRate);
        peakDecay = 1.0 - 1.0 / static_cast<double>(peakDecayInSamples);

        rmsDecayInSamples = static_cast<int>(rmsDecayInSeconds * sampleRate);
        rmsDecay = 1.0 - 1.0 / static_cast<double>(rmsDecayInSamples);

        double attackTimeInSeconds = 0.0; //Time it takes to reach 1-1/e = 0.63
        alphaAttack = 0.0;// exp(-1.0 / (sampleRate * attackTimeInSeconds)); //aA = e^(-1/TA*fs)

        double releaseTimeInSeconds = peakDecayInSeconds; //Time it takes to reach 1 - (1-1/e) = 0.37
        alphaRelease = exp(-1.0 / (sampleRate * releaseTimeInSeconds)); //aR = e^(-1/TR*fs)

        currMaxPeak[0] = DecibelConverter::ToDecibel(0.0);
        currMaxPeak[1] = DecibelConverter::ToDecibel(0.0);
        currMaxRMS[0] = DecibelConverter::ToDecibel(0.0);
        currMaxRMS[1] = DecibelConverter::ToDecibel(0.0);
        state[0] = DecibelConverter::ToDecibel(0.0);
        state[1] = DecibelConverter::ToDecibel(0.0);
    };

    // Set peak decay
    void setPeakDecay(double dc)
    {
        peakDecayInSeconds = dc;
        prepare(sampleRate);
    };

    // Set rms decay
    void setRmsDecay(double dc) {
        rmsDecayInSeconds = dc;
        prepare(sampleRate);
    };

    // Updates peak envelope follower from given audio buffer
    template <typename SampleType>
    void updatePeak(const SampleType* const* channelData, int numChannels, int numSamples)
    {
        if (numChannels > 0 && numSamples > 0) {
            for (int ch = 0; ch < numChannels; ch++) {
                for (int i = 0; i < numSamples; i++) {
                    double in = DecibelConverter::ToDecibel(std::abs(channelData[ch][i]));

                    if (in > state[ch])
                        state[ch] = alphaAttack * state[ch] + (1 - alphaAttack) * in;
                    else
                        state[ch] = alphaRelease * state[ch] + (1 - alphaRelease) * in;

                    currMaxPeak[ch] = (state[ch]); //y_L
                }
                currMaxPeak[ch] = std::min(6.0, (currMaxPeak[ch]));
            }
        }
    };

    // Updates rms envelope follower from given audio buffer
    void updateRMS(const double* const* channelData, int numChannels, int numSamples)
    {
        if (numChannels > 0 && numSamples > 0) {
            for (int ch = 0; ch < numChannels; ch++) {
                for (int i = 0; i < numSamples; i++) {
                    double sum = std::abs(channelData[ch][i]);
                    sum *= sum;

                    if (sum > currMaxRMS[ch])
                        currMaxRMS[ch] = sum * rmsDecay;
                    else if (currMaxRMS[ch] > 0.001f)
                        currMaxRMS[ch] *= peakDecay;
                    else currMaxRMS[ch] = 0.0f;
                }
            }
        }
    };

    // Gets current peak, call after updatePeak
    double getPeak(int channel) {
        return currMaxPeak[channel];
    };

    // Gets current rms, vall after updateRMS
    double getRMS(int channel) {
        return sqrt(currMaxRMS[channel]);
    };

private:
    double currMaxPeak[2] = { 0.0, 0.0 };
    double currMaxRMS[2] = { 0.0, 0.0 };
    double peakDecay{ 0.99992 };
    double rmsDecay{ 0.95 };
    double peakDecayInSeconds{ 0.5 };
    double rmsDecayInSeconds{ 0.0 };

    int peakDecayInSamples{ 0 };
    int rmsDecayInSamples{ 0 };

    double state[2] = { 0.0, 0.0 };
    double alphaAttack = 0.0;
    double alphaRelease = 0.0;
    double sampleRate{ 0.0 };
};

class simpleSVF {
public:
    enum f_type
    {
        HighPass,
        BandPass,
        LowPass,
        AllPass,
        Bell,
        LowShelf,
        HighShelf,
        
        kFltNum
    };

    typedef struct _dataset{
        bool   In = true;
        double dB = 0.0;
        double Hz = 1000.0;
        double Q = 1.0;
        f_type Type = Bell;
        double Fs = 48000.0;

        double w = Hz * M_PI / Fs;;
        double g = tan(w);
        double k = 2.0 / Q;
        double gt0 = 1 / (1 + g * (g + k));
        double gk0 = (g + k) * gt0;

        double m0 = 1.0, m1 = 1.0, m2 = 1.0;
        double v0 = 0.0, v1 = 0.0, v2 = 0.0;
        double t0 = 0.0, t1 = 0.0, t2 = 0.0;
        double ic1eq = 0.0;
        double ic2eq = 0.0;
    } dataset;

    simpleSVF() {initSVF();}
    
    /*
    simpleSVF& operator=(const simpleSVF& obj) {
        memcpy(&flt, &obj.flt, sizeof(flt));
        return *this;
    }
     */
    
    void setSVF(bool In, double Hz, double Q, double dB, f_type type, double Fs)
    {
        flt.In    = In;
        flt.Fs    = Fs;
        flt.Hz    = Hz;
        flt.Q     = Q;
        flt.dB    = dB;
        flt.Type  = type;
        
        makeSVF();
    }
    
    void setIn(bool In) { flt.In = In; }

    void initSVF() {
        flt.ic1eq = 0.0; flt.ic2eq = 0.0;
    }

    void makeSVF()
    {
        if (flt.Hz > flt.Fs / 2.0) flt.Hz = flt.Fs / 2.0;
        flt.w = flt.Hz * M_PI / flt.Fs;
        double g0 = tan(flt.w);
        double k0 = 1.0 / flt.Q;
        double A = pow(10.0, flt.dB / 40.0);

        double kdA = k0 / A;
        double kmA = k0 * A;
        double gdSA = g0 / sqrt(A);
        double gmSA = g0 * sqrt(A);
        double AmA = A * A;

        switch (flt.Type)
        {
            case HighPass:  flt.m0 = 1;   flt.m1 = 0;   flt.m2 = 0;   flt.g = g0;   flt.k = k0;   break;
            case BandPass:  flt.m0 = 0;   flt.m1 = 1;   flt.m2 = 0;   flt.g = g0;   flt.k = k0;   break;
            case LowPass:   flt.m0 = 0;   flt.m1 = 0;   flt.m2 = 1;   flt.g = g0;   flt.k = k0;   break;
            case AllPass:   flt.m0 = 1;   flt.m1 = -k0; flt.m2 = 1;   flt.g = g0;   flt.k = k0;   break;
            case Bell:      flt.m0 = 1;   flt.m1 = kmA; flt.m2 = 1;   flt.g = g0;   flt.k = kdA;  break;
            case LowShelf:  flt.m0 = 1;   flt.m1 = kmA; flt.m2 = AmA; flt.g = gdSA; flt.k = k0;   break;
            case HighShelf: flt.m0 = AmA; flt.m1 = kmA; flt.m2 = 1;   flt.g = gmSA; flt.k = k0;   break;
            default: break;
        }

        flt.gt0 = 1.0 / (1.0 + flt.g * (flt.g + flt.k));
        flt.gk0 = (flt.g + flt.k) * flt.gt0;
        return;
    }
    
    inline double computeSVF (double vin)
    {
        // tick serial(possibly quicker on cpus with low latencies)
        flt.t0 = vin - flt.ic2eq;
        flt.v0 = flt.gt0 * flt.t0 - flt.gk0 * flt.ic1eq; // high
        flt.t1 = flt.g * flt.v0;
        flt.v1 = flt.ic1eq + flt.t1; // band
        flt.t2 = flt.g * flt.v1;
        flt.v2 = flt.ic2eq + flt.t2; // low
        flt.ic1eq += 2.0 * flt.t1;
        flt.ic2eq += 2.0 * flt.t2;

        if (!flt.In) return vin; // NOT early return, to keep filter running while off for continuity
        return flt.m0 * flt.v0 + flt.m1 * flt.v1 + flt.m2 * flt.v2;
    }

    inline double mag_response (double freq) {
        if (!flt.In) return 1.0;

        double ONE_OVER_SAMPLE_RATE = 1.0 / flt.Fs;

        // exp(complex(0.0, -2.0 * pi) * frequency / sampleRate)
        double _zr = (0.0) * freq * ONE_OVER_SAMPLE_RATE;
        double _zi = (-2.0 * M_PI) * freq * ONE_OVER_SAMPLE_RATE;

        // z = zr + zi;
        double zr = exp(_zr) * cos(_zi);
        double zi = exp(_zr) * sin(_zi);

        double nr = 0, ni = 0;
        double dr = 0, di = 0;

        // z * z
        double zsq_r = zr * zr - zi * zi;
        double zsq_i = zi * zr + zr * zi;
        double gsq = flt.g * flt.g;

        // Numerator complex
        double c_nzsq = (flt.m0 + flt.m1 * flt.g + flt.m2 * gsq);
        double c_nz = (flt.m0 * -2.0 + flt.m2 * 2.0 * gsq);
        double c_n = (flt.m0 + flt.m1 * -flt.g + flt.m2 * gsq);
        nr = zsq_r * c_nzsq + zr * c_nz + c_n;
        ni = zsq_i * c_nzsq + zi * c_nz;

        // Denominator complex
        double c_dzsq = (1.0 + flt.k * flt.g + gsq);
        double c_dz = (-2.0 + 2.0 * gsq);
        double c_d = (1.0 + flt.k * -flt.g + gsq);
        dr = zsq_r * c_dzsq + zr * c_dz + c_d;
        di = zsq_i * c_dzsq + zi * c_dz;


        // Numerator / Denominator
        double norm = dr * dr + di * di;
        double ddr = (nr * dr + ni * di) / norm;
        double ddi = (ni * dr - nr * di) / norm;

        return sqrt(ddr * ddr + ddi * ddi);
    }

private:
    dataset flt;
};






class SVF_FirstOrder {
public:
    enum filter_Type
    {
        kLowPass,
        kHighPass
    };
    
    typedef struct ds{
        int    In = 0;
        double Hz = 1000.0;
        filter_Type Type = kLowPass;
        double Fs = 48000.0;

        double w = Hz * M_PI / Fs;;
        double g = tan(w);

        double m0 = 1.0, m2 = 1.0;
        double v0 = 0.0, v2 = 0.0;
        double t0 = 0.0, t2 = 0.0;
        double iceq = 0.0;
    } dataset;

    SVF_FirstOrder()
    {
        initSVF();
    };

    void initSVF() {
        flt.iceq = 0.0;
    };

    void setSVF(double fParamIn, double fParamHz, filter_Type fParamtype, double fParamFs)
    {
        flt.In    = fParamIn ? 1 : 0;
        flt.Fs    = fParamFs;
        flt.Hz    = fParamHz;
        flt.Type  = fParamtype;

        makeSVF();
    }

    void makeSVF()
    {
        if (flt.Hz > flt.Fs / 2.0) flt.Hz = flt.Fs / 2.0;
        flt.w = flt.Hz * M_PI / flt.Fs;
        flt.g = tan(flt.w);

        switch (flt.Type)
        {
            case kLowPass:      flt.m0 = 0; flt.m2 = 1;   break;
            case kHighPass:     flt.m0 = 1; flt.m2 = 0;   break;
            default: break;
        }

        return;
    };

    double computeSVF (double vin)
    {
        // disable v1 stage
        flt.t0 = vin - flt.iceq;
        flt.v0 = flt.t0 / (1.0 + flt.g);// gt0 * t0;
        flt.t2 = flt.g * flt.v0;
        flt.v2 = flt.iceq + flt.t2;
        flt.iceq += 2.0 * flt.t2;

        if (flt.In != 1) return vin;
        return flt.m0 * flt.v0 + flt.m2 * flt.v2;
    };

private:
    dataset flt;
};

//------------------------------------------------------------------------






class delayLine
{
public:
    delayLine (int maxDelayInSamples = minTotalSize)
    {
        setMaximumDelayInSamples (maxDelayInSamples);
    }
    
    void setMaximumDelayInSamples (int maxDelayInSamples)
    {
        if (maxDelayInSamples < minTotalSize) maxDelayInSamples = minTotalSize;
        
        totalSize = maxDelayInSamples;
        
        for (auto& iter : bufferData)
            iter.resize(totalSize);
        
        reset();
    }
    int getMaximumDelayInSamples() const noexcept       { return totalSize; }
    
    
    void setDelay (int newDelayInSamples)
    {
        if (newDelayInSamples < 0) newDelayInSamples = 0;
        if (newDelayInSamples > totalSize) newDelayInSamples = totalSize;

        delay = newDelayInSamples;
    }
    int getDelay() const
    {
        return delay;
    }
    
    void prepare (int numChannels)
    {
        if (numChannels <= 0) numChannels = 2;

        bufferData.resize(numChannels);
        for (auto& iter : bufferData)
            iter.resize(totalSize);

        writePos.resize (numChannels);
        readPos.resize  (numChannels);

        reset();
    }

    void reset()
    {
        for (auto vec : { &writePos, &readPos })
            std::fill (vec->begin(), vec->end(), 0);

        for (auto& iter : bufferData)
            std::fill (iter.begin(), iter.end(), 0.0);
    }

    void pushSample (int channel, double sample)
    {
        bufferData[channel][writePos[(size_t) channel]] = sample;
        
        writePos[(size_t) channel] = (writePos[(size_t) channel] + totalSize - 1) % totalSize;
    }

    double popSample (int channel)
    {
        auto index = (readPos[(size_t) channel] + delay) % totalSize;

        double result = bufferData[channel][index];

        readPos[(size_t) channel] = (readPos[(size_t) channel] + totalSize - 1) % totalSize;

        return result;
    }
    
    double getSample (int channel, int pos)
    {
        auto index = (readPos[(size_t) channel] + pos) % totalSize;

        return bufferData[channel][index];
    }


private:
    //==============================================================================
    static constexpr int minTotalSize = 1;
    
    //==============================================================================
    std::vector<std::vector<double>> bufferData;
    std::vector<int> writePos, readPos;
    int delay = 0;
    int totalSize = minTotalSize;
};

} // namespace yg331
