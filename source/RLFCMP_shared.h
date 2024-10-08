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
#include <deque>
#include <complex>
#include <array>

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

enum detectorType
{
    Bold,
    Smooth,
    Clean,
    detectorNum = 2
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

    inline double mag_response (double freq)
    {
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

// for sake of inline
// private:
    dataset flt;
};




class PassShelfFilter {
public:
    enum type : int
    {
        Pass,
        Shelf,
        FltNum = 1
    };
    enum range : int
    {
        Low,
        High
    };
    
    typedef struct ds
    {
        bool   In = true;
        double Freq = 100.0;
        double dB = 0.0;
        type   Type = Pass;
        range  Range = Low;
        double Fs = 48000.0;

        double w = Freq * M_PI / Fs;;
        double g = tan(w);

        double m0 = 1.0, m2 = 1.0;
        double v0 = 0.0, v2 = 0.0;
        double t0 = 0.0, t2 = 0.0;
        double iceq = 0.0;
    } dataset;

    PassShelfFilter()
    {
        initSVF();
    }

    void initSVF() {
        flt.iceq = 0.0;
    }
    
    void   setIn(bool v) { flt.In = v; }
    bool   getIn() const { return flt.In; }
    
    void   setFreq(double v) { flt.Freq = v; }
    double getFreq() const { return flt.Freq; }
    
    void   setGain(double v) { flt.dB = v; }
    double getGain() const { return flt.dB; }
    
    void   setType(int v) { flt.Type = static_cast<type>(v); }
    int    getType() const { return static_cast<int>(flt.Type); }
    
    void   setRange(int v) { flt.Range = static_cast<range>(v); }
    int    getRange() const { return static_cast<int>(flt.Range); }

    void   setFs(double v) { flt.Fs = v; }
    double getFs() const { return flt.Fs; }
    
    void setSVF(double In, double Freq, double dB, int _type, double Fs)
    {
        flt.In    = In ? 1 : 0;
        flt.Freq  = Freq;
        flt.dB    = dB;
        flt.Type  = static_cast<type>(_type);
        flt.Fs    = Fs;

        makeSVF();
    }

    void makeSVF()
    {
        if (flt.Freq > flt.Fs / 2.0) flt.Freq = flt.Fs / 2.0;
        flt.w = flt.Freq * M_PI / flt.Fs;
        flt.g = tan(flt.w);
        double A = pow(10.0, flt.dB / 40.0);
        double AmA = A * A;
        
        if (flt.Range == Low)
        {
            switch (flt.Type)
            {
                case Pass:  flt.m0 = 1;   flt.m2 = 0;   break;
                case Shelf: flt.m0 = 1;   flt.m2 = AmA; break;
                default: break;
            }
        }
        else // if (flt.Range == High)
        {
            switch (flt.Type)
            {
                case Pass:  flt.m0 = 0;   flt.m2 = 1;   break; // HighCut, LowPass
                case Shelf: flt.m0 = AmA; flt.m2 = 1;   break;
                default: break;
            }
        }

        return;
    };

    inline double computeSVF (double vin)
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
    
    inline double mag_response (double freq)
    {
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

        // Numerator complex
        nr = zr * (-flt.m0 /* + m1 * (g - 1) */ + flt.m2 * flt.g) + (flt.m0 /* + m1 * (g + 1) */ + flt.m2 * flt.g);
        ni = zi * (-flt.m0 /* + m1 * (g - 1) */ + flt.m2 * flt.g);

        // Denominator complex
        dr = zr * (flt.g - 1) + (flt.g + 1);
        di = zi * (flt.g - 1);

        // Numerator / Denominator
        double norm = dr * dr + di * di;
        double ddr = (nr * dr + ni * di) / norm;
        double ddi = (ni * dr - nr * di) / norm;

        return sqrt(ddr * ddr + ddi * ddi);
    }

// private:
    dataset flt;
};

//------------------------------------------------------------------------

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

//------------------------------------------------------------------------
//  Min, Max, Default of Parameters
//------------------------------------------------------------------------
static constexpr ParamValue dftBypass          = 0.0;
static constexpr ParamValue dftSoftBypass      = 0.0;
static constexpr ParamValue dftLookaheadEnable = 1.0;
static constexpr ParamValue dftScLfIn          = 1.0;
static constexpr ParamValue dftScHfIn          = 1.0;
static constexpr ParamValue dftType            = detectorType::Smooth;

static constexpr ParamValue minScLfFreq  = 20.0;
static constexpr ParamValue maxScLfFreq  = 18000.0;
static constexpr ParamValue dftScLfFreq  = 60.0;

static constexpr ParamValue minScLfGain  = -20.0;
static constexpr ParamValue maxScLfGain  = 20.0;
static constexpr ParamValue dftScLfGain  = 0.0;

static constexpr ParamValue minScHfFreq  = 20.0;
static constexpr ParamValue maxScHfFreq  = 18000.0;
static constexpr ParamValue dftScHfFreq  = 1880.0;

static constexpr ParamValue minScHfGain  = -20.0;
static constexpr ParamValue maxScHfGain  = 20.0;
static constexpr ParamValue dftScHfGain  = 4.0;

static constexpr ParamValue minAttack    = 0.5;
static constexpr ParamValue maxAttack    = 70.0;
static constexpr ParamValue dftAttack    = 8.0;

static constexpr ParamValue minRelease   = 20.0;
static constexpr ParamValue maxRelease   = 2000.0;
static constexpr ParamValue dftRelease   = 100.0;

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

static constexpr ParamValue minInput     = -12.0;
static constexpr ParamValue maxInput     = 12.0;
static constexpr ParamValue dftInput     = 0.0;

static constexpr ParamValue minOutput    = -12.0;
static constexpr ParamValue maxOutput    = 12.0;
static constexpr ParamValue dftOutput    = 0.0;

static const ParameterConverter paramScLfType (0,  0,  ParameterConverter::paramType::list, PassShelfFilter::type::FltNum);
static const ParameterConverter paramScLfFreq (minScLfFreq,  maxScLfFreq,  ParameterConverter::paramType::log);
static const ParameterConverter paramScLfGain (minScLfGain,  maxScLfGain,  ParameterConverter::paramType::range);
static const ParameterConverter paramScHfType (0,  0,  ParameterConverter::paramType::list, PassShelfFilter::type::FltNum);
static const ParameterConverter paramScHfFreq (minScHfFreq,  maxScHfFreq,  ParameterConverter::paramType::log);
static const ParameterConverter paramScHfGain (minScHfGain,  maxScHfGain,  ParameterConverter::paramType::range);
static const ParameterConverter paramType     (0,  0,  ParameterConverter::paramType::list, detectorType::detectorNum);
static const ParameterConverter paramAttack   (minAttack,    maxAttack,    ParameterConverter::paramType::log);
static const ParameterConverter paramRelease  (minRelease,   maxRelease,   ParameterConverter::paramType::log);
static const ParameterConverter paramThreshold(minThreshold, maxThreshold, ParameterConverter::paramType::range);
static const ParameterConverter paramRatio    (minRatio,     maxRatio,     ParameterConverter::paramType::range);
static const ParameterConverter paramKnee     (minKnee,      maxKnee,      ParameterConverter::paramType::range);
static const ParameterConverter paramMakeup   (minMakeup,    maxMakeup,    ParameterConverter::paramType::range);
static const ParameterConverter paramMix      (minMix,       maxMix,       ParameterConverter::paramType::range);
static const ParameterConverter paramInput    (minInput,     maxInput,     ParameterConverter::paramType::range);
static const ParameterConverter paramOutput   (minOutput,    maxOutput,    ParameterConverter::paramType::range);

static constexpr ParamValue PF_FREQ  = 1500.0;
static constexpr ParamValue PF_Q     = M_SQRT1_2; // == 1.0 ~= 1.007
static constexpr ParamValue PF_dB    = 4.0;
static constexpr ParamValue RLB_FREQ = 38.134;
static constexpr ParamValue RLB_Q    = 0.70758;

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
        
        bufferData.resize(totalSize);
        
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

    void reset()
    {
        writePos = 0;
        readPos = 0;
        std::fill (bufferData.begin(), bufferData.end(), 0.0);
    }

    void pushSample (int channel, double sample)
    {
        bufferData[writePos] = sample;
        
        writePos = (writePos + totalSize - 1) % totalSize;
    }

    double popSample (int channel)
    {
        auto index = (readPos + delay) % totalSize;

        double result = bufferData[index];

        readPos = (readPos + totalSize - 1) % totalSize;

        return result;
    }
    
    double getSample (int channel, int pos)
    {
        auto index = (readPos + pos) % totalSize;

        return bufferData[index];
    }

private:
    //==============================================================================
    static constexpr int minTotalSize = 1;
    
    //==============================================================================
    std::vector<double> bufferData;
    int writePos = 0, readPos = 0;
    int delay = 0;
    int totalSize = minTotalSize;
};












static constexpr int _fftOrder = 12;
static constexpr int _fftSize = 1 << _fftOrder;      // 4096 samples
static constexpr int _numBins = _fftSize / 2 + 1;    // 2049 bins

/*
   PFFFT : a Pretty Fast FFT.

   This is basically an adaptation of the single precision fftpack
   (v4) as found on netlib taking advantage of SIMD instruction found
   on cpus such as intel x86 (SSE1), powerpc (Altivec), and arm (NEON).

   For architectures where no SIMD instruction is available, the code
   falls back to a scalar version.

   Restrictions:

   - 1D transforms only, with 32-bit single precision.

   - supports only transforms for inputs of length N of the form
   N=(2^a)*(3^b)*(5^c), a >= 5, b >=0, c >= 0 (32, 48, 64, 96, 128,
   144, 160, etc are all acceptable lengths). Performance is best for
   128<=N<=8192.

   - all (float*) pointers in the functions below are expected to
   have an "simd-compatible" alignment, that is 16 bytes on x86 and
   powerpc CPUs.

   You can allocate such buffers with the functions
   pffft_aligned_malloc / pffft_aligned_free (or with stuff like
   posix_memalign..)

*/

#ifndef PFFFT_H
#define PFFFT_H

#include <stddef.h> // for size_t

#ifdef __cplusplus
extern "C" {
#endif

    /* opaque struct holding internal stuff (precomputed twiddle factors)
       this struct can be shared by many threads as it contains only
       read-only data.
    */
    typedef struct PFFFT_Setup PFFFT_Setup;

    /* direction of the transform */
    typedef enum { PFFFT_FORWARD, PFFFT_BACKWARD } pffft_direction_t;

    /* type of transform */
    typedef enum { PFFFT_REAL, PFFFT_COMPLEX } pffft_transform_t;

    /*
      prepare for performing transforms of size N -- the returned
      PFFFT_Setup structure is read-only so it can safely be shared by
      multiple concurrent threads.
    */
    PFFFT_Setup* pffft_new_setup(int N, pffft_transform_t transform);
    void pffft_destroy_setup(PFFFT_Setup*);

    /*
       Perform a Fourier transform , The z-domain data is stored in the
       most efficient order for transforming it back, or using it for
       convolution. If you need to have its content sorted in the
       "usual" way, that is as an array of interleaved complex numbers,
       either use pffft_transform_ordered , or call pffft_zreorder after
       the forward fft, and before the backward fft.

       Transforms are not scaled: PFFFT_BACKWARD(PFFFT_FORWARD(x)) = N*x.
       Typically you will want to scale the backward transform by 1/N.

       The 'work' pointer should point to an area of N (2*N for complex
       fft) floats, properly aligned. If 'work' is NULL, then stack will
       be used instead (this is probably the best strategy for small
       FFTs, say for N < 16384).

       input and output may alias.
    */
    void pffft_transform(PFFFT_Setup* setup, const float* input, float* output, float* work, pffft_direction_t direction);

    /*
       Similar to pffft_transform, but makes sure that the output is
       ordered as expected (interleaved complex numbers).  This is
       similar to calling pffft_transform and then pffft_zreorder.

       input and output may alias.
    */
    void pffft_transform_ordered(PFFFT_Setup* setup, const float* input, float* output, float* work, pffft_direction_t direction);

    /*
       call pffft_zreorder(.., PFFFT_FORWARD) after pffft_transform(...,
       PFFFT_FORWARD) if you want to have the frequency components in
       the correct "canonical" order, as interleaved complex numbers.

       (for real transforms, both 0-frequency and half frequency
       components, which are real, are assembled in the first entry as
       F(0)+i*F(n/2+1). Note that the original fftpack did place
       F(n/2+1) at the end of the arrays).

       input and output should not alias.
    */
    void pffft_zreorder(PFFFT_Setup* setup, const float* input, float* output, pffft_direction_t direction);

    /*
       Perform a multiplication of the frequency components of dft_a and
       dft_b and accumulate them into dft_ab. The arrays should have
       been obtained with pffft_transform(.., PFFFT_FORWARD) and should
       *not* have been reordered with pffft_zreorder (otherwise just
       perform the operation yourself as the dft coefs are stored as
       interleaved complex numbers).

       the operation performed is: dft_ab += (dft_a * fdt_b)*scaling

       The dft_a, dft_b and dft_ab pointers may alias.
    */
    void pffft_zconvolve_accumulate(PFFFT_Setup* setup, const float* dft_a, const float* dft_b, float* dft_ab, float scaling);

    /*
      the float buffers must have the correct alignment (16-byte boundary
      on intel and powerpc). This function may be used to obtain such
      correctly aligned buffers.
    */
    void* pffft_aligned_malloc(size_t nb_bytes);
    void  pffft_aligned_free(void*);

    /* return 4 or 1 wether support SSE/Altivec instructions was enable when building pffft.c */
    int pffft_simd_size(void);

#ifndef PFFFT_SIMD_DISABLE
    void validate_pffft_simd(); // a small function inside pffft.c that will detect compiler bugs with respect to simd instruction
#endif

#ifdef __cplusplus
}
#endif

#endif // PFFFT_H

/**
 * C++ Wrapper for pffft, a reasonably fast FFT library.
 *  The class here reflects closely the Juce FFT class and is a drop
 *  in replacement.
 *  See: https://bitbucket.org/jpommier/pffft/src/master/
 */
class PFFFT
{
public:
    PFFFT(int order)
    {
        size_ = 1 << order;
        scale_ = 1.f / size_;
        pSetup_ = pffft_new_setup(size_, PFFFT_REAL);
    }

    ~PFFFT()
    {
        pffft_destroy_setup(pSetup_);
    }

    void performRealOnlyForwardTransform(float* pBuffer, bool onlyCalculateNonNegativeFrequencies = false)
    {
        pffft_transform_ordered(pSetup_, pBuffer, pBuffer, NULL, PFFFT_FORWARD);
    }

    void performRealOnlyInverseTransform(float* pBuffer)
    {
        pffft_transform_ordered(pSetup_, pBuffer, pBuffer, NULL, PFFFT_BACKWARD);

        for (int i = 0; i < size_; ++i)
        {
            pBuffer[i] *= scale_;
        }
    }

    void performFrequencyOnlyForwardTransform(float* inputOutputData, bool ignoreNegativeFreqs = false) const noexcept
    {
        if (size_ == 1) return;

        pffft_transform_ordered(pSetup_, inputOutputData, inputOutputData, NULL, PFFFT_FORWARD);

        auto* out = reinterpret_cast<std::complex<float> *>(inputOutputData);

        const auto limit = ignoreNegativeFreqs ? (size_ / 2) + 1 : size_;

        for (int i = 0; i < limit; ++i)
        {
            inputOutputData[i] = std::abs(out[i]);
        }

        std::fill(inputOutputData + limit, inputOutputData + size_ * 2, 0.f);
    }

    [[nodiscard]] int getSize() const noexcept { return size_; }

private:
    int size_;
    float scale_;

    PFFFT_Setup* pSetup_;
};

/**
STFT analysis and resynthesis of audio data.

Each channel should have its own FFTProcessor.
*/
class FFTProcessor
{
public:
    FFTProcessor();

    int getLatencyInSamples() const { return fftSize; }

    void reset();
    float processSample(float sample, bool bypassed);
    void processBlock(float* data, int numSamples, bool bypassed);
    void processBlock(double* data, int numSamples, bool bypassed);
    int getData(float* out) {
        if (!data_avail)
            return 0;

        auto* cdata = reinterpret_cast<std::complex<float>*>(&fftData);
        for (int i = 0; i < numBins; ++i) {
            float magnitude = std::abs(cdata[i]);
            out[i] = magnitude;
        }
        data_avail = 0;
        return 1;
    };
    static void hannWindow(float* window, int length)
    {
        float delta = (M_PI * 2) / float(length);
        float phase = 0.0f;
        for (int i = 0; i < length; ++i) {
            window[i] = 0.5f * (1.0f - std::cos(phase));
            phase += delta;
            // phase = i * delta
        }
        float pwr = 0.0;
        for (int i = 0; i < length; ++i) {
            pwr += window[i];
        }
        pwr = 1.0 / pwr;
        for (int i = 0; i < length; ++i) {
            window[i] *= pwr;
        }
    }
    static void bkhsWindow(float* window, int length)
    {
        double dwindowpos = (M_PI * 2) / length;
        double pwr = 0.0;
        for (int i = 0; i < 0.5 * length + 1; i++) {
            double windowpos = i * dwindowpos;
            window[i] = (0.35875 - 0.48829 * cos(windowpos) + 0.14128 * cos(2.0 * windowpos) - 0.01168 * cos(3.0 * windowpos));
            pwr += window[i];
        }
        pwr = 0.5 / (pwr * 2 - window[(int)(0.5 * length)]);
        for (int i = 0; i < 0.5 * length + 1; i++) {
            window[i] *= pwr;
        }
        for (int i = 0; i < 0.5 * length; i++) {
            window[length - i - 1] = window[i];
        }
    }
    static inline float Ino(float x)
    {
        double d = 0, ds = 1, s = 1;
        do
        {
            d += 2;
            ds *= x * x / (d * d);
            s += ds;
        } while (ds > s * 1e-6);
        return s;
    };

    static void ksblWindow(float* window, int length)
    {
        int Np = (length - 1) / 2;
        float Alpha;
        float Inoalpha;

        Alpha = 3.0 * M_PI;

        Inoalpha = Ino(Alpha);

        for (int j = 0; j <= Np; j++)
        {
            window[Np + j] = Ino(Alpha * std::sqrt(1.0 - ((float)(j * j) / (float)(Np * Np)))) / Inoalpha;
        }
        for (int j = 0; j < Np; j++)
        {
            window[j] = window[length - 1 - j];
        }

        float pwr = 0.0;
        for (int i = 0; i < length; i++) {
            pwr += window[i];
        }
        pwr = 1.0 / pwr;
        for (int i = 0; i < length; i++) {
            window[i] *= pwr * 1.86; // Normalization & Amplitube correction
        }
    };

private:
    void processFrame(bool bypassed);
    void processSpectrum(float* data, int numBins);

    // The FFT has 2^order points and fftSize/2 + 1 bins.
    static constexpr int fftOrder = _fftOrder;
    static constexpr int fftSize = 1 << fftOrder;      // 4096 samples
    static constexpr int numBins = fftSize / 2 + 1;    // 2049 bins
    static constexpr int overlap = 4;                  // 75% overlap
    static constexpr int hopSize = fftSize / overlap;  // //256 samples

    // Gain correction for using Hann window with 75% overlap.
    static constexpr float windowCorrection = 2.0f / 3.0f;

    PFFFT fft;
    std::array<float, fftSize> window = { 0.0, };

    // Counts up until the next hop.
    int count = 0;

    int data_avail = 0;

    // Write position in input FIFO and read position in output FIFO.
    int pos = 0;

    // Circular buffers for incoming and outgoing audio data.
    /* SSE and co like 16-bytes aligned pointers */
    alignas(16) std::array<float, fftSize>  inputFifo;
    alignas(16) std::array<float, fftSize> outputFifo;

    // The FFT working space. Contains interleaved complex numbers.
    alignas(16) std::array<float, fftSize * 2>  fftData;
};

} // namespace yg331


// License for PFFFT
/* Copyright (c) 2013  Julien Pommier ( pommier@modartt.com )

   Based on original fortran 77 code from FFTPACKv4 from NETLIB,
   authored by Dr Paul Swarztrauber of NCAR, in 1985.

   As confirmed by the NCAR fftpack software curators, the following
   FFTPACKv5 license applies to FFTPACKv4 sources. My changes are
   released under the same terms.

   FFTPACK license:

   http://www.cisl.ucar.edu/css/software/fftpack5/ftpk.html

   Copyright (c) 2004 the University Corporation for Atmospheric
   Research ("UCAR"). All rights reserved. Developed by NCAR's
   Computational and Information Systems Laboratory, UCAR,
   www.cisl.ucar.edu.

   Redistribution and use of the Software in source and binary forms,
   with or without modification, is permitted provided that the
   following conditions are met:

   - Neither the names of NCAR's Computational and Information Systems
   Laboratory, the University Corporation for Atmospheric Research,
   nor the names of its sponsors or contributors may be used to
   endorse or promote products derived from this Software without
   specific prior written permission.

   - Redistributions of source code must retain the above copyright
   notices, this list of conditions, and the disclaimer below.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions, and the disclaimer below in the
   documentation and/or other materials provided with the
   distribution.

   THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT
   HOLDERS BE LIABLE FOR ANY CLAIM, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE
   SOFTWARE.
*/

// License for PFFFT c++ wrapper
// https://www.kvraudio.com/forum/viewtopic.php?p=8726913#p8726913
/*
    ==============================================================================

    PFFT.h
    Created: 29 Aug 2022 1:08:12pm
    Author:  Justin Johnson

    ==============================================================================

    MIT License

    Copyright (c) 2021 Justin Johnson

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

// License for FFTProcessor
/*  Copyright (c) 2023 Matthijs Hollemans
 
    MIT License

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/
