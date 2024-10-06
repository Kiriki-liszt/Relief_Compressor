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
        double Hz = 100.0;
        double dB = 0.0;
        type   Type = Pass;
        range  Range = Low;
        double Fs = 48000.0;

        double w = Hz * M_PI / Fs;;
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
    
    void setRange (int _range)
    {
        flt.Range = static_cast<range>(_range);
    }
    
    void setSVF(double In, double Hz, double dB, int _type, double Fs)
    {
        flt.In    = In ? 1 : 0;
        flt.Hz    = Hz;
        flt.dB    = dB;
        flt.Type  = static_cast<type>(_type);
        flt.Fs    = Fs;

        makeSVF();
    }

    void makeSVF()
    {
        if (flt.Hz > flt.Fs / 2.0) flt.Hz = flt.Fs / 2.0;
        flt.w = flt.Hz * M_PI / flt.Fs;
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

static constexpr ParamValue minOutput    = -12.0;
static constexpr ParamValue maxOutput    = 12.0;
static constexpr ParamValue dftOutput    = 0.0;

static const ParameterConverter paramScLfType (0,  0,  ParameterConverter::paramType::list, PassShelfFilter::type::FltNum);
static const ParameterConverter paramScLfFreq (minScLfFreq,  maxScLfFreq,  ParameterConverter::paramType::log);
static const ParameterConverter paramScLfGain (minScLfGain,  maxScLfGain,  ParameterConverter::paramType::range);
static const ParameterConverter paramScHfType (0,  0,  ParameterConverter::paramType::list, PassShelfFilter::type::FltNum);
static const ParameterConverter paramScHfFreq (minScHfFreq,  maxScHfFreq,  ParameterConverter::paramType::log);
static const ParameterConverter paramScHfGain (minScHfGain,  maxScHfGain,  ParameterConverter::paramType::range);
static const ParameterConverter paramAttack   (minAttack,    maxAttack,    ParameterConverter::paramType::log);
static const ParameterConverter paramRelease  (minRelease,   maxRelease,   ParameterConverter::paramType::log);
static const ParameterConverter paramThreshold(minThreshold, maxThreshold, ParameterConverter::paramType::range);
static const ParameterConverter paramRatio    (minRatio,     maxRatio,     ParameterConverter::paramType::range);
static const ParameterConverter paramKnee     (minKnee,      maxKnee,      ParameterConverter::paramType::range);
static const ParameterConverter paramMakeup   (minMakeup,    maxMakeup,    ParameterConverter::paramType::range);
static const ParameterConverter paramMix      (minMix,       maxMix,       ParameterConverter::paramType::range);
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
