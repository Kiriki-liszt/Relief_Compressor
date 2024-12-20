//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#include "RLFCMP_controller.h"
#include "RLFCMP_cids.h"
#include "vstgui/plugin-bindings/vst3editor.h"
#include "base/source/fstreamer.h"
#include "pluginterfaces/base/ustring.h"
#include "vstgui/vstgui_uidescription.h"
#include "vstgui/uidescription/detail/uiviewcreatorattributes.h"


using namespace Steinberg;

namespace VSTGUI
{
//------------------------------------------------------------------------
// EQCurveView
//------------------------------------------------------------------------
EQCurveView::EQCurveView(
    const CRect& size,
    IControlListener* listener,
    int32_t tag,
    CBitmap* background
)
    : CControl(size, listener, tag, background)
{
    BackColor = kWhiteCColor;
    LineColor = kBlackCColor;
    BorderColor = kBlackCColor;
    FFTLineColor = kBlackCColor;
    FFTFillColor = kBlackCColor;
    idleRate = 60;
    
    LF_SVF.setQ(M_SQRT1_2);
    HF_SVF.setQ(M_SQRT1_2);
    LF_SVF.setFs(48000.0);
    HF_SVF.setFs(48000.0);
    
    setWantsIdle(true);
}
EQCurveView::EQCurveView(const EQCurveView& v)
    : CControl(v)
    , BackColor(v.BackColor)
    , LineColor(v.LineColor)
    , BorderColor(v.BorderColor)
    , FFTLineColor(v.FFTLineColor)
    , FFTFillColor(v.FFTFillColor)
{
    setWantsIdle(true);
}

#define safe_bin(bin, x)   std::max(std::min((int)bin  + (int)x, (int)numBins   - 1) , 0)
#define safe_band(band, x) std::max(std::min((int)band + (int)x, (int)MAX_BANDS - 1) , 0)

void EQCurveView::setFFTArray(float* fftArray, int sampleBlockSize, double sampleRate)
{
    // Unit frequency per bin, with sample rate
    double freqBin_width = sampleRate / fftSize;
    
    double SR = sampleRate / (double)sampleBlockSize;
    double coeff = exp(-1.0 / (100.0 * 0.001 * SR));

    for (int i = 0; i < numBins; ++i)
    {
        fft_freq[i] = (i + 0.5) * freqBin_width;
        fft_RMS[i]  = coeff * fft_RMS[i] + (1.0 - coeff) * fftArray[i] * fftArray[i];
        fft_linear[i] = std::sqrt(fft_RMS[i]);
    }
}

void EQCurveView::draw(CDrawContext* pContext)
{
    pContext->setLineWidth(1);
    pContext->setFillColor(getBackColor());
    pContext->setFrameColor(getBorderColor());
    pContext->drawRect(getViewSize(), VSTGUI::kDrawFilledAndStroked);

    double MAX_FREQ = 22000.0;
    double MIN_FREQ = 10.0;
    double FREQ_LOG_MAX = log(MAX_FREQ / MIN_FREQ);
    double ceiling = 0.0;
    double noise_floor = -72.0;
    double DB_EQ_RANGE = 15.0;

// Given frequency, return screen x position
#define freq_to_x(view, freq) \
(view.getWidth() * log(freq / MIN_FREQ) / FREQ_LOG_MAX)

// Given screen x position, return frequency
#define x_to_freq(view, x) \
std::max(std::min(MIN_FREQ * exp(FREQ_LOG_MAX * x / view.getWidth()), MAX_FREQ), MIN_FREQ)

// Given a magnitude, return y screen position as 0..1 with applied tilt
#define mag_to_01(m, freq) \
1.0 - (( ((20 * log10(m)) + (4.5 * ((log(freq) / log(2)) - (log(1024) / log(2))))) - ceiling) / (noise_floor - ceiling));

// Given a magnitude (1.0 .... very small number), return y screen position
#define mag_to_y(view, m) \
((((20.0 * log10(m)) - ceiling) / (noise_floor - ceiling)) * (view.getHeight()))

// Given decibels, return screen y position
#define db_to_y(view, db) \
(((db - ceiling) / (noise_floor - ceiling)) * view.getHeight())

// Given screen y position, return decibels
#define y_to_db(view, y) \
ceiling + ((y / view.getHeight()) * (noise_floor - ceiling));

#define dB_to_y_EQ(view, dB) \
    view.getHeight() * (1.0 - (((dB / DB_EQ_RANGE) / 2) + 0.5));

    auto border = getBorderColor();
    border.setNormAlpha(0.5);

    {
        VSTGUI::CRect r(getViewSize());
        pContext->setFrameColor(border);
        for (int x = 2; x < 10; x++) {
            VSTGUI::CCoord Hz_10 = freq_to_x(r, 10.0 * x);
            const VSTGUI::CPoint _p1(r.left + Hz_10, r.bottom);
            const VSTGUI::CPoint _p2(r.left + Hz_10, r.top);
            pContext->drawLine(_p1, _p2);
        }
        for (int x = 1; x < 10; x++) {
            VSTGUI::CCoord Hz_100 = freq_to_x(r, 100.0 * x);
            const VSTGUI::CPoint _p1(r.left + Hz_100, r.bottom);
            const VSTGUI::CPoint _p2(r.left + Hz_100, r.top);
            pContext->drawLine(_p1, _p2);
        }
        for (int x = 1; x < 10; x++) {
            VSTGUI::CCoord Hz_1000 = freq_to_x(r, 1000.0 * x);
            const VSTGUI::CPoint _p1(r.left + Hz_1000, r.bottom);
            const VSTGUI::CPoint _p2(r.left + Hz_1000, r.top);
            pContext->drawLine(_p1, _p2);
        }

        for (int x = 1; x < 3; x++) {
            VSTGUI::CCoord Hz_10000 = freq_to_x(r, 10000.0 * x);
            const VSTGUI::CPoint _p1(r.left + Hz_10000, r.bottom);
            const VSTGUI::CPoint _p2(r.left + Hz_10000, r.top);
            pContext->drawLine(_p1, _p2);
        }
    }

    {
        VSTGUI::CRect r(getViewSize());
        pContext->setFrameColor(border);

        for (int cnt = -(int)DB_EQ_RANGE; cnt < (int)DB_EQ_RANGE; cnt += 5)
        {
            VSTGUI::CCoord dB = dB_to_y_EQ(r, cnt);
            const VSTGUI::CPoint _p1(r.left,  r.bottom - dB);
            const VSTGUI::CPoint _p2(r.right, r.bottom - dB);
            pContext->drawLine(_p1, _p2);
        }
    }


    VSTGUI::CGraphicsPath* FFT_curve = pContext->createGraphicsPath();
    
    if (FFT_curve)
    {
        VSTGUI::CRect r(getViewSize());

        double y_start = mag_to_01(fft_linear[0], fft_freq[0]);
        y_start = (std::max)((std::min)(y_start, 1.0), 0.0);
        y_start *= r.getHeight();
        FFT_curve->beginSubpath(VSTGUI::CPoint(r.left - 1, r.bottom - y_start));
        double x_last = 0.0;
        // RAW
        for (int bin = 0; bin < numBins; ++bin) {
            double x = freq_to_x(r, fft_freq[bin]);
            x = (std::max)((std::min)(x, r.getWidth()), 0.0);
            double y = mag_to_01(fft_linear[bin], fft_freq[bin]);
            y = (std::max)((std::min)(y, 1.0), 0.0);
            y *= r.getHeight();
            if (x - x_last > 0.1)
            {
                x_last = x;
                FFT_curve->addLine(VSTGUI::CPoint(r.left + x, r.bottom - y));
            }
        }

        FFT_curve->addLine(VSTGUI::CPoint(r.right + 1, r.bottom + 1));
        FFT_curve->addLine(VSTGUI::CPoint(r.left - 1, r.bottom + 1));
        FFT_curve->closeSubpath();


        VSTGUI::CColor ff = getFFTFillColor();
        ff.setNormAlpha(0.8);
        pContext->setFrameColor(VSTGUI::kTransparentCColor);
        pContext->setFillColor(ff);
        pContext->setDrawMode(VSTGUI::kAntiAliasing);
        pContext->setLineWidth(0.0);
        pContext->setLineStyle(VSTGUI::kLineSolid);
        pContext->drawGraphicsPath(FFT_curve, VSTGUI::CDrawContext::kPathFilled);


        pContext->setFrameColor(getFFTLineColor());
        pContext->setDrawMode(VSTGUI::kAntiAliasing);
        pContext->setLineWidth(1.0);
        pContext->setLineStyle(VSTGUI::kLineSolid);
        pContext->drawGraphicsPath(FFT_curve, VSTGUI::CDrawContext::kPathStroked);

        FFT_curve->forget();
    }

    VSTGUI::CGraphicsPath* EQ_curve = pContext->createGraphicsPath();
    if (EQ_curve)
    {
        VSTGUI::CRect r(getViewSize());

        VSTGUI::CCoord y_mid = r.bottom - (r.getHeight() / 2.0);
        EQ_curve->beginSubpath(VSTGUI::CPoint(r.left - 1, y_mid));
        static constexpr double inc = 0.05;
        for (double x = -inc; x <= r.getWidth() + inc; x += inc)
        {
            double tmp = MIN_FREQ * exp(FREQ_LOG_MAX * x / r.getWidth());
            double freq = (std::max)((std::min)(tmp, MAX_FREQ), MIN_FREQ);

            double dB_LF = 20 * log10(LF_SVF.mag_response(freq));
            double dB_HF = 20 * log10(HF_SVF.mag_response(freq));
            double dB = dB_LF + dB_HF;

            double m = 1.0 - (((dB / DB_EQ_RANGE) / 2) + 0.5);
            double scy = m * r.getHeight();

            EQ_curve->addLine(VSTGUI::CPoint(r.left + x, r.top + scy));
        }
        EQ_curve->addLine(VSTGUI::CPoint(r.right + 1, r.bottom + 1));
        EQ_curve->addLine(VSTGUI::CPoint(r.left - 1, r.bottom + 1));
        EQ_curve->closeSubpath();

        pContext->setFrameColor(getLineColor());
        pContext->setDrawMode(VSTGUI::kAntiAliasing);
        pContext->setLineWidth(1.5);
        pContext->setLineStyle(VSTGUI::kLineSolid);
        pContext->drawGraphicsPath(EQ_curve, VSTGUI::CDrawContext::kPathStroked);
        EQ_curve->forget();
    }

    // box outline
    pContext->setLineWidth(1);
    pContext->setFrameColor(getBorderColor());
    pContext->drawRect(getViewSize(), VSTGUI::kDrawStroked);

    setDirty(false);
};

bool EQCurveView::sizeToFit() {
    if (getDrawBackground())
    {
        CRect vs(getViewSize());
        vs.setWidth(getDrawBackground()->getWidth());
        vs.setHeight(getDrawBackground()->getHeight());
        setViewSize(vs);
        setMouseableArea(vs);
        return true;
    }
    return false;
};

//------------------------------------------------------------------------
//  Transfer Curve View
//------------------------------------------------------------------------
void TransferCurveView::draw(CDrawContext* _pContext)
{
    CDrawContext* pContext = _pContext;

    // draw border
    pContext->setLineWidth(1);
    pContext->setFillColor(BackColor);
    pContext->setFrameColor(VSTGUI::CColor(223, 233, 233, 255));
    pContext->drawRect(getViewSize(), VSTGUI::kDrawFilledAndStroked);

    // draw db lines
    
    //double FREQ_LOG_MAX = log(MAX_FREQ / MIN_FREQ);
    //double DB_EQ_RANGE = 15.0;
    {
        VSTGUI::CRect r(getViewSize());
        auto width = r.getWidth();
        auto height = r.getHeight();

        pContext->setFrameColor(VSTGUI::CColor(255, 255, 255, 55));
        for (int x = MIN_dB; x < MAX_dB; x += 10) {
            VSTGUI::CCoord ver = width * -x * Inv_DB_R;
            const VSTGUI::CPoint _p1(r.left + ver, r.bottom);
            const VSTGUI::CPoint _p2(r.left + ver, r.top);
            pContext->drawLine(_p1, _p2);
        }

        for (int x = MIN_dB; x < MAX_dB; x += 10) {
            VSTGUI::CCoord ver = height * -x * Inv_DB_R;
            const VSTGUI::CPoint _p1(r.left, r.bottom - ver);
            const VSTGUI::CPoint _p2(r.right, r.bottom - ver);
            pContext->drawLine(_p1, _p2);
        }
    }


    // draw curve
    VSTGUI::CGraphicsPath* pathFull = pContext->createGraphicsPath();
    VSTGUI::CGraphicsPath* pathInput = pContext->createGraphicsPath();
    if (pathFull && pathInput)
    {
        VSTGUI::CRect r(getViewSize());
        // VSTGUI::CCoord inset = 30;
        // r.inset(inset, 0);

        double width = r.getWidth();
        double inv_width = 1.0 / r.getWidth();
        double height = r.getHeight();

        pathFull->beginSubpath(VSTGUI::CPoint(r.left - 1, r.bottom));
        pathInput->beginSubpath(VSTGUI::CPoint(r.left - 1, r.bottom));
        std::vector<double> pp_x, pp_y;
        for (int x = -1; x <= width + 1; x++)
        {
            // x * inv_width = [0 ~ 1]
            // 1.0 - x * inv_width = [1 ~ 0]
            // -DB_Range == -60
            // -DB_Range * (1.0 - x * inv_width) = [-60 ~ 0]
            double x_dB = -DB_Range * (1.0 - x * inv_width);

            double overshoot = x_dB - threshold;
            
            double gain = 0.0;
            if (overshoot <= -kneeHalf)
                gain = 0.0;
            else if (overshoot > -kneeHalf && overshoot <= kneeHalf)
                gain = 0.5 * slope * ((overshoot + kneeHalf) * (overshoot + kneeHalf)) / knee;
            else
                gain = slope * overshoot;

            double y_dB = x_dB + gain + makeup;
            double a = yg331::DecibelConverter::ToGain(y_dB);
            double b = yg331::DecibelConverter::ToGain(x_dB);
            y_dB = mix * a + (1.0 - mix) * b; // -60 ~ 0
            y_dB = yg331::DecibelConverter::ToDecibel(y_dB);
            
            double y = -y_dB * Inv_DB_R; // 1 ~ 0
            y = (1.0 - y) * height;

            pathFull->addLine(VSTGUI::CPoint(r.left + x, r.bottom - y));
            if (x_dB < vuInMono)
            {
                pp_x.push_back(x);
                pp_y.push_back(y);
                pathInput->addLine(VSTGUI::CPoint(r.left + x, r.bottom - y));
            }
        }
        while (!pp_x.empty())
        {
            double x = *(pp_x.end()-1);
            double y = *(pp_y.end()-1);
            pathInput->addLine(VSTGUI::CPoint(r.left + x, r.bottom - y));
            pp_x.pop_back();
            pp_y.pop_back();
        }
        pathFull->addLine(VSTGUI::CPoint(r.right + 1, r.bottom + 1));
        pathFull->addLine(VSTGUI::CPoint(r.left - 1, r.bottom + 1));
        pathFull->closeSubpath();
        pathInput->closeSubpath();
        
        pContext->setFrameColor(LineOffColor);
        pContext->setDrawMode(VSTGUI::kAntiAliasing);
        pContext->setLineWidth(2.0);
        pContext->setLineStyle(VSTGUI::kLineSolid);
        pContext->drawGraphicsPath(pathFull, VSTGUI::CDrawContext::kPathStroked);
        pathFull->forget();
        
        pContext->setFrameColor(LineColor);
        pContext->setDrawMode(VSTGUI::kAntiAliasing);
        pContext->setLineWidth(2.0);
        pContext->setLineStyle(VSTGUI::kLineSolid);
        pContext->drawGraphicsPath(pathInput, VSTGUI::CDrawContext::kPathStroked);
        pathInput->forget();
    }
    
    // draw border
    pContext->setLineWidth(1);
    pContext->setFillColor(kTransparentCColor);
    pContext->setFrameColor(VSTGUI::CColor(223, 233, 233, 255)); // black borders
    pContext->drawRect(getViewSize(), VSTGUI::kDrawStroked);


    setDirty(false);
};
bool TransferCurveView::sizeToFit()
{
    if (getDrawBackground())
    {
        CRect vs(getViewSize());
        vs.setWidth(getDrawBackground()->getWidth());
        vs.setHeight(getDrawBackground()->getHeight());
        setViewSize(vs);
        setMouseableArea(vs);
        return true;
    }
    return false;
};

//------------------------------------------------------------------------
// MyVuMeter
//------------------------------------------------------------------------
MyVuMeter::MyVuMeter(const CRect& size, int32_t _style = kVertical)
    : CControl(size, nullptr, 0)
    , style(_style)
{
    vuOnColor = kWhiteCColor;
    vuOffColor = kBlackCColor;

    rectOn(size.left, size.top, size.right, size.bottom);
    rectOff(size.left, size.top, size.right, size.bottom);

    setWantsIdle(true);
};
MyVuMeter::MyVuMeter(const MyVuMeter& vuMeter)
    : CControl(vuMeter)
    , style(vuMeter.style)
    , vuOnColor(vuMeter.vuOnColor)
    , vuOffColor(vuMeter.vuOffColor)
    , rectOn(vuMeter.rectOn)
    , rectOff(vuMeter.rectOff)
{
    setWantsIdle(true);
};

void MyVuMeter::draw(CDrawContext* _pContext) 
{
    CRect _rectOn(rectOn);
    CRect _rectOff(rectOff);
    CPoint pointOn;
    CPoint pointOff;
    CDrawContext* pContext = _pContext;

    bounceValue();

    // float newValue = getValueNormalized(); // normalize
    float newValue = (plainValue - getMin ()) / static_cast<float> (getMax () - getMin ());

    if (style & kHorizontal)
    {
        auto tmp = (CCoord)((int32_t)(newValue * getViewSize().getWidth()));
        // pointOff(tmp, 0); // x, y
        // pointOn(tmp, 0); // x, y

        // _rectOff.left += tmp;
        // _rectOn.right = tmp + rectOn.left;
        _rectOn.right = _rectOn.left + tmp;
    }
    else
    {
        auto tmp = (CCoord)((int32_t)(newValue * getViewSize().getHeight()));
        pointOn(0, tmp); // x, y

        //_rectOff.bottom = tmp + rectOff.top;
        //_rectOn.top += tmp;
        _rectOn.top = _rectOff.bottom - tmp;
    }

    pContext->setFillColor(vuOffColor);
    pContext->drawRect(rectOff, kDrawFilled);

    pContext->setFillColor(vuOnColor);
    pContext->drawRect(_rectOn, kDrawFilled);

    setDirty(false);
};
void MyVuMeter::setViewSize(const CRect& newSize, bool invalid = true) 
{
    CControl::setViewSize(newSize, invalid);
    rectOn = getViewSize();
    rectOff = getViewSize();
};
bool MyVuMeter::sizeToFit() 
{
    if (getDrawBackground())
    {
        CRect vs(getViewSize());
        vs.setWidth(getDrawBackground()->getWidth());
        vs.setHeight(getDrawBackground()->getHeight());
        setViewSize(vs);
        setMouseableArea(vs);
        return true;
    }
    return false;
};



static const std::string kAttrBackColor    = "back-color";
static const std::string kAttrBorderColor  = "border-color";
static const std::string kAttrLineColor    = "line-color";
static const std::string kAttrLineOffColor = "line-off-color";
static const std::string kAttrFFTLineColor = "FFT-line-color";
static const std::string kAttrFFTFillColor = "FFT-fill-color";

static const std::string kAttrVuOnColor  = "vu-on-color";
static const std::string kAttrVuOffColor = "vu-off-color";
static const std::string kAttrPDclick    = "click-behave";
static const std::string kAttrPDMin      = "update-min";
static const std::string kAttrPDMax      = "update-max";

//------------------------------------------------------------------------
// EQ Curve View Factory
//------------------------------------------------------------------------
class MyEQCurveViewFactory : public ViewCreatorAdapter
{
public:
    //register this class with the view factory
    MyEQCurveViewFactory() { UIViewFactory::registerViewCreator(*this); }

    //return an unique name here
    IdStringPtr getViewName() const override { return "EQ Curve View"; }

    //return the name here from where your custom view inherites.
    //    Your view automatically supports the attributes from it.
    IdStringPtr getBaseViewName() const override { return UIViewCreator::kCView; }

    //create your view here.
    //    Note you don't need to apply attributes here as
    //    the apply method will be called with this new view
    CView* create(const UIAttributes& attributes, const IUIDescription* description) const override
    {
        return new EQCurveView(CRect(0, 0, 100, 20), nullptr, -1, nullptr);
    }
    bool apply(
        CView* view,
        const UIAttributes& attributes,
        const IUIDescription* description) const override
    {
        auto* v = dynamic_cast<EQCurveView*> (view);

        if (!v)
            return false;

        CColor color;
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrBackColor), color, description))
            v->setBackColor(color);
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrBorderColor), color, description))
            v->setBorderColor(color);
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrLineColor), color, description))
            v->setLineColor(color);
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrFFTLineColor), color, description))
            v->setFFTLineColor(color);
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrFFTFillColor), color, description))
            v->setFFTFillColor(color);

        return true;
    }

    bool getAttributeNames(StringList& attributeNames) const override
    {
        attributeNames.emplace_back(kAttrBackColor);
        attributeNames.emplace_back(kAttrBorderColor);
        attributeNames.emplace_back(kAttrLineColor);
        attributeNames.emplace_back(kAttrFFTLineColor);
        attributeNames.emplace_back(kAttrFFTFillColor);
        return true;
    }

    AttrType getAttributeType(const std::string& attributeName) const override
    {
        if (attributeName == kAttrBackColor)
            return kColorType;
        if (attributeName == kAttrBorderColor)
            return kColorType;
        if (attributeName == kAttrLineColor)
            return kColorType;
        if (attributeName == kAttrFFTLineColor)
            return kColorType;
        if (attributeName == kAttrFFTFillColor)
            return kColorType;
        return kUnknownType;
    }

    //------------------------------------------------------------------------
    bool getAttributeValue(
        CView* view,
        const string& attributeName,
        string& stringValue,
        const IUIDescription* desc) const override
    {
        auto* v = dynamic_cast<EQCurveView*> (view);

        if (!v)
            return false;

        if (attributeName == kAttrBackColor)
        {
            UIViewCreator::colorToString(v->getBackColor(), stringValue, desc);
            return true;
        }
        else if (attributeName == kAttrBorderColor)
        {
            UIViewCreator::colorToString(v->getBorderColor(), stringValue, desc);
            return true;
        }
        else if (attributeName == kAttrLineColor)
        {
            UIViewCreator::colorToString(v->getLineColor(), stringValue, desc);
            return true;
        }
        else if (attributeName == kAttrFFTLineColor)
        {
            UIViewCreator::colorToString(v->getFFTLineColor(), stringValue, desc);
            return true;
        }
        else if (attributeName == kAttrFFTFillColor)
        {
            UIViewCreator::colorToString(v->getFFTFillColor(), stringValue, desc);
            return true;
        }

        return false;
    }
};

//------------------------------------------------------------------------
// Transfer Curve View Factory
//------------------------------------------------------------------------
class MyTransferCurveControlFactory : public ViewCreatorAdapter
{
public:
    MyTransferCurveControlFactory() { UIViewFactory::registerViewCreator(*this); }
    IdStringPtr getViewName() const override { return "Transfer Curve Viewer"; }
    IdStringPtr getBaseViewName() const override { return UIViewCreator::kCControl; }
    CView* create(const UIAttributes& attributes, const IUIDescription* description) const override
    {
        CRect size(CPoint(45, 45), CPoint(400, 150));
        return new TransferCurveView(size);
    }
    bool apply(
        CView* view,
        const UIAttributes& attributes,
        const IUIDescription* description) const override
    {
        auto* transferCurveView = dynamic_cast<TransferCurveView*> (view);

        if (!transferCurveView)
            return false;

        CColor color;
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrBackColor), color, description))
            transferCurveView->setBackColor(color);
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrLineColor), color, description))
            transferCurveView->setLineColor(color);
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrLineOffColor), color, description))
            transferCurveView->setLineOffColor(color);

        return true;
    }

    bool getAttributeNames(StringList& attributeNames) const override
    {
        attributeNames.emplace_back(kAttrBackColor);
        attributeNames.emplace_back(kAttrLineColor);
        attributeNames.emplace_back(kAttrLineOffColor);
        return true;
    }

    AttrType getAttributeType(const std::string& attributeName) const override
    {
        if (attributeName == kAttrBackColor)
            return kColorType;
        if (attributeName == kAttrLineColor)
            return kColorType;
        if (attributeName == kAttrLineOffColor)
            return kColorType;
        return kUnknownType;
    }

    //------------------------------------------------------------------------
    bool getAttributeValue(
        CView* view,
        const string& attributeName,
        string& stringValue,
        const IUIDescription* desc) const override
    {
        auto* transferCurveView = dynamic_cast<TransferCurveView*> (view);

        if (!transferCurveView)
            return false;

        if (attributeName == kAttrBackColor)
        {
            UIViewCreator::colorToString(transferCurveView->getBackColor(), stringValue, desc);
            return true;
        }
        else if (attributeName == kAttrLineColor)
        {
            UIViewCreator::colorToString(transferCurveView->getLineColor(), stringValue, desc);
            return true;
        }
        else if (attributeName == kAttrLineOffColor)
        {
            UIViewCreator::colorToString(transferCurveView->getLineOffColor(), stringValue, desc);
            return true;
        }
        return false;
    }
};

//------------------------------------------------------------------------
// Click Reset Parameter Display View Factory
//------------------------------------------------------------------------
class MyClickResetParamDisplayFactory : public ViewCreatorAdapter
{
public:
    MyClickResetParamDisplayFactory() { UIViewFactory::registerViewCreator(*this); }
    IdStringPtr getViewName() const override { return "Click to Reset Parameter Display"; }
    IdStringPtr getBaseViewName() const override { return UIViewCreator::kCParamDisplay; }
    CView* create(const UIAttributes& attributes, const IUIDescription* description) const override
    {
        CRect ss(0, 0, 100, 20);
        return new ClickResetParamDisplay(ss, nullptr, 0);
    }
    bool apply(
        CView* view,
        const UIAttributes& attributes,
        const IUIDescription* description) const SMTG_OVERRIDE
    {
        auto* vv = dynamic_cast<ClickResetParamDisplay*> (view);

        if (!vv)
            return false;

        const auto* attr = attributes.getAttributeValue(kAttrPDclick);
        if (attr)
            vv->setStyle_(*attr == kAttrPDMin ? ClickResetParamDisplay::kUpdateMin : ClickResetParamDisplay::kUpdateMax);

        return true;
    }

    bool getAttributeNames(StringList& attributeNames) const SMTG_OVERRIDE
    {
        attributeNames.emplace_back(kAttrPDclick);
        return true;
    }

    AttrType getAttributeType(const std::string& attributeName) const SMTG_OVERRIDE
    {
        if (attributeName == kAttrPDclick)
            return kListType;
        return kUnknownType;
    }

    //------------------------------------------------------------------------
    bool getAttributeValue(
        CView* view,
        const string& attributeName,
        string& stringValue,
        const IUIDescription* desc) const SMTG_OVERRIDE
    {
        auto* vv = dynamic_cast<ClickResetParamDisplay*> (view);

        if (!vv)
            return false;

        if (attributeName == kAttrPDclick)
        {
            if (vv->getStyle_() & ClickResetParamDisplay::kUpdateMin)
                stringValue = kAttrPDMin;
            else
                stringValue = kAttrPDMax;
            return true;
        }

        return false;
    }

    //------------------------------------------------------------------------
    bool getPossibleListValues(
        const string& attributeName,
        ConstStringPtrList& values) const SMTG_OVERRIDE
    {
        if (attributeName == kAttrPDclick)
        {
            values.emplace_back(&kAttrPDMin);
            values.emplace_back(&kAttrPDMax);
            return true;
        }
        return false;
    }
};

//------------------------------------------------------------------------
// VU Meter View Factory
//------------------------------------------------------------------------
class MyVUMeterFactory : public ViewCreatorAdapter
{
public:
    //register this class with the view factory
    MyVUMeterFactory() { UIViewFactory::registerViewCreator(*this); }
    
    //return an unique name here
    IdStringPtr getViewName() const override { return "My Vu Meter"; }
    
    //return the name here from where your custom view inherites.
    //    Your view automatically supports the attributes from it.
    IdStringPtr getBaseViewName() const override { return UIViewCreator::kCControl; }
    
    //create your view here.
    //    Note you don't need to apply attributes here as
    //    the apply method will be called with this new view
    CView* create(const UIAttributes & attributes, const IUIDescription * description) const override
    {
        CRect size(CPoint(45, 45), CPoint(400, 150));
        return new MyVuMeter(size, 2);
    }
    
    // apply custom attributes to your view
    bool apply(CView* view, const UIAttributes& attributes, const IUIDescription* description) const SMTG_OVERRIDE
    {
        auto* vuMeter = dynamic_cast<MyVuMeter*> (view);
        
        if (!vuMeter)
            return false;
        
        const auto* attr = attributes.getAttributeValue(UIViewCreator::kAttrOrientation);
        if (attr)
            vuMeter->setStyle(*attr == UIViewCreator::strVertical ? MyVuMeter::kVertical : MyVuMeter::kHorizontal);
        
        CColor color;
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrVuOnColor), color, description))
            vuMeter->setVuOnColor(color);
        if (UIViewCreator::stringToColor(attributes.getAttributeValue(kAttrVuOffColor), color, description))
            vuMeter->setVuOffColor(color);
        
        return true;
    }
    
    // add your custom attributes to the list
    bool getAttributeNames(StringList& attributeNames) const SMTG_OVERRIDE
    {
        attributeNames.emplace_back(UIViewCreator::kAttrOrientation);
        attributeNames.emplace_back(kAttrVuOnColor);
        attributeNames.emplace_back(kAttrVuOffColor);
        return true;
    }
    
    // return the type of your custom attributes
    AttrType getAttributeType(const std::string& attributeName) const SMTG_OVERRIDE
    {
        if (attributeName == UIViewCreator::kAttrOrientation)
            return kListType;
        if (attributeName == kAttrVuOnColor)
            return kColorType;
        if (attributeName == kAttrVuOffColor)
            return kColorType;
        return kUnknownType;
    }
    
    // return the string value of the custom attributes of the view
    bool getAttributeValue(
                           CView* view,
                           const string& attributeName,
                           string& stringValue,
                           const IUIDescription* desc) const SMTG_OVERRIDE
    {
        auto* vuMeter = dynamic_cast<MyVuMeter*> (view);
        
        if (!vuMeter)
            return false;
        
        if (attributeName == UIViewCreator::kAttrOrientation)
        {
            if (vuMeter->getStyle() & MyVuMeter::kVertical)
                stringValue = UIViewCreator::strVertical;
            else
                stringValue = UIViewCreator::strHorizontal;
            return true;
        }
        else if (attributeName == kAttrVuOnColor)
        {
            UIViewCreator::colorToString(vuMeter->getVuOnColor(), stringValue, desc);
            return true;
        }
        else if (attributeName == kAttrVuOffColor)
        {
            UIViewCreator::colorToString(vuMeter->getVuOffColor(), stringValue, desc);
            return true;
        }
        return false;
    }
    
    //------------------------------------------------------------------------
    bool getPossibleListValues(
                               const string& attributeName,
                               ConstStringPtrList& values) const SMTG_OVERRIDE
    {
        if (attributeName == UIViewCreator::kAttrOrientation)
        {
            return UIViewCreator::getStandardAttributeListValues(UIViewCreator::kAttrOrientation, values);
        }
        return false;
    }
};

//create a static instance so that it registers itself with the view factory
MyEQCurveViewFactory            __gMyEQCurveViewFactory;
MyTransferCurveControlFactory   __gMyTransferCurveControlFactory;
MyClickResetParamDisplayFactory __gMyClickResetParamDisplayFactory;
MyVUMeterFactory                __gMyVUMeterFactory;
} // namespace VSTGUI


namespace yg331 {
//------------------------------------------------------------------------
// VuMeterController
//------------------------------------------------------------------------
VSTGUI::CView* VuMeterController::verifyView(CView* view,
                  const UIAttributes&   /*attributes*/,
                  const IUIDescription* /*description*/)
{
#define minVU   -60.0
#define minGRVU -20.0
#define maxVU     0.0
    if (MyVuMeter* control = dynamic_cast<MyVuMeter*>(view); control)
    {
        if (control->getTag() == kInLRMS || control->getTag() == kInLPeak)  {
            vuMeterInL  = control;
            vuMeterInL->setMin(minVU);
            vuMeterInL->setMax(maxVU);
            vuMeterInL->setDefaultValue(minVU);
            //vuMeterInL->registerViewListener(this);
        }
        if (control->getTag() == kInRRMS || control->getTag() == kInRPeak)  {
            vuMeterInR  = control;
            vuMeterInR->setMin(minVU);
            vuMeterInR->setMax(maxVU);
            vuMeterInR->setDefaultValue(minVU);
            //vuMeterInR->registerViewListener(this);
        }
        if (control->getTag() == kOutLRMS || control->getTag() == kOutLPeak) {
            vuMeterOutL = control;
            vuMeterOutL->setMin(minVU);
            vuMeterOutL->setMax(maxVU);
            vuMeterOutL->setDefaultValue(minVU);
            //vuMeterOutL->registerViewListener(this);
        }
        if (control->getTag() == kOutRRMS || control->getTag() == kOutRPeak) {
            vuMeterOutR = control;
            vuMeterOutR->setMin(minVU);
            vuMeterOutR->setMax(maxVU);
            vuMeterOutR->setDefaultValue(minVU);
            //vuMeterOutR->registerViewListener(this);
        }
        if (control->getTag() == kGainReduction)   {
            vuMeterGR = control;
            vuMeterGR->setMin(minGRVU);
            vuMeterGR->setMax(maxVU);
            vuMeterGR->setDefaultValue(maxVU);
            //vuMeterGR->registerViewListener(this);
        }
    }
    if (ClickResetParamDisplay* control = dynamic_cast<ClickResetParamDisplay*>(view); control)
    {
        if (control->getTag() == kGainReduction)   {
            pdGainReduction = control;
            pdGainReduction->setMin(minGRVU);
            pdGainReduction->setMax(maxVU);
            pdGainReduction->setDefaultValue(maxVU);
            pdGainReduction->setPrecision(1);
            // pdGainReduction->registerViewListener(this);
        }
    }
    

    return view;
};

//------------------------------------------------------------------------
// LogRangeParameter Declaration
//------------------------------------------------------------------------
class LogRangeParameter : public Vst::RangeParameter
{
public:
    using RangeParameter::RangeParameter;
    
    LogRangeParameter (const Vst::TChar* title, Vst::ParamID tag, const Vst::TChar* units = nullptr,
                       Vst::ParamValue minPlain = 0., Vst::ParamValue maxPlain = 1.,
                       Vst::ParamValue defaultValuePlain = 0., int32 stepCount = 0,
                       int32 flags = Steinberg::Vst::ParameterInfo::kCanAutomate, Vst::UnitID unitID = Steinberg::Vst::kRootUnitId,
                       const Vst::TChar* shortTitle = nullptr)
    : Vst::RangeParameter(title, tag, units, minPlain, maxPlain, defaultValuePlain, stepCount, flags, unitID, shortTitle)
    {
        UString (info.title, str16BufferSize (Vst::String128)).assign (title);
        if (units)
            UString (info.units, str16BufferSize (Vst::String128)).assign (units);
        if (shortTitle)
            UString (info.shortTitle, str16BufferSize (Vst::String128)).assign (shortTitle);

        info.stepCount = stepCount;
        info.defaultNormalizedValue = valueNormalized = toNormalized (defaultValuePlain);
        info.flags = flags;
        info.id = tag;
        info.unitId = unitID;
    }
    
    /** Converts a normalized value to plain value (e.g. 0.5 to 10000.0Hz). */
    Vst::ParamValue toPlain(Vst::ParamValue _valueNormalized) const SMTG_OVERRIDE;
    
    /** Converts a plain value to a normalized value (e.g. 10000 to 0.5). */
    Vst::ParamValue toNormalized(Vst::ParamValue plainValue) const SMTG_OVERRIDE;
    
    /** Converts a normalized value to a string. */
    void toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const SMTG_OVERRIDE;
    
    OBJ_METHODS (LogRangeParameter, RangeParameter)
};
//------------------------------------------------------------------------
// LogRangeParameter Implementation
//------------------------------------------------------------------------
Vst::ParamValue LogRangeParameter::toPlain(Vst::ParamValue _valueNormalized) const
{
    double FREQ_LOG_MAX = std::log(getMax() / getMin());
    double tmp = getMin() * std::exp(FREQ_LOG_MAX * _valueNormalized);
    double freq = (std::max)((std::min)(tmp, getMax()), getMin());
    return freq;
    //return _valueNormalized * (getMax() - getMin()) + getMin();
}

//------------------------------------------------------------------------
Vst::ParamValue LogRangeParameter::toNormalized(Vst::ParamValue plainValue) const
{
    SMTG_ASSERT(getMax() - getMin() != 0);
    double FREQ_LOG_MAX = std::log(getMax() / getMin());
    return std::log(plainValue / getMin()) / FREQ_LOG_MAX;
    //return (plainValue - getMin()) / (getMax() - getMin());
}

void LogRangeParameter::toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const
{
    {
        //Parameter::toString(toPlain(_valueNormalized), string);
        UString wrapper(string, str16BufferSize(Vst::String128));
        {
            if (!wrapper.printFloat(toPlain(_valueNormalized), precision))
                string[0] = 0;
            wrapper.append(STR16(" "));
            wrapper.append(getInfo().units);
        }
    }
}

//------------------------------------------------------------------------
// LinRangeParameter Declaration
//------------------------------------------------------------------------
class LinRangeParameter : public Vst::RangeParameter
{
public:
    using RangeParameter::RangeParameter;
    void toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const SMTG_OVERRIDE;
    
    OBJ_METHODS (LinRangeParameter, RangeParameter)
};
//------------------------------------------------------------------------
// LinRangeParameter Implementation
//------------------------------------------------------------------------
void LinRangeParameter::toString(Vst::ParamValue _valueNormalized, Vst::String128 string) const
{
    {
        //Parameter::toString(toPlain(_valueNormalized), string);
        UString wrapper(string, str16BufferSize(Vst::String128));
        {
            if (!wrapper.printFloat(toPlain(_valueNormalized), precision))
                string[0] = 0;
            wrapper.append(STR16(" "));
            wrapper.append(getInfo().units);
        }
    }
}

//------------------------------------------------------------------------
// RLFCMP_Controller Implementation
//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Controller::initialize (FUnknown* context)
{
    // Here the Plug-in will be instantiated

    //---do not forget to call parent ------
    tresult result = EditControllerEx1::initialize (context);
    if (result != kResultOk)
    {
        return result;
    }

    // Here you could register some parameters
    int32 stepCount;
    int32 flags;
    int32 tag;
    Vst::ParamValue defaultVal;
    Vst::ParamValue minPlain;
    Vst::ParamValue maxPlain;
    Vst::ParamValue defaultPlain;

    tag          = kParamBypass;
    stepCount    = 1;
    defaultVal   = 0;
    flags        = Vst::ParameterInfo::kIsBypass | Vst::ParameterInfo::kCanAutomate;
    parameters.addParameter(STR16("Bypass"), nullptr, stepCount, defaultVal, flags, tag);
    
    tag          = kParamScLfIn;
    auto* ParamScLfIn = new Vst::StringListParameter(STR16("SC LF In"), tag);
    ParamScLfIn->appendString (STR16("OFF"));
    ParamScLfIn->appendString (STR16("ON"));
    ParamScLfIn->getInfo().defaultNormalizedValue = 1.0;
    ParamScLfIn->setNormalized(1.0);
    parameters.addParameter (ParamScLfIn);
    
    tag          = kParamScLfType;
    auto* ParamScLfType = new Vst::StringListParameter(STR16("SC LF Type"), tag);
    ParamScLfType->appendString (STR16("CUT"));
    ParamScLfType->appendString (STR16("SHELF"));
    ParamScLfType->getInfo().defaultNormalizedValue = 0.0;
    parameters.addParameter (ParamScLfType);
    
    tag          = kParamScLfFreq;
    flags        = Vst::ParameterInfo::kCanAutomate;
    minPlain     = minScLfFreq;
    maxPlain     = maxScLfFreq;
    defaultPlain = dftScLfFreq;
    stepCount    = 0;
    auto* ParamScLfFreq = new LogRangeParameter(STR16("SC LF Freq"), tag, STR16("Hz"), minPlain, maxPlain, defaultPlain, stepCount, flags);
    ParamScLfFreq->setPrecision(0);
    parameters.addParameter(ParamScLfFreq);
    
    tag          = kParamScLfGain;
    flags        = Vst::ParameterInfo::kCanAutomate;
    minPlain     = minScLfGain;
    maxPlain     = maxScLfGain;
    defaultPlain = dftScLfGain;
    stepCount    = 0;
    auto* ParamScLfGain = new LinRangeParameter(STR16("SC LF Gain"), tag, STR16("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
    ParamScLfGain->setPrecision(1);
    parameters.addParameter(ParamScLfGain);
    
    tag          = kParamScHfIn;
    auto* ParamScHfIn = new Vst::StringListParameter(STR16("SC HF In"), tag);
    ParamScHfIn->appendString (STR16("OFF"));
    ParamScHfIn->appendString (STR16("ON"));
    ParamScHfIn->getInfo().defaultNormalizedValue = 1.0;
    ParamScHfIn->setNormalized(1.0);
    parameters.addParameter (ParamScHfIn);
    
    tag          = kParamScHfType;
    auto* ParamScHfType = new Vst::StringListParameter(STR16("SC HF Type"), tag);
    ParamScHfType->appendString (STR16("CUT"));
    ParamScHfType->appendString (STR16("SHELF"));
    ParamScHfType->getInfo().defaultNormalizedValue = 1.0;
    ParamScHfType->setNormalized(1.0);
    parameters.addParameter (ParamScHfType);
    
    tag          = kParamScHfFreq;
    flags        = Vst::ParameterInfo::kCanAutomate;
    minPlain     = minScHfFreq;
    maxPlain     = maxScHfFreq;
    defaultPlain = dftScHfFreq;
    stepCount    = 0;
    auto* ParamScHfFreq = new LogRangeParameter(STR16("SC HF Freq"), tag, STR16("Hz"), minPlain, maxPlain, defaultPlain, stepCount, flags);
    ParamScHfFreq->setPrecision(0);
    parameters.addParameter(ParamScHfFreq);
    
    tag          = kParamScHfGain;
    flags        = Vst::ParameterInfo::kCanAutomate;
    minPlain     = minScHfGain;
    maxPlain     = maxScHfGain;
    defaultPlain = dftScHfGain;
    stepCount    = 0;
    auto* ParamScHfGain = new LinRangeParameter(STR16("SC HF Gain"), tag, STR16("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
    ParamScHfGain->setPrecision(1);
    parameters.addParameter(ParamScHfGain);
    
    tag          = kParamScListen;
    auto* ParamScListen = new Vst::StringListParameter(STR16("SC Listen"), tag);
    ParamScListen->appendString (STR16("OFF"));
    ParamScListen->appendString (STR16("ON"));
    parameters.addParameter (ParamScListen);
    
    tag          = kParamDetectorType;
    auto* ParamDetectorType = new Vst::StringListParameter(STR16("Detector Type"), tag);
    ParamDetectorType->appendString (STR16("PEAK"));
    ParamDetectorType->appendString (STR16("RMS"));
    parameters.addParameter (ParamDetectorType);
    
    tag          = kParamSidechainTopology;
    auto* ParamSidechainTopology = new Vst::StringListParameter(STR16("Sidechain Topology"), tag);
    ParamSidechainTopology->appendString (STR16("LIN"));
    ParamSidechainTopology->appendString (STR16("LOG"));
    ParamSidechainTopology->getInfo().defaultNormalizedValue = ParamSidechainTopology->toNormalized(1.0);
    ParamSidechainTopology->setNormalized(ParamSidechainTopology->toNormalized(1.0));
    parameters.addParameter (ParamSidechainTopology);
    
    tag          = kParamHilbertEnable;
    auto* ParamHilbertEnable = new Vst::StringListParameter(STR16("Hilbert Enable"), tag);
    ParamHilbertEnable->appendString (STR16("OFF"));
    ParamHilbertEnable->appendString (STR16("ON"));
    parameters.addParameter (ParamHilbertEnable);
    
    tag          = kParamLookaheadEnable;
    auto* ParamLookaheadEnable = new Vst::StringListParameter(STR16("Lookahead Enable"), tag);
    ParamLookaheadEnable->appendString (STR16("OFF"));
    ParamLookaheadEnable->appendString (STR16("ON"));
    ParamLookaheadEnable->getInfo().defaultNormalizedValue = ParamLookaheadEnable->toNormalized(1.0);
    ParamLookaheadEnable->setNormalized(ParamLookaheadEnable->toNormalized(1.0));
    parameters.addParameter (ParamLookaheadEnable);

    tag          = kParamAttack;
    flags        = Vst::ParameterInfo::kCanAutomate;
    minPlain     = minAttack;
    maxPlain     = maxAttack;
    defaultPlain = dftAttack;
    stepCount    = 0;
    auto* ParamAttack = new LogRangeParameter(STR16("Attack"), tag, STR16("ms"), minPlain, maxPlain, defaultPlain, stepCount, flags);
    ParamAttack->setPrecision(1);
    parameters.addParameter(ParamAttack);

    tag          = kParamRelease;
    flags        = Vst::ParameterInfo::kCanAutomate;
    minPlain     = minRelease;
    maxPlain     = maxRelease;
    defaultPlain = dftRelease;
    stepCount    = 0;
    auto* ParamRelease = new LogRangeParameter(STR16("Release"), tag, STR16("ms"), minPlain, maxPlain, defaultPlain, stepCount, flags);
    ParamRelease->setPrecision(1);
    parameters.addParameter(ParamRelease);

    tag          = kParamThreshold;
    flags        = Vst::ParameterInfo::kCanAutomate;
    minPlain     = minThreshold;
    maxPlain     = maxThreshold;
    defaultPlain = dftThreshold;
    stepCount    = 0;
    auto* ParamThreshold = new LinRangeParameter(STR16("Threshold"), tag, STR16("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
    ParamThreshold->setPrecision(1);
    parameters.addParameter(ParamThreshold);

    tag          = kParamRatio;
    flags        = Vst::ParameterInfo::kCanAutomate;
    minPlain     = minRatio;
    maxPlain     = maxRatio;
    defaultPlain = dftRatio;
    stepCount    = 0;
    auto* ParamRatio = new LogRangeParameter(STR16("Ratio"), tag, STR16(""), minPlain, maxPlain, defaultPlain, stepCount, flags);
    ParamRatio->setPrecision(1);
    parameters.addParameter(ParamRatio);

    tag          = kParamKnee;
    flags        = Vst::ParameterInfo::kCanAutomate;
    minPlain     = minKnee;
    maxPlain     = maxKnee;
    defaultPlain = dftKnee;
    stepCount    = 0;
    auto* ParamKnee = new LinRangeParameter(STR16("Knee"), tag, STR16("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
    ParamKnee->setPrecision(1);
    parameters.addParameter(ParamKnee);

    tag          = kParamMakeup;
    flags        = Vst::ParameterInfo::kCanAutomate;
    minPlain     = minMakeup;
    maxPlain     = maxMakeup;
    defaultPlain = dftMakeup;
    stepCount    = 0;
    auto* ParamMakeup = new LinRangeParameter(STR16("Makeup"), tag, STR16("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
    ParamMakeup->setPrecision(1);
    parameters.addParameter(ParamMakeup);

    tag          = kParamMix;
    flags        = Vst::ParameterInfo::kCanAutomate;
    minPlain     = minMix;
    maxPlain     = maxMix;
    defaultPlain = dftMix;
    stepCount    = 0;
    auto* ParamMix = new LinRangeParameter(STR16("Mix"), tag, STR16("%"), minPlain, maxPlain, defaultPlain, stepCount, flags);
    ParamMix->setPrecision(1);
    parameters.addParameter(ParamMix);
    
    tag          = kParamInput;
    flags        = Vst::ParameterInfo::kCanAutomate;
    minPlain     = minInput;
    maxPlain     = maxInput;
    defaultPlain = dftInput;
    stepCount    = 0;
    auto* ParamInput = new LinRangeParameter(STR16("Input"), tag, STR16("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
    ParamInput->setPrecision(1);
    parameters.addParameter(ParamInput);
    
    tag          = kParamOutput;
    flags        = Vst::ParameterInfo::kCanAutomate;
    minPlain     = minOutput;
    maxPlain     = maxOutput;
    defaultPlain = dftOutput;
    stepCount    = 0;
    auto* ParamOut = new LinRangeParameter(STR16("Output"), tag, STR16("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
    ParamOut->setPrecision(1);
    parameters.addParameter(ParamOut);

    tag          = kParamSoftBypass;
    stepCount    = 1;
    defaultVal   = 0;
    flags        = Vst::ParameterInfo::kCanAutomate | Vst::ParameterInfo::kIsList;
    parameters.addParameter(STR16("SoftBypass"), nullptr, stepCount, defaultVal, flags, tag);

    // GUI only parameter
    Vst::StringListParameter* zoomParameter = new Vst::StringListParameter(STR("Zoom"), kParamZoom);
    for (ZoomFactorVector::const_iterator it = zoomFactors.begin(), end = zoomFactors.end(); it != end; ++it)
    {
        zoomParameter->appendString(it->title);
    }
    zoomParameter->setNormalized(zoomParameter->toNormalized(dftZoom));
    zoomParameter->getInfo().defaultNormalizedValue = zoomParameter->toNormalized(dftZoom);
    zoomParameter->addDependent(this);
    uiParameters.addParameter(zoomParameter);
    

    tag          = paramControllerPage;
    auto* paramControllerPage = new Vst::StringListParameter(STR16("controllerPage"), tag);
    paramControllerPage->appendString (STR16("COMPRESSOR")); // first of StringList is default
    paramControllerPage->appendString (STR16("SIDECHAIN EQ"));
    paramControllerPage->addDependent(this);
    uiParameters.addParameter (paramControllerPage);
    
    return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Controller::terminate ()
{
    // Here the Plug-in will be de-instantiated, last possibility to remove some memory!
    getParameterObject(kParamZoom)->removeDependent(this);
    getParameterObject(paramControllerPage)->removeDependent(this);

    //---do not forget to call parent ------
    return EditControllerEx1::terminate ();
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Controller::setComponentState (IBStream* state)
{
    // Here you get the state of the component (Processor part)
    if (!state)
        return kResultFalse;
    
    IBStreamer streamer(state, kLittleEndian);

    bool    savedBypass     = dftBypass;
    // int32   savedZoom       = dftZoom;
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
    // if (streamer.readInt32 (savedZoom)          == false) savedZoom       = dftZoom;
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

    setParamNormalized(kParamBypass,    savedBypass ? 1 : 0);
    // setParamNormalized(kParamZoom,      savedZoom);
    // setParamNormalized(kParamOS,        savedOS);

    setParamNormalized(kParamScLfIn,    savedScLfIn ? 1 : 0);
    setParamNormalized(kParamScLfType,  paramScLfType.ToNormalized(savedScLfType));
    setParamNormalized(kParamScLfFreq,  paramScLfFreq.ToNormalized(savedScLfFreq));
    setParamNormalized(kParamScLfGain,  paramScLfGain.ToNormalized(savedScLfGain));
    setParamNormalized(kParamScHfIn,    savedScHfIn ? 1 : 0);
    setParamNormalized(kParamScHfType,  paramScHfType.ToNormalized(savedScHfType));
    setParamNormalized(kParamScHfFreq,  paramScHfFreq.ToNormalized(savedScHfFreq));
    setParamNormalized(kParamScHfGain,  paramScHfGain.ToNormalized(savedScHfGain));
    setParamNormalized(kParamScListen,  savedScListen ? 1 : 0);
    
    setParamNormalized(kParamDetectorType,      paramDetectorType.ToNormalized(savedDType));
    setParamNormalized(kParamSidechainTopology, paramSidechainTopology.ToNormalized(savedScTopology));
    setParamNormalized(kParamHilbertEnable,     savedHilbertEnable ? 1 : 0);
    setParamNormalized(kParamLookaheadEnable,   savedLookaheadEnable ? 1 : 0);
    setParamNormalized(kParamAttack,            paramAttack.ToNormalized(savedAttack));
    setParamNormalized(kParamRelease,           paramRelease.ToNormalized(savedRelease));
    
    setParamNormalized(kParamThreshold,     paramThreshold.ToNormalized(savedThreshold));
    setParamNormalized(kParamRatio,         paramRatio.ToNormalized(savedRatio));
    setParamNormalized(kParamKnee,          paramKnee.ToNormalized(savedKnee));
    setParamNormalized(kParamMakeup,        paramMakeup.ToNormalized(DecibelConverter::ToDecibel(savedMakeup)));
    
    setParamNormalized(kParamMix,           savedMix);
    setParamNormalized(kParamInput,         paramInput.ToNormalized(DecibelConverter::ToDecibel(savedInputGain)));
    setParamNormalized(kParamOutput,        paramOutput.ToNormalized(DecibelConverter::ToDecibel(savedOutputGain)));
    setParamNormalized(kParamSoftBypass,    savedSoftBypass ? 1 : 0);
    
    bypass       = savedBypass;
    // pZoom;
    OS           = savedOS;
    
    scLfIn       = savedScLfIn;
    scLfType     = savedScLfType;
    scLfFreq     = savedScLfFreq;
    scLfGain     = savedScLfGain;
    scHfIn       = savedScHfIn;
    scHfType     = savedScHfType;
    scHfFreq     = savedScHfFreq;
    scHfGain     = savedScHfGain;
    scListen     = savedScListen;
    
    dType        = savedDType;
    scTopology   = savedScTopology;
    hilbertEnable    = savedHilbertEnable;
    lookaheadEnable  = savedLookaheadEnable;
    attack       = savedAttack;
    _release      = savedRelease;

    threshold    = savedThreshold;
    ratio        = savedRatio;
    knee         = savedKnee;
    makeup       = DecibelConverter::ToDecibel(savedMakeup);

    mix          = savedMix;
    inputGain    = DecibelConverter::ToDecibel(savedInputGain);
    outputGain   = DecibelConverter::ToDecibel(savedOutputGain);
    softBypass   = savedSoftBypass;
    
    return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Controller::setState (IBStream* state)
{
    // Here you get the state of the controller
    if (!state)
        return kResultFalse;
    
    IBStreamer streamer(state, kLittleEndian);

    // bool    savedBypass     = dftBypass;
    int32   savedZoom       = dftZoom;
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

    // if (streamer.readBool  (savedBypass)        == false) savedBypass     = dftBypass;
    if (streamer.readInt32 (savedZoom)          == false) savedZoom       = dftZoom;
    if (streamer.readInt32 (savedOS)            == false) savedOS         = overSample_1x;
    
    if (streamer.readBool  (savedScLfIn)        == false) savedScLfIn     = scLfIn;
    if (streamer.readInt32 (savedScLfType)      == false) savedScLfType   = scLfType;
    if (streamer.readDouble(savedScLfFreq)      == false) savedScLfFreq   = scLfFreq;
    if (streamer.readDouble(savedScLfGain)      == false) savedScLfGain   = scLfGain;
    if (streamer.readBool  (savedScHfIn)        == false) savedScHfIn     = scHfIn;
    if (streamer.readInt32 (savedScHfType)      == false) savedScHfType   = scHfType;
    if (streamer.readDouble(savedScHfFreq)      == false) savedScHfFreq   = scHfFreq;
    if (streamer.readDouble(savedScHfGain)      == false) savedScHfGain   = scHfGain;
    if (streamer.readBool  (savedScListen)      == false) savedScListen   = scListen;
    
    if (streamer.readInt32 (savedDType)             == false) savedDType           = dType;
    if (streamer.readInt32 (savedScTopology)        == false) savedScTopology      = scTopology;
    if (streamer.readBool  (savedHilbertEnable)     == false) savedHilbertEnable   = hilbertEnable;
    if (streamer.readBool  (savedLookaheadEnable)   == false) savedLookaheadEnable = lookaheadEnable;
    if (streamer.readDouble(savedAttack)            == false) savedAttack          = attack;
    if (streamer.readDouble(savedRelease)           == false) savedRelease         = _release;
    
    if (streamer.readDouble(savedThreshold)     == false) savedThreshold       = threshold;
    if (streamer.readDouble(savedRatio)         == false) savedRatio           = ratio;
    if (streamer.readDouble(savedKnee)          == false) savedKnee            = knee;
    if (streamer.readDouble(savedMakeup)        == false) savedMakeup          = DecibelConverter::ToGain(makeup);
    
    if (streamer.readDouble(savedMix)           == false) savedMix             = mix;
    if (streamer.readDouble(savedInputGain)     == false) savedInputGain       = DecibelConverter::ToGain(inputGain);
    if (streamer.readDouble(savedOutputGain)    == false) savedOutputGain      = DecibelConverter::ToGain(outputGain);
    if (streamer.readBool  (savedSoftBypass)    == false) savedSoftBypass      = softBypass;

    // setParamNormalized(kParamBypass,    savedBypass ? 1 : 0);
    setParamNormalized(kParamZoom,      paramZoom.ToNormalized(savedZoom));
    // setParamNormalized(kParamOS,        savedOS); // UNUSED

    setParamNormalized(kParamScLfIn,    savedScLfIn ? 1 : 0);
    setParamNormalized(kParamScLfType,  paramScLfType.ToNormalized(savedScLfType));
    setParamNormalized(kParamScLfFreq,  paramScLfFreq.ToNormalized(savedScLfFreq));
    setParamNormalized(kParamScLfGain,  paramScLfGain.ToNormalized(savedScLfGain));
    setParamNormalized(kParamScHfIn,    savedScHfIn ? 1 : 0);
    setParamNormalized(kParamScHfType,  paramScHfType.ToNormalized(savedScHfType));
    setParamNormalized(kParamScHfFreq,  paramScHfFreq.ToNormalized(savedScHfFreq));
    setParamNormalized(kParamScHfGain,  paramScHfGain.ToNormalized(savedScHfGain));
    setParamNormalized(kParamScListen,  savedScListen ? 1 : 0);
    
    setParamNormalized(kParamDetectorType,      paramDetectorType.ToNormalized(savedDType));
    setParamNormalized(kParamSidechainTopology, paramSidechainTopology.ToNormalized(savedScTopology));
    setParamNormalized(kParamHilbertEnable,     savedHilbertEnable ? 1 : 0);
    setParamNormalized(kParamLookaheadEnable,   savedLookaheadEnable ? 1 : 0);
    setParamNormalized(kParamAttack,            paramAttack.ToNormalized(savedAttack));
    setParamNormalized(kParamRelease,           paramRelease.ToNormalized(savedRelease));
    
    setParamNormalized(kParamThreshold,     paramThreshold.ToNormalized(savedThreshold));
    setParamNormalized(kParamRatio,         paramRatio.ToNormalized(savedRatio));
    setParamNormalized(kParamKnee,          paramKnee.ToNormalized(savedKnee));
    setParamNormalized(kParamMakeup,        paramMakeup.ToNormalized(DecibelConverter::ToDecibel(savedMakeup)));
    
    setParamNormalized(kParamMix,           savedMix);
    setParamNormalized(kParamInput,         paramInput.ToNormalized(DecibelConverter::ToDecibel(savedInputGain)));
    setParamNormalized(kParamOutput,        paramOutput.ToNormalized(DecibelConverter::ToDecibel(savedOutputGain)));
    setParamNormalized(kParamSoftBypass,    savedSoftBypass ? 1 : 0);
    
    // bypass       = savedBypass;
    zoom         = savedZoom;
    // OS           = savedOS;  // UNUSED
    
    scLfIn       = savedScLfIn;
    scLfType     = savedScLfType;
    scLfFreq     = savedScLfFreq;
    scLfGain     = savedScLfGain;
    scHfIn       = savedScHfIn;
    scHfType     = savedScHfType;
    scHfFreq     = savedScHfFreq;
    scHfGain     = savedScHfGain;
    scListen     = savedScListen;
    
    dType        = savedDType;
    scTopology   = savedScTopology;
    hilbertEnable    = savedHilbertEnable;
    lookaheadEnable  = savedLookaheadEnable;
    attack       = savedAttack;
    _release      = savedRelease;

    threshold    = savedThreshold;
    ratio        = savedRatio;
    knee         = savedKnee;
    makeup       = DecibelConverter::ToDecibel(savedMakeup);

    mix          = savedMix;
    inputGain    = DecibelConverter::ToDecibel(savedInputGain);
    outputGain   = DecibelConverter::ToDecibel(savedOutputGain);
    softBypass   = savedSoftBypass;
    
    return kResultTrue;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Controller::getState (IBStream* state)
{
    // Here you are asked to deliver the state of the controller (if needed)
    // Note: the real state of your plug-in is saved in the processor

    if (!state)
        return kResultFalse;
    
    IBStreamer streamer(state, kLittleEndian);

    // bypass       = getParamNormalized(kParamBypass) > 0 ? true : false;
    zoom         = paramZoom.ToPlainList(getParamNormalized(kParamZoom));
    // OS           = getParamNormalized(kParamOS); // UNUSED, don't update
    
    scLfIn       = getParamNormalized(kParamScLfIn) > 0 ? true : false;
    scLfType     = paramScLfType.ToPlainList(getParamNormalized(kParamScLfType));
    scLfFreq     = paramScLfFreq.ToPlain(getParamNormalized(kParamScLfFreq));
    scLfGain     = paramScLfGain.ToPlain(getParamNormalized(kParamScLfGain));
    scHfIn       = getParamNormalized(kParamScHfIn) > 0 ? true : false;
    scHfType     = paramScHfType.ToPlainList(getParamNormalized(kParamScHfType));
    scHfFreq     = paramScHfFreq.ToPlain(getParamNormalized(kParamScHfFreq));
    scHfGain     = paramScHfGain.ToPlain(getParamNormalized(kParamScHfGain));
    scListen     = getParamNormalized(kParamScListen) > 0 ? true : false;
    
    dType        = paramDetectorType.ToPlainList(getParamNormalized(kParamDetectorType));
    scTopology   = paramSidechainTopology.ToPlainList(getParamNormalized(kParamSidechainTopology));
    hilbertEnable    = getParamNormalized(kParamHilbertEnable) > 0 ? true : false;
    lookaheadEnable  = getParamNormalized(kParamLookaheadEnable) > 0 ? true : false;
    attack       = paramAttack.ToPlain(getParamNormalized(kParamAttack));
    _release      = paramRelease.ToPlain(getParamNormalized(kParamRelease));

    threshold    = paramThreshold.ToPlain(getParamNormalized(kParamThreshold));
    ratio        = paramRatio.ToPlain(getParamNormalized(kParamRatio));
    knee         = paramKnee.ToPlain(getParamNormalized(kParamKnee));
    makeup       = DecibelConverter::ToGain(paramMakeup.ToPlain(getParamNormalized(kParamMakeup)));

    mix          = getParamNormalized(kParamMix);
    inputGain    = DecibelConverter::ToGain(paramInput.ToPlain(getParamNormalized(kParamInput)));
    outputGain   = DecibelConverter::ToGain(paramOutput.ToPlain(getParamNormalized(kParamOutput)));
    softBypass   = getParamNormalized(kParamSoftBypass) > 0 ? true : false;
    
    // streamer.writeBool(bypass);
    streamer.writeInt32(zoom);
    streamer.writeInt32(OS);            // UNUSED
    
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
    streamer.writeDouble(_release);
    
    streamer.writeDouble(threshold);
    streamer.writeDouble(ratio);
    streamer.writeDouble(knee);
    streamer.writeDouble(makeup);     // saved in lin gain, not dB
    
    streamer.writeDouble(mix);
    streamer.writeDouble(inputGain);  // saved in lin gain, not dB
    streamer.writeDouble(outputGain); // saved in lin gain, not dB
    streamer.writeBool(softBypass);
    
    return kResultTrue;
}

//------------------------------------------------------------------------
VSTGUI::IController* RLFCMP_Controller::createSubController (VSTGUI::UTF8StringPtr name,
                                                             const VSTGUI::IUIDescription* description,
                                                             VSTGUI::VST3Editor* editor)
{
    if (VSTGUI::UTF8StringView(name) == "VuMeterController")
    {
        auto* controller = new VuMeterController(editor, this);
        addUIVuMeterController(controller);
        return controller;
    }
    if (VSTGUI::UTF8StringView(name) == "eqCurveController")
    {
        Steinberg::Vst::Parameter* ScLfInParam = getParameterObject(kParamScLfIn);
        Steinberg::Vst::Parameter* ScLfTypeParam = getParameterObject(kParamScLfType);
        Steinberg::Vst::Parameter* ScLfFreqParam = getParameterObject(kParamScLfFreq);
        Steinberg::Vst::Parameter* ScLfGainParam = getParameterObject(kParamScLfGain);
        Steinberg::Vst::Parameter* ScHfInParam = getParameterObject(kParamScHfIn);
        Steinberg::Vst::Parameter* ScHfTypeParam = getParameterObject(kParamScHfType);
        Steinberg::Vst::Parameter* ScHfFreqParam = getParameterObject(kParamScHfFreq);
        Steinberg::Vst::Parameter* ScHfGainParam = getParameterObject(kParamScHfGain);
        auto* controller = new EQCurveViewController(editor,
                                                     this,
                                                     ScLfInParam,
                                                     ScLfTypeParam,
                                                     ScLfFreqParam,
                                                     ScLfGainParam,
                                                     ScHfInParam,
                                                     ScHfTypeParam,
                                                     ScHfFreqParam,
                                                     ScHfGainParam);
        addUIEQCurveViewController(controller);
        return controller;
    }
    if (VSTGUI::UTF8StringView(name) == "transferCurveController")
    {
        Steinberg::Vst::Parameter* ThresholdParam = getParameterObject(kParamThreshold);
        Steinberg::Vst::Parameter* KneeParam = getParameterObject(kParamKnee);
        Steinberg::Vst::Parameter* RatioParam = getParameterObject(kParamRatio);
        Steinberg::Vst::Parameter* MakeupParam = getParameterObject(kParamMakeup);
        Steinberg::Vst::Parameter* MixParam = getParameterObject(kParamMix);
        auto* controller = new TransferCurveViewController(editor, this, ThresholdParam, KneeParam, RatioParam, MakeupParam, MixParam);
        addUITransferCurveViewController(controller);
        return controller;
    }
    return nullptr;
}

//------------------------------------------------------------------------
IPlugView* PLUGIN_API RLFCMP_Controller::createView (FIDString name)
{
    // Here the Host wants to open your editor (if you have one)
    if (FIDStringsEqual (name, Vst::ViewType::kEditor))
    {
        // create your editor here and return a IPlugView ptr of it
        auto* view = new VSTGUI::VST3Editor (this, "view", "RLFCMP_editor.uidesc");
        
        std::vector<double> _zoomFactors;
        _zoomFactors.push_back(0.50);
        _zoomFactors.push_back(0.75);
        _zoomFactors.push_back(1.00);
        _zoomFactors.push_back(1.25);
        _zoomFactors.push_back(1.50);
        _zoomFactors.push_back(1.75);
        _zoomFactors.push_back(2.00);
        view->setAllowedZoomFactors(_zoomFactors);
        view->setZoomFactor(1.0);
        view->setIdleRate(1000.0/60.0);
        
        setKnobMode(Steinberg::Vst::KnobModes::kLinearMode);

        return view;
    }
    return nullptr;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Controller::setParamNormalized (Vst::ParamID tag, Vst::ParamValue value)
{
    // called by host to update your parameters
    tresult result = EditControllerEx1::setParamNormalized (tag, value);
    return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Controller::getParamStringByValue (Vst::ParamID tag, Vst::ParamValue valueNormalized, Vst::String128 string)
{
    // called by host to get a string for given normalized value of a specific parameter
    // (without having to set the value!)
    return EditControllerEx1::getParamStringByValue (tag, valueNormalized, string);
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Controller::getParamValueByString (Vst::ParamID tag, Vst::TChar* string, Vst::ParamValue& valueNormalized)
{
    // called by host to get a normalized value from a string representation of a specific parameter
    // (without having to set the value!)
    return EditControllerEx1::getParamValueByString (tag, string, valueNormalized);
}

//------------------------------------------------------------------------
void RLFCMP_Controller::editorAttached(Vst::EditorView* editor)
{
    editors.push_back(editor);
    if (!editors.empty())
    {
        if (IPtr<Vst::IMessage> message = owned (allocateMessage ()))
        {
           message->setMessageID ("GUI");
           message->getAttributes()->setInt ("start", 1);
           sendMessage (message);
        }
    }
}

//------------------------------------------------------------------------
void RLFCMP_Controller::editorRemoved(Vst::EditorView* editor)
{
    editors.erase(std::find(editors.begin(), editors.end(), editor));
    if (editors.empty())
    {
        if (IPtr<Vst::IMessage> message = owned (allocateMessage ()))
        {
           message->setMessageID ("GUI");
           message->getAttributes()->setInt ("start", 0);
           sendMessage (message);
        }
    }
}

//------------------------------------------------------------------------
void PLUGIN_API RLFCMP_Controller::update(FUnknown* changedUnknown, int32 message)
{
    EditControllerEx1::update(changedUnknown, message);

    // GUI Resizing
    // check 'zoomtest' code at
    // https://github.com/steinbergmedia/vstgui/tree/vstgui4_10/vstgui/tests/uidescription%20vst3/source

    Vst::Parameter* param = FCast<Vst::Parameter>(changedUnknown);
    if (!param)
        return;

    if (param->getInfo().id == kParamZoom)
    {
        size_t index = static_cast<size_t> (param->toPlain(param->getNormalized()));

        if (index >= zoomFactors.size())
            return;

        for (EditorVector::const_iterator it = editors.begin(), end = editors.end(); it != end; ++it)
        {
            VSTGUI::VST3Editor* editor = dynamic_cast<VSTGUI::VST3Editor*>(*it);
            if (editor)
                editor->setZoomFactor(zoomFactors[index].factor);
        }
    }
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Controller::notify(Vst::IMessage* message)
{
    if (!message)
        return kInvalidArgument;
    
    if (strcmp (message->getMessageID (), "GUI") == 0)
    {
        ParamValue getValue = 0.0;

        if (message->getAttributes ()->getFloat (msgInputPeakL, getValue) == kResultTrue) vuInLPeak  = getValue;
        if (message->getAttributes ()->getFloat (msgInputPeakR, getValue) == kResultTrue) vuInRPeak  = getValue;
        if (message->getAttributes ()->getFloat (msgInputRMSL, getValue) == kResultTrue) vuInLRMS  = getValue;
        if (message->getAttributes ()->getFloat (msgInputRMSR, getValue) == kResultTrue) vuInRRMS  = getValue;
        if (message->getAttributes ()->getFloat (msgOutputPeakL, getValue) == kResultTrue) vuOutLPeak  = getValue;
        if (message->getAttributes ()->getFloat (msgOutputPeakR, getValue) == kResultTrue) vuOutRPeak  = getValue;
        if (message->getAttributes ()->getFloat (msgOutputRMSL, getValue) == kResultTrue) vuOutLRMS  = getValue;
        if (message->getAttributes ()->getFloat (msgOutputRMSR, getValue) == kResultTrue) vuOutRRMS  = getValue;
        if (message->getAttributes ()->getFloat (msgDetectorLevel, getValue) == kResultTrue)
        {
            if (!transferCurveViewControllerControllers.empty())
            {
                for (auto &iter : transferCurveViewControllerControllers)
                {
                    iter->setVuInMono(getValue);
                }
            }
        }
        if (message->getAttributes ()->getFloat (msgGainReduction, getValue) == kResultTrue)
        {
            vuGainReduction  = getValue;
            
            if (!vuMeterControllers.empty())
            {
                for (auto &iter : vuMeterControllers)
                {
                    iter->setVuMeterByTag(vuInLPeak, kInLPeak);
                    iter->setVuMeterByTag(vuInRPeak, kInRPeak);
                    iter->setVuMeterByTag(vuInLRMS, kInLRMS);
                    iter->setVuMeterByTag(vuInRRMS, kInRRMS);
                    iter->setVuMeterByTag(vuOutLPeak, kOutLPeak);
                    iter->setVuMeterByTag(vuOutRPeak, kOutRPeak);
                    iter->setVuMeterByTag(vuOutLRMS, kOutLRMS);
                    iter->setVuMeterByTag(vuOutRRMS, kOutRRMS);
                    iter->setVuMeterByTag(vuGainReduction, kGainReduction);
                    iter->updateVuMeterValue();
                }
            }
        }
        return kResultOk;
    }
    return EditControllerEx1::notify(message);
}
//------------------------------------------------------------------------
} // namespace yg331
