//------------------------------------------------------------------------
// Copyright(c) 2024 yg331.
//------------------------------------------------------------------------

#pragma once

#include "RLFCMP_shared.h"
#include "public.sdk/source/vst/vsteditcontroller.h"
#include "vstgui/plugin-bindings/vst3editor.h"
#include "vstgui/uidescription/uiviewswitchcontainer.h"
#include "vstgui/uidescription/delegationcontroller.h"

namespace VSTGUI {
//------------------------------------------------------------------------
// EQ Curve Display
//------------------------------------------------------------------------
class EQCurveView : public CControl
{
public:
    EQCurveView(const CRect& size, IControlListener* listener, int32_t tag, CBitmap* background);
    EQCurveView(const EQCurveView& v);

    void setFFTArray(float* array, int sampleBlockSize, double sampleRate);

    // get/set Attributes
    void setBackColor(CColor color) { if (BackColor != color) { BackColor = color; setDirty(true); } }
    CColor getBackColor() const { return BackColor; }

    void setBorderColor(CColor color) { if (BorderColor != color) { BorderColor = color; setDirty(true); } }
    CColor getBorderColor() const { return BorderColor; }

    void setLineColor(CColor color) { if (LineColor != color) { LineColor = color; setDirty(true); } }
    CColor getLineColor() const { return LineColor; }

    void setFFTLineColor(CColor color) { if (FFTLineColor != color) { FFTLineColor = color; setDirty(true); } }
    CColor getFFTLineColor() const { return FFTLineColor; }

    void setFFTFillColor(CColor color) { if (FFTFillColor != color) { FFTFillColor = color; setDirty(true); } }
    CColor getFFTFillColor() const { return FFTFillColor; }
    
    // get/set Parameters
    void   setIn(bool value, int type) { if (type == yg331::PassShelfFilter::rLow) {LF_SVF.setIn(value);}   else {HF_SVF.setIn(value);} setDirty(true); }
    bool   getIn(int type) const       { if (type == yg331::PassShelfFilter::rLow) {return LF_SVF.getIn();} else {return HF_SVF.getIn();} }
    
    void   setType(int value, int type) { if (type == yg331::PassShelfFilter::rLow) {LF_SVF.setType(value);} else {HF_SVF.setType(value);} setDirty(true); }
    int    getType(int type) const { if (type == yg331::PassShelfFilter::rLow) {return LF_SVF.getType();} else {return HF_SVF.getType();} }
    
    void   setFreq(double value, int type) { if (type == yg331::PassShelfFilter::rLow) {LF_SVF.setFreq(value);} else {HF_SVF.setFreq(value);} setDirty(true); }
    int    getFreq(int type) const { if (type == yg331::PassShelfFilter::rLow) {return LF_SVF.getFreq();} else {return HF_SVF.getFreq();} }
    
    void   setGain(double value, int type) { if (type == yg331::PassShelfFilter::rLow) {LF_SVF.setGain(value);} else {HF_SVF.setGain(value);} setDirty(true); }
    int    getGain(int type) const { if (type == yg331::PassShelfFilter::rLow) {return LF_SVF.getGain();} else {return HF_SVF.getGain();} }
    
    void   makeSVF() { LF_SVF.makeSVF(); HF_SVF.makeSVF(); setDirty(true); }

    // overrides
    void setDirty(bool state) override { CView::setDirty(state); };
    void draw(CDrawContext* pContext) override;
    void setViewSize(const CRect& newSize, bool invalid = true) override { CControl::setViewSize(newSize, invalid); };
    bool sizeToFit() override;

    /** called on idle when view wants idle */
    void onIdle() override { invalid(); };

    CLASS_METHODS(EQCurveView, CControl)

    //------------------------------------------------------------------------

protected:
    ~EQCurveView() noexcept override
    {
    };

    CColor BackColor;
    CColor LineColor;
    CColor BorderColor;
    CColor FFTLineColor;
    CColor FFTFillColor;

    yg331::PassShelfFilter LF_SVF;
    yg331::PassShelfFilter HF_SVF;

    static constexpr int fftOrder = yg331::_fftOrder;
    static constexpr int fftSize = 1 << fftOrder;      // 4096 samples
    static constexpr int numBins = fftSize / 2 + 1;    // 2049 bins

    float fft_linear[numBins] = { 0.0 };
    float fft_RMS[numBins]    = { 0.0 };
    float fft_freq[numBins]   = { 0.0 };
};

//------------------------------------------------------------------------
//  Transfer Curve View
//------------------------------------------------------------------------
class TransferCurveView : public CControl {
public:
    TransferCurveView(
        const VSTGUI::CRect& size,
        Steinberg::Vst::EditController* _editController = nullptr)
        : CControl(size, nullptr, 0)
    {
        BackColor = kWhiteCColor;
        LineColor = kBlackCColor;
        setWantsIdle(true);
    }
    TransferCurveView(const TransferCurveView& v)
        : CControl(v)
        , BackColor(v.BackColor)
        , LineColor(v.LineColor)
    {
        setWantsIdle(true);
    }

    // get/set Parameters
    void   setThreshold(double value) { if (threshold != value) { threshold = value; setDirty(true); } }
    double getThreshold() const { return threshold; }

    void   setKnee(double value) { if (knee != value) { knee = value; kneeHalf  = knee / 2.0; setDirty(true); } }
    double getKnee() const { return knee; }

    void   setRatio(double value) { if (ratio != value) { ratio = value; slope = 1.0 / ratio - 1.0; setDirty(true); } }
    double getRatio() const { return ratio; }
    
    void   setMakeup(double value) { if (makeup != value) { makeup = value; setDirty(true); } }
    double getMakeup() const { return makeup; }

    void   setMix(double value) { if (mix != value) { mix = value; setDirty(true); } }
    double getMix() const { return mix; }

    // get/set Attributes
    void   setBackColor(CColor color) { if (BackColor != color) { BackColor = color; setDirty(true); } }
    CColor getBackColor() const { return BackColor; }

    void   setLineColor(CColor color) { if (LineColor != color) { LineColor = color; setDirty(true); } }
    CColor getLineColor() const { return LineColor; }

    // overrides
    void setDirty(bool state) override { CView::setDirty(state); };
    void draw(CDrawContext* _pContext) override;
    void setViewSize(const CRect& newSize, bool invalid = true) override { CControl::setViewSize(newSize, invalid); };
    bool sizeToFit() override ;
    void onIdle() override { invalid(); };

    CLASS_METHODS(TransferCurveView, CControl)

protected:
    ~TransferCurveView() noexcept override
    {
    };

    CColor        BackColor;
    CColor        LineColor;

    double threshold = yg331::dftThreshold;
    double knee      = yg331::dftKnee;
    double ratio     = yg331::dftRatio;
    double makeup    = yg331::dftMakeup;
    double mix       = yg331::dftMix/100.0;
    
    double slope     = 1.0 / ratio - 1.0;
    double kneeHalf  = knee / 2.0;

    static constexpr double MAX_dB = 0.0;
    static constexpr double MIN_dB = -60.0;
    static constexpr double DB_Range = MAX_dB - MIN_dB;
    static constexpr double Inv_DB_R = 1.0 / DB_Range;
};

//------------------------------------------------------------------------
//  Click to Reset Param Display View
//------------------------------------------------------------------------
class ClickResetParamDisplay : public CParamDisplay {
public:
    enum updateStyle
    {
        kUpdateMin = 1 << 0,
        kUpdateMax = 1 << 1,
    };
    ClickResetParamDisplay
        (const CRect& size, CBitmap* background = nullptr, int32_t style = 0)
        : CParamDisplay(size, background, style)
    {
        originalBack = getBackColor();
        setWantsIdle(true);
    };
    ClickResetParamDisplay
        (const CParamDisplay& paramDisplay)
        : CParamDisplay(paramDisplay)
    {
        originalBack = getBackColor();
        setWantsIdle(true);
    };
    
    // get/set Attributes
    void    setStyle_(int32_t newStyle) { _style = newStyle; }
    int32_t getStyle_() const { return _style; }
    
    // overrides
    void setValue(float val) override
    {
        directValue = val;
        if (_style == kUpdateMax)
            CParamDisplay::setValue(std::max(getValue(), val));
        else
            CParamDisplay::setValue(std::min(getValue(), val));
        
        if (_style == kUpdateMax)
            if (getValue()>0.0)
                if (!over) {
                    originalBack = getBackColor();
                    setBackColor(VSTGUI::kRedCColor);
                    over = true;
                }
    };
    
    void onMouseDownEvent(MouseDownEvent& event) override {
        if (over) {
            setBackColor(originalBack);
            over = false;
        }
        CParamDisplay::setValue(directValue);
        CParamDisplay::onMouseDownEvent(event);
    };
    
    void onIdle() override { invalid(); };
 
protected:
    int32_t _style;
    float   directValue = 0.0;
    bool    over = false;
    CColor  originalBack;
};

//------------------------------------------------------------------------
//  Metering reset container
//------------------------------------------------------------------------
class MeterViewContainer : public CViewContainer
{
public:
    MeterViewContainer(const CRect& size) : CViewContainer(size) {};
    void onMouseDownEvent(MouseDownEvent& event) override {
        for (auto& child : getChildren())
        {
            child->onMouseDownEvent(event);
        }
        CViewContainer::onMouseDownEvent(event);
    };
};

//------------------------------------------------------------------------
//  VU meter view
//------------------------------------------------------------------------
class MyVuMeter : public CControl {
private:
    enum StyleEnum
    {
        StyleHorizontal = 0,
        StyleVertical,
    };
public:
    enum Style
    {
        kHorizontal = 1 << StyleHorizontal,
        kVertical = 1 << StyleVertical,
    };

    MyVuMeter(const CRect& size, int32_t style);
    MyVuMeter(const MyVuMeter& vuMeter);

    // get/set Attributes
    void    setStyle(int32_t newStyle) { style = newStyle; invalid(); }
    int32_t getStyle() const { return style; }

    CColor getVuOnColor() const { return vuOnColor; }
    void   setVuOnColor(CColor color) { if (vuOnColor != color) { vuOnColor = color; setDirty(true); } }
    
    CColor getVuOffColor() const { return vuOffColor; }
    void   setVuOffColor(CColor color) { if (vuOffColor != color) { vuOffColor = color; setDirty(true); } }

    // overrides
    void setDirty(bool state) override { CView::setDirty(state); };
    void draw(CDrawContext* _pContext) override;
    void setViewSize(const CRect& newSize, bool invalid) override;
    bool sizeToFit() override;
    
    /** called on idle when view wants idle */
    void onIdle() override { invalid(); };
    
    CLASS_METHODS(MyVuMeter, CControl)

protected:
    ~MyVuMeter() noexcept override
    {
    };

    int32_t  style;

    CColor   vuOnColor;
    CColor   vuOffColor;

    CRect    rectOn;
    CRect    rectOff;
};
}

namespace yg331 {

class EQCurveViewController;
class TransferCurveViewController;
class VuMeterController;

//------------------------------------------------------------------------
//  RLFCMP_Controller
//------------------------------------------------------------------------
class RLFCMP_Controller : public Steinberg::Vst::EditControllerEx1, public VSTGUI::VST3EditorDelegate
{
public:
//------------------------------------------------------------------------
    RLFCMP_Controller () = default;
    ~RLFCMP_Controller () SMTG_OVERRIDE = default;

    // Create function
    static Steinberg::FUnknown* createInstance (void* /*context*/)
    {
        return (Steinberg::Vst::IEditController*)new RLFCMP_Controller;
    }

    // IPluginBase
    Steinberg::tresult PLUGIN_API initialize (Steinberg::FUnknown* context) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API terminate () SMTG_OVERRIDE;

    // EditController
    Steinberg::tresult PLUGIN_API setComponentState (Steinberg::IBStream* state) SMTG_OVERRIDE;
    Steinberg::IPlugView* PLUGIN_API createView (Steinberg::FIDString name) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API setState (Steinberg::IBStream* state) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API getState (Steinberg::IBStream* state) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API setParamNormalized (Steinberg::Vst::ParamID tag,
                                                      Steinberg::Vst::ParamValue value) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API getParamStringByValue (Steinberg::Vst::ParamID tag,
                                                         Steinberg::Vst::ParamValue valueNormalized,
                                                         Steinberg::Vst::String128 string) SMTG_OVERRIDE;
    Steinberg::tresult PLUGIN_API getParamValueByString (Steinberg::Vst::ParamID tag,
                                                         Steinberg::Vst::TChar* string,
                                                         Steinberg::Vst::ParamValue& valueNormalized) SMTG_OVERRIDE;
    
    //---from VST3EditorDelegate-----------
    VSTGUI::IController* createSubController (VSTGUI::UTF8StringPtr name,
                                              const VSTGUI::IUIDescription* description,
                                              VSTGUI::VST3Editor* editor) SMTG_OVERRIDE;
    //---from ComponentBase-----
    // EditController
    Steinberg::tresult PLUGIN_API notify(Steinberg::Vst::IMessage* message) SMTG_OVERRIDE;
    // Steinberg::tresult PLUGIN_API receiveText(const char* text) SMTG_OVERRIDE;
    void PLUGIN_API update(Steinberg::FUnknown* changedUnknown, Steinberg::int32 message) SMTG_OVERRIDE;
    void editorAttached(Steinberg::Vst::EditorView* editor) SMTG_OVERRIDE; ///< called from EditorView if it was attached to a parent
    void editorRemoved (Steinberg::Vst::EditorView* editor) SMTG_OVERRIDE; ///< called from EditorView if it was removed from a parent

    //------------------------------------------------------------------------
    Steinberg::Vst::Parameter* getParameterObject(Steinberg::Vst::ParamID tag) SMTG_OVERRIDE
    {
        Steinberg::Vst::Parameter* param = EditControllerEx1::getParameterObject(tag);
        if (param == 0)
        {
            param = uiParameters.getParameter(tag);
        }
        return param;
    }
    bool isPrivateParameter(const Steinberg::Vst::ParamID paramID) SMTG_OVERRIDE
    {
        return uiParameters.getParameter(paramID) != 0 ? true : false;
    }

    // make sure that our UI only parameters doesn't call the following three EditController methods: beginEdit, endEdit, performEdit
    //------------------------------------------------------------------------
    Steinberg::tresult beginEdit(Steinberg::Vst::ParamID tag) SMTG_OVERRIDE
    {
        if (EditControllerEx1::getParameterObject(tag))
            return EditControllerEx1::beginEdit(tag);
        return Steinberg::kResultFalse;
    }

    //------------------------------------------------------------------------
    Steinberg::tresult performEdit(Steinberg::Vst::ParamID tag, Steinberg::Vst::ParamValue valueNormalized) SMTG_OVERRIDE
    {
        if (EditControllerEx1::getParameterObject(tag))
            return EditControllerEx1::performEdit(tag, valueNormalized);
        return Steinberg::kResultFalse;
    }

    //------------------------------------------------------------------------
    Steinberg::tresult endEdit(Steinberg::Vst::ParamID tag) SMTG_OVERRIDE
    {
        if (EditControllerEx1::getParameterObject(tag))
            return EditControllerEx1::endEdit(tag);
        return Steinberg::kResultFalse;
    }
    
    //---Internal functions-------
    void addUIVuMeterController(VuMeterController* controller)
    {
        vuMeterControllers.push_back(controller);
    }
    void removeUIVuMeterController(VuMeterController* controller)
    {
        auto it = std::find(vuMeterControllers.begin(), vuMeterControllers.end(), controller);
        if (it != vuMeterControllers.end())
            vuMeterControllers.erase(it);
    }
    
    void addUIEQCurveViewController(EQCurveViewController* controller)
    {
        eqCurveViewControllers.push_back(controller);
    }
    void removeUIEQCurveViewController(EQCurveViewController* controller)
    {
        auto it = std::find(eqCurveViewControllers.begin(), eqCurveViewControllers.end(), controller);
        if (it != eqCurveViewControllers.end())
            eqCurveViewControllers.erase(it);
    }
    
    void addUITransferCurveViewController(TransferCurveViewController* controller)
    {
        transferCurveViewControllerControllers.push_back(controller);
    }
    void removeUITransferCurveViewController(TransferCurveViewController* controller)
    {
        auto it = std::find(transferCurveViewControllerControllers.begin(), transferCurveViewControllerControllers.end(), controller);
        if (it != transferCurveViewControllerControllers.end())
            transferCurveViewControllerControllers.erase(it);
    }
    
    Steinberg::Vst::ParamValue getVuMeterByTag(Steinberg::Vst::ParamID tag);

    //---Interface---------
    DEFINE_INTERFACES
        // Here you can add more supported VST3 interfaces
        // DEF_INTERFACE (Vst::IXXX)
    END_DEFINE_INTERFACES (EditController)
    DELEGATE_REFCOUNT (EditController)

//------------------------------------------------------------------------
protected:
    // UI only parameter list
    Steinberg::Vst::ParameterContainer uiParameters;

    // editor list
    typedef std::vector<Steinberg::Vst::EditorView*> EditorVector;
    EditorVector editors;
    
    // zoom title-value struct
    struct ZoomFactor
    {
        const Steinberg::tchar* title;
        double factor;

        ZoomFactor(const Steinberg::tchar* title, double factor) : title(title), factor(factor) {}
    };
    typedef std::vector<ZoomFactor> ZoomFactorVector;
    ZoomFactorVector zoomFactors;
    
    // sub-controller list
    using UIVuMeterControllerList = std::vector<VuMeterController*>;
    using UIEQCurveViewControllerList = std::vector<EQCurveViewController*>;
    using UITransferCurveViewControllerList = std::vector<TransferCurveViewController*>;
    
    UIVuMeterControllerList vuMeterControllers;
    UIEQCurveViewControllerList eqCurveViewControllers;
    UITransferCurveViewControllerList transferCurveViewControllerControllers;
    
    ParamValue pZoom;
    Steinberg::Vst::ParamValue vuInLRMS = 0.0, vuInRRMS = 0.0;
    Steinberg::Vst::ParamValue vuInLPeak = 0.0, vuInRPeak = 0.0;
    Steinberg::Vst::ParamValue vuOutLRMS = 0.0, vuOutRRMS = 0.0;
    Steinberg::Vst::ParamValue vuOutLPeak = 0.0, vuOutRPeak = 0.0;
    Steinberg::Vst::ParamValue vuGainReduction = 0.0;
    Steinberg::Vst::ParamValue vuInMono = 0.0, vuOutMono = 0.0;
};

//------------------------------------------------------------------------
// EQCurveViewController
//------------------------------------------------------------------------
class EQCurveViewController
    : public Steinberg::FObject
    , public VSTGUI::DelegationController
    , public VSTGUI::CBaseObject
{
public:
    EQCurveViewController (IController* baseController,
                           RLFCMP_Controller* _mainController,
                           Steinberg::Vst::Parameter* ParamScLfIn,
                           Steinberg::Vst::Parameter* ParamScLfType,
                           Steinberg::Vst::Parameter* ParamScLfFreq,
                           Steinberg::Vst::Parameter* ParamScLfGain,
                           Steinberg::Vst::Parameter* ParamScHfIn,
                           Steinberg::Vst::Parameter* ParamScHfType,
                           Steinberg::Vst::Parameter* ParamScHfFreq,
                           Steinberg::Vst::Parameter* ParamScHfGain)
        : DelegationController(baseController)
        , mainController(_mainController)
        , ParamScLfIn  (ParamScLfIn)
        , ParamScLfType(ParamScLfType)
        , ParamScLfFreq(ParamScLfFreq)
        , ParamScLfGain(ParamScLfGain)
        , ParamScHfIn  (ParamScHfIn)
        , ParamScHfType(ParamScHfType)
        , ParamScHfFreq(ParamScHfFreq)
        , ParamScHfGain(ParamScHfGain)
        , eqCurveView(nullptr)
    {
        if (ParamScLfIn  ) ParamScLfIn  ->addDependent(this);
        if (ParamScLfType) ParamScLfType->addDependent(this);
        if (ParamScLfFreq) ParamScLfFreq->addDependent(this);
        if (ParamScLfGain) ParamScLfGain->addDependent(this);
        if (ParamScHfIn  ) ParamScHfIn  ->addDependent(this);
        if (ParamScHfType) ParamScHfType->addDependent(this);
        if (ParamScHfFreq) ParamScHfFreq->addDependent(this);
        if (ParamScHfGain) ParamScHfGain->addDependent(this);
    }

    ~EQCurveViewController() override
    {
        if (ParamScLfIn  ) ParamScLfIn  ->removeDependent(this);
        if (ParamScLfType) ParamScLfType->removeDependent(this);
        if (ParamScLfFreq) ParamScLfFreq->removeDependent(this);
        if (ParamScLfGain) ParamScLfGain->removeDependent(this);
        if (ParamScHfIn  ) ParamScHfIn  ->removeDependent(this);
        if (ParamScHfType) ParamScHfType->removeDependent(this);
        if (ParamScHfFreq) ParamScHfFreq->removeDependent(this);
        if (ParamScHfGain) ParamScHfGain->removeDependent(this);
        
        if (eqCurveView)
        {
            eqCurveView->unregisterControlListener (this);
            eqCurveView->forget ();
        }

        mainController->removeUIEQCurveViewController(this);
    }
    
    void setFFTArray(float* array, int sampleBlockSize, double sampleRate)
    { eqCurveView->setFFTArray(array, sampleBlockSize, sampleRate); }
    
private:
    using CControl       = VSTGUI::CControl;
    using CView          = VSTGUI::CView;
    using UTF8String     = VSTGUI::UTF8String;
    using UIAttributes   = VSTGUI::UIAttributes;
    using IUIDescription = VSTGUI::IUIDescription;
    using EQCurveView    = VSTGUI::EQCurveView;
    
    // FObject
    void PLUGIN_API update( Steinberg::FUnknown* changedUnknown, Steinberg::int32 message) SMTG_OVERRIDE
    {
        if (eqCurveView)
        {
            if (auto* p = Steinberg::FCast<Steinberg::Vst::Parameter>(changedUnknown); p)
            {
                if (message == kChanged)
                {
                    
                    if (p == ParamScLfIn   && ParamScLfIn  ) eqCurveView->setIn  (p->getNormalized(), yg331::PassShelfFilter::rLow);
                    if (p == ParamScLfType && ParamScLfType) eqCurveView->setType(paramScLfType.ToPlain(p->getNormalized()), yg331::PassShelfFilter::rLow);
                    if (p == ParamScLfFreq && ParamScLfFreq) eqCurveView->setFreq(paramScLfFreq.ToPlain(p->getNormalized()), yg331::PassShelfFilter::rLow);
                    if (p == ParamScLfGain && ParamScLfGain) eqCurveView->setGain(paramScLfGain.ToPlain(p->getNormalized()), yg331::PassShelfFilter::rLow);
                    if (p == ParamScHfIn   && ParamScHfIn  ) eqCurveView->setIn  (p->getNormalized(), yg331::PassShelfFilter::rHigh);
                    if (p == ParamScHfType && ParamScHfType) eqCurveView->setType(paramScHfType.ToPlain(p->getNormalized()), yg331::PassShelfFilter::rHigh);
                    if (p == ParamScHfFreq && ParamScHfFreq) eqCurveView->setFreq(paramScHfFreq.ToPlain(p->getNormalized()), yg331::PassShelfFilter::rHigh);
                    if (p == ParamScHfGain && ParamScHfGain) eqCurveView->setGain(paramScHfGain.ToPlain(p->getNormalized()), yg331::PassShelfFilter::rHigh);
                    // curveView->invalid();
                    eqCurveView->makeSVF();
                }
                else if (message == kWillDestroy)
                {
                    if (p == ParamScLfIn   && ParamScLfIn  ) { ParamScLfIn  ->removeDependent(this); ParamScLfIn   = nullptr; }
                    if (p == ParamScLfType && ParamScLfType) { ParamScLfType->removeDependent(this); ParamScLfType = nullptr; }
                    if (p == ParamScLfFreq && ParamScLfFreq) { ParamScLfFreq->removeDependent(this); ParamScLfFreq = nullptr; }
                    if (p == ParamScLfGain && ParamScLfGain) { ParamScLfGain->removeDependent(this); ParamScLfGain = nullptr; }
                    if (p == ParamScHfIn   && ParamScHfIn  ) { ParamScHfIn  ->removeDependent(this); ParamScHfIn   = nullptr; }
                    if (p == ParamScHfType && ParamScHfType) { ParamScHfType->removeDependent(this); ParamScHfType = nullptr; }
                    if (p == ParamScHfFreq && ParamScHfFreq) { ParamScHfFreq->removeDependent(this); ParamScHfFreq = nullptr; }
                    if (p == ParamScHfGain && ParamScHfGain) { ParamScHfGain->removeDependent(this); ParamScHfGain = nullptr; }
                }
            }
        }
    }

    //--- is called when a view is created -----
    CView* verifyView(
        CView* view,
        const UIAttributes& attributes,
        const IUIDescription* description) override
    {
        if (EQCurveView* control = dynamic_cast<EQCurveView*>(view); control)
        {
            eqCurveView = control;
            eqCurveView->registerControlListener(this);
            eqCurveView->remember();
            
            update(ParamScLfIn,   kChanged);
            update(ParamScLfType, kChanged);
            update(ParamScLfFreq, kChanged);
            update(ParamScLfGain, kChanged);
            update(ParamScHfIn,   kChanged);
            update(ParamScHfType, kChanged);
            update(ParamScHfFreq, kChanged);
            update(ParamScHfGain, kChanged);
        }
        return view;
    }

    RLFCMP_Controller* mainController;
    Steinberg::Vst::Parameter* ParamScLfIn;
    Steinberg::Vst::Parameter* ParamScLfType;
    Steinberg::Vst::Parameter* ParamScLfFreq;
    Steinberg::Vst::Parameter* ParamScLfGain;
    Steinberg::Vst::Parameter* ParamScHfIn;
    Steinberg::Vst::Parameter* ParamScHfType;
    Steinberg::Vst::Parameter* ParamScHfFreq;
    Steinberg::Vst::Parameter* ParamScHfGain;
    EQCurveView* eqCurveView;
};

//------------------------------------------------------------------------
// Transfer Curve View Controller
//------------------------------------------------------------------------
class TransferCurveViewController
    : public Steinberg::FObject
    , public VSTGUI::DelegationController
    , public VSTGUI::CBaseObject
{
public:
    TransferCurveViewController(
        IController* baseController,
        RLFCMP_Controller* _mainController,
        Steinberg::Vst::Parameter* ParamThreshold,
        Steinberg::Vst::Parameter* ParamKnee,
        Steinberg::Vst::Parameter* ParamRatio,
        Steinberg::Vst::Parameter* ParamMakeup,
        Steinberg::Vst::Parameter* ParamMix
    )
        : DelegationController(baseController)
        , mainController(_mainController)
        , ParamThreshold(ParamThreshold)
        , ParamKnee     (ParamKnee)
        , ParamRatio    (ParamRatio)
        , ParamMakeup   (ParamMakeup)
        , ParamMix      (ParamMix)
        , transferCurveView(nullptr)
    {
        if (ParamThreshold) ParamThreshold->addDependent(this);
        if (ParamKnee     ) ParamKnee     ->addDependent(this);
        if (ParamRatio    ) ParamRatio    ->addDependent(this);
        if (ParamMakeup   ) ParamMakeup   ->addDependent(this);
        if (ParamMix      ) ParamMix      ->addDependent(this);
    }

    ~TransferCurveViewController() override
    {
        if (ParamThreshold) ParamThreshold->removeDependent(this);
        if (ParamKnee     ) ParamKnee     ->removeDependent(this);
        if (ParamRatio    ) ParamRatio    ->removeDependent(this);
        if (ParamMakeup   ) ParamMakeup   ->removeDependent(this);
        if (ParamMix      ) ParamMix      ->removeDependent(this);

        if (transferCurveView)
        {
            transferCurveView->unregisterControlListener (this);
            transferCurveView->forget ();
        }

        mainController->removeUITransferCurveViewController(this);
    }

private:
    using CControl       = VSTGUI::CControl;
    using CView          = VSTGUI::CView;
    using UTF8String     = VSTGUI::UTF8String;
    using UIAttributes   = VSTGUI::UIAttributes;
    using IUIDescription = VSTGUI::IUIDescription;
    using TransferCurveView = VSTGUI::TransferCurveView;

    // FObject
    void PLUGIN_API update( Steinberg::FUnknown* changedUnknown, Steinberg::int32 message) SMTG_OVERRIDE
    {
        if (transferCurveView)
        {
            if (auto* p = Steinberg::FCast<Steinberg::Vst::Parameter>(changedUnknown); p)
            {
                if (message == kChanged)
                {
                    if (p == ParamThreshold && ParamThreshold) transferCurveView->setThreshold(paramThreshold.ToPlain(p->getNormalized()));
                    if (p == ParamKnee      && ParamKnee     ) transferCurveView->setKnee     (paramKnee     .ToPlain(p->getNormalized()));
                    if (p == ParamRatio     && ParamRatio    ) transferCurveView->setRatio    (paramRatio    .ToPlain(p->getNormalized()));
                    if (p == ParamMakeup    && ParamMakeup   ) transferCurveView->setMakeup   (paramMakeup   .ToPlain(p->getNormalized()));
                    if (p == ParamMix       && ParamMix      ) transferCurveView->setMix      (p->getNormalized()); // mix is in 0~1
                    // curveView->invalid();
                }
                else if (message == kWillDestroy)
                {
                    if (p == ParamThreshold && ParamThreshold) { ParamThreshold->removeDependent(this); ParamThreshold = nullptr; }
                    if (p == ParamKnee      && ParamKnee     ) { ParamKnee     ->removeDependent(this); ParamKnee      = nullptr; }
                    if (p == ParamRatio     && ParamRatio    ) { ParamRatio    ->removeDependent(this); ParamRatio     = nullptr; }
                    if (p == ParamMakeup    && ParamMakeup   ) { ParamMakeup   ->removeDependent(this); ParamMakeup    = nullptr; }
                    if (p == ParamMix       && ParamMix      ) { ParamMix      ->removeDependent(this); ParamMix       = nullptr; }
                }
            }
        }
    }
    
    //--- is called when a view is created -----
    CView* verifyView(
        CView* view,
        const UIAttributes& attributes,
        const IUIDescription* description) override
    {
        if (TransferCurveView* control = dynamic_cast<TransferCurveView*>(view); control)
        {
            transferCurveView = control;
            transferCurveView->registerControlListener(this);
            transferCurveView->remember();
            //curveView->setValue(0.0);

            update(ParamKnee,      kChanged);
            update(ParamThreshold, kChanged);
            update(ParamRatio,     kChanged);
            update(ParamMakeup,    kChanged);
            update(ParamMix,       kChanged);
        }
        return view;
    }

    RLFCMP_Controller*         mainController;
    Steinberg::Vst::Parameter* ParamKnee;
    Steinberg::Vst::Parameter* ParamThreshold;
    Steinberg::Vst::Parameter* ParamRatio;
    Steinberg::Vst::Parameter* ParamMakeup;
    Steinberg::Vst::Parameter* ParamMix;
    TransferCurveView*         transferCurveView;
};

//------------------------------------------------------------------------
// VuMeterController
//------------------------------------------------------------------------
class VuMeterController : public VSTGUI::IController, public VSTGUI::ViewListenerAdapter
{
public:
    VuMeterController(RLFCMP_Controller* _mainController) :
        mainController(_mainController),
        vuMeterInL(nullptr),
        vuMeterInR(nullptr),
        vuMeterOutL(nullptr),
        vuMeterOutR(nullptr),
        vuMeterGR(nullptr)
    {
    }
    ~VuMeterController() override
    {
        if (vuMeterInL)  viewWillDelete(vuMeterInL);
        if (vuMeterInR)  viewWillDelete(vuMeterInR);
        if (vuMeterOutL) viewWillDelete(vuMeterOutL);
        if (vuMeterOutR) viewWillDelete(vuMeterOutR);
        if (vuMeterGR)   viewWillDelete(vuMeterGR);

        mainController->removeUIVuMeterController(this);
    }

    void updateVuMeterValue()
    {
        if (mainController) {
            if (vuMeterInL)  vuMeterInL-> setValue(mainController->getVuMeterByTag(vuMeterInL->getTag()));
            if (vuMeterInR)  vuMeterInR-> setValue(mainController->getVuMeterByTag(vuMeterInR->getTag()));
            if (vuMeterOutL) vuMeterOutL->setValue(mainController->getVuMeterByTag(vuMeterOutL->getTag()));
            if (vuMeterOutR) vuMeterOutR->setValue(mainController->getVuMeterByTag(vuMeterOutR->getTag()));
            if (vuMeterGR)   vuMeterGR->  setValue(mainController->getVuMeterByTag(vuMeterGR->getTag()));
        }
    }

private:
    using CControl       = VSTGUI::CControl;
    using CView          = VSTGUI::CView;
    using MyVuMeter      = VSTGUI::MyVuMeter;
    using UTF8String     = VSTGUI::UTF8String;
    using UIAttributes   = VSTGUI::UIAttributes;
    using IUIDescription = VSTGUI::IUIDescription;

    //--- from IControlListener ----------------------
    void valueChanged(CControl* /*pControl*/) SMTG_OVERRIDE {}
    void controlBeginEdit(CControl* /*pControl*/) SMTG_OVERRIDE {}
    void controlEndEdit(CControl* pControl) SMTG_OVERRIDE {}
    //--- is called when a view is created -----
    CView* verifyView(CView* view,
                      const UIAttributes&   /*attributes*/,
                      const IUIDescription* /*description*/) SMTG_OVERRIDE;
    
    //--- from IViewListenerAdapter ----------------------
    //--- is called when a view will be deleted: the editor is closed -----
    void viewWillDelete(CView* view) SMTG_OVERRIDE
    {
        if (dynamic_cast<MyVuMeter*>(view) == vuMeterInL  && vuMeterInL)  { vuMeterInL-> unregisterViewListener(this); vuMeterInL  = nullptr; }
        if (dynamic_cast<MyVuMeter*>(view) == vuMeterInR  && vuMeterInR)  { vuMeterInR-> unregisterViewListener(this); vuMeterInR  = nullptr; }
        if (dynamic_cast<MyVuMeter*>(view) == vuMeterOutL && vuMeterOutL) { vuMeterOutL->unregisterViewListener(this); vuMeterOutL = nullptr; }
        if (dynamic_cast<MyVuMeter*>(view) == vuMeterOutR && vuMeterOutR) { vuMeterOutR->unregisterViewListener(this); vuMeterOutR = nullptr; }
        if (dynamic_cast<MyVuMeter*>(view) == vuMeterGR   && vuMeterGR)   { vuMeterGR->  unregisterViewListener(this); vuMeterGR   = nullptr; }
    }

    RLFCMP_Controller* mainController;
    MyVuMeter* vuMeterInL;
    MyVuMeter* vuMeterInR;
    MyVuMeter* vuMeterOutL;
    MyVuMeter* vuMeterOutR;
    MyVuMeter* vuMeterGR;
};

//------------------------------------------------------------------------
} // namespace yg331
