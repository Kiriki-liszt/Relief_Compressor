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
    void   setLfIn  (bool value)   { LF_SVF.setIn(value); }
    bool   getLfIn  () const       { return LF_SVF.getIn(); }
    void   setLfType(int value)    { LF_SVF.setType(value); }
    int    getLfType() const       { return LF_SVF.getType(); }
    void   setLfFreq(double value) { LF_SVF.setFreq(value); }
    int    getLfFreq() const       { return LF_SVF.getFreq(); }
    void   setLfGain(double value) { LF_SVF.setGain(value); }
    int    getLfGain() const       { return LF_SVF.getGain(); }

    void   setHfIn  (bool value)   { HF_SVF.setIn(value); }
    bool   getHfIn  () const       { return HF_SVF.getIn(); }
    void   setHfType(int value)    { HF_SVF.setType(value); }
    int    getHfType() const       { return HF_SVF.getType(); }
    void   setHfFreq(double value) { HF_SVF.setFreq(value); }
    int    getHfFreq() const       { return HF_SVF.getFreq(); }
    void   setHfGain(double value) { HF_SVF.setGain(value); }
    int    getHfGain() const       { return HF_SVF.getGain(); }
    
    void   makeSVF() { LF_SVF.makeSVF(); HF_SVF.makeSVF(); }

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

    yg331::SVF_12 LF_SVF;
    yg331::SVF_12 HF_SVF;

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
    void   setThreshold(double value) { threshold = value; }
    double getThreshold() const { return threshold; }

    void   setKnee(double value) { knee = value; kneeHalf  = knee / 2.0; }
    double getKnee() const { return knee; }

    void   setRatio(double value) { ratio = value; slope = 1.0 / ratio - 1.0; }
    double getRatio() const { return ratio; }
    
    void   setMakeup(double value) { makeup = value; }
    double getMakeup() const { return makeup; }

    void   setMix(double value) { mix = value; }
    double getMix() const { return mix; }
    
    void   setVuInMono(double value) { vuInMono = value; }
    double getVuInMono() const { return vuInMono; }

    // get/set Attributes
    void   setBackColor(CColor color) { if (BackColor != color) { BackColor = color; setDirty(true); } }
    CColor getBackColor() const { return BackColor; }

    void   setLineColor(CColor color) { if (LineColor != color) { LineColor = color; setDirty(true); } }
    CColor getLineColor() const { return LineColor; }

    void   setLineOffColor(CColor color) { if (LineOffColor != color) { LineOffColor = color; setDirty(true); } }
    CColor getLineOffColor() const { return LineOffColor; }

    // overrides
    void setDirty(bool state) override { CView::setDirty(state); };
    void draw(CDrawContext* _pContext) override;
    void setViewSize(const CRect& newSize, bool invalid = true) override { CControl::setViewSize(newSize, invalid); };
    bool sizeToFit() override ;
    void onIdle() override { invalid(); };

    CLASS_METHODS(TransferCurveView, CControl)

protected:
    ~TransferCurveView() noexcept override {};

    CColor        BackColor;
    CColor        LineColor;
    CColor        LineOffColor;
    
    double vuInMono = -60.0;

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
class ClickResetParamDisplay : public CParamDisplay
{
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
    
    CLASS_METHODS(ClickResetParamDisplay, CParamDisplay)
 
protected:
    ~ClickResetParamDisplay() noexcept override {};
    
    int32_t _style;
    float   directValue = 0.0;
    bool    over = false;
    CColor  originalBack;
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
    
    // don't mess with it's real value, else it need some listener and gets messy...
    void setValue_(float v) {plainValue = v;}

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
    
    double plainValue = 0.0;
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
    
    // sub-controller list
    using UIVuMeterControllerList = std::vector<VuMeterController*>;
    using UIEQCurveViewControllerList = std::vector<EQCurveViewController*>;
    using UITransferCurveViewControllerList = std::vector<TransferCurveViewController*>;
    
    UIVuMeterControllerList vuMeterControllers;
    UIEQCurveViewControllerList eqCurveViewControllers;
    UITransferCurveViewControllerList transferCurveViewControllerControllers;

    bool       bypass       = dftBypass;
    int32      zoom         = dftZoom;
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
    ParamValue _release      = dftRelease;

    ParamValue threshold    = dftThreshold;
    ParamValue ratio        = dftRatio;
    ParamValue knee         = dftKnee;
    ParamValue makeup       = DecibelConverter::ToGain(dftMakeup);

    ParamValue mix          = dftMix/maxMix;
    ParamValue inputGain    = DecibelConverter::ToGain(dftInput);
    ParamValue outputGain   = DecibelConverter::ToGain(dftOutput);
    bool       softBypass   = dftSoftBypass;
    
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
class EQCurveViewController :
    public Steinberg::FObject,
    public VSTGUI::DelegationController
    //public VSTGUI::CBaseObject
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
                           Steinberg::Vst::Parameter* ParamScHfGain) :
        DelegationController(baseController),
        mainController(_mainController),
        ParamScLfIn  (ParamScLfIn),
        ParamScLfType(ParamScLfType),
        ParamScLfFreq(ParamScLfFreq),
        ParamScLfGain(ParamScLfGain),
        ParamScHfIn  (ParamScHfIn),
        ParamScHfType(ParamScHfType),
        ParamScHfFreq(ParamScHfFreq),
        ParamScHfGain(ParamScHfGain),
        eqCurveView(nullptr)
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
        
        //if (eqCurveView)
        //{
        //    eqCurveView->unregisterControlListener (this);
        //    eqCurveView->forget ();
        //}

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
    void PLUGIN_API update (Steinberg::FUnknown* changedUnknown, Steinberg::int32 message) SMTG_OVERRIDE
    {
        if (eqCurveView)
        {
            if (auto* p = Steinberg::FCast<Steinberg::Vst::Parameter>(changedUnknown); p)
            {
                if (message == kChanged)
                {
                    if (p == ParamScLfIn   && ParamScLfIn  ) eqCurveView->setLfIn  (p->getNormalized());
                    if (p == ParamScLfType && ParamScLfType) eqCurveView->setLfType((paramScLfType.ToPlain(p->getNormalized())==0)?SVF_12::tHighPass:SVF_12::tLowShelf);
                    if (p == ParamScLfFreq && ParamScLfFreq) eqCurveView->setLfFreq(paramScLfFreq.ToPlain(p->getNormalized()));
                    if (p == ParamScLfGain && ParamScLfGain) eqCurveView->setLfGain(paramScLfGain.ToPlain(p->getNormalized()));
                    if (p == ParamScHfIn   && ParamScHfIn  ) eqCurveView->setHfIn  (p->getNormalized());
                    if (p == ParamScHfType && ParamScHfType) eqCurveView->setHfType((paramScHfType.ToPlain(p->getNormalized())==0)?SVF_12::tLowPass:SVF_12::tHighShelf);
                    if (p == ParamScHfFreq && ParamScHfFreq) eqCurveView->setHfFreq(paramScHfFreq.ToPlain(p->getNormalized()));
                    if (p == ParamScHfGain && ParamScHfGain) eqCurveView->setHfGain(paramScHfGain.ToPlain(p->getNormalized()));
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
            //eqCurveView->registerControlListener(this);
            //eqCurveView->remember();
            
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
class TransferCurveViewController :
    public Steinberg::FObject,
    public VSTGUI::DelegationController // IControlListener + IController
    //public VSTGUI::CBaseObject
{
public:
    TransferCurveViewController(
        IController* baseController,
        RLFCMP_Controller* mainController,
        Steinberg::Vst::Parameter* ParamThreshold,
        Steinberg::Vst::Parameter* ParamKnee,
        Steinberg::Vst::Parameter* ParamRatio,
        Steinberg::Vst::Parameter* ParamMakeup,
        Steinberg::Vst::Parameter* ParamMix
    ) :
        DelegationController(baseController),
        mainController(mainController),
        ParamThreshold(ParamThreshold),
        ParamKnee     (ParamKnee),
        ParamRatio    (ParamRatio),
        ParamMakeup   (ParamMakeup),
        ParamMix      (ParamMix),
        transferCurveView(nullptr)
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

        //if (transferCurveView)
        //{
        //    transferCurveView->unregisterControlListener (this);
        //    transferCurveView->forget ();
        //}

        mainController->removeUITransferCurveViewController(this);
    }
    
    void setVuInMono(double val)
    {
        if (transferCurveView) transferCurveView->setVuInMono(val);
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
    //------------------------------------------------------------------------
    // void valueChanged (CControl* pControl)
    /* do sth like these
        {
            CXYPad::calculateXY (pControl->getValue (), x, y);

            auto xId = xParam->getInfo ().id;
            if (editController->setParamNormalized (xId, x) == Steinberg::kResultTrue)
                editController->performEdit (xId, editController->getParamNormalized (xId));
        }
     */
    
    //--- is called when a view is created -----
    CView* verifyView(
        CView* view,
        const UIAttributes& attributes,
        const IUIDescription* description) override
    {
        if (TransferCurveView* control = dynamic_cast<TransferCurveView*>(view); control)
        {
            transferCurveView = control;
            // CControl has main IControlListener* listener, and sub ViewEventListenerAdapter as impl
            // propagtes subListener's (beginEdit, endEdit, valueChanged) methods
            // so, when CControl::beginEdit () is called, it iterates subListeners and call controlBeginEdit (this)
            // However, this TransferCurveView does not change it's value, so listener is not needed.
            //transferCurveView->registerControlListener(this);
            //transferCurveView->remember();
            //curveView->setValue(0.0);

            update(ParamKnee,      kChanged);
            update(ParamThreshold, kChanged);
            update(ParamRatio,     kChanged);
            update(ParamMakeup,    kChanged);
            update(ParamMix,       kChanged);
        }
        return view;
    }

    IController*                baseController;
    RLFCMP_Controller*          mainController;
    Steinberg::Vst::Parameter*  ParamKnee;
    Steinberg::Vst::Parameter*  ParamThreshold;
    Steinberg::Vst::Parameter*  ParamRatio;
    Steinberg::Vst::Parameter*  ParamMakeup;
    Steinberg::Vst::Parameter*  ParamMix;
    TransferCurveView*          transferCurveView;
};

//------------------------------------------------------------------------
// VuMeterController
//------------------------------------------------------------------------
class VuMeterController :
    public VSTGUI::CBaseObject,
    public VSTGUI::DelegationController // IControlListener + IController
    // public VSTGUI::ViewListenerAdapter   // IViewListener
{
public:
    VuMeterController(IController* baseController, RLFCMP_Controller* mainController) :
        DelegationController(baseController),
        mainController(mainController),
        vuMeterInL(nullptr),
        vuMeterInR(nullptr),
        vuMeterOutL(nullptr),
        vuMeterOutR(nullptr),
        vuMeterGR(nullptr),
        pdGainReduction(nullptr)
    {}
    ~VuMeterController() override
    {
        vuMeterInL      = nullptr;
        vuMeterInR      = nullptr;
        vuMeterOutL     = nullptr;
        vuMeterOutR     = nullptr;
        vuMeterGR       = nullptr;
        pdGainReduction = nullptr;

        mainController->removeUIVuMeterController(this);
    }
    enum
    {
        kInMono = 100,
        kInLRMS,
        kInRRMS,
        kInLPeak,
        kInRPeak,
        kOutMono,
        kOutLRMS,
        kOutRRMS,
        kOutLPeak,
        kOutRPeak,
        kGainReduction
    };
    
    double getVuMeterByTag(Steinberg::Vst::ParamID tag)
    {
       switch (tag) {
           case kInMono:        return vuInMono;        break;
           case kInLRMS:        return vuInLRMS;        break;
           case kInRRMS:        return vuInRRMS;        break;
           case kInLPeak:       return vuInLPeak;       break;
           case kInRPeak:       return vuInRPeak;       break;
           case kOutMono:       return vuOutMono;       break;
           case kOutLRMS:       return vuOutLRMS;       break;
           case kOutRRMS:       return vuOutRRMS;       break;
           case kOutLPeak:      return vuOutLPeak;      break;
           case kOutRPeak:      return vuOutRPeak;      break;
           case kGainReduction: return vuGainReduction; break;
           default: break;
       }
       return 0;
    }
    
    void setVuMeterByTag(double v, Steinberg::Vst::ParamID tag)
    {
       switch (tag) {
           case kInMono:        vuInMono = v;        break;
           case kInLRMS:        vuInLRMS = v;        break;
           case kInRRMS:        vuInRRMS = v;        break;
           case kInLPeak:       vuInLPeak = v;       break;
           case kInRPeak:       vuInRPeak = v;       break;
           case kOutMono:       vuOutMono = v;       break;
           case kOutLRMS:       vuOutLRMS = v;       break;
           case kOutRRMS:       vuOutRRMS = v;       break;
           case kOutLPeak:      vuOutLPeak = v;      break;
           case kOutRPeak:      vuOutRPeak = v;      break;
           case kGainReduction: vuGainReduction = v; break;
           default: break;
       }
    }
    
    void updateVuMeterValue()
    {
        if (mainController != nullptr) {
            if (vuMeterInL)         vuMeterInL->    setValue_(getVuMeterByTag(vuMeterInL->getTag()));
            if (vuMeterInR)         vuMeterInR->    setValue_(getVuMeterByTag(vuMeterInR->getTag()));
            if (vuMeterOutL)        vuMeterOutL->   setValue_(getVuMeterByTag(vuMeterOutL->getTag()));
            if (vuMeterOutR)        vuMeterOutR->   setValue_(getVuMeterByTag(vuMeterOutR->getTag()));
            if (vuMeterGR)          vuMeterGR->     setValue_(getVuMeterByTag(vuMeterGR->getTag()));
            if (pdGainReduction)    pdGainReduction->setValue(getVuMeterByTag(vuMeterGR->getTag()));
        }
    }
    
private:
    using CControl                  = VSTGUI::CControl;
    using CView                     = VSTGUI::CView;
    using MyVuMeter                 = VSTGUI::MyVuMeter;
    using ClickResetParamDisplay    = VSTGUI::ClickResetParamDisplay;
    using UTF8String                = VSTGUI::UTF8String;
    using UIAttributes              = VSTGUI::UIAttributes;
    using IUIDescription            = VSTGUI::IUIDescription;

    //--- from IControlListener ----------------------
    //void valueChanged(CControl* /*pControl*/) SMTG_OVERRIDE {}
    //void controlBeginEdit(CControl* /*pControl*/) SMTG_OVERRIDE {}
    //void controlEndEdit(CControl* pControl) SMTG_OVERRIDE {}
    //--- is called when a view is created -----
    CView* verifyView(CView* view,
                      const UIAttributes&   /*attributes*/,
                      const IUIDescription* /*description*/) SMTG_OVERRIDE;

    RLFCMP_Controller* mainController = nullptr;
    MyVuMeter* vuMeterInL   = nullptr;
    MyVuMeter* vuMeterInR   = nullptr;
    MyVuMeter* vuMeterOutL  = nullptr;
    MyVuMeter* vuMeterOutR  = nullptr;
    MyVuMeter* vuMeterGR    = nullptr;
    ClickResetParamDisplay* pdGainReduction = nullptr;
    
    double vuInMono   = -120.0;
    double vuInLRMS   = -120.0;
    double vuInRRMS   = -120.0;
    double vuInLPeak  = -120.0;
    double vuInRPeak  = -120.0;
    double vuOutMono  = -120.0;
    double vuOutLRMS  = -120.0;
    double vuOutRRMS  = -120.0;
    double vuOutLPeak = -120.0;
    double vuOutRPeak = -120.0;
    double vuGainReduction = 0.0;
};

//------------------------------------------------------------------------
} // namespace yg331
