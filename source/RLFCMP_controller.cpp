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


static const std::string kAttrVuOnColor  = "vu-on-color";
static const std::string kAttrVuOffColor = "vu-off-color";
static const std::string kAttrPDclick    = "click-behave";
static const std::string kAttrPDMin      = "update-min";
static const std::string kAttrPDMax      = "update-max";
namespace VSTGUI {
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
MyVUMeterFactory          __gMyVUMeterFactory;
} // namespace VSTGUI


namespace yg331 {
//------------------------------------------------------------------------
// VuMeterController
//------------------------------------------------------------------------
VSTGUI::CView* DetectorIndicatorController::verifyView(
                                            VSTGUI::CView* view,
                                            const VSTGUI::UIAttributes&   attributes,
                                            const VSTGUI::IUIDescription* /*description*/
)
{
    if (UIViewSwitchContainer* control = dynamic_cast<UIViewSwitchContainer*>(view); control)
    {
        viewSwitch = control;
        viewSwitch->registerViewListener(this);
    }

    return view;
}

void DetectorIndicatorController::viewWillDelete(VSTGUI::CView* view)
{
    if (dynamic_cast<UIViewSwitchContainer*>(view) == viewSwitch  && viewSwitch)    { viewSwitch-> unregisterViewListener(this); viewSwitch  = nullptr; }
}

VSTGUI::CView* VuMeterController::verifyView(CView* view,
                  const UIAttributes&   /*attributes*/,
                  const IUIDescription* /*description*/)
{
#define minVU -30.0 // need that margin at bottom
#define maxVU 0.0
    if (MyVuMeter* control = dynamic_cast<MyVuMeter*>(view); control) {
        if (control->getTag() == kInLRMS || control->getTag() == kInLPeak)  {
            vuMeterInL  = control;
            vuMeterInL->setMin(minVU);
            vuMeterInL->setMax(maxVU);
            vuMeterInL->setDefaultValue(minVU);
            vuMeterInL->registerViewListener(this);
        }
        if (control->getTag() == kInRRMS || control->getTag() == kInRPeak)  {
            vuMeterInR  = control;
            vuMeterInR->setMin(minVU);
            vuMeterInR->setMax(maxVU);
            vuMeterInR->setDefaultValue(minVU);
            vuMeterInR-> registerViewListener(this);
        }
        if (control->getTag() == kOutLRMS || control->getTag() == kOutLPeak) {
            vuMeterOutL = control;
            vuMeterOutL->setMin(minVU);
            vuMeterOutL->setMax(maxVU);
            vuMeterOutL->setDefaultValue(minVU);
            vuMeterOutL->registerViewListener(this);
        }
        if (control->getTag() == kOutRRMS || control->getTag() == kOutRPeak) {
            vuMeterOutR = control;
            vuMeterOutR->setMin(minVU);
            vuMeterOutR->setMax(maxVU);
            vuMeterOutR->setDefaultValue(minVU);
            vuMeterOutR->registerViewListener(this);
        }
        if (control->getTag() == kGainReduction)   {
            vuMeterGR   = control;
            vuMeterGR->setMin(minVU);
            vuMeterGR->setMax(maxVU);
            vuMeterGR->setDefaultValue(0.0);
            vuMeterGR->  registerViewListener(this);
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
    
    tag          = kParamSidechainFilter;
    stepCount    = 1;
    defaultVal   = 0;
    flags        = Vst::ParameterInfo::kCanAutomate | Vst::ParameterInfo::kIsList;
    parameters.addParameter(STR16("SidechainFilter"), nullptr, stepCount, defaultVal, flags, tag);

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
    
    tag          = kParamBias;
    flags        = Vst::ParameterInfo::kCanAutomate;
    minPlain     = minBias;
    maxPlain     = maxBias;
    defaultPlain = dftBias;
    stepCount    = 0;
    auto* ParamBias = new LinRangeParameter(STR16("Bias"), tag, STR16("dB"), minPlain, maxPlain, defaultPlain, stepCount, flags);
    ParamBias->setPrecision(1);
    parameters.addParameter(ParamBias);

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
    auto* ParamRatio = new Vst::RangeParameter(STR16("Ratio"), tag, STR16(""), minPlain, maxPlain, defaultPlain, stepCount, flags);
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
    if (zoomFactors.empty())
    {
        zoomFactors.push_back(ZoomFactor(STR("50%"),  0.50)); // 0/6
        zoomFactors.push_back(ZoomFactor(STR("75%"),  0.75)); // 1/6
        zoomFactors.push_back(ZoomFactor(STR("100%"), 1.00)); // 2/6
        zoomFactors.push_back(ZoomFactor(STR("125%"), 1.25)); // 3/6
        zoomFactors.push_back(ZoomFactor(STR("150%"), 1.50)); // 4/6
        zoomFactors.push_back(ZoomFactor(STR("175%"), 1.75)); // 5/6
        zoomFactors.push_back(ZoomFactor(STR("200%"), 2.00)); // 6/6
    }

    Vst::StringListParameter* zoomParameter = new Vst::StringListParameter(STR("Zoom"), kParamZoom);
    for (ZoomFactorVector::const_iterator it = zoomFactors.begin(), end = zoomFactors.end(); it != end; ++it)
    {
        zoomParameter->appendString(it->title);
    }
    zoomParameter->setNormalized(zoomParameter->toNormalized(2));
    zoomParameter->addDependent(this);
    uiParameters.addParameter(zoomParameter);
    
    return result;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Controller::terminate ()
{
    // Here the Plug-in will be de-instantiated, last possibility to remove some memory!
    getParameterObject(kParamZoom)->removeDependent(this);

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

    int32           savedBypass     = 0;
    // Vst::ParamValue savedZoom       = 0.0;
    // Vst::ParamValue savedOS         = 0.0;
    int32           savedSidechainFilter = 0;
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
    // if (streamer.readDouble(savedZoom)       == false) savedZoom       = 2.0 / 6.0;
    // if (streamer.readDouble(savedOS)         == false) savedOS         = 0.0;
    if (streamer.readInt32 (savedSidechainFilter) == false) savedSidechainFilter = 0;
    if (streamer.readDouble(savedAttack)     == false) savedAttack     = paramAttack.ToNormalized(dftAttack);
    if (streamer.readDouble(savedRelease)    == false) savedRelease    = paramRelease.ToNormalized(dftRelease);
    if (streamer.readDouble(savedBias)       == false) savedBias       = paramBias.ToNormalized(dftBias);
    if (streamer.readDouble(savedThreshold)  == false) savedThreshold  = paramThreshold.ToNormalized(dftThreshold);
    if (streamer.readDouble(savedRatio)      == false) savedRatio      = paramRatio.ToNormalized(dftRatio);
    if (streamer.readDouble(savedKnee)       == false) savedKnee       = paramKnee.ToNormalized(dftKnee);
    if (streamer.readDouble(savedMakeup)     == false) savedMakeup     = paramMakeup.ToNormalized(dftMakeup);
    if (streamer.readDouble(savedMix)        == false) savedMix        = paramMix.ToNormalized(dftMix);
    if (streamer.readDouble(savedOutput)     == false) savedOutput     = paramOutput.ToNormalized(dftOutput);
    if (streamer.readInt32 (savedSoftBypass) == false) savedSoftBypass = 0;

    setParamNormalized(kParamBypass,     savedBypass ? 1 : 0);
    // setParamNormalized(kParamZoom,       savedZoom);
    // setParamNormalized(kParamOS,         savedOS);
    setParamNormalized(kParamSidechainFilter,     savedSidechainFilter ? 1 : 0);
    setParamNormalized(kParamAttack,     savedAttack);
    setParamNormalized(kParamRelease,    savedRelease);
    setParamNormalized(kParamBias,       savedBias);
    setParamNormalized(kParamThreshold,  savedThreshold);
    setParamNormalized(kParamRatio,      savedRatio);
    setParamNormalized(kParamKnee,       savedKnee);
    setParamNormalized(kParamMakeup,     savedMakeup);
    setParamNormalized(kParamMix,        savedMix);
    setParamNormalized(kParamOutput,     savedOutput);
    setParamNormalized(kParamSoftBypass, savedSoftBypass ? 1 : 0);
    
    return kResultOk;
}

//------------------------------------------------------------------------
tresult PLUGIN_API RLFCMP_Controller::setState (IBStream* state)
{
    // Here you get the state of the controller
    if (!state)
        return kResultFalse;

    IBStreamer streamer(state, kLittleEndian);
    
    ParamValue savedZoom = 2.0 / 6.0;
    
    if (streamer.readDouble(savedZoom)          == false) return kResultFalse;
    
    pZoom          = savedZoom;
    
    setParamNormalized(kParamZoom,       savedZoom);
    
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
    
    pZoom = getParamNormalized(kParamZoom);
    
    if (streamer.writeDouble(pZoom) == false) return kResultFalse;
    
    return kResultTrue;
}

//------------------------------------------------------------------------
VSTGUI::IController* RLFCMP_Controller::createSubController (VSTGUI::UTF8StringPtr name,
                                                             const VSTGUI::IUIDescription* /*description*/,
                                                             VSTGUI::VST3Editor* /*editor*/)
{
    if (VSTGUI::UTF8StringView (name) == "UIDetectorIndicatorController")
    {
        auto* controller = new DetectorIndicatorController (this);
        addDetectorIndicatorController (controller);
        return controller;
    }
    if (VSTGUI::UTF8StringView(name) == "VuMeterController")
    {
        auto* controller = new VuMeterController(this);
        addUIVuMeterController(controller);
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
}

//------------------------------------------------------------------------
void RLFCMP_Controller::editorRemoved(Vst::EditorView* editor)
{
    editors.erase(std::find(editors.begin(), editors.end(), editor));
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

        if (message->getAttributes ()->getFloat ("Input L", getValue) == kResultTrue) vuInLPeak  = getValue;
        if (message->getAttributes ()->getFloat ("Input R", getValue) == kResultTrue) vuInRPeak  = getValue;
        if (message->getAttributes ()->getFloat ("Output L", getValue) == kResultTrue) vuOutLPeak  = getValue;
        if (message->getAttributes ()->getFloat ("Output R", getValue) == kResultTrue) vuOutRPeak  = getValue;
        if (message->getAttributes ()->getFloat ("Gain Reduction", getValue) == kResultTrue)
        {
            vuGainReduction  = getValue;
            
            if (!vuMeterControllers.empty())
            {
                for (auto iter = vuMeterControllers.begin(); iter != vuMeterControllers.end(); iter++)
                {
                    (*iter)->updateVuMeterValue();
                }
            }
        }
        
        return kResultOk;
    }
    return EditControllerEx1::notify(message);
}

Steinberg::Vst::ParamValue RLFCMP_Controller::getVuMeterByTag(Steinberg::Vst::ParamID tag)
{
   switch (tag) {
       case kInMono:     return vuInMono;    break;
       case kInLRMS:    return vuInLRMS;    break;
       case kInRRMS:    return vuInRRMS;    break;
       case kInLPeak:    return vuInLPeak;    break;
       case kInRPeak:    return vuInRPeak;    break;
       case kOutMono:    return vuOutMono;   break;
       case kOutLRMS:   return vuOutLRMS;   break;
       case kOutRRMS:   return vuOutRRMS;   break;
       case kOutLPeak:   return vuOutLPeak;   break;
       case kOutRPeak:   return vuOutRPeak;   break;
       case kGainReduction:     return vuGainReduction;     break;
       default: break;
   }
   return 0;
}
//------------------------------------------------------------------------
} // namespace yg331
