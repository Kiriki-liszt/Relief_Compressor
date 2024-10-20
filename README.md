# Relief_Compressor  

Simple compressor for typical use.  
Comes in three type - Bold, Smooth, Clean.  

|  Type  | Detecter |
|:------:|:--------:|
|  Bold  |   Peak   |
| Smooth |    RMS   |
|  Clean |  Hilbert |

LookAhead about 0.4ms.  
It uses FIR as smoother described in Waves Audio patent(US6535846B1, Expired).  

Peak Linear - Metric Halo ChannelStrip MIO Comp  
Peak Logarithmic - Sonnox Oxford Dynamics  
RMS Linear - AMEK Mastering Compressor  
RMS Logarithmic - SSL Native Bus Compressor 2  

## Types of Detectors  

### Bold  

It's a typical Peak detector.  

- Metric Halo ChannelStrip MIO compressor
- Weiss DS1-MK3 (yes, it is)
- Sononx Oxford Dynamics
- bx SSL, Waves SSL EV2
- UADx SSL G Bus Compressor
- Waves Audio R-comp

### Smooth  

It's a semi-true-RMS detector, using 1-pole filter.  
A true RMS would use square-normalize-sum-sqrt.  

- AMEK Mastering Compressor
- SSL Native Bus Compressor 2
- bx Shadow Hills Mastering Comp
- Ozone 11 Dynamics (RMS)

### Clean  

It's a Hilbert detector used in RMS style.  
It uses two path with about 90-degree phase difference.  
If both path's signal are square-add-sqrted, it gets ideal level envelope.  

- FabFilter Pro-C 2 Mastering
- TDR Kotelvnikov RMS
- Ozone 11 Dynamics (Peak)

## Sidechain Topology  

### Linear  

It detects level in linear gain scale, and then calculates gain reduction.  
Not typical, but certainly favorable.  

- Weiss DS1-MK3
- Metric Halo ChannelStrip MIO compressor
- AMEK Mastering Compressor
- bx SSL, Waves SSL EV2
- bx Shadow Hills Mastering Comp

### Logarithmic  

It calculates gain reduction in log dB, and applies level envelope to gain reduction.  
More common method to use.  

- FabFilter Pro-C 2
- Sonnox Oxford Dynamics
- TDR Kotelvnikov
- SSL Native Channel Strip 2
- UADx SSL G Bus Comp, E Channel Strip

## Why Attack and Release affects Threshold?  

It happens if ripples of detector are not controlled, but it's in sync with signal path so that ripple does not become a problem.  
It can be solved by appling zero-attack moderate-release to rectified signal.  
In case of 'Clean' type, it uses Hilbert transform so this ripple is non-existent.  
However, it overshoots so it need a mild smoother.  

Need to know is, that it changes harmonic pattern generated from compressor.  
Reduction rate of harmonics are more steep.  
Since I don't have a full size HW compressor(what a shame), I have to use MXR studio comp.  
MXR studio compressor is copy of UAD 1176 LN version.  
That compressor had harmonic pattern of one-stage level detection.  
BUT this is analog gear, it does not worry about aliasing, so it's not a problem...  

Have to check out how other HW comprssor's harmonic pattern looks like.  

### List of timing independent threshold(two-stage level detection) compressors  

- Sonnox Oxford Dynamics
- Waves Audio R-Comp
- Weiss DS1-MK3
- Ozone 11 Dynamics
- Pro-C 2 Mastering
- AMEK Mastering Compressor
- Maag MAGNUM-K
- TBTECH Cenozoix Compressor
- TBProAudio Impress 3
- and many analog compressors

### List of timing dependent threshold(one-stage level detection) compressors  

- TDR Kotelnikov
- Metric Halo ChannelStrip
- Pro-C 2 Clean
- most of Waves Audio Compressors
- elysia alpha
- BX Shadow Hills Mastering Comp
- Kiive XTComp
- SPL IRON
- SSL Native X-Comp
- SSL Native Bus Compressor 2
- UADx SSL G Bus Compressor
- and many other digital compressors

## Why Clean type sounds different?  

When compressing, harmonic distortion happens and fills up a removed loudness.  
In Clean type, that harmonic distortion is surpressed, leading to a sence of feel that it grabs more then it should.  
Our ears are used to hear that volume drop compensation by harmonic distortion, so when it doesn't, it feels somewhat unexpected.  

UADx dbx 160 has some kind of Linear RMS, Hilbert detection going on...  
