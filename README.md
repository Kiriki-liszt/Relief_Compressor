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

| SampleRate[Hz] | LookAhead[smpl] |
|:--------------:|:---------------:|
|      44100     |        16       |
|      48000     |        18       |
|      88200     |        34       |
|      96000     |        38       |
|     192000     |        78       |
|     384000     |       158       |

## Why Attack and Release affects Threshold?  

It happens as ripple of detector are not controlled, but it's in sync with signal path so that ripple does not become problem.  
It can be solved by appling zero-attack moderate-release to rectified signal.  
In case of 'Clean' type, it uses Hilbert transform so this ripple is non-existent.  

However, it changes harmonic pattern generated from compressor.  
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

## Bold type  

It resembles Metric Halo Channel Strip - MIO compressor.  
Without LookAhead, it's almost identical.  
Attack and Release affects effective threshold.  

## Why Clean type sounds different?  

When compressing, harmonic distortion happens and fills up a removed loudness.  
In Clean type, that harmonic distortion is surpressed, leading to a sence of feel that it grabs more then it should.  
Our ears are used to hear that volume drop compensation by harmonic distortion, so when it doesn't, it feels somewhat unexpected.  
