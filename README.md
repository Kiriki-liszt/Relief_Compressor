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

## Why Clean type sounds different?  

When compressing, harmonic distortion happens and fills up a removed loudness.  
In Clean type, that harmonic distortion is surpressed, leading to a sence of feel that it grabs more then it should.  
Our ears are used to hear that volume drop compensation by harmonic distortion, so when it doesn't, it feels somewhat unexpected.  
