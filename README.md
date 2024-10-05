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
