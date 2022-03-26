#define N 45                     //This is the number of samples taken in one ms
int PinOutR = 19, PinOutL=17;      // Output Pins
int PinInR = 20, PinInL=18;        // Input Pins
int sample_rate=44100;
float valR[N],valL[N];           // variable to store the input value
float outR[N],outL[N];           //variable to store output values
float a0=1.0000, a1=-4.7456, a2=9.0140, a3=-8.5659, a4=4.0723, a5=-0.7748;     //Denominator of transfer function (Chebyshev II)
float b0=0.0175, b1=-0.0520, b2=0.0346, b3=0.0346, b4=-0.0520, b5=0.0175;      // Numerator of transfer function  (Chebyshev II)

void setup()
{
  pinMode(PinOutR, OUTPUT);                       // sets the output pins
  pinMode(PinOutL, OUTPUT);
  memset(valR,0,N);                               // initializes all array values to 0
  memset(valL,0,N);
  memset(outR,0,N);
  memset(outL,0,N);
}

void loop()
{ 
for(int i=0;i<N;i++)                                                                            //This loop reads values from input pins
{
 valR[i] = analogRead(PinInR); 
 if (valR[i]>0.2*1023)
  valR[i]=0;
 outR[i]= a1/a0*outR[i-1] + a2/a0*outR[i-2] +a3/a0*outR[i-3] + a4/a0*outR[i-4] + a5/a0*outR[i-5];       // low-pass filter
 outR[i]=outR[i]-b0/a0*valR[i] +b1/a0*valR[i-1] + b2/a0*valR[i-2] +b3/a0*valR[i-3] + b4/a0*valR[i-4] + b5/a0*valR[i-5];
 analogWrite(PinOutR, -outR[i]/4);                                                               // analogRead values go from 0 to 1023, analogWrite values from 0 to 255
 
 valL[i] = analogRead(PinInL);
  if (valL[i]>0.2*1023)
  valL[i]=0;
 outL[i] = a0*outL[i] +a1*outL[i-1] + a2*outL[i-2] +a3*outL[i-3] + a4*outL[i-4] + a5*outL[i-5];      // low-pass filter
 outL[i]=outL[i]-b0*valL[i] +b1*valL[i-1] + b2*valL[i-2] +b3*valL[i-3] + b4*valL[i-4] + b5*valL[i-5];
 analogWrite(PinOutL, -outL[i]/4);
}
}
