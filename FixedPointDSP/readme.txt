Transmitter and Receiver Design and PART II
for ENEE 630 Project

There are two main executable scripts.

'mainFixed.m' - works on Fixed Point. If you run it, it generates and save 
data file (.mat) and figures (.png)

 'main.m' - works on Double, it use more accurate interpolation for frequency 
estimation and can calculate phase. If you run it, it generates figures (.fig)


Common simulation parameters:
sim_set:
    When it is 0, 
    SNR changes from -3 dB to 15 dB with 0.5 dB increments. Frequenct shift is
    f0 and time shift is t0.

    When it is 1, 
    Time delay changes from -2.5 ms to 2.5 ms with 1/16 ms increments. 
    Frequenct shift is f0 and SNR is 100dB.

    When it is 2, 
    Frequency changes from -1500 Hz to 1500 Hz with 125 Hz increments. 
    Time shift is t0 and SNR is 100dB.


    When it is 3, (only in main.m)
    Phase changes from -0.5 pi to 0.5 pi . 
    Time shift is t0, Frequency shift is f0 and SNR is 100dB.

f0: Frequency shift in Hz when sim_set is 0 or 1 (default 0 Hz)

t0: Time shift in ms when sim_set is 0 or 2 (default 0 ms)

p0=0: Phase shift takes value from -0.5 to 0.5, in Radyan default (0 Rad)

SNR0: SNR value when sim_set is 1 or 2 default (100 dB)

FERlim: Number of the frames with error, simulation stops if it is exceeded.

frame_number: Maximim number of the frames, simulation stops if it is exceeded.

intpol_par: (only in main.m)
	When it is 0, two point linear interpolation is used.
	When it is 1, threshold aided N-points linear interpolation is used.

strsim: (only in mainFixed.m)
	String name to save files.


Some additional functions:

ffttable(size,column_start,column_end,row_start,row_end): Calculate a matrix
W, multipication of W*A gives the FFT of the A where between indices from 
row_start to row_end, where size and column_end is the length of A,
column_start is 1.

mySqr(V): Calculates the magnitude square of the column vector V, which is
V^HV. (^H is the Hermitian operation)

MatrixVecMult(M,V): Performs Matrix Vector multipication for fixed points

myConv16(V1,V2): Apply convoultion for fixed points.


Other functions are defined as comments and/or in the report.