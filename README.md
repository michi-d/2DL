# 2DL algorithm

This is the C++ implementation of the 2-dimensional matrix similarity matching algorithm which I developed in order to assess similarity among neurograms to predict sound or speech quality from computational models of hearing.

The C++ function is wrapped as a MEX-file for usage with MATLAB.

The algorithm is an extension of the 1-dimensional Victor-Purpura spike metric, generalized for the 2-dimensional case.

For further information on the method see the publications:

* Drews, M., Nicoletti, M., Hemmert, W., & Rini, S. (2013, May). The neurogram matching similarity index (NMSI) for the assessment of similarities among neurograms. In 2013 IEEE International Conference on Acoustics, Speech and Signal Processing (pp. 1162-1166). IEEE.
Chicago	


* Drews, M., Schlapak, S., Rini, S., Nicoletti, M., & Hemmert, W. (2013). A Neurogram Matching Similarity Index (NMSI) for the Assessment of Audio Quality. In Proc. 4th International Workshop on Perceptual Quality of Systems (PQS 2013) (pp. 96-101).
Chicago	

## Usage

 ´spike_metric(A, B, approx_level, qt, qch, qa, c_add, c_del, one);'

**INPUT**      : A, B, approx_steps, qt, qch, qa, c_add, c_del, one


**A,B**         : MATLAB arrays (max. 2-dimensional) which contain the amplitudes of a spike signal at a certain time and in a certain channel. The time is given by the horizontal position in the array (second dimension: j). The channel is given by the vertical position in the array (first dimension: i)

**approx_steps**: Level of approximation, number of time sections to build from the arrays A and B. For approx_steps == 1 the result is exact.

**qt**          : Cost of shifting an element by an index 1 in time direction. Cost in relation to real time is then dependent from the sampling rate.

**qc**          : Cost of shifting an element by an index 1 in channel direction. Cost in relation to frequency is then dependent from the frequency resolution.

**qa**          : Cost of changing the amplitude of an element by a value of 1.

**c_add**       : Cost of adding an element at an arbitrary position.

**c_del**       : Cost of deleting an element at an arbitrary position (normally equals c_add).

**one**         : Neutral element/value which indicates the absence of a spike signal (normally 0)


**OUTPUT**      : out

**out**         : An 1D-array which contains the distances of the corresponding time sections. Length is naturally given by approx_steps. For approx_steps == 1 the result is exact.
