# SNR Estimation

This repository contains the Matlab code for estimating the Signal-to-noise ratio (SNR) in linear systems with Gaussian matrices. This function is the official implementation of the paper [SNR Estimation in Linear Systems With Gaussian Matrices](https://ieeexplore.ieee.org/abstract/document/8052123), IEEE Signal Processing Letters, v. 24, issue 12, December 2017.  


## Summary 

Given a linear system in the form:

$y = Wx + n$  

where:
- $W$  =  $\Psi^{1/2} \bar{W}$.   
- $\bar{W}$ : is a Gaussian matrix with i.i.d. entries of zero mean and unit variance.
- The entries of $x$, and those of the noise $n$, are drawn from any two distributions and they are i.i.d. with zero mean and unknown variances.

the function snr_estimate estimates the snr of the system. 


### Inputs:

- $y$: The bservations vector of length $M$.
- $\bar{W}$:  Channel/Data matrix of dimension $M \times N$.
- $\Psi$ : Hermitian non negative left correlation matrix of dimension $M \times M$ (can be set to the identity matrix).

### Outputs:

- SNR       :   The SNR of the linear system.
- signal_var:   The varince of $x$ (i.e., $\sigma_{x}^{2})$.
- noise_var :   The noise varince (i.e., $\sigma_{n}^{2})$.




## Citation

Please cite this work if you found it useful:

[SNR Estimation in Linear Systems With Gaussian Matrices](https://ieeexplore.ieee.org/abstract/document/8052123)

```
@article{suliman2017snr,
  title={SNR estimation in linear systems with Gaussian matrices},
  author={Suliman, Mohamed A and Alrashdi, Ayed M and Ballal, Tarig and Al-Naffouri, Tareq Y},
  journal={IEEE Signal Processing Letters},
  volume={24},
  number={12},
  pages={1867--1871},
  year={2017},
  publisher={IEEE}
}
