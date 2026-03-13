# FFT-DRT
This repository contains some of the source code used for the paper titled Fast Fourier Transform-Based Distribution of Relaxation Times Analysis for Efficient and Flexible Time-Domain Electrochemical Impedance Characterization in the Journal of the Electrochemical Society, 172 (2025) 026504. https://iopscience.iop.org/article/10.1149/1945-7111/adb33f/meta. It is available online and in the docs folder.

Electrochemical impedance spectroscopy (EIS) has been widely utilized to characterize a plethora of electrochemical devices, including batteries and fuel cells [1-2] due to its simplicity, non-invasiveness, and ability to probe frequencies from several MHz to a few mHz [3-5]. To overcome the limitations of the standard equivalent circuits for modeling EIS spectra, the distribution of relaxation times (DRT) has arisen as a valuable tool to identify and investigate the properties of physico-chemical processes [1-2]. Typically, the DRT is recovered using regularized regression using the data measured in the frequency domain [6-8]. Alternatively, and to reduce the measurement time at low frequencies, the DRT can also be deconvolved using time data [2]. In this article, we leverage the Fourier transform to advance DRT deconvolution from time data, both in terms of accuracy and flexibility, which we illustrate with artificial battery data [9]. 


# Dependencies
numpy 

scipy

matplotlib

time

cxpy

bayes_opt

pandas


# Tutorials
1. **PRBS_2ZARC.ipynb**: this notebook shows how to deconvolve the DRT from time data generated with the pseudo-random binary sequence and the series association of two ZARC circuits;
2. **logarithmic-chirp_PWC.ipynb**: this notebook shows how to deconvolve the DRT from time data generated with the logarithmic chirp and the piecewise-constant element;
3. **PRBS_battery.ipynb** : this notebook shows how to deconvolve the DRT from time data generated with the pseudo-random binary sequence and the series association of a ZARC circuit and a Warburg element.

# Citation

```
@article{py2025fast,
  title={Fast Fourier transform-based distribution of relaxation times analysis for efficient and flexible time-domain electrochemical impedance characterization},
  author={Py, Baptiste and Ciucci, Francesco},
  journal={Journal of The Electrochemical Society},
  volume={172},
  number={2},
  pages={026504},
  year={2025},
  publisher={IOP Publishing}
}
```

# References

[1] Z. Wang, Y. Wang, F. Ciucci, Distribution of relaxation times: Foundations, methods, diagnostics, and prognosis for electrochemical systems, Curr. Opin. Electrochem. (2025) 101789. https://doi.org/10.1016/j.coelec.2025.101789.

[2] A. Maradesa, B. Py, F. Ciucci et al., Advancing electrochemical impedance analysis through innovations in the distribution of relaxation times method, Joule.  8 (2024) 1958-1981. https://www.cell.com/joule/fulltext/S2542-4351(24)00236-8.

[3] B. Py, A. Maradesa, F. Ciucci, Gaussian processes for the analysis of electrochemical impedance spectroscopy data: Prediction, filtering, and active learning. Electrochim. Acta. 439 (2023) 141688. https://doi.org/10.1016/j.electacta.2022.141688.

[4] B. Py, C. Zhao, F. Ciucci, Q. Meyer, Gaussian processes for fast and accurate measurements of the polarization resistance of hydrogen fuel cells from impedance spectroscopy, J. Electrochem. Society. 172 (2025) 074502.https://iopscience.iop.org/article/10.1149/1945-7111/ade82c/meta.

[5] B. Py, Z. Wang, Y. Wang, F. Ciucci, Entropy-based regularized regression for advanced distribution of relaxation times deconvolution, J. Power Sources. 644 (2025) 236910. https://doi.org/10.1016/j.jpowsour.2025.236910.

[6] T.H. Wan, M. Saccoccio, C. Chen, F. Ciucci, Influence of the discretization methods on the distribution of relaxation times deconvolution: implementing radial basis functions with DRTtools, Electrochim. Acta. 184 (2015) 483-499. https://doi.org/10.1016/j.electacta.2015.09.097.

[7] A. Maradesa, B. Py, T.H. Wan, M.B. Effat, F. Ciucci, Selecting the regularization parameter in the distribution of relaxation times, J. Electrochem. Society. 170-3 (2023) 030502. https://doi.org/10.1149/1945-7111/acbca4.

[8] B. Py, F. Ciucci, Beyond ridge regression: Enhancing distribution of relaxation times deconvolution, J. Electrochem. Society. 171 (2024) 060529. https://iopscience.iop.org/article/10.1149/1945-7111/ad576a/meta.

[9] B. Py, F. Ciucci, Fast Fourier transform-based distribution of relaxation times analysis for efficient and flexible time-domain electrochemical impedance characterization, J. Electrochem. Society. 172 (2025) 026504. https://iopscience.iop.org/article/10.1149/1945-7111/adb33f/meta.
