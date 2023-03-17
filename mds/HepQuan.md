**HepQuant** section provieds a friendly way to interactively search the corresponding structure of  low-molecular-weight heparin (LMWH) or heparin/heparan sulfate oligosaccharides by  provided of LC-MS m/z and charge.


**Step1**: Deconvolution by DeconTools. The parameter setting XML file can be download at [here](https://raw.githubusercontent.com/XieShaoshuai/MsPHep/main/Example%20data/HeparinParameterFile.xml)

**Step2**: Selective database and adducted form.

**Step3**: Parameters setting. The parameters include ppm, degree of polymerization range, and Min Scan numbers. Users can define the degree of polymerization range to reduce searching time. The value of Min Scan numbers is set for the minimum scan number (MSN) of each peak. Peaks for which the scan number is less than the MSN will be removed before further searching.
**example dataset** Examples for deconvolution Enoxaparin dataset can download at [here] (https://github.com/XieShaoshuai/MsPHep/tree/main/Example%20data/Enoxaparin%20dataset)

