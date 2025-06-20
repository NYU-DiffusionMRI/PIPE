# PIPE (Protocol Independent Parameter Estimation) toolbox
[This MATLAB toolbox](https://github.com/NYU-DiffusionMRI/PIPE/blob/master/PIPE.m) contains all necessary functions for parameter estimation of the [Standard Model (SM) of diffusion](https://github.com/NYU-DiffusionMRI/SMI) in white matter[^note]. Check [our paper](https://doi.org/10.1016/j.neuroimage.2022.119290) for details on this implementation and on the Standard Model in general. Below we provide instructions on how to run the PIPE toolbox. See the '[PIPE_SM_estimation_human_data.m](https://github.com/NYU-DiffusionMRI/SMI/blob/master/example.m)' script that performs the parameter estimation in an [example dataset](https://cai2r.net/resources/standard-model-of-diffusion-in-white-matter-the-smi-toolbox/).

<br>

## PIPE Authors
- [Santiago Coelho](https://santiagocoelho.github.io/)
- [Els Fieremans](https://www.diffusion-mri.com/who-we-are/els-fieremans/)
- [Dmitry Novikov](https://www.diffusion-mri.com/who-we-are/dmitry-novikov/)

Do not hesitate to reach out to Santiago.Coelho@nyulangone.org (or [@santicoelho](https://x.com/santicoelho) in X/Twitter) for feedback, suggestions, questions, or comments[^note].

## LICENSE

A [US and PCT patent application](https://patents.google.com/patent/WO2024073122A1/en?q=(santiago+coelho%2c+System%2c+method+and+computer-accessible+medium+for+diffusion+mri+shells)&oq=santiago+coelho%2c+System%2c+method+and+computer-accessible+medium+for+diffusion+mri+without+shells+) contains some of the related developments. 

```
%  Authors: Santiago Coelho (santiago.coelho@nyulangone.org), Els Fieremans, Dmitry Novikov
%  Copyright (c) 2025 New York University
%              
%   Permission is hereby granted, free of charge, to any non-commercial entity ('Recipient') obtaining a 
%   copy of this software and associated documentation files (the 'Software'), to the Software solely for
%   non-commercial research, including the rights to use, copy and modify the Software, subject to the 
%   following conditions: 
% 
%     1. The above copyright notice and this permission notice shall be included by Recipient in all copies
%     or substantial portions of the Software. 
% 
%     2. THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
%     NOT LIMITED TO THE WARRANTIESOF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
%     IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BELIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
%     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF ORIN CONNECTION WITH THE
%     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
% 
%     3. In no event shall NYU be liable for direct, indirect, special, incidental or consequential damages
%     in connection with the Software. Recipient will defend, indemnify and hold NYU harmless from any 
%     claims or liability resulting from the use of the Software by recipient. 
% 
%     4. Neither anything contained herein nor the delivery of the Software to recipient shall be deemed to
%     grant the Recipient any right or licenses under any patents or patent application owned by NYU. 
% 
%     5. The Software may only be used for non-commercial research and may not be used for clinical care. 
% 
%     6. Any publication by Recipient of research involving the Software shall cite the references listed
%     below.
%
% REFERENCES:
% - Coelho, S., ..., 2025. ArXiv
% - Coelho, S., Baete, S., Lemberskiy, G., Ades-Aron, B., Barrol, G., Veraart, J., Novikov, D.S., Fieremans, E., 2022. Reproducibility of the Standard Model of diffusion in white matter on clinical MRI systems. Neuroimage. doi: 10.1016/j.neuroimage.2022.119290. Epub 2022 May 8. PMID: 35545197; PMCID: PMC9248353.


[^note]:
    Please cite these works if you use the PIPE toolbox in your publication:
    - Coelho, S., ..., 2025. ArXiv.


