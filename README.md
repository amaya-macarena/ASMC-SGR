
# ASMC-SGR
## Description
Adaptive Sequential Monte Carlo inversion combined with Sequential Geostatistical Resampling for posterior inference and evidence computation

These codes correspond to the article by Amaya et al. (2022). It is a Python 3.7 implementation of the Adaptive Sequential Monte Carlo (ASMC) method (Zhou et al., 2016; algorithm 4) to estimate 
the posterior probability density function (PDF) and the evidence (marginal likelihood) trough a Bayesian inversion. ASMC is a particle approach that relies on importance sampling over a sequence 
of intermediate distributions (power posteriors) that link the prior and the posterior PDF. Each power posterior is approximated by updating the particle importance weights and states using a small 
pre-defined number MCMC proposal steps. ASMC method adaptively tunes the sequence of power posteriors and performs resampling of particles when the variance of their importance weights 
becomes too large.

The test case is a synthetic groundwater transport problem from Laloy et al. (2016) with channelized categorical models representing the hydraulic conductivity spatial distribution. 

This ASMC implementation (referred to as ASMC-SGR, algorithm 1 in Amaya et al. (2022)) includes:

-an adaptation of the Sequential Geostatistical Resampling algorithm from Laloy et al. (2016) to generate new model proposals troughout the inversion using the direct sampling
approach (DS, Mariethoz et al., 2010). DS perfomrs multiple-point statistics (MPS) conditioned siumlations based on a training image, and it is implemented using the
DeeSse algorithm (patented by the University of Neuchâtel, http://www.randlab.org/research/deesse/),

-an adaptation of the DREAMzs algorithm to perform MCMC steps approximating each power posterior (ter Braak and Vrugt, 2008; Vrugt, 2009; Laloy and Vrugt, 2012),

-the finite-volume open-source code MaFloT for transport simulations in porous media (Kunze & Lunati, 2012).

Clarification: In this codes beta indicates the inverse temperature of the power posterior, which in Amaya et al. (2021, 2022) is indicated as alpha. 

## Files

run_asmc_sgr.py : control the user-defined parameters and run the inversion.

asmc_sgr.py : main ASMC-SGR code.

asmc_sgr_func.py : auxiliar functions called by asmc_sgr.py.

deesse: algorithm for MPS simulations (version 2016).

lic: folder to place the DeeSse license to be requested to University of Neuchâtel.

data_03.mat: tracer concentration and data noise.

ti_rotated.gslib: training image. 

run_DS.py: run DeeSse for the initial simulation (conly conditioned to cond_init.dat)

cond_init.dat: initial conditioning points corresponding to the facies at the pumping wells. '

test.unc: DeeSse parameters input file for the initial simulation. 
test.con: DeeSse parameters input file.

Sequential Geostatistical Resampling (SGR) scripts:
SetInput.py
InitializeInput.py
GenCond.py

Matlab scripts to run hydrological simulations with MaFloT:
Diffusion.m
Dispersion.m
DisplayVariable.m
Initialize.m
InputFile.m
InputFile_Original.m
MaFloT.m
Output.m
PresMat.m
runMaFloT.m
Transport.m
UpMat.m
Velocity.m



## Citation 

Amaya, M., Linde, N., Laloy, E. (2022). Advances in Water Resources (under review). 

## License

See LICENSE.txt

## Contact

Macarena Amaya (macarena.amaya@unil.ch)

## References

Amaya, M., Linde, N., Laloy, E. (2022). Advances in Water Resources (under review). 

Amaya, M., Linde, N., & Laloy, E. (2021). Adaptive sequential Monte Carlo for posterior inference and model selection among complex geological priors. Geophysical Journal International, 226(2), 1220-1238.

Künze, R., & Lunati, I. (2012). An adaptive multiscale method for density-driven instabilities. Journal of Computational Physics, 231(17), 5557-5570.

Laloy, E., & Vrugt, J. A. (2012). High‐dimensional posterior exploration of hydrologic models using multiple‐try DREAM (ZS) and high‐performance computing. Water Resources Research, 48(1).

Laloy, E., Linde, N., Jacques, D., & Mariethoz, G. (2016). Merging parallel tempering with sequential geostatistical resampling for improved posterior exploration of high-dimensional subsurface categorical fields. Advances in water resources, 90, 57-69.

Mariethoz, G., Renard, P., & Straubhaar, J. (2010). The direct sampling method to perform multiple‐point geostatistical simulations. Water Resources Research, 46(11).

Ter Braak, C. J., & Vrugt, J. A. (2008). Differential evolution Markov chain with snooker updater and fewer chains. Statistics and Computing, 18(4), 435-446.

Vrugt, J. A., ter Braak, C., Diks, C., Robinson, B. A., Hyman, J. M., & Higdon, D. (2009). Accelerating Markov chain Monte Carlo simulation by differential evolution with self-adaptive randomized subspace sampling. International Journal of Nonlin ear Sciences and Numerical Simu- lation, 10(3), 273–290.

Zhou, Y., Johansen, A. M., & Aston, J. A., (2016). Toward automatic model comparison: an adaptive sequential Monte Carlo approach, Journal of Computational and Graphical Statistics,69925(3), 701–726.



