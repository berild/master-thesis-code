<h2>Integrated Nested Laplace Approximations within Monte Carlo methods</h2>

------

* Master Thesis delivered June 27, 2020
* Title: Integrated Nested Laplace Approximations within Monte Carlo methods
* Author: Martin Outzen Berild

------

This repository contains the code used in the master thesis. The different methodologies was applied to four examples; biviariate linear model, spatial autoregressive combined model, model-aware Bayesian parametric quantile regression, and a gamma frailty model. The work done here is a continuation of my project thesis, and the related code can be found <a href="https://github.com/berild/project-thesis">here</a>.

<h3 id="abstract">Abstract</h3>

The Integrated Nested Laplace Approximation (INLA) is a deterministic approach to Bayesian inference on latent Gaussian models (LGMs) and focuses on fast and accurate approximation of posterior marginals for the parameters in the models. In practice, applications of INLA are limited to the class of models implemented in the R package R-INLA. Recently, methods have been developed to extend this class of models to those that can be expressed as conditional LGMs by fixing some of the parameters in the models to descriptive values. These methods differ in the manner descriptive values are chosen. This thesis considers the three following INLA within Monte Carlo methods: Markov chain Monte Carlo (MCMC) with INLA, importance sampling (IS) with INLA, and a novel approach that combines INLA and the adaptive multiple importance sampling (AMIS) algorithm. This thesis compares the INLA within Monte Carlo methods on a series of applications with simulated and observed datasets and evaluates their performance based on accuracy, efficiency, and robustness. The implementation of the methods are validated by exact posteriors in a simple bivariate linear model and tested on a spatial autoregressive combined model. Then, it presents a new approach to Bayesian quantile regression using AMIS with INLA, which is verified in a simulation study and applied to two observed datasets. Also, this thesis attempts to approximate the posteriors in a Gamma frailty model using AMIS with INLA. The examples show that the AMIS with INLA approach, in general, outperformed the other methods on more complex models, but the IS with INLA algorithm could be considered for faster inference when good proposals are available. Also, the Bayesian quantile regression approach produced promising quantile curves in the simulation study, and the applications present a small portion of the large class of models that are facilitated through INLA for this type of quantile regression. In addition, the AMIS with INLA algorithm produced accurate posteriors in the Gamma frailty model with few clusters but attained a slight bias for a higher number of dimensions in the AMIS algorithm.


<h3 id="examples">Examples</h3>

<details>
  <summary id="blm" style ="cursor: pointer; font-size: 1.5em;">Bivariate linear model (click to view)</summary>
  
To apply the combined methods on the bivariate linear model, run the <a href="https://github.com/berild/master-thesis-code/blob/master/toy/toy.R">toy/toy.R</a> script. The functions for each algorithm is presented in <a href="https://github.com/berild/master-thesis-code/blob/master/toy/amis_w_inla.R">toy/amis_w_inla.R</a>, <a href="https://github.com/berild/master-thesis-code/blob/master/toy/is_w_inla.R">toy/is_w_inla.R</a>, and <a href="https://github.com/berild/master-thesis-code/blob/master/toy/mcmc_w_inla.R">toy/mcmc_w_inla.R</a>. General functions used in the algorithms, plotting, and evaluation is given in <a href="https://github.com/berild/master-thesis-code/blob/master/toy/general_functions.R">toy/general_functions.R</a>. Result of our simulation is given with <a href="https://git-lfs.github.com/">git-lfs</a> in <a href="https://github.com/berild/master-thesis-code/tree/master/sims/toy">sims</a> and use <a href="https://github.com/berild/master-thesis-code/blob/master/toy/plot_toy.R">plot_toy.R</a> to replicate our plots.

  <h4>Result</h4>
  <h5>Univariate posterior marginals </h5>
  <img src="https://imgur.com/jkNOmGi.png"
       alt="univariate posterior marginals bivariate linear model"
       style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
  <h5>Bivariate posterior marginals</h5>
  <img src="https://imgur.com/L3M5qkU.png"
       alt="bivariate posterior marginals bivariate linear model"
       style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
  <h5>Adaptation with importance sampling algorithms</h5>
  <img src="https://i.imgur.com/lF5zuHX.png"
       alt="AMIS w/ INLA and IS w/ INLA adaptation"
       style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
  <h5>Trace plots MCMC</h5>
  <img src="https://imgur.com/bDNpEz0.png"
       alt="MCMC w/ INLA trace plot"
       style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
  </details>


<details>
  <summary id="bl" style ="cursor: pointer; font-size: 1.5em;">Spatial autoregressive combined model (click to view)</summary>
  
To apply the combined methods on the Spatial autoregressive combined model, run the <a href="https://github.com/berild/master-thesis-code/blob/master/sem/sem.R">sem/sem.R</a> script. The functions for each algorithm is presented in <a href="https://github.com/berild/master-thesis-code/blob/master/sem/amis_w_inla.R">sem/amis_w_inla.R</a>, <a href="https://github.com/berild/master-thesis-code/blob/master/sem/is_w_inla.R">sem/is_w_inla.R</a>, and <a href="https://github.com/berild/master-thesis-code/blob/master/sem/mcmc_w_inla.R">sem/mcmc_w_inla.R</a>. General functions used in the algorithms, plotting, and evaluation is given in <a href="https://github.com/berild/master-thesis-code/blob/master/sem/general_functions.R">sem/general_functions.R</a>. Result of our simulation is given with <a href="https://git-lfs.github.com/">git-lfs</a> in <a href="https://github.com/berild/master-thesis-code/tree/master/sims/sem">sims/sem</a> and use <a href="https://github.com/berild/master-thesis-code/blob/master/sem/plot_sem.R">sem/sem/plot_sem.R</a> to replicate our plots.

<h5>Election turnover in Italy 2001</h5>
<img src="https://imgur.com/1gbfjKI.png"
     alt="Election turnover Italy 2001"
     style="width: 40%; display: block; margin-left: auto; margin-right: auto;" /> 

<h5>GDP per capita Italy 1997</h5>
<img src="https://imgur.com/HsuMk6o.png"
     alt="GDP per capita Italy 1997"
     style="width: 40%; display: block; margin-left: auto; margin-right: auto;" />
     
<h4>Result</h4>
  
  
| Parameter |    MCMC    |  IS w/INLA | AMIS w/INLA|
|:---------:|:----------:|:----------:|------------|
| Intercept | 5.76(2.34) | 6.17(2.46) | 6.11(2.42) |
| GDPPCAP   | 1.75(0.59) | 1.84(0.61) | 1.83(0.61) |
|   &rho;   | 0.86(0.04) | 0.84(0.07) | 0.84(0.07) |
|  &lambda; | 0.21(0.11) | 0.25(0.13) | 0.24(0.13) |
|    &tau;  | 0.26(0.02) | 0.26(0.02) | 0.26(0.02) |

<h5>Posterior marginals of parameters in SAC model</h5>
<img src="https://imgur.com/iy7XWiR.png"
     alt="SAC posterior marginals IS w/ INLA, AMIS w/ INLA and MCMC"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
     
<img src="https://imgur.com/KYjqWMO.png"
     alt="SAC posterior marginals MCMC w/ INLA and AMIS w/ INLA"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 

</details>


<details>
  <summary id="mc" style ="cursor: pointer; font-size: 1.5em;">Model-aware Bayesian parametric quantile regression (click to view)</summary>
    
 To apply the combined methods for Bayesian parametric quantile regression, run the <a href="https://github.com/berild/master-thesis-code/blob/master/PQR/pqr.R">PQR/pqr.R</a> script. The functions for each algorithm is presented in <a href="https://github.com/berild/master-thesis-code/blob/master/PQR/amis_w_inla.R">PQR/amis_w_inla.R</a>, <a href="https://github.com/berild/master-thesis-code/blob/master/PQR/is_w_inla.R">PQR/is_w_inla.R</a>, and <a href="https://github.com/berild/master-thesis-code/blob/master/PQR/mcmc_w_inla.R">PQR/mcmc_w_inla.R</a>. General functions used in the algorithms, plotting, and evaluation is given in <a href="https://github.com/berild/master-thesis-code/blob/master/PQR/general_functions.R">PQR/general_functions.R</a>. Result of our simulation is given with <a href="https://git-lfs.github.com/">git-lfs</a> in <a href="https://github.com/berild/master-thesis-code/tree/master/sims/pqr">sims/pqr</a> and use <a href="https://github.com/berild/master-thesis-code/blob/master/PQR/plot_pqr.R">PQR/plot_pqr.R</a> to replicate our plots.   
  
<h4>Result</h4>

<h5> Simulated datasets Bayesian PQR </h5>
<img src="https://imgur.com/BoBSs3s.png"
     alt="Simulated dataset Bayesian PQR"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" />
<h5> RW2 model Bayesian PQR</h5>     
<img src="https://imgur.com/Tqy3wa5.png"
     alt="RW2 model Bayesian PQR"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
<h5>Gamma Bayesian PQR</h5>
<img src="https://imgur.com/84qFSfP.png"
     alt="Gamma Bayesian PQR"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
</details>


<details>
  <summary id="sem" style ="cursor: pointer; font-size: 1.5em;">Gamma Frailty model (click to view)</summary>
  
   To apply the AMIS with INLA method on Gamma frailty models, run the <a href="https://github.com/berild/master-thesis-code/blob/master/survival/frailty.R">survival/frailty.R</a> script. The functions in the AMIS w/ INLA algorithm is presented in <a href="https://github.com/berild/project-thesis/blob/master/sem/sem_amis_w_inla.R">survival/frailty_amis_w_inla.R</a>. href="https://github.com/berild/master-thesis-code/blob/master/survival/frailty_general_functions.R">survival/frailty_general_functions.R</a>. Result of our simulation is given with <a href="https://git-lfs.github.com/">git-lfs</a> in <a href="https://github.com/berild/master-thesis-code/tree/master/sims/frailty">sims/frailty</a> and use <a href="https://github.com/berild/master-thesis-code/blob/master/survival/plot_frailty.R">survival/plot_frailty.R</a> to replicate our plots.   

<h4>Result</h4>
<h5>Posterior marginals Gamma frailty model 4 clusters</h5>
<img src="https://imgur.com/VE04gLY.png"
     alt="Spatial Econometric Model"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
<h5>Posterior mean and the 0.025 and 0.975 quantiles of log frailty with 4 clusters</h5>     
<img src="https://imgur.com/fN56z1d.png"
     alt="Spatial Econometric Model"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
<h5>Posterior marginals Gamma frailty model 20 clusters</h5>
<img src="https://imgur.com/4s7hcLN.png"
     alt="Spatial Econometric Model"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
<h5>Posterior mean and the 0.025 and 0.975 quantiles of log frailty with 20 clusters</h5>  
<img src="https://imgur.com/8wtbwXc.png"
     alt="Spatial Econometric Model"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
</details>

<h3> Code Explanation </h3>
... In progress...
<details>
  <summary id="blm" style ="cursor: pointer; font-size: 1.5em;">AMIS with INLA (click to view)</summary>
  
</details>

<details>
  <summary id="blm" style ="cursor: pointer; font-size: 1.5em;">IS with INLA (click to view)</summary>
  
</details>

<details>
  <summary id="blm" style ="cursor: pointer; font-size: 1.5em;">MCMC with INLA (click to view)</summary>
  
</details>

<details>
  <summary id="blm" style ="cursor: pointer; font-size: 1.5em;">General (click to view)</summary>
  
</details>