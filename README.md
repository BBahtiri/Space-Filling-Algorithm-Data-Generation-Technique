# Data generation from a nonlinear constitutive model for training a Deep Learning model
A space-filling procedure to generate data from a constitutive model (viscoelastic-viscoplastic-damage) including moisture, strain rate, and nanoparticle volume fraction dependency. 

The driving force for the generation of loading paths in our finite deformation model is the deformation gradient F. Each Fi j (i, j = 1, . . . , 3) component leads to a different loading scenario. Therefore, we generate data using the deformation gradient in a spatiotemporal space. The sampling process starts from an undeformed configuration in which the diagonal elements of the deformation gradient are set to 1.0, and the upper/lower triangular components are set to 0.0. This ensures a bounded spatial space within a realistic viscoelastic-viscoplastic regime for the nanocomposites. To cover the ninedimensional spatial space with sufficient random points, quasi-random numbers are produced using the Halton sequence generation algorithm. To generate loading paths for training the DL model, we first produce the uniform data points for each of the nine components of the total deformation gradient. It allows capturing loading scenarios related to uniaxial tension or compression, triaxial loading, biaxial loading, and pure or simple shear loading. We utilize the uniform data points for each component of the deformation gradient to generate 5% of the training data. Secondly, we implement an algorithm to randomly visit the quasi-random data points in the nine-dimensional spatial space of the deformation gradient components. The created sequence then serves as an input to integrate the constitutive model using the Euler backward algorithm and create the training sequence. The loading paths are created using different time and deformation increments. This ensures a realistic time and deformation step within FE simulations and allows us to constrain the strain rate. 
# References

For more informations, refer to our paper:

[Paper]([https://doi.org/10.1016/j.cma.2023.116293)])

```
article{bahtiri2023machine,
  title={A machine learning-based viscoelastic--viscoplastic model for epoxy nanocomposites with moisture content},
  author={Bahtiri, Betim and Arash, Behrouz and Scheffler, Sven and Jux, Maximilian and Rolfes, Raimund},
  journal={Computer Methods in Applied Mechanics and Engineering},
  volume={415},
  pages={116293},
  year={2023},
  publisher={Elsevier}
}
```
