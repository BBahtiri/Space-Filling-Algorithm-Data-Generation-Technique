# PHYSICS-AI Data Generation: Constitutive Model-Based Training Data

[![MATLAB](https://img.shields.io/badge/MATLAB-R2019b+-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Paper](https://img.shields.io/badge/Paper-CMA%202023-blue.svg)](https://doi.org/10.1016/j.cma.2023.116293)

A **space-filling data generation framework** for creating synthetic training datasets from nonlinear constitutive models. This repository generates comprehensive datasets for training physics-informed deep learning models on complex material behavior, specifically targeting viscoelastic-viscoplastic-damage materials with environmental dependencies.

![Data Generation Workflow](rheological_model_workflow.png)

## 🎯 Overview

This framework addresses the critical challenge of generating sufficient, diverse training data for physics-informed neural networks in materials science. Traditional experimental data is often limited in scope and expensive to obtain. Our approach generates synthetic but physically realistic datasets using established constitutive models, enabling robust training of deep learning models for material behavior prediction.

### Key Innovation

The system generates training data by:
1. **Space-Filling Design** → Systematic exploration of deformation gradient space using Halton sequences
2. **Realistic Loading Paths** → Multi-dimensional deformation scenarios including uniaxial, biaxial, triaxial, and shear loading
3. **Environmental Coupling** → Integration of temperature, moisture, and nanoparticle effects
4. **Physics Consistency** → All generated data respects fundamental thermodynamic principles

## 🚀 Features

### ✨ Core Capabilities
- **Multi-Scale Data Generation**: Covers molecular to macro-scale material responses
- **Space-Filling Algorithms**: Halton sequence-based quasi-random sampling for comprehensive coverage
- **Complex Loading Scenarios**: Uniaxial, biaxial, triaxial, and shear deformation modes
- **Environmental Sensitivity**: Temperature, moisture content, and nanoparticle volume fraction effects
- **Strain Rate Control**: Realistic time and deformation increments matching FE simulations

### 🔬 Scientific Features
- **Deformation Gradient Sampling**: 9-dimensional space exploration of F tensor components
- **Constitutive Model Integration**: Euler backward integration for robust numerical stability
- **Material Model Validation**: Comparison between traditional and ML-accelerated predictions
- **Thermodynamic Consistency**: All generated data satisfies energy conservation principles

### 🛠️ Technical Features
- **MATLAB Implementation**: Optimized for computational materials science workflows
- **Modular Architecture**: Easy adaptation to different material systems
- **Comprehensive Validation**: Built-in comparison tools for model verification
- **Batch Processing**: Efficient generation of large datasets

## 📋 Requirements

### System Requirements
- **MATLAB**: R2019b or higher
- **Toolboxes**: Statistics and Machine Learning Toolbox (recommended)
- **Memory**: Minimum 8GB RAM, 16GB+ recommended for large datasets
- **Storage**: 10GB+ free space for generated datasets

### Dependencies
```matlab
% Core MATLAB toolboxes
Statistics and Machine Learning Toolbox
Parallel Computing Toolbox (optional, for acceleration)
```

## 🔧 Installation

### Option 1: Clone Repository
```bash
# Clone the repository
git clone [https://github.com/yourusername/PHYSICS-AI-DataGen.git](https://github.com/BBahtiri/Space-Filling-Algorithm-Data-Generation-Technique)
cd PHYSICS-AI-DataGen

# Open MATLAB and add to path
addpath(genpath(pwd))
```

### Option 2: Direct Download
1. Download the repository as ZIP
2. Extract to your desired location
3. Add the folder to your MATLAB path

## 📁 Project Structure

```
PHYSICS-AI-DataGen/
├── 📄 Main_rheological_model_random.m  # Main data generation script
├── 🔧 lstm_forward.m                   # LSTM forward propagation utilities
├── 🔧 maeRegressionLayer.m             # Custom regression layer for validation
├── 📁 lib/                             # Supporting library functions
├── 📊 simulations_const_random_*/       # Generated dataset directories
├── 📁 validation/                      # Model validation results
├── 📄 requirements.txt                 # MATLAB toolbox requirements
├── 📄 README.md                        # This file
└── 🖼️ rheological_model_workflow.png   # Workflow diagram
```

## 🏃‍♂️ Quick Start

### 1. Basic Data Generation
```matlab
% Open MATLAB and navigate to the repository directory
cd('path/to/PHYSICS-AI-DataGen')

% Run the main data generation script
Main_rheological_model_random
```

### 2. Configure Generation Parameters
Edit the key parameters in `Main_rheological_model_random.m`:

```matlab
% Material Properties
wnp = 0.0;              % Nanoparticle weight fraction (0-0.2)
wgf = 0.0;              % Glass fiber weight fraction (0-0.6)
Temper = 296;           % Temperature (K)
zita = 0.0;             % Moisture content (0-1)

% Generation Parameters
peaks = 5;              % Number of deformation peaks per path
sims = 1000;            % Number of simulation paths to generate

% Deformation Control
et = 1e-3;              % Strain rate (1/s)
de1 = 1e-5;             % Primary deformation increment
de2 = 1e-6;             % Secondary deformation increment
```

### 3. Run Validation (Optional)
```matlab
% Set validation mode
validate = 1;           % Enable validation mode

% Load trained ML model for comparison
load('lstm_3_seqtoseq_2150_64.mat');

% Run validation simulations
Main_rheological_model_random
```

## 📊 Data Generation Process

### Space-Filling Design

The framework uses a sophisticated approach to sample the 9-dimensional deformation gradient space:

1. **Halton Sequence Generation**:
```matlab
% Generate quasi-random points
p = haltonset(9);
p = scramble(p,'RR2');
X0 = p(1:1000,:);

% Normalize to physical ranges
M1 = normalize(X0(:,1),'range',[0.9,1.1]);  % F11 diagonal terms
M4 = normalize(X0(:,4),'range',[-0.05,0.05]); % F12 off-diagonal terms
```

2. **Deformation Path Creation**:
```matlab
% Generate loading paths between random states
for i = 1:peaks
    prev_peak = lpd(i,j-1);
    current_peak = states(j-1,i);
    
    if current_peak > prev_peak
        current_strain = prev_peak:de:current_peak;
    else
        current_strain = prev_peak:-de:current_peak;
    end
end
```

### Loading Scenarios

The system generates diverse loading conditions by varying deformation gradient components:

| Combination | Loading Type | F Matrix Structure |
|-------------|--------------|-------------------|
| `1` | Uniaxial X | `[1+ε 0 0; 0 1 0; 0 0 1]` |
| `123` | Triaxial | `[1+ε 0 0; 0 1+ε 0; 0 0 1+ε]` |
| `9` | General Shear | `[1+ε γ γ; γ 1+ε γ; γ γ 1+ε]` |
| `112` | Uniaxial + Shear | `[1+ε γ 0; 0 1 0; 0 0 1]` |

### Material Response Calculation

For each deformation path, the system:

1. **Integrates Constitutive Model**:
```matlab
[stress,F,F_ve,F_v,F_e,F_vp,J_final,t_tot,sigma,sigma_inf,sigma_t,sigma_v,d,...
 strain,tau,trueStrain,b_almansi,C_pert,rdn] = model_vevpd_fmix(params,...
 de_mono,d0,dt_mono,vnp,params_F,nnetwork,0,et_load,et_unload,...
 trueStrain11,trueStrain22,trueStrain33,trueStrain12,trueStrain13,...
 trueStrain21,trueStrain23,trueStrain31,trueStrain32,0);
```

2. **Extracts Training Features**:
```matlab
% Create input-output pairs for ML training
inputs{1,k} = [current_b; dt; wnp; zita; Temper];
outputs{1,k} = [current_stress_1];
```

## ⚙️ Advanced Usage

### Custom Material Systems

To adapt for different materials, modify the constitutive model parameters:

```matlab
% Material parameters (example for different polymer)
mu1 = 800;          % First shear modulus
mu2 = 850;          % Second shear modulus  
nu1 = 0.25;         % First Poisson's ratio
m = 0.7;            % Strain rate sensitivity
gamma_dot_0 = 1e12; % Reference strain rate
dG = 2e-19;         % Activation energy
```

### Environmental Conditions

Systematically vary environmental parameters:

```matlab
% Temperature sweep
temps = [273, 296, 333, 363];  % Kelvin

% Moisture content variation
zita_values = [0, 0.005, 0.011, 0.02];  % Fractional moisture

% Nanoparticle content
wnp_values = [0, 0.05, 0.1, 0.15];  % Weight fraction
```

### Custom Deformation Paths

Define specific loading scenarios:

```matlab
% Custom deformation gradient sequence
custom_F11 = [1.0, 1.02, 1.05, 1.03, 1.0];  % Loading-unloading
custom_F22 = ones(size(custom_F11));         % Constrained
custom_F33 = ones(size(custom_F11));

% Generate corresponding stress-strain data
[stress_custom, ~] = generate_custom_path(custom_F11, custom_F22, custom_F33);
```

### Parallel Processing

For large dataset generation:

```matlab
% Enable parallel processing
parpool('local', 4);  % Use 4 cores

% Parallel data generation
parfor sim_id = 1:num_simulations
    [stress_data{sim_id}, strain_data{sim_id}] = generate_single_path(sim_id);
end
```

## 📈 Output Data Structure

### Generated Files

The framework creates structured output directories:

```
simulations_const_random_5/
├── input_target_ml_et_1_1_1_1.mat    # Input-output pairs
├── space.fig                          # Deformation space visualization
├── load_path.fig                      # Loading path plots
└── results.fig                        # Stress-strain results
```

### Data Format

Each `.mat` file contains:

```matlab
% Input features (x): [n_features × n_timesteps]
% - Row 1-6: Symmetric stress tensor components (S11, S12, S13, S22, S23, S33)
% - Row 7: Time increment (dt)
% - Row 8: Nanoparticle weight fraction (wnp)
% - Row 9: Moisture content (zita)
% - Row 10: Temperature (Temper)

% Output targets (t): [n_stress_components × n_timesteps]
% - Stress tensor components corresponding to input deformation
```

### Validation Results

When validation mode is enabled:

```
validation/
├── sim_1.fig                          # Individual simulation comparison
├── sim_2.fig
├── ...
├── dist_rmse.fig                      # RMSE distribution
└── dist_mae.fig                       # MAE distribution with statistics
```

## 🔬 Scientific Background

### Theoretical Foundation

The data generation is based on a comprehensive constitutive model incorporating:

1. **Hyperelasticity**: Large deformation elastic response
2. **Viscoelasticity**: Time-dependent material behavior
3. **Viscoplasticity**: Rate-dependent permanent deformation
4. **Damage Mechanics**: Progressive material degradation
5. **Environmental Coupling**: Temperature and moisture effects

### Material Model Equations

The framework implements a sophisticated material model:

```
σ = σ_∞ + σ_v + σ_d
```

Where:
- `σ_∞`: Equilibrium stress
- `σ_v`: Viscous overstress  
- `σ_d`: Damage-related stress

### Key Advantages

- **Physics Consistency**: All generated data respects thermodynamic principles
- **Comprehensive Coverage**: Systematic exploration of material response space
- **Computational Efficiency**: Orders of magnitude faster than full FE simulations
- **Validation Ready**: Built-in comparison with ML model predictions

## 📚 Applications

### Primary Use Cases

1. **ML Model Training**: Generate large datasets for physics-informed neural networks
2. **Material Characterization**: Explore material response under diverse conditions
3. **Model Validation**: Compare constitutive models with experimental data
4. **Design Optimization**: Screen material compositions and processing conditions

### Integration with ML Frameworks

The generated data is designed for seamless integration with:

- **TensorFlow/Keras**: Direct loading of `.mat` files
- **PyTorch**: Convert using `scipy.io.loadmat`
- **MATLAB Deep Learning Toolbox**: Native compatibility


### Areas for Contribution
- **New Material Models**: Implement additional constitutive laws
- **Enhanced Sampling**: Improve space-filling algorithms
- **Validation Tools**: Expand model comparison capabilities
- **Documentation**: Tutorials and example workflows
- **Performance**: Optimization for large-scale generation

### Development Setup
```bash
# Fork and clone your fork
git clone https://github.com/yourusername/PHYSICS-AI-DataGen.git
cd PHYSICS-AI-DataGen

# Create feature branch
git checkout -b feature/new-material-model

# Make changes and test
% Run validation suite in MATLAB
run_validation_tests

# Submit pull request
```

## 🐛 Troubleshooting

### Common Issues

**Memory Errors**
```matlab
% Reduce simulation parameters
peaks = 3;          % Reduce from 5
sims = 500;         % Reduce from 1000
```

**Convergence Problems**
```matlab
% Adjust numerical parameters
de1 = 5e-6;         % Smaller deformation increments
dt_max = 1.0;       % Larger time steps
```

**Invalid Stress Values**
```matlab
% Check stress magnitude filter
if any(t > abs(400), 'all')
    % Skip simulation
    continue
end
```

### Performance Optimization

**For Large Datasets**:
```matlab
% Use parallel processing
parpool('local', maxNumCompThreads);

% Optimize memory usage
clear unnecessary_variables;
pack;  % Defragment memory
```

**For Speed**:
```matlab
% Reduce figure generation
set(0,'DefaultFigureVisible','off');

% Use faster integration schemes
% (modify model_vevpd_fmix parameters)
```


## 📚 Citation

If you use this code in your research, please cite:

```bibtex
@article{bahtiri2023machine,
  title={A machine learning-based viscoelastic--viscoplastic model for epoxy nanocomposites with moisture content},
  author={Bahtiri, Betim and Arash, Behrouz and Scheffler, Sven and Jux, Maximilian and Rolfes, Raimund},
  journal={Computer Methods in Applied Mechanics and Engineering},
  volume={415},
  pages={116293},
  year={2023},
  publisher={Elsevier},
  doi={https://doi.org/10.1016/j.cma.2023.116293}
}
```



---

**⭐ Star this repository if you find it useful!**

