# Recon3D

A 3-Dimensional View of Human Metabolism and Disease

****************************************************

## GEM-PRO (Recon3D_GP)

GEM-PRO for the *H. sapiens* Recon3D metabolic model.

### Instructions

**Requirements**
- Jupyter notebook installation and Python 2/3
- [ssbio](https://github.com/SBRG/ssbio)

**Installing nglview and ssbio**
```bash
pip install nglview
pip install ssbio
```

**Cloning**
```bash
git clone https://github.com/SBRG/Recon3D
```

**Opening**
- Launch the notebook `Recon3D - Loading and Exploring the GEM-PRO.ipynb` from a Jupyter instance

**Updating**
- Launch and run the notebook `Recon3D - Updating the GEM-PRO.ipynb` from a Jupyter instance

### Citation

TBD


## Single gene deletion simulations (sgd)

The description of these simulations is made available through the Supplementary Information.

### Instructions
 
**Requirements**
- MATLAB R2015B+
- Mosek v.7+
- [RAVEN v.1.08+](https://github.com/SysBioChalmers/RAVEN)
- R 3.4.0 and R libraries gplots (3.0.1) and RColorBrewer (1.1-2)

**Workflow**
- Set /sgd as working directory 
- Convert Recond3D model into RAVEN-format and save it to models/Recon3d_toRaven.mat) by running convertRecon3toRaven.m
- Simulate in silico library of sgds for HMR2, Recon3d, and HMR2-derived GBM models under Ham's medium composition (described in media/) by running main.m. Medium-constrained wild-type models are saved in models/constrained/.
- Analyze in silico sgd results by running main.R. Plots are saved in /plots