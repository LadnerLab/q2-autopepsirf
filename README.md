# q2-autopepsirf

Qiime2 plugin used for the automation of q2-pepsirf and q2-ps-plot.

https://ladnerlab.github.io/pepsirf-q2-plugin-docs/

## Installation

### Qiime2 Installation:

Visit: https://docs.qiime2.org/2021.8/install/ for intallation documentation on Qiime2

### Pepsirf Installation:

Visit: https://github.com/LadnerLab/PepSIRF for installation documentation on PepSIRF

### q2-pepsirf1 Installation:

Visit: https://github.com/LadnerLab/q2-pepsirf1 for installation documentation on q2-pepsirf1

### q2-ps-plot Installation:

Visit: https://github.com/LadnerLab/q2-ps-plot for installation documentation on q2-ps-plot

### q2-pepsirf installation:
#### Dependencies:
- `qiime2-2021.11 +`
- `PepSIRF`
- `q2-pepsirf1`
- `q2-ps-plot`

#### Directions:
Make sure your Qiime2 conda environment is activated by running the command:
  
```
conda activate qiime2-2021.11
```

You can replace `qiime2-2021.11` above with whichever version of QIIME 2 you have currently installed.

Now you are ready to install q2-autopepsirf. Run the following commands:

```
pip install git+https://github.com/LadnerLab/q2-autopepsirf.git
```

Run `qiime info` to check for a successful installation. If installation was successful, you should see `autopepsirf: version` in the list of installed plugins.
