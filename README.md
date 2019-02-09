# HPV-immunization-for-MSM-multimodel-approach

# Summary
R code and data used to fit various models for penile-anal HPV16 transmission among MSM. The models vary in sexual contact structure and natural history of infection, and are used to provide model-averaged predictions about the effectiveness of targeted vaccination. See the corresponding publication with supplementary annex: "Potential effectiveness of prophylactic HPV immunization for men who have sex with men in the Netherlands: a multi-model approach", PLOS Medicine 2019 (doi: ...)

# Data
The data from the HIV & HPV in MSM (H2M) study that was used for modeling is provided in .csv files within [H2M-data](H2M-data). In addition, [model-data](model-data) contains the output of HPV16 transmission models, written as external representations of R objects stored in .zip archives. Download the data and unzip the R objects into a directory, from where the code provided in [R-scripts](R-scripts) can be run.

# Model
R code with corresponding .dll is given for the SISPRS model in [SISPRS-model](SISPRS-model). The model can be run from within the directory where the data are unzipped. The output from this particular model is also provided, but need not be loaded for running the model. The code is given as an example to illustrate the structure of HPV16 transmission models and their output, source code is not provided. The .dll was compiled under R x64 3.4.1 (2017-06-30), performance under other R versions is not guaranteed.

# Scripts
The scripts provided can be used to reproduce figures of the corresponding publication. Each script is named after a specific figure. Whether reproductions have identical appearance as the published figures may depend on the R version and associated packages, and on the screen resolution. The published figures were created under R x64 3.5.1 (2018-07-02) with 1680 x 1050 screen resolution.
