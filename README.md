# HPV-immunization-for-MSM-multimodel-approach

# Summary
R code and data used to fit various models for penile-anal HPV16 transmission among MSM. The models vary in sexual contact structure and natural history of infection, and are used to provide model-averaged predictions about the effectiveness of targeted vaccination. See the corresponding publication with supplementary annex: "Potential effectiveness of prophylactic HPV immunization for men who have sex with men in the Netherlands: a multi-model approach", PLOS Medicine 2019 (doi: ...)

# Data
The data from the HIV & HPV in MSM (H2M) study that was used for modeling is provided in .csv files within [H2M-data](H2M-data). In addition, [model-data](model-data) contains the output of HPV16 transmission models, written as external representations of R objects stored in .zip archives. Download the data and unzip the R objects into a directory, from where the code provided in [R-scripts](R-scripts) can be run.

# Scripts
The code provided can be used to reproduce figures of the corresponding publication. Each script is named after a specific figure. Whether reproductions have identical appearance as the published figures may depend on the R version and on the screen resolution. The published figures were created under R version 3.5.1 (2018-07-02) -- "Feather Spray" -- on a x86_64-w64-mingw32/x64 (64-bit) platform with 1680 x 1050 resolution.
