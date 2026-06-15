.
├── code/                              Main analysis scripts (Figures 1–5)
│   ├── figure1.R                      Data integration, batch correction, DEGs, model construction
│   ├── figure2.R                      Risk-score modelling and validation plots
│   ├── figure3.R                      Single-cell processing, subclustering, trajectory, CNV
│   ├── figure4.R                      CAF / myeloid analysis and cell–cell communication
│   └── figure5.R                      Drug-sensitivity prediction
│
├── data/                              Processed input matrices, organised by figure
│   ├── Figure1/                       Bulk expression, clinical, DEGs, train/test partitions
│   ├── Figure2/                       Inputs for prognostic-factor / nomogram analysis
│   ├── Figure3/                       Single-cell-derived inputs (CNV positions, annotation)
│   ├── Figure4/                       Microenvironment / CellChat outputs
│   └── Figure5/                       Drug-prediction inputs and outputs
│
└── External validation and CIBERSORT/ Revision analyses (external validation, immune, drug)
