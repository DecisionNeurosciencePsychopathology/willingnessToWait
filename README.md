# Willingness to wait

Code to analyze the willingness to wait task administrated by the dnpl lab

## Getting Started

You will need to clone this repo and have willingness to wait data

### Prerequisites

MATLAB, R and VBA toolbox (written in MATLAB)

## Main scripts

wtw_analysis - main parent level script to loop over different models & subejcts in VBA analysis
wtw_sceptic_vba - Main model function that will run a specific flavor of the sceptic model using the VBA toolbox

Evolution functions  - functions with prefix 'h'

```
h_wtwsceptic_fixed.m
```

Observation functions - functions with prefix 'g'

```
g_wtwsceptic.m
```

diagnose_mux.m - will make a verbose graph of the hidden states and reward rate for a single subject