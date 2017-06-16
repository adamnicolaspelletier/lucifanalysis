# lucifanalysis.py Version 2.0
Script to analyze high-throughput luciferase experiments in an automated fashion

### Author: Adam-Nicolas Pelletier
### Github: https://github.com/adamnicolaspelletier/lucifanalysis.git

## Table of Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Information](#information)
- [Usage](#usage)
- [In Development](#in-development)
- [Known Bugs](#known-bugs)
- [Support](#support)
- [Contributing](#contributing)



## Requirements

 * Python 2.7+
 * Have MAFFT command line tools installed . Can be obtained at http://mafft.cbrc.jp/alignment/software/

 * Python Libraries :

	Pandas (+ Numpy) 0.20.1+
	BioPython 1.69+
	Matplotlib 2.0.2
	outlier_utils 0.0.3+
	
 * Luciferase Data. Example data [is included](lucifanalysis/docs/Input/template_AP170428.txt)


## Installation

Download to your project directory from GitHub, or clone it on the command line using the following command:

```sh
git clone git://github.com/adamnicolaspelletier/lucifanalysis.git

```

## Information

This script was made to facilitate analysis of high-throughput dual-luciferase experimental data. In such experiments, the Firefly luciferase signal is driven by a transcription factor, measured in a given number of replicates, and compared to other experimental conditions to derive functional differences. This signal is normalized on Renilla Luciferase, which is meant to be a constant accross conditions, allowing to compensate for variations in transfection efficiency, cell number or lysis efficiency. 

The ratio of Firefly/Renilla gives Relative Luciferase Units (RLU), a standardized measurment of the luciferase signal for each sample. Computing Fold changes is usually the more "human brain" friendly way to detect variations: one condition is set as the control against which other conditions will be compared to __(RLU/Control_RLU)__. Thus, the control has 100% Fold Change, or 1, and experiemental conditions can either be similar, higher or lower than the control. 

Experiments are usually done in the 96 well plate format: and biologists analyze these using simple template spreadsheets. 
**The advantage of this script is its scalability**: template spreadsheets adapt very poorly to divergent patterns, such as variations in number of replicates, or different positions of controls, etc. Figure generation inherits the same problem, requiring the biologist to manually assign values for **__EACH EXPERIMENT__**. This is time- consuming AND error-prone in a high-throughput setting. Finally, spreadsheets struggle to generate a standardized report for each experiment in a non 96 well plate format, making it harder to integrate all data from a collection of experiments together in databases. *This scripts does it all for you, and strives to do more* (see [In Development](#in development)). 

**IN SUMMARY**: while this script can definately perform extremely well for luciferase analysis of a single experiment, it really shines when multiple experiments need to be analyzed. 


## Usage

Download the [template file](lucifanalysis/docs/Input/template_AP170428.txt) to fill in luciferase experimental values, as well as the plate plan and legend. It can be easily edited in any spreadsheet software, such as Microsoft Excel or Open Office. 

The script uses optional flags for personalized analysis: input file (-lu) followed by path to luciferase results and plan will analyze this file instead of the default template file. Accepts both relative and absolute paths

```sh
python lucifanalysis.py -lu path/to/luciferase/file.txt  ## Location based on current working directory
python lucifanalysis.py -lu ~/absolute/path/to/luciferase/file.txt  ## Absolute path, current working directory is irrelevant 

```

Grubbs tests for detection of outliers can be activated with the -g flag, eliminating any outlier replicates automatically.
**WARNING**:Grubbs test is not known to be reliable in small n, meaning it is not recommended. More practical way is being implemented in future release. 

```sh
python lucifanalysis.py -g

```

A PDF report can also be generated from the text report using the -p flag, which has figures. 

```sh
python lucifanalysis.py -p

```

Use the -h flag for further options and documentation



## In Development

- Add an outlier detection function. 
- Remove Grubbs test
- Add support to add reports to a MySQL database
- Add parameter for output PATH instead of in file. 
- Add support for personalized figure layout instead of automatic


## Known Bugs

IPYTHON has problems with the argparse module, other Python distributions do not. I suggest using Official Python in the meantime. 
Script has not been tested in Python 3.


## Support

Please [open an issue](https://github.com/adamnicolaspelletier/lucifanalysis.git/issues/new) for support.


## Contributing

Please contribute using [Github Flow](https://guides.github.com/introduction/flow/). Create a branch, add commits, and [open a pull request](https://github.com/adamnicolaspelletier/lucifanalysis/compare/).
