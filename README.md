# Description
These files are to be used with R and did run in version 4.3.0. 

## Data files (sources are described in the paper)
The file 'data_raw.rda' includes the GDP series and financial stress indicator (CISS). These series were downloaded from the OECD MEI database and the ECB statistical data warehouse.

## Code files
The individual R-source files are the following:

- 'designmat.R' produces the design matrices depending on model specification and forecast horizon
- 'qf-bart.R' contains the main estimation function
- 'forecast-ex.R' sets up a grid for all required estimations, defines algorithm settings and was run in a parallel computing environment
- 'utils.R' contains several helper functions

Some auxiliary functions may have been copied from other respositories. This code comes without technical support of any kind. Please report any typos
or errors to: mpfarrho@gmail.com. The code is free to use, provided that the paper is cited properly.
