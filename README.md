# Sequential Specification Tests to Choose a Model: A Change-Point Approach

Replication materials for "Sequential Specification Tests to Choose a Model: A
Change-Point Approach"

Forthcoming in _Communications in Statistics--Theory and Methods_

https://arxiv.org/abs/1708.00907


To replicate:

1. Download Academic Probation data from https://www.aeaweb.org/aej/app/data/2008-0202_data.zip, unzip, and save as 'data/AEJApp2008-0202_data/data_for_analysis.dta'
2. Download unemployment data from https://www.bls.gov/cps/cpsaat01.xlsx and save as .csv as 'data/unemployment2.csv'
3. In `R` run
```
rmarkdown::render('sst.Rnw')
```

Details on the packages loaded and the `R` version can be found in the file sessionInfo.txt
