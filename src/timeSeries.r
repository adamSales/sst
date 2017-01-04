require(vars)
#### do the VAR analysis
#### modeled after Bernanke and Gertler 1995 "Inside the Black Box..."
#### Fig 1 in B&G
#### data downloaded 1/2/17 from FRED

### assemble data
gdp <- read.csv('data/GDP.csv')
gdpdef <- read.csv('data/GDPDEF.csv')
priceIndex <- read.csv('data/PPIACO.csv')
fedfund <- read.csv('data/FEDFUNDS.csv')

tsdat <- merge(gdp,gdpdef)
tsdat <- merge(tsdat,priceIndex)
names(tsdat)[4] <- 'priceIndex'
tsdat <- merge(tsdat,fedfund)
names(tsdat)[5] <- 'fedfund'

mod <- VAR(tsdat[,-1],100)