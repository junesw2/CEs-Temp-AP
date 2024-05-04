########################################################
### This Code Includes a WARNING, See "README" Below ###
########################################################

### TO REMOVE ###



################################################################################
###
### Trends in Population Exposure to Compound Extreme-risk Temperature and Air Pollution in Europe
### February 2024
### Joan Ballester, Hicham Achebak, Zhao Yue chen
### 
################################################################################

rm( list = ls() );
cat("\014");
setwd("/PROJECTES/ADAPTATION/proj/zchen/P20230815_ZC_heat_AP_CEs/20230612_joan_final_version_figures/")
# Required Libraries
suppressMessages( library(lubridate) ); # Required Functions: wday
suppressMessages( library(tidyverse) ); # Required Functions: read_delim
suppressMessages( library(plyr) ); # Required Functions: join
suppressMessages( library(ISOweek) ); # Required Functions: ISOweek2date
suppressMessages( library(dlnm) ); # Required Functions: crossbasis
suppressMessages( library(splines) ); # Required Functions: ns, bs
suppressMessages( library(mixmeta) ); # Required Functions: mixmeta
suppressMessages( library(DistributionUtils) ); # Required Functions: is.wholenumber
suppressMessages( library(tsModel) ); # Required Functions: Lag
suppressMessages( library(MASS) ); # Required Functions: mvrnorm
suppressMessages( library(giscoR) ); # Required Functions: gisco_get_countries
suppressMessages( library(sf) ); # Required Functions: st_read
suppressMessages( library(grid) ); # Required Functions: viewport
#suppressMessages( library(writexl) ); # Required Functions: write_xlsx
suppressMessages( library(scales) ); # Required Functions: alpha

# Required Codes
source("QAIC.R"); # Akaike's Information Criterion for Overdispersed Count Data
source("mean_annual_cycle_temp_mort.R"); # Time Series of the Mean Annual Cycle of a Weekly Time Series of Temperature or Mortality Counts
source("FWALD.R"); # Wald Test
# source("map_generator.R"); # Map Drawer
source("map_generator_nature_medicine.R"); # Map Drawer



################################################################################
### Definition of the Sensitivity Parameters
################################################################################

# Boolean of Factual (TRUE) or Counterfactual (FALSE) Temperatures
bFACTUAL =  TRUE; sFACTUAL =        "FACTUAL";
# bFACTUAL = FALSE; sFACTUAL = "COUNTERFACTUAL";

# Sex and Age Groups
iSEXAGE = 11;
if      ( iSEXAGE == 11 ){ sSEX = "allsex"; sAGE =  "allage";
}else if( iSEXAGE == 21 ){ sSEX =  "women"; sAGE =  "allage";
}else if( iSEXAGE == 31 ){ sSEX =    "men"; sAGE =  "allage";
}else if( iSEXAGE == 12 ){ sSEX = "allsex"; sAGE =  "age80p";
}else if( iSEXAGE == 22 ){ sSEX =  "women"; sAGE =  "age80p";
}else if( iSEXAGE == 32 ){ sSEX =    "men"; sAGE =  "age80p";
}else if( iSEXAGE == 13 ){ sSEX = "allsex"; sAGE = "age6579";
}else if( iSEXAGE == 23 ){ sSEX =  "women"; sAGE = "age6579";
}else if( iSEXAGE == 33 ){ sSEX =    "men"; sAGE = "age6579";
}else if( iSEXAGE == 14 ){ sSEX = "allsex"; sAGE = "age0064";
}else if( iSEXAGE == 24 ){ sSEX =  "women"; sAGE = "age0064";
}else if( iSEXAGE == 34 ){ sSEX =    "men"; sAGE = "age0064";
}else                    { stop("ERROR: Invalid Sex and Age Groups !!!"); }
rm(iSEXAGE);

# Formula for the Best Linear Unbiased Predictions and Pooled Coefficients
FORMULA_META = COEF_MODEL ~ TEMP_AVG + TEMP_IQR + AGE_80MAX;

# Boolean for the Country Random Effects
bCOU_RAND_EFF = TRUE;

# Period for the Baseline Temperature
DAY1_TEMP =  1; MON1_TEMP =  1; YEA1_TEMP = 1991; DATE1_TEMP = as.Date( paste( YEA1_TEMP, MON1_TEMP, DAY1_TEMP, sep = "-" ) );
DAY2_TEMP = 31; MON2_TEMP = 12; YEA2_TEMP = 2020; DATE2_TEMP = as.Date( paste( YEA2_TEMP, MON2_TEMP, DAY2_TEMP, sep = "-" ) ); # All Data in 1991-2020
if( DAY1_TEMP !=  1 | MON1_TEMP !=  1 |
    DAY2_TEMP != 31 | MON2_TEMP != 12 ){ stop("ERROR: Invalid Baseline Temperature Dates !!!"); }

# Period for the Calibration of the Epidemiological Associations
DAY1_CALI =  2; MON1_CALI =  1; YEA1_CALI = 2003; DATE1_CALI = as.Date( paste( YEA1_CALI, MON1_CALI, DAY1_CALI, sep = "-" ) );
DAY2_CALI = 26; MON2_CALI = 12; YEA2_CALI = 2019; DATE2_CALI = as.Date( paste( YEA2_CALI, MON2_CALI, DAY2_CALI, sep = "-" ) ); # All Data in 2003-2019
if( lubridate::wday( DATE1_CALI, week_start = 1 ) != 4 | lubridate::wday( DATE2_CALI, week_start = 1 ) != 4 ){ stop("ERROR: Invalid Calibration Dates !!!"); }

# Period for the Predictions
DAY1_PRED =  2; MON1_PRED =  1; YEA1_PRED = 2003; DATE1_PRED = as.Date( paste( YEA1_PRED, MON1_PRED, DAY1_PRED, sep = "-" ) );
DAY2_PRED =  3; MON2_PRED = 11; YEA2_PRED = 2022; DATE2_PRED = as.Date( paste( YEA2_PRED, MON2_PRED, DAY2_PRED, sep = "-" ) ); # All Data in 2003-2022
if( lubridate::wday( DATE1_PRED, week_start = 1 ) != 4 | lubridate::wday( DATE2_PRED, week_start = 1 ) != 4 ){ stop("ERROR: Invalid Prediction Dates !!!"); }

# Period for Year 2022
DAY1_YE22 =  6; MON1_YE22 =  1; YEA1_YE22 = 2022; DATE1_YE22 = as.Date( paste( YEA1_YE22, MON1_YE22, DAY1_YE22, sep = "-" ) );
DAY2_YE22 =  3; MON2_YE22 = 11; YEA2_YE22 = 2022; DATE2_YE22 = as.Date( paste( YEA2_YE22, MON2_YE22, DAY2_YE22, sep = "-" ) ); # Year 2022
if( lubridate::wday( DATE1_YE22, week_start = 1 ) != 4 | lubridate::wday( DATE2_YE22, week_start = 1 ) != 4 ){ stop("ERROR: Invalid Year 2022 Dates !!!"); }

# Period for Summer 2022
DAY1_SU22 =  2; MON1_SU22 =  6; YEA1_SU22 = 2022; DATE1_SU22 = as.Date( paste( YEA1_SU22, MON1_SU22, DAY1_SU22, sep = "-" ) );
DAY2_SU22 =  1; MON2_SU22 =  9; YEA2_SU22 = 2022; DATE2_SU22 = as.Date( paste( YEA2_SU22, MON2_SU22, DAY2_SU22, sep = "-" ) ); # Summer 2022
if( lubridate::wday( DATE1_SU22, week_start = 1 ) != 4 | lubridate::wday( DATE2_SU22, week_start = 1 ) != 4 ){ stop("ERROR: Invalid Summer 2022 Dates !!!"); }

# Period for Weeks 28-32, 2022
DAY1_PEAK = 14; MON1_PEAK =  7; YEA1_PEAK = 2022; DATE1_PEAK = as.Date( paste( YEA1_PEAK, MON1_PEAK, DAY1_PEAK, sep = "-" ) );
DAY2_PEAK = 11; MON2_PEAK =  8; YEA2_PEAK = 2022; DATE2_PEAK = as.Date( paste( YEA2_PEAK, MON2_PEAK, DAY2_PEAK, sep = "-" ) ); # Weeks 28-32, 2022
if( lubridate::wday( DATE1_PEAK, week_start = 1 ) != 4 | lubridate::wday( DATE2_PEAK, week_start = 1 ) != 4 ){ stop("ERROR: Invalid Summer Weeks 28-32 of 2022 Dates !!!"); }

# Exposure-Response: Type and Degree of the Spline
VAR_FUN = "ns"; VAR_DEG = NA;

# Exposure-Response: Percentiles of the Temperature Knots
VAR_PRC = c(10,50,90) / 100;

# Lag-Response: Minimum Lag (in Weeks)
MIN_LAG = 0; if( MIN_LAG <     0    ){ stop( "ERROR: Invalid MIN_LAG !!!" ); }
MAX_LAG = 3; if( MAX_LAG <= MIN_LAG ){ stop( "ERROR: Invalid MAX_LAG !!!" ); }

# Degrees of Freedom per Year for the Seasonal and Long-Term Trends for the Periods Pre-Covid (A, before DATE1_C19), Covid (B, DATE1_C19-DATE2_C19) and Post-Covid (C, after DATE2_C19)
DF_SEAS_C19A = DF_SEAS_C19B = DF_SEAS_C19C = 8;
if( DF_SEAS_C19A <= 0 | DF_SEAS_C19B <= 0 | DF_SEAS_C19C <= 0 ){ stop( "ERROR: Invalid DF_SEAS !!!" ); }

#iSENSITIVITY_OPTION = 7 ###split by regions
iSENSITIVITY_OPTION = 1;####using avg and TQR
if      ( iSENSITIVITY_OPTION ==  0 ){ # Baseline Model
}else if( iSENSITIVITY_OPTION ==  1 ){ FORMULA_META = COEF_MODEL ~ TEMP_AVG + TEMP_IQR            ; bCOU_RAND_EFF =  TRUE;
}else if( iSENSITIVITY_OPTION ==  2 ){ FORMULA_META = COEF_MODEL ~ TEMP_AVG +            AGE_80MAX; bCOU_RAND_EFF =  TRUE;
}else if( iSENSITIVITY_OPTION ==  3 ){ FORMULA_META = COEF_MODEL ~            TEMP_IQR + AGE_80MAX; bCOU_RAND_EFF =  TRUE;
}else if( iSENSITIVITY_OPTION ==  4 ){ FORMULA_META = COEF_MODEL ~ TEMP_AVG                       ; bCOU_RAND_EFF =  TRUE;
}else if( iSENSITIVITY_OPTION ==  5 ){ FORMULA_META = COEF_MODEL ~            TEMP_IQR            ; bCOU_RAND_EFF =  TRUE;
}else if( iSENSITIVITY_OPTION ==  6 ){ FORMULA_META = COEF_MODEL ~                       AGE_80MAX; bCOU_RAND_EFF =  TRUE;
}else if( iSENSITIVITY_OPTION ==  7 ){ FORMULA_META = COEF_MODEL ~                REG0               ; bCOU_RAND_EFF =  TRUE;
}else if( iSENSITIVITY_OPTION ==  8 ){ FORMULA_META = COEF_MODEL ~ TEMP_AVG + TEMP_IQR + AGE_80MAX; bCOU_RAND_EFF = FALSE;
}else if( iSENSITIVITY_OPTION ==  9 ){
  
  DAY1_CALI =  2; MON1_CALI =  1; YEA1_CALI = 2003; DATE1_CALI = as.Date( paste( YEA1_CALI, MON1_CALI, DAY1_CALI, sep = "-" ) );
  DAY2_CALI =  3; MON2_CALI = 11; YEA2_CALI = 2022; DATE2_CALI = as.Date( paste( YEA2_CALI, MON2_CALI, DAY2_CALI, sep = "-" ) ); # All Data in 2003-2022
  if( lubridate::wday( DATE1_CALI, week_start = 1 ) != 4 | lubridate::wday( DATE2_CALI, week_start = 1 ) != 4 ){ stop("ERROR: Invalid Calibration Dates !!!"); }
  
}else if( iSENSITIVITY_OPTION == 10 ){
  
  DAY1_CALI =  2; MON1_CALI =  1; YEA1_CALI = 2020; DATE1_CALI = as.Date( paste( YEA1_CALI, MON1_CALI, DAY1_CALI, sep = "-" ) );
  DAY2_CALI =  3; MON2_CALI = 11; YEA2_CALI = 2022; DATE2_CALI = as.Date( paste( YEA2_CALI, MON2_CALI, DAY2_CALI, sep = "-" ) ); # All Data in 2020-2022
  if( lubridate::wday( DATE1_CALI, week_start = 1 ) != 4 | lubridate::wday( DATE2_CALI, week_start = 1 ) != 4 ){ stop("ERROR: Invalid Calibration Dates !!!"); }
  
}else if( iSENSITIVITY_OPTION == 11 ){ VAR_FUN = "bs"; VAR_DEG = 2;
}else if( iSENSITIVITY_OPTION == 12 ){ VAR_PRC = c(10,25,50,75,90) / 100;
}else if( iSENSITIVITY_OPTION == 13 ){ VAR_PRC = c(10,25,75,90) / 100;
}else if( iSENSITIVITY_OPTION == 14 ){ VAR_PRC = c(10,75,90) / 100;
}else if( iSENSITIVITY_OPTION == 15 ){ VAR_PRC = c(10,90) / 100;
}else if( iSENSITIVITY_OPTION == 16 ){ VAR_PRC = c(25,75) / 100;
}else if( iSENSITIVITY_OPTION == 17 ){ MAX_LAG =  1; if( MAX_LAG <= MIN_LAG ){ stop( "ERROR: Invalid MAX_LAG !!!" ); }
}else if( iSENSITIVITY_OPTION == 18 ){ MAX_LAG =  2; if( MAX_LAG <= MIN_LAG ){ stop( "ERROR: Invalid MAX_LAG !!!" ); }
}else if( iSENSITIVITY_OPTION == 19 ){ MAX_LAG =  4; if( MAX_LAG <= MIN_LAG ){ stop( "ERROR: Invalid MAX_LAG !!!" ); }
}else if( iSENSITIVITY_OPTION == 20 ){ DF_SEAS_C19A = DF_SEAS_C19B = DF_SEAS_C19C =  6; if( DF_SEAS_C19A <= 0 | DF_SEAS_C19B <= 0 | DF_SEAS_C19C <= 0 ){ stop( "ERROR: Invalid DF_SEAS !!!" ); }
}else if( iSENSITIVITY_OPTION == 21 ){ DF_SEAS_C19A = DF_SEAS_C19B = DF_SEAS_C19C =  7; if( DF_SEAS_C19A <= 0 | DF_SEAS_C19B <= 0 | DF_SEAS_C19C <= 0 ){ stop( "ERROR: Invalid DF_SEAS !!!" ); }
}else if( iSENSITIVITY_OPTION == 22 ){ DF_SEAS_C19A = DF_SEAS_C19B = DF_SEAS_C19C =  9; if( DF_SEAS_C19A <= 0 | DF_SEAS_C19B <= 0 | DF_SEAS_C19C <= 0 ){ stop( "ERROR: Invalid DF_SEAS !!!" ); }
}else if( iSENSITIVITY_OPTION == 23 ){ DF_SEAS_C19A = DF_SEAS_C19B = DF_SEAS_C19C = 10; if( DF_SEAS_C19A <= 0 | DF_SEAS_C19B <= 0 | DF_SEAS_C19C <= 0 ){ stop( "ERROR: Invalid DF_SEAS !!!" ); }
}else if( iSENSITIVITY_OPTION == 24 ){
  
  DF_SEAS_C19B = 24; if( DF_SEAS_C19A <= 0 | DF_SEAS_C19B <= 0 | DF_SEAS_C19C <= 0 ){ stop( "ERROR: Invalid DF_SEAS !!!" ); }
  
}else if( iSENSITIVITY_OPTION == 25 ){
  
  DF_SEAS_C19B = 24; if( DF_SEAS_C19A <= 0 | DF_SEAS_C19B <= 0 | DF_SEAS_C19C <= 0 ){ stop( "ERROR: Invalid DF_SEAS !!!" ); }
  
  DAY1_CALI =  2; MON1_CALI =  1; YEA1_CALI = 2003; DATE1_CALI = as.Date( paste( YEA1_CALI, MON1_CALI, DAY1_CALI, sep = "-" ) );
  DAY2_CALI =  3; MON2_CALI = 11; YEA2_CALI = 2022; DATE2_CALI = as.Date( paste( YEA2_CALI, MON2_CALI, DAY2_CALI, sep = "-" ) ); # All Data in 2003-2022
  if( lubridate::wday( DATE1_CALI, week_start = 1 ) != 4 | lubridate::wday( DATE2_CALI, week_start = 1 ) != 4 ){ stop("ERROR: Invalid Calibration Dates !!!"); }
  
}else if( iSENSITIVITY_OPTION == 26 ){
  
  DF_SEAS_C19B = 24; if( DF_SEAS_C19A <= 0 | DF_SEAS_C19B <= 0 | DF_SEAS_C19C <= 0 ){ stop( "ERROR: Invalid DF_SEAS !!!" ); }
  
  DAY1_CALI =  2; MON1_CALI =  1; YEA1_CALI = 2020; DATE1_CALI = as.Date( paste( YEA1_CALI, MON1_CALI, DAY1_CALI, sep = "-" ) );
  DAY2_CALI =  3; MON2_CALI = 11; YEA2_CALI = 2022; DATE2_CALI = as.Date( paste( YEA2_CALI, MON2_CALI, DAY2_CALI, sep = "-" ) ); # All Data in 2020-2022
  if( lubridate::wday( DATE1_CALI, week_start = 1 ) != 4 | lubridate::wday( DATE2_CALI, week_start = 1 ) != 4 ){ stop("ERROR: Invalid Calibration Dates !!!"); }
  
}else                                { stop("ERROR: Invalid Sensitivity Option !!!");
}

# Period for the Initial Covid-19 Waves
#    This Is Only Used if DF_SEAS_C19A, DF_SEAS_C19B and DF_SEAS_C19C Are Not All Equal
DAY1_C19 =  6; MON1_C19 =  2; YEA1_C19 = 2020; DATE1_C19 = as.Date( paste( YEA1_C19, MON1_C19, DAY1_C19, sep = "-" ) );
DAY2_C19 = 27; MON2_C19 =  5; YEA2_C19 = 2021; DATE2_C19 = as.Date( paste( YEA2_C19, MON2_C19, DAY2_C19, sep = "-" ) ); # Initial Covid-19 Waves
if( lubridate::wday( DATE1_C19, week_start = 1 ) != 4 | lubridate::wday( DATE2_C19, week_start = 1 ) != 4 ){ stop("ERROR: Invalid Initial Covid-19 Waves Dates !!!"); }

# Cumulative Exposure-Response: Percentiles for the Predictions
PRED_PRC = sort( unique( c( seq(  0.0,   1.0, 0.1 ),
                            seq(  1.5,   5.0, 0.5 ),
                            seq(  6.0,  94.0, 1.0 ),
                            seq( 95.0,  98.5, 0.5 ),
                            seq( 99.0, 100.0, 0.1 ) ) / 100 ) );
if( any( 0 > PRED_PRC | PRED_PRC > 1 ) ){ stop("ERROR: Invalid Percentile Vector for the Predictions of the Cumulative Exposure-Response !!!"); }

# Boolean Indicating if the Calculation of the Minimum Mortality Temperature is Based on the Local Minima of the Cumulative Exposure-Response
bLOCAL_MINIMUM_MMT = FALSE;

# Otherwise, or if there Is No Local Minimum, Lower and Upper Percentiles for the Calculation of the Minimum Mortality Temperature
MIN_PMMT =   5 / 100;
MAX_PMMT = 100 / 100;

# Start of Redirection of All Console Output
foldout = paste0( "./dataout/", sSEX, "_", sAGE, "_iSENSITIVITY.", iSENSITIVITY_OPTION, "_", sFACTUAL, "/" );
if( !file_test( "-d", foldout ) ){ dir.create( foldout, recursive = TRUE ); }
sink( paste0( foldout, "log_", sSEX, "_", sAGE, "_iSENSITIVITY.", iSENSITIVITY_OPTION, "_", sFACTUAL, ".txt" ) );



################################################################################
### Preparation of the Data
################################################################################

print("");
print("= Preparation of the Data =");
print("");

# Reading the Weekly Temperature Table
if( sSEX == "allsex" & sAGE != "allage" ){
  DATATABLE_TEMP = read_delim( paste0( "./datain/weekly_table_temp_women_", sAGE, "_1950_2022.csv" ), show_col_types = FALSE );
}else{
  DATATABLE_TEMP = read_delim( paste0( "./datain/weekly_table_temp_", sSEX, "_", sAGE, "_1950_2022.csv" ), show_col_types = FALSE );
}
DATATABLE_TEMP$year = as.numeric( DATATABLE_TEMP$year );
DATATABLE_TEMP$woy  = as.numeric( DATATABLE_TEMP$woy  );
DATATABLE_TEMP$temp = as.numeric( DATATABLE_TEMP$temp );

# Transforming Temperatures from Kelvin to Celsius
DATATABLE_TEMP$temp = DATATABLE_TEMP$temp - 273.15;

# Reading the Annual Population Table
if( sSEX == "allsex" & sAGE != "allage" ){
  DATATABLE_POPU_WOM = read_delim( paste0( "./datain/weekly_table_popu_women_", sAGE, "_1950_2022_demo_r_pjanaggr3.csv" ), show_col_types = FALSE );
  DATATABLE_POPU_MEN = read_delim( paste0( "./datain/weekly_table_popu_men_",   sAGE, "_1950_2022_demo_r_pjanaggr3.csv" ), show_col_types = FALSE );
  if( any( dim( DATATABLE_POPU_WOM          ) != dim( DATATABLE_POPU_MEN          ) ) ){ stop( "ERROR: Invalid DATATABLE_POPU Dimensions !!!" ); }
  if( any(      DATATABLE_POPU_WOM$location   !=      DATATABLE_POPU_MEN$location   ) ){ stop( "ERROR: Invalid DATATABLE_POPU$location !!!" ); }
  if( any(      DATATABLE_POPU_WOM$year       !=      DATATABLE_POPU_MEN$year       ) ){ stop( "ERROR: Invalid DATATABLE_POPU$year !!!" ); }
  if( any(      DATATABLE_POPU_WOM$age        !=      DATATABLE_POPU_MEN$age        ) ){ stop( "ERROR: Invalid DATATABLE_POPU$age !!!" ); }
  DATATABLE_POPU = DATATABLE_POPU_WOM;
  DATATABLE_POPU$sex = "allsex";
  DATATABLE_POPU$popu = DATATABLE_POPU_WOM$popu + DATATABLE_POPU_MEN$popu;
  rm(DATATABLE_POPU_WOM,DATATABLE_POPU_MEN);
}else{
  DATATABLE_POPU = read_delim( paste0( "./datain/weekly_table_popu_", sSEX, "_", sAGE, "_1950_2022_demo_r_pjanaggr3.csv" ), show_col_types = FALSE );
}
DATATABLE_POPU$year = as.numeric( DATATABLE_POPU$year );
DATATABLE_POPU$popu = as.numeric( DATATABLE_POPU$popu );

# Checking that the Population Data Corresponds to the Correct Sex and Age Categories
if( any( DATATABLE_POPU$sex != sSEX ) ){ stop( "ERROR: Invalid DATATABLE_POPU$sex Vector !!!" ); }
if( any( DATATABLE_POPU$age != sAGE ) ){ stop( "ERROR: Invalid DATATABLE_POPU$age Vector !!!" ); }

# Deleting the Age and Sex Columns from the Population Data
DATATABLE_POPU$sex = DATATABLE_POPU$age = NULL;

# Merging the Weekly Temperature and Annual Population Tables
#    Missing Data is Assigned to NA
#    Annual Population Is Transformed into Weekly Population
DATATABLE_DATA = join( DATATABLE_TEMP, DATATABLE_POPU, by = c( "location", "year" ), type = "left" );
rm(DATATABLE_TEMP,DATATABLE_POPU);

# Reading the Weekly Mortality Table
if( sSEX == "allsex" & sAGE != "allage" ){
  DATATABLE_MORT_WOM = read_delim( paste0( "./datain/weekly_table_mort_women_", sAGE, "_2000_2023.csv" ), delim = " ", show_col_types = FALSE );
  DATATABLE_MORT_MEN = read_delim( paste0( "./datain/weekly_table_mort_men_",   sAGE, "_2000_2023.csv" ), delim = " ", show_col_types = FALSE );
  if( any( dim( DATATABLE_MORT_WOM          ) != dim( DATATABLE_MORT_MEN          ) ) ){ stop( "ERROR: Invalid DATATABLE_MORT Dimensions !!!" ); }
  if( any(      DATATABLE_MORT_WOM$location   !=      DATATABLE_MORT_MEN$location   ) ){ stop( "ERROR: Invalid DATATABLE_MORT$location !!!" ); }
  if( any(      DATATABLE_MORT_WOM$year       !=      DATATABLE_MORT_MEN$year       ) ){ stop( "ERROR: Invalid DATATABLE_MORT$year !!!" ); }
  if( any(      DATATABLE_MORT_WOM$woy        !=      DATATABLE_MORT_MEN$woy        ) ){ stop( "ERROR: Invalid DATATABLE_MORT$woy !!!" ); }
  if( any(      DATATABLE_MORT_WOM$cause      !=      DATATABLE_MORT_MEN$cause      ) ){ stop( "ERROR: Invalid DATATABLE_MORT$cause !!!" ); }
  if( any(      DATATABLE_MORT_WOM$age        !=      DATATABLE_MORT_MEN$age        ) ){ stop( "ERROR: Invalid DATATABLE_MORT$age !!!" ); }
  DATATABLE_MORT = DATATABLE_MORT_WOM;
  DATATABLE_MORT$sex = "allsex";
  DATATABLE_MORT$mort = DATATABLE_MORT_WOM$mort + DATATABLE_MORT_MEN$mort;
  rm(DATATABLE_MORT_WOM,DATATABLE_MORT_MEN);
}else{
  DATATABLE_MORT = read_delim( paste0( "./datain/weekly_table_mort_", sSEX, "_", sAGE, "_2000_2023.csv" ), delim = " ", show_col_types = FALSE );
}
DATATABLE_MORT$year = as.numeric( DATATABLE_MORT$year );
DATATABLE_MORT$woy  = as.numeric( DATATABLE_MORT$woy  );
DATATABLE_MORT$mort = as.numeric( DATATABLE_MORT$mort );

# Checking that the Mortality Data Corresponds to the Correct Cause, Sex and Age Categories
if( any( DATATABLE_MORT$cause != "allcause" ) ){ stop( "ERROR: Invalid DATATABLE_MORT$cause Vector !!!" ); }
if( any( DATATABLE_MORT$sex   !=    sSEX    ) ){ stop( "ERROR: Invalid DATATABLE_MORT$sex Vector !!!" ); }
if( any( DATATABLE_MORT$age   !=    sAGE    ) ){ stop( "ERROR: Invalid DATATABLE_MORT$age Vector !!!" ); }

# Deleting the Cause, Age and Sex Columns from the Mortality Data
DATATABLE_MORT$cause = DATATABLE_MORT$sex = DATATABLE_MORT$age = NULL;

# Merging the Weekly Temperature, Population and Mortality Tables
#    Missing Data is Assigned to NA
DATATABLE_DATA = join( DATATABLE_DATA, DATATABLE_MORT, by = c( "location", "year", "woy" ), type = "left" );
rm(DATATABLE_MORT);

# Adding the Date of the Thursday of the Week
DATATABLE_DATA$date = ISOweek2date( paste0( DATATABLE_DATA$year, "-W", sprintf( "%02d", DATATABLE_DATA$woy ), "-", 4 ) );
print( paste0( "Total Mortality (2003-): ", sum( DATATABLE_DATA$mort[ which( 2003 <= DATATABLE_DATA$year ) ], na.rm = TRUE ) ) );

# Reading the Weekly Metatable
if( sSEX == "allsex" & sAGE != "allage" ){
  METATABLE = read_delim( paste0( "./datain/weekly_table_meta_women_", sAGE, ".csv" ), show_col_types = FALSE );
}else{
  METATABLE = read_delim( paste0( "./datain/weekly_table_meta_", sSEX, "_", sAGE, ".csv" ), show_col_types = FALSE );
}

# Removing the Missing Countries from the Analysis
if( sAGE != "allage" ){
  paste0( "=== WARNING: Removing the United Kingdom, Ireland and Germany from the Analysis ===" );
  METATABLE = METATABLE[ which( METATABLE$nuts0 != "UK" ), ];
  METATABLE = METATABLE[ which( METATABLE$nuts0 != "IE" ), ];
  METATABLE = METATABLE[ which( METATABLE$nuts0 != "DE" ), ];
}else if( sSEX != "allsex" & sAGE == "allage" ){
  paste0( "=== WARNING: Removing the United Kingdom from the Analysis ===" );
  METATABLE = METATABLE[ which( METATABLE$nuts0 != "UK" ), ];
}else if( sSEX == "allsex" & sAGE == "allage" ){
  # No Country to Remove, Do Nothing
}else{
  stop( "ERROR: This Case Should Not Exist !!!" );
}

### TO REMOVE ###
# Removing the Few Regions in Which the Fitting Does Not Converge from the Analysis
if      ( sSEX == "women" & sAGE == "age0064" ){ iREG = which( METATABLE$nuts3 == "EL643" ); METATABLE = METATABLE[ c( 1 : ( iREG - 1 ), ( iREG + 1 ) : dim(METATABLE)[1] ), ]; rm(iREG);
}else if( sSEX == "women" & sAGE == "age6579" ){ iREG = which( METATABLE$nuts3 == "CH054" ); METATABLE = METATABLE[ c( 1 : ( iREG - 1 ), ( iREG + 1 ) : dim(METATABLE)[1] ), ]; rm(iREG);
}

# Creating the Vector of NUTS Region Codes
vREG = METATABLE$location;
nREG = length(vREG);

# Creating the Vector of English Region Names
vREG_NAME = METATABLE$name_eng;

# Creating the Vector with the Countries of Each Region
vCOU_of_REG = METATABLE$nuts0;
print( "Number of Regions in Each Country" );
print( count( vCOU_of_REG ) );

# Creating the Vector of Countries
vCOU = unique(METATABLE$nuts0);
nCOU = length(vCOU);

# Creating the Vector of English Country Names
vCOU_NAME = vCOU;
vCOU_NAME[ which( vCOU_NAME == "AL" ) ] = "Albania";
vCOU_NAME[ which( vCOU_NAME == "AT" ) ] = "Austria";
vCOU_NAME[ which( vCOU_NAME == "BE" ) ] = "Belgium";
vCOU_NAME[ which( vCOU_NAME == "BG" ) ] = "Bulgaria";
vCOU_NAME[ which( vCOU_NAME == "CH" ) ] = "Switzerland";
vCOU_NAME[ which( vCOU_NAME == "CY" ) ] = "Cyprus";
vCOU_NAME[ which( vCOU_NAME == "CZ" ) ] = "Czechia";
vCOU_NAME[ which( vCOU_NAME == "DE" ) ] = "Germany";
vCOU_NAME[ which( vCOU_NAME == "DK" ) ] = "Denmark";
vCOU_NAME[ which( vCOU_NAME == "EE" ) ] = "Estonia";
vCOU_NAME[ which( vCOU_NAME == "EL" ) ] = "Greece";
vCOU_NAME[ which( vCOU_NAME == "ES" ) ] = "Spain";
vCOU_NAME[ which( vCOU_NAME == "FI" ) ] = "Finland";
vCOU_NAME[ which( vCOU_NAME == "FR" ) ] = "France";
vCOU_NAME[ which( vCOU_NAME == "HR" ) ] = "Croatia";
vCOU_NAME[ which( vCOU_NAME == "HU" ) ] = "Hungary";
vCOU_NAME[ which( vCOU_NAME == "IE" ) ] = "Ireland";
vCOU_NAME[ which( vCOU_NAME == "IS" ) ] = "Iceland";
vCOU_NAME[ which( vCOU_NAME == "IT" ) ] = "Italy";
vCOU_NAME[ which( vCOU_NAME == "LI" ) ] = "Liechtenstein";
vCOU_NAME[ which( vCOU_NAME == "LT" ) ] = "Lithuania";
vCOU_NAME[ which( vCOU_NAME == "LU" ) ] = "Luxembourg";
vCOU_NAME[ which( vCOU_NAME == "LV" ) ] = "Latvia";
vCOU_NAME[ which( vCOU_NAME == "ME" ) ] = "Montenegro";
vCOU_NAME[ which( vCOU_NAME == "MT" ) ] = "Malta";
vCOU_NAME[ which( vCOU_NAME == "NL" ) ] = "Netherlands";
vCOU_NAME[ which( vCOU_NAME == "NO" ) ] = "Norway";
vCOU_NAME[ which( vCOU_NAME == "PL" ) ] = "Poland";
vCOU_NAME[ which( vCOU_NAME == "PT" ) ] = "Portugal";
vCOU_NAME[ which( vCOU_NAME == "RO" ) ] = "Romania";
vCOU_NAME[ which( vCOU_NAME == "RS" ) ] = "Serbia";
vCOU_NAME[ which( vCOU_NAME == "SE" ) ] = "Sweden";
vCOU_NAME[ which( vCOU_NAME == "SI" ) ] = "Slovenia";
vCOU_NAME[ which( vCOU_NAME == "SK" ) ] = "Slovakia";
vCOU_NAME[ which( vCOU_NAME == "UK" ) ] = "United Kingdom";

# Restricting the Data to the Specified Periods
#    The Datatables for the Calibration of the Epidemiological Associations and the Predictions Include MAX_LAG Additional Weeks to Calculate the Attributable Mortality
DATATABLE_TEMP = DATATABLE_DATA[ which( DATE1_TEMP <= DATATABLE_DATA$date & DATATABLE_DATA$date <= DATE2_TEMP                         ), ];
DATATABLE_CALI = DATATABLE_DATA[ which( DATE1_CALI <= DATATABLE_DATA$date & DATATABLE_DATA$date <= DATE2_CALI + max( 7 * MAX_LAG, 0 ) ), ];
DATATABLE_PRED = DATATABLE_DATA[ which( DATE1_PRED <= DATATABLE_DATA$date & DATATABLE_DATA$date <= DATE2_PRED + max( 7 * MAX_LAG, 0 ) ), ];
rm(DATATABLE_DATA);

# Creating the Datalists for the Specified Periods
DATALIST_TEMP = lapply( vREG, function(x) DATATABLE_TEMP[ DATATABLE_TEMP$location == x, ] ); names(DATALIST_TEMP) = vREG;
DATALIST_CALI = lapply( vREG, function(x) DATATABLE_CALI[ DATATABLE_CALI$location == x, ] ); names(DATALIST_CALI) = vREG;
DATALIST_PRED = lapply( vREG, function(x) DATATABLE_PRED[ DATATABLE_PRED$location == x, ] ); names(DATALIST_PRED) = vREG;
rm(DATATABLE_TEMP,DATATABLE_CALI,DATATABLE_PRED,METATABLE);

TOTAL_POPU = 0;
for( iREG in 1:nREG ){
  TOTAL_POPU = TOTAL_POPU + rev( DATALIST_PRED[[iREG]]$popu )[1];
}
print( paste0( "Total Population (2022): ", TOTAL_POPU ) );
print( paste0( "Number of Regions: ", nREG ) );
print( paste0( "Number of Countries: ", nCOU ) );
rm(TOTAL_POPU,iREG);

# Adding the Mean Annual Cycle of Temperature and Mortality to the Period for the Predictions
for( iREG in 1:nREG ){
  DATALIST_PRED[[iREG]]$temp_mac = mean_annual_cycle_temp_mort( TRUE,
                                                                DATALIST_TEMP[[iREG]]$year,
                                                                DATALIST_TEMP[[iREG]]$woy,
                                                                DATALIST_TEMP[[iREG]]$temp,
                                                                DATALIST_PRED[[iREG]]$year,
                                                                DATALIST_PRED[[iREG]]$woy );
  DATALIST_PRED[[iREG]]$mort_mac = mean_annual_cycle_temp_mort( FALSE,
                                                                DATALIST_CALI[[iREG]]$year,
                                                                DATALIST_CALI[[iREG]]$woy,
                                                                DATALIST_CALI[[iREG]]$mort,
                                                                DATALIST_PRED[[iREG]]$year,
                                                                DATALIST_PRED[[iREG]]$woy );
}
rm(iREG);

# Filling the Missing Mortality Data in the Prediction Period by Using the Mean Annual Cycle during the Calibration Period
print( "=== WARNING: Missing Mortality Values in the Prediction Period Are Replaced by its Mean Annual Cycle during the Calibration Period ===" );
for( iREG in 1:nREG ){
  isNA = !is.finite( DATALIST_PRED[[iREG]]$mort );
  if( any( isNA ) ){
    if( any( !is.finite( DATALIST_PRED[[iREG]]$mort_mac[ which(isNA) ] ) ) ){ stop( "ERROR: Impossible to Fill Missing Mortality Data with its Mean Annual Cycle !!!" ); }
    if( sum(isNA) > 0 ){ print( paste0( "      ", sum(isNA), " Missing Mortality Values in Region ", vREG[iREG], ": from ", DATALIST_PRED[[iREG]]$year[      which( isNA )  [1] ], "-W" , DATALIST_PRED[[iREG]]$woy[      which( isNA )  [1] ],
                                        " to ", DATALIST_PRED[[iREG]]$year[ rev( which( isNA ) )[1] ], "-W" , DATALIST_PRED[[iREG]]$woy[ rev( which( isNA ) )[1] ] ) ); }
    DATALIST_PRED[[iREG]]$mort[ which(isNA) ] = DATALIST_PRED[[iREG]]$mort_mac[ which(isNA) ];
  }
  rm(isNA);
}
rm(iREG);

# Adding the Variable Week-of-Period to the Specified Periods
#    This Variable is Used for the Spline of Time, See FORMULA Below
for( iREG in 1:nREG ){
  DATALIST_CALI[[iREG]]$wop = 1:length(DATALIST_CALI[[iREG]]$woy);
}
rm(iREG);

# Checking that the Datalist for the Baseline Temperature Has the Right Format
for( iREG in 1:nREG ){
  if( 0 >      DATALIST_TEMP[[iREG]]$date  [1] - DATE1_TEMP |      DATALIST_TEMP[[iREG]]$date  [1] - DATE1_TEMP > +6 ){ stop( "ERROR: Invalid Initial Week in DATALIST_TEMP !!!" ); }
  if( 0 < rev( DATALIST_TEMP[[iREG]]$date )[1] - DATE2_TEMP | rev( DATALIST_TEMP[[iREG]]$date )[1] - DATE2_TEMP < -6 ){ stop( "ERROR: Invalid Final Week in DATALIST_TEMP !!!" ); }
  if( any( diff( DATALIST_TEMP[[iREG]]$date ) != 7 ) ){ stop( "ERROR: Missing or Disordered Weeks in DATALIST_TEMP !!!" ); }
  if( any( !is.finite( DATALIST_TEMP[[iREG]]$temp ) | DATALIST_TEMP[[iREG]]$temp <  -50 | DATALIST_TEMP[[iREG]]$temp >  50 ) ){ stop( "ERROR: Missing or Invalid Temperature Data in DATALIST_TEMP !!!" ); }
  if( any( !is.finite( DATALIST_TEMP[[iREG]]$popu ) | DATALIST_TEMP[[iREG]]$popu <=   0                                    ) ){ stop( "ERROR: Missing or Invalid Population Data in DATALIST_TEMP !!!" ); }
}
rm(iREG);

# Checking that the Datalist for the Calibration of the Epidemiological Associations Has the Right Format
for( iREG in 1:nREG ){
  if( 0 >      DATALIST_CALI[[iREG]]$date  [1] - DATE1_CALI                         |      DATALIST_CALI[[iREG]]$date  [1] - DATE1_CALI                         > +6 ){ stop( "ERROR: Invalid Initial Week in DATALIST_CALI !!!" ); }
  if( 0 < rev( DATALIST_CALI[[iREG]]$date )[1] - DATE2_CALI - max( 7 * MAX_LAG, 0 ) | rev( DATALIST_CALI[[iREG]]$date )[1] - DATE2_CALI - max( 7 * MAX_LAG, 0 ) < -6 ){ stop( "ERROR: Invalid Final Week in DATALIST_CALI !!!" ); }
  if( any( diff( DATALIST_CALI[[iREG]]$date ) != 7 ) ){ stop( "ERROR: Missing or Disordered Weeks in DATALIST_CALI !!!" ); }
  if( any( !is.finite( DATALIST_CALI[[iREG]]$temp ) | DATALIST_CALI[[iREG]]$temp <  -50 | DATALIST_CALI[[iREG]]$temp >  50 ) ){ stop( "ERROR: Missing or Invalid Temperature Data in DATALIST_CALI !!!" ); }
  if( any( !is.finite( DATALIST_CALI[[iREG]]$popu ) | DATALIST_CALI[[iREG]]$popu <=   0                                    ) ){ stop( "ERROR: Missing or Invalid Population Data in DATALIST_CALI !!!" ); }
  ### README: THE CODE IS DESIGNED SO THAT THE CALIBRATION PERIOD ENDS AT 2019. AFTER 2019, THERE ARE MISSING MORTALITY VALUES. THIS CASE IS ONLY CONSIDERED IN SENSITIVITY ANALYSES ###
  ### TO REMOVE ###
  if( YEA2_CALI <= 2019 & !( vREG[iREG] == "DK014" & sAGE == "age0064" & sSEX ==  "women" |
                             vREG[iREG] == "DK014" & sAGE == "age0064" & sSEX == "allsex" |
                             vREG[iREG] == "DK014" & sAGE == "age80p"  & sSEX ==   "men"  |
                             vREG[iREG] == "DK014" & sAGE == "age80p"  & sSEX == "allsex" ) ){
    if( any( !is.finite( DATALIST_CALI[[iREG]]$mort ) | DATALIST_CALI[[iREG]]$mort <    0                                    ) ){ DATALIST_CALI[[iREG]]=DATALIST_CALI[[iREG]][!is.na(DATALIST_CALI[[iREG]]$mort),] #stop( "ERROR: Missing or Invalid Mortality Data in DATALIST_CALI !!!" ); 
    }
    if( any( !is.finite( DATALIST_CALI[[iREG]]$wop  ) | DATALIST_CALI[[iREG]]$wop  <=   0                                    ) ){ stop( "ERROR: Missing or Invalid Week-of-the-Period Data in DATALIST_CALI !!!" ); }
  }}
rm(iREG);

# Checking that the Datalist for the Predictions Has the Right Format
for( iREG in 1:nREG ){
  if( 0 >      DATALIST_PRED[[iREG]]$date  [1] - DATE1_PRED                         |      DATALIST_PRED[[iREG]]$date  [1] - DATE1_PRED                         > +6 ){ stop( "ERROR: Invalid Initial Week in DATALIST_PRED !!!" ); }
  if( 0 < rev( DATALIST_PRED[[iREG]]$date )[1] - DATE2_PRED - max( 7 * MAX_LAG, 0 ) | rev( DATALIST_PRED[[iREG]]$date )[1] - DATE2_PRED - max( 7 * MAX_LAG, 0 ) < -6 ){ stop( "ERROR: Invalid Final Week in DATALIST_PRED !!!" ); }
  if( any( diff( DATALIST_PRED[[iREG]]$date ) != 7 ) ){ stop( "ERROR: Missing or Disordered Weeks in DATALIST_PRED !!!" ); }
  if( any( !is.finite( DATALIST_PRED[[iREG]]$temp     ) | DATALIST_PRED[[iREG]]$temp     <  -50 | DATALIST_PRED[[iREG]]$temp     >  50 ) ){ stop( "ERROR: Missing or Invalid Temperature Data in DATALIST_PRED !!!" ); }
  if( any( !is.finite( DATALIST_PRED[[iREG]]$popu     ) | DATALIST_PRED[[iREG]]$popu     <=   0                                        ) ){ stop( "ERROR: Missing or Invalid Population Data in DATALIST_PRED !!!" ); }
  if( any( !is.finite( DATALIST_PRED[[iREG]]$mort     ) | DATALIST_PRED[[iREG]]$mort     <    0                                        ) ){ stop( "ERROR: Missing or Invalid Mortality Data in DATALIST_PRED !!!" ); }
  if( any( !is.finite( DATALIST_PRED[[iREG]]$temp_mac ) | DATALIST_PRED[[iREG]]$temp_mac <  -50 | DATALIST_PRED[[iREG]]$temp_mac >  50 ) ){ stop( "ERROR: Missing or Invalid Mean Annual Cycle of Temperature Data in DATALIST_PRED !!!" ); }
  if( any( !is.finite( DATALIST_PRED[[iREG]]$mort_mac ) | DATALIST_PRED[[iREG]]$mort_mac <    0                                        ) ){ stop( "ERROR: Missing or Invalid Mean Annual Cycle of Mortality Data in DATALIST_PRED !!!" ); }
}
rm(iREG);



################################################################################
### Calculation of the Scenario Temperatures
################################################################################

print("");
print("= Calculation of the Scenario Temperatures =");
print("");

if( bFACTUAL ){
  
  # Creating the Factual Temperatures
  for( iREG in 1:nREG ){
    DATALIST_TEMP[[iREG]]$temp_sce = DATALIST_TEMP[[iREG]]$temp;
    DATALIST_PRED[[iREG]]$temp_sce = DATALIST_PRED[[iREG]]$temp;
  }
  rm(iREG);
  
}else{
  
  # Loading the Counterfactual Temperature Differences
  load( file = paste0( "./datagmst/ANTHROPO_WARMING_", sSEX, "_", sAGE, "_REG.RData" ) );
  
  # Creating the Counterfactual Temperatures
  for( iREG in 1:nREG ){
    iiREG = which( ANTHROPO_WARMING_REG$location == vREG[iREG] );
    if( length(iiREG) != 1 ){ stop("ERROR: Invalid iiREG !!!"); }
    
    DATALIST_TEMP[[iREG]]$temp_sce = DATALIST_TEMP[[iREG]]$temp - ANTHROPO_WARMING_REG$temp_dif[iiREG];
    DATALIST_PRED[[iREG]]$temp_sce = DATALIST_PRED[[iREG]]$temp - ANTHROPO_WARMING_REG$temp_dif[iiREG];
    
    rm(iiREG);
  }
  rm(iREG);
  
  # Deleting the Counterfactual Temperature Differences
  rm(ANTHROPO_WARMING_REG);
  
}

# Adding the Mean Annual Cycle of Temperature to the Periods of Baseline Temperatures and Predictions
for( iREG in 1:nREG ){
  DATALIST_TEMP[[iREG]]$temp_sce_mac = mean_annual_cycle_temp_mort( TRUE,
                                                                    DATALIST_TEMP[[iREG]]$year,
                                                                    DATALIST_TEMP[[iREG]]$woy,
                                                                    DATALIST_TEMP[[iREG]]$temp_sce,
                                                                    DATALIST_TEMP[[iREG]]$year,
                                                                    DATALIST_TEMP[[iREG]]$woy );
  DATALIST_PRED[[iREG]]$temp_sce_mac = mean_annual_cycle_temp_mort( TRUE,
                                                                    DATALIST_TEMP[[iREG]]$year,
                                                                    DATALIST_TEMP[[iREG]]$woy,
                                                                    DATALIST_TEMP[[iREG]]$temp_sce,
                                                                    DATALIST_PRED[[iREG]]$year,
                                                                    DATALIST_PRED[[iREG]]$woy );
}
rm(iREG);



################################################################################
### Calculation of the Location-Specific Associations
################################################################################

print("");
print("= Calculation of the Location-Specific Associations =");
print("");

# Reduced Coefficients
if      ( VAR_FUN == "ns" ){ COEF_MODEL = matrix( data = NA, nREG, length(VAR_PRC) +    1   , dimnames = list(vREG) );
}else if( VAR_FUN == "bs" ){ COEF_MODEL = matrix( data = NA, nREG, length(VAR_PRC) + VAR_DEG, dimnames = list(vREG) );
}else                      { stop("ERROR: Invalid VAR_FUN !!!"); }

# Reduced Covariance Matrices
VCOV_MODEL = vector( "list", nREG );
names(VCOV_MODEL) = vREG;

# Cumulative Exposure-Response Before the Meta-Analysis
CROSS_PRED_REG_NOMETA = vector( "list", nREG );
names(CROSS_PRED_REG_NOMETA) = vREG;

# Akaike's Information Criterion for Overdispersed Count Data
vQAIC = array( NA, dim = c( nREG ), dimnames = list( vREG ) );

for( iREG in 1:nREG ){
  print( paste0( "   Region ", iREG, ": ", vREG_NAME[iREG], " (", vREG[iREG], ")" ) );
  
  # If the Degrees of Freedom per Year for the Seasonal and Long-Term Trends Are Equal for the Periods Pre-Covid (A), Covid (B) and Post-Covid (C)
  if( DF_SEAS_C19A == DF_SEAS_C19B & DF_SEAS_C19A == DF_SEAS_C19C ){
    
    # Formula of the Models
    #   length(wop)              :: Number of Weeks of the Period
    #   length(wop) * 7          :: Number of Days of the Period
    #   length(wop) * 7 / 365.25 :: Number of Years of the Period
    FORMULA_SEA = mort ~ ns( wop, df = round( DF_SEAS_C19A * length(wop) * 7 / 365.25 ) );
    FORMULA_CRB = mort ~ ns( wop, df = round( DF_SEAS_C19A * length(wop) * 7 / 365.25 ) ) + CROSS_BASIS;
    
    # Otherwise
  }else{
    
    # Vectors of All Dates of the Periods Pre-Covid, Covid and Post-Covid
    vDATEA = DATALIST_CALI[[iREG]]$date[ which(      DATALIST_CALI[[iREG]]$date  [1] <= DATALIST_CALI[[iREG]]$date & DATALIST_CALI[[iREG]]$date <=               DATE1_C19              ) ];
    vDATEB = DATALIST_CALI[[iREG]]$date[ which(               DATE1_C19              <= DATALIST_CALI[[iREG]]$date & DATALIST_CALI[[iREG]]$date <=               DATE2_C19              ) ];
    vDATEC = DATALIST_CALI[[iREG]]$date[ which(               DATE2_C19              <= DATALIST_CALI[[iREG]]$date & DATALIST_CALI[[iREG]]$date <= rev( DATALIST_CALI[[iREG]]$date )[1] ) ];
    
    # Vectors of Dates of the Knots of the Periods Pre-Covid, Covid and Post-Covid
    #   length(wop)              :: Number of Weeks of the Period
    #   length(wop) * 7          :: Number of Days of the Period
    #   length(wop) * 7 / 365.25 :: Number of Years of the Period
    vDATEA = vDATEA[ seq( from = 1, to = length(vDATEA), length.out = round( DF_SEAS_C19A * length(vDATEA) * 7 / 365.25 ) ) ];
    vDATEB = vDATEB[ seq( from = 1, to = length(vDATEB), length.out = round( DF_SEAS_C19B * length(vDATEB) * 7 / 365.25 ) ) ];
    vDATEC = vDATEC[ seq( from = 1, to = length(vDATEC), length.out = round( DF_SEAS_C19C * length(vDATEC) * 7 / 365.25 ) ) ];
    
    # Vector of Week-of-Period of the Knots of the Periods Pre-Covid, Covid and Post-Covid
    vKNOTS = DATALIST_CALI[[iREG]]$wop[ DATALIST_CALI[[iREG]]$date %in% unique( c( vDATEA, vDATEB, vDATEC ) ) ];
    rm(vDATEA,vDATEB,vDATEC);
    
    # Removing the First and Last Knot at the Boundaries
    vKNOTS = vKNOTS[ 2 : ( length(vKNOTS) - 1 ) ];
    
    # Formula of the Models
    FORMULA_SEA = mort ~ ns( wop, knots = vKNOTS );
    FORMULA_CRB = mort ~ ns( wop, knots = vKNOTS ) + CROSS_BASIS;
    
  }
  
  # Fitting of the Seasonality Model
  GLM_MODEL_SEA = glm( formula = FORMULA_SEA, DATALIST_CALI[[iREG]], family = quasipoisson, na.action = "na.exclude" );
  
  # Prediction of the Seasonality Model
  DATALIST_CALI[[iREG]]$mort_pred_seas = predict( GLM_MODEL_SEA, type = "response" );
  rm(GLM_MODEL_SEA);
  
  # Cross-Basis of Temperature
  if( VAR_FUN == "ns" ){
    CROSS_BASIS = crossbasis( DATALIST_CALI[[iREG]]$temp,
                              c( MIN_LAG, MAX_LAG ),
                              argvar = list( fun = VAR_FUN,                   knots = quantile( DATALIST_CALI[[iREG]]$temp, VAR_PRC, na.rm = TRUE ), Boundary.knots = range( DATALIST_CALI[[iREG]]$temp, na.rm = TRUE ) ),
                              arglag = list( fun = "integer" ) );
  }else if( VAR_FUN == "bs" ){
    CROSS_BASIS = crossbasis( DATALIST_CALI[[iREG]]$temp,
                              c( MIN_LAG, MAX_LAG ),
                              argvar = list( fun = VAR_FUN, degree = VAR_DEG, knots = quantile( DATALIST_CALI[[iREG]]$temp, VAR_PRC, na.rm = TRUE ) ),
                              arglag = list( fun = "integer" ) );
  }else{
    stop("ERROR: Invalid VAR_FUN !!!");
  }
  
  # Fitting of the Cross-Basis Model
  GLM_MODEL_CRB = glm( formula = FORMULA_CRB, DATALIST_CALI[[iREG]], family = quasipoisson, na.action = "na.exclude" );
  vQAIC[iREG] = QAIC( GLM_MODEL_CRB );
  
  # Cumulative Exposure-Response without Centring
  suppressMessages( CROSS_PRED_REG_NOMETA[[iREG]] <- crosspred( CROSS_BASIS, GLM_MODEL_CRB, at = quantile( DATALIST_CALI[[iREG]]$temp, PRED_PRC, na.rm = TRUE ) ) );
  
  # Vector of Local Minima of the Cumulative Exposure-Response
  iMMT = 1 + which( diff( sign( diff( CROSS_PRED_REG_NOMETA[[iREG]]$allRRfit ) ) ) == 2 );
  
  # Minimum Mortality Temperature
  if( bLOCAL_MINIMUM_MMT & length(iMMT) > 0 ){
    # If There Is at Least One Local Minimum, the Minimum Mortality Temperature Is the Local Minimum with the Lowest Relative Risk
    MMT = CROSS_PRED_REG_NOMETA[[iREG]]$predvar[ iMMT[ which.min( CROSS_PRED_REG_NOMETA[[iREG]]$allRRfit[iMMT] ) ] ];
  }else{
    # If There Are No Local Minima, the Minimum Mortality Temperature Is the Lowest Relative Risk within the Predefined Temperature Percentile Range [MIN_PMMT,MAX_PMMT]
    MMT = CROSS_PRED_REG_NOMETA[[iREG]]$predvar[ which( PRED_PRC == MIN_PMMT ) - 1 + which.min( CROSS_PRED_REG_NOMETA[[iREG]]$allRRfit[ which( PRED_PRC == MIN_PMMT ) : which( PRED_PRC == MAX_PMMT ) ] ) ];
  }
  rm(iMMT);
  
  # Cumulative Exposure-Response with Centring
  CROSS_PRED_REG_NOMETA[[iREG]] = crosspred( CROSS_BASIS, GLM_MODEL_CRB, at = quantile( DATALIST_CALI[[iREG]]$temp, PRED_PRC, na.rm = TRUE ), bylag = 1, cen = MMT );
  
  # Reduced Coefficients and Covariance Matrices
  REDUCED = crossreduce( CROSS_BASIS, GLM_MODEL_CRB, cen = MMT );
  COEF_MODEL[iREG,] = coef( REDUCED );
  VCOV_MODEL[[iREG]] = vcov( REDUCED );
  rm(REDUCED);
  
  rm(FORMULA_SEA,FORMULA_CRB, CROSS_BASIS,GLM_MODEL_CRB, MMT);
  if( !( DF_SEAS_C19A == DF_SEAS_C19B & DF_SEAS_C19A == DF_SEAS_C19C ) ){ rm(vKNOTS); }
  
}
rm(iREG);

print( paste0( "Akaike's Information Criterion for Overdispersed Count Data: ", round( sum( vQAIC ) ) ) );



################################################################################
### Calculation of the Best Linear Unbiased Predictions
################################################################################

print("");
print("= Calculation of the Best Linear Unbiased Predictions =");
print("");

# Temperature Meta-Predictors: Temperature Average (TEMP_AVG), Temperature Inter-Quartile Range (TEMP_IQR)
TEMP_AVG = sapply( DATALIST_CALI, function(x) mean( x$temp, na.rm = TRUE ) );
TEMP_IQR = sapply( DATALIST_CALI, function(x)  IQR( x$temp, na.rm = TRUE ) );

# Societal Meta-Predictor: Percentage of People Aged 80 or More (AGE_80MAX)
DATA_AGE_80MAX = read.csv( "./datasocietal/pc_y80_max.csv" );
AGE_80MAX = array( NA, dim = c( nREG ), dimnames = list( vREG ) );
for( iREG in 1:nREG ){
  if     ( vREG[iREG] == "BE224" ){ vREG1 = "BE221"; vREG2 =    NA  ; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "BE225" ){ vREG1 = "BE222"; vREG2 =    NA  ; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "BE328" ){ vREG1 = "BE324"; vREG2 = "BE327"; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "BE329" ){ vREG1 = "BE321"; vREG2 = "BE322"; vREG3 = "BE325"; vREG4 = "BE326"; }
  else if( vREG[iREG] == "BE32A" ){ vREG1 = "BE321"; vREG2 =    NA  ; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "BE32B" ){ vREG1 = "BE322"; vREG2 =    NA  ; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "BE32C" ){ vREG1 = "BE325"; vREG2 =    NA  ; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "BE32D" ){ vREG1 = "BE326"; vREG2 =    NA  ; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "EE009" ){ vREG1 = "EE006"; vREG2 =    NA  ; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "EE00A" ){ vREG1 = "EE007"; vREG2 =    NA  ; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "NO020" ){ vREG1 = "NO021"; vREG2 = "NO022"; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "NO074" ){ vREG1 = "NO072"; vREG2 = "NO073"; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "NO081" ){ vREG1 = "NO011"; vREG2 =    NA  ; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "NO082" ){ vREG1 = "NO012"; vREG2 = "NO031"; vREG3 = "NO032"; vREG4 =    NA  ; }
  else if( vREG[iREG] == "NO091" ){ vREG1 = "NO033"; vREG2 = "NO034"; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "NO092" ){ vREG1 = "NO041"; vREG2 = "NO042"; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "NO0A1" ){ vREG1 = "NO043"; vREG2 =    NA  ; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "NO0A2" ){ vREG1 = "NO051"; vREG2 = "NO052"; vREG3 =    NA  ; vREG4 =    NA  ; }
  else if( vREG[iREG] == "NO0A3" ){ vREG1 = "NO053"; vREG2 =    NA  ; vREG3 =    NA  ; vREG4 =    NA  ; }
  else                            { vREG1 = vREG[iREG]; vREG2 =    NA  ; vREG3 =    NA  ; vREG4 =    NA  ; }
  
  iROW = rep(NA,4);
  if( is.character(vREG1) ){ iROW[1] = which( DATA_AGE_80MAX$indic_de == "PC_Y80_MAX" & DATA_AGE_80MAX$unit == "PC" & DATA_AGE_80MAX$geo == vREG1 & DATA_AGE_80MAX$TIME_PERIOD == 2019 ); }
  if( is.character(vREG2) ){ iROW[2] = which( DATA_AGE_80MAX$indic_de == "PC_Y80_MAX" & DATA_AGE_80MAX$unit == "PC" & DATA_AGE_80MAX$geo == vREG2 & DATA_AGE_80MAX$TIME_PERIOD == 2019 ); }
  if( is.character(vREG3) ){ iROW[3] = which( DATA_AGE_80MAX$indic_de == "PC_Y80_MAX" & DATA_AGE_80MAX$unit == "PC" & DATA_AGE_80MAX$geo == vREG3 & DATA_AGE_80MAX$TIME_PERIOD == 2019 ); }
  if( is.character(vREG4) ){ iROW[4] = which( DATA_AGE_80MAX$indic_de == "PC_Y80_MAX" & DATA_AGE_80MAX$unit == "PC" & DATA_AGE_80MAX$geo == vREG4 & DATA_AGE_80MAX$TIME_PERIOD == 2019 ); }
  rm(vREG1,vREG2,vREG3,vREG4);
  
  if( all( !is.finite(iROW) ) ){ stop("ERROR: Elderly Percentage Not Found !!!"); }
  else                         { AGE_80MAX[iREG] = sum( DATA_AGE_80MAX[ iROW[ which( is.finite(iROW) ) ], "OBS_VALUE" ] ); }
  rm(iROW);
}
rm(DATA_AGE_80MAX, iREG);

# Minimum Mortality Temperature
MMT_REG = array( NA, dim = c( nREG ), dimnames = list( vREG ) );

# Cumulative Exposure-Response
CROSS_PRED_REG_META = vector( "list", nREG );
names(CROSS_PRED_REG_META) = vREG;

# Multivariate Meta-Analysis of the Reduced Coefficients
if( bCOU_RAND_EFF ){
  fCOU = factor( vCOU_of_REG );
  fREG = factor( vREG );
  MULTIVAR = mixmeta( FORMULA_META,
                      VCOV_MODEL,
                      data = data.frame( vREG = vREG ),
                      control = list( showiter = TRUE, igls.inititer = 10 ),
                      method = "reml",
                      random =~ 1 | fCOU/fREG );
  rm(fCOU,fREG);
}else{
  MULTIVAR = mixmeta( FORMULA_META,
                      VCOV_MODEL,
                      data = data.frame( vREG = vREG ),
                      control = list( showiter = TRUE, igls.inititer = 10 ),
                      method = "reml" );
}
print( summary( MULTIVAR ) );

# Wald Test of the Meta-Predictors
if( length( summary(MULTIVAR)$lab$p ) > 1 ){
  for( iMETAPRED in 1:length( summary(MULTIVAR)$lab$p ) ){
    print( paste0( "Wald Test of ", summary(MULTIVAR)$lab$p[iMETAPRED], ": p = ", sprintf( "%.10f", FWALD( MULTIVAR, summary(MULTIVAR)$lab$p[iMETAPRED] ) ) ) );
  }
  rm(iMETAPRED);
}

# Best Linear Unbiased Predictions
BLUP = blup( MULTIVAR, vcov = TRUE );

for( iREG in 1:nREG ){
  
  # Basis in the Temperature Domain
  if( VAR_FUN == "ns" ){
    BASIS_VAR = onebasis( quantile( DATALIST_CALI[[iREG]]$temp, PRED_PRC, na.rm = TRUE ),
                          fun = VAR_FUN,
                          knots = quantile( DATALIST_CALI[[iREG]]$temp, VAR_PRC, na.rm = TRUE ),
                          Boundary.knots = range( DATALIST_CALI[[iREG]]$temp, na.rm = TRUE ) );
  }else if( VAR_FUN == "bs" ){
    BASIS_VAR = onebasis( quantile( DATALIST_CALI[[iREG]]$temp, PRED_PRC, na.rm = TRUE ),
                          fun = VAR_FUN,
                          degree = VAR_DEG,
                          knots = quantile( DATALIST_CALI[[iREG]]$temp, VAR_PRC, na.rm = TRUE ) );
  }else{
    stop("ERROR: Invalid VAR_FUN !!!");
  }
  
  # Vector of Percentile Temperatures for the Predictions
  PRC_VAR = quantile( DATALIST_CALI[[iREG]]$temp, PRED_PRC, na.rm = TRUE );
  
  # Cumulative Exposure-Response without Centring
  suppressMessages( PRED_MORT <- crosspred( BASIS_VAR, coef = BLUP[[iREG]]$blup, vcov = BLUP[[iREG]]$vcov, model.link = "log", at = PRC_VAR ) );
  
  # Vector of Local Minima of the Cumulative Exposure-Response
  iMMT = 1 + which( diff( sign( diff( PRED_MORT$allRRfit ) ) ) == 2 );
  
  # Minimum Mortality Temperature
  if( bLOCAL_MINIMUM_MMT & length(iMMT) > 0 ){
    # If There Is at Least One Local Minimum, the Minimum Mortality Temperature Is the Local Minimum with the Lowest Relative Risk
    MMT_REG[iREG] = PRC_VAR[ iMMT[ which.min( PRED_MORT$allRRfit[iMMT] ) ] ];
  }else{
    # If There Are No Local Minima, the Minimum Mortality Temperature Is the Lowest Relative Risk within the Predefined Temperature Percentile Range [MIN_PMMT,MAX_PMMT]
    MMT_REG[iREG] = PRC_VAR[ which( PRED_PRC == MIN_PMMT ) - 1 + which.min( PRED_MORT$allRRfit[ which( PRED_PRC == MIN_PMMT ) : which( PRED_PRC == MAX_PMMT ) ] ) ];
  }
  rm(PRED_MORT,iMMT);
  
  # Cumulative Exposure-Response
  CROSS_PRED_REG_META[[iREG]] = crosspred( BASIS_VAR, coef = BLUP[[iREG]]$blup, vcov = BLUP[[iREG]]$vcov, model.link = "log", at = PRC_VAR, cen = MMT_REG[iREG] );
  rm(BASIS_VAR,PRC_VAR);
  
}
rm(iREG);

##################### read the temperature data for NUTS-3 regions in 35 countries

library(data.table)
####start to extract 5% and 95% for each nuts3
foldout = paste0( "/PROJECTES/ADAPTATION/proj/zchen/p20230314_ZC_weekly_comparison_epi/indata/popw_temp/" );
gpw_temp=fread( paste0(foldout,"poldata_gpw_temp_ma.csv"))


# Group by location and transform temperature from K to degree Celsius
df=na.omit(gpw_temp)
vREG=levels(factor(df$location))
# Split the data frame into a list based on the location column
DATALIST_CALI1 <- split(df, df$location)
# Keep only the data frames corresponding to the regions in vREG
DATALIST_CALI1 <- DATALIST_CALI1[vREG]
# Rename the list elements with the names from vREG
names(DATALIST_CALI1) <- vREG

################################################################################
### Calculation of the Pooled Coefficients
################################################################################
# Temperature Meta-Predictors: Temperature Average (TEMP_AVG), Temperature Inter-Quartile Range (TEMP_IQR)
TEMP_AVG = sapply( DATALIST_CALI1, function(x) mean( x$temp, na.rm = TRUE ) );
TEMP_IQR = sapply( DATALIST_CALI1, function(x)  IQR( x$temp, na.rm = TRUE ) );


print("");
print("= Calculation of the Pooled Coefficients =");
print("");

# Cumulative Exposure-Response
CROSS_PRED_NUT3_META = vector( "list", length(vREG));
names(CROSS_PRED_NUT3_META) = vREG;

# New Data for the Predictions 
NEW_DATA = data.frame( TEMP_AVG ,   
                       TEMP_IQR )  

# Pooled Coefficients of the Cumulative Exposure-Response
MULTIVAR_PRED = predict( MULTIVAR, NEW_DATA, vcov = TRUE, format = "list" )

# Multi-Location Temperature for Each Percentile
set.seed(13041975);
POOLED_TEMP_AVG = sapply( DATALIST_CALI1, function(x) quantile( jitter( x$temp ), PRED_PRC, na.rm = TRUE ) ) ;

crosspred_nuts_curve=function( POOLED_TEMP_AVG,VAR_FUN,MULTIVAR_PRED,REG){
  # Basis in the Variable Domain
  if( VAR_FUN == "ns" ){
    BASIS_VAR = onebasis( POOLED_TEMP_AVG,
                          fun = VAR_FUN,
                          knots = POOLED_TEMP_AVG[ paste0( 100 * VAR_PRC, ".0%" ) ],
                          Boundary.knots = range( POOLED_TEMP_AVG, na.rm = TRUE ) );
  }else if( VAR_FUN == "bs" ){
    BASIS_VAR = onebasis( POOLED_TEMP_AVG,
                          fun = VAR_FUN,
                          degree = VAR_DEG,
                          knots = POOLED_TEMP_AVG[ paste0( 100 * VAR_PRC, ".0%" ) ] );
  }else{
    stop("ERROR: Invalid VAR_FUN !!!");
  }
  
  # Cumulative Exposure-Response for the Calculation of the Minimum Mortality Temperature
  PRED_MORT = BASIS_VAR %*% MULTIVAR_PRED$fit;
  
  # Vector of Local Minima of the Cumulative Exposure-Response
  iMMT = 1 + which( diff( sign( diff( PRED_MORT ) ) ) == 2 );
  
  # Index of the Minimum Mortality Temperature
  if( bLOCAL_MINIMUM_MMT & length(iMMT) > 0 ){
    # If There Is at Least One Local Minimum, the Minimum Mortality Temperature Is the Local Minimum with the Lowest Relative Risk
    POOLED_CENT_IND = iMMT[ which.min( PRED_MORT[iMMT] ) ];
  }else{
    # If There Are No Local Minima, the Minimum Mortality Temperature Is the Lowest Relative Risk within the Predefined Temperature Percentile Range [MIN_PMMT,MAX_PMMT]
    POOLED_CENT_IND = which( PRED_PRC == MIN_PMMT ) - 1 + which.min( PRED_MORT[ which( PRED_PRC == MIN_PMMT ) : which( PRED_PRC == MAX_PMMT ) ] );
  }
  rm(PRED_MORT,iMMT);
  
  # Percentile of the Minimum Mortality Temperature
  POOLED_CENT_PRC = pmin( pmax( 100 * PRED_PRC[ POOLED_CENT_IND ], 0 ), 100 );
  rm(POOLED_CENT_IND);
  
  # Value of the Minimum Mortality Temperature
  if( is.wholenumber( POOLED_CENT_PRC ) ){
    TEMP_AVG_CEN = POOLED_TEMP_AVG[ paste0( POOLED_CENT_PRC, ".0%" ) ];
  }else{
    TEMP_AVG_CEN = POOLED_TEMP_AVG[ paste0( POOLED_CENT_PRC,   "%" ) ];
  }
  rm(POOLED_CENT_PRC);
  
  # Cumulative Exposure-Response
  CROSS_PRED_TOT_META= crosspred( BASIS_VAR, coef = MULTIVAR_PRED$fit, vcov = MULTIVAR_PRED$vcov, model.link = "log", at = seq(-20, 40, by = 0.5), cen = TEMP_AVG_CEN ); # <= ZHAO
  rm(MULTIVAR_PRED, POOLED_TEMP_AVG, BASIS_VAR, TEMP_AVG_CEN);
  
  # Vector of Local Minima of the Cumulative Exposure-Response
  iMMT = 1 + which( diff( sign( diff( CROSS_PRED_TOT_META$allRRfit ) ) ) == 2 );
  
  # Minimum Mortality Temperature
  if( bLOCAL_MINIMUM_MMT & length(iMMT) > 0 ){
    # If There Is at Least One Local Minimum, the Minimum Mortality Temperature Is the Local Minimum with the Lowest Relative Risk
    MMT_TOT = CROSS_PRED_TOT_META$predvar[                                             iMMT[ which.min( CROSS_PRED_TOT_META$allRRfit[iMMT] ) ]                                            ];
  }else{
    # If There Are No Local Minima, the Minimum Mortality Temperature Is the Lowest Relative Risk within the Predefined Temperature Percentile Range [MIN_PMMT,MAX_PMMT]
    MMT_TOT = CROSS_PRED_TOT_META$predvar[ which( PRED_PRC == MIN_PMMT ) - 1 + which.min( CROSS_PRED_TOT_META$allRRfit[ which( PRED_PRC == MIN_PMMT ) : which( PRED_PRC == MAX_PMMT ) ] ) ];
  }
  CROSS_PRED_TOT_META$REG=REG
  CROSS_PRED_TOT_META$MMT=MMT_TOT
  return(CROSS_PRED_TOT_META)}
#rm(COEF_MODEL,VCOV_MODEL, iMMT);

CROSS_PRED_LIS = vector( "list", length(vREG) );
names(CROSS_PRED_LIS) = vREG;
for( REG in vREG ){
  CROSS_PRED_LIS[[REG]] =crosspred_nuts_curve( POOLED_TEMP_AVG[,REG],VAR_FUN,MULTIVAR_PRED[[REG]],REG)}

####save the cross basis

foldout = paste0( "/PROJECTES/ADAPTATION/proj/zchen/P20230815_ZC_heat_AP_CEs/dataout/", sSEX, "_", sAGE, "_cross_pred/" );
if( !file_test( "-d", foldout ) ){ dir.create( foldout, recursive = TRUE ); }


saveRDS(CROSS_PRED_LIS,file=paste0(foldout,sSEX, "_", sAGE, "_cross_pred2003-2019_weekly.Rdata"))
