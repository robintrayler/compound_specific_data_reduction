# Readme

This repository contains the data reduction code used by the [Stable Isotope Ecosystem Lab of UC Merced](https://isotopes.ucmerced.edu).
## Setup
1. Open a terminal window and navigate to the appropriate folder. 

```bash
# navigate to the GCC folder
cd ~/Box\ Sync/Data\ Repository/GCC
```

2. Define the run ID for the data you are going to correct. 

```bash
# define run ID
run_ID="GCC19850910" 
```

3. clone the GitHub repository containing the code into a new folder with the `run_ID` as the name. 

```bash
# clone the data reduction repository
git clone https://github.com/robintrayler/compound_specific_data_reduction $run_ID

# navigate into the new folder
cd $run_ID

# rename the data reduction script
mv data_reduction.R $run_ID\_data_reduction.R
```

4. Paste the raw isodat csv file into your newly created folder. The file should be named with the `run_ID`. 
## Use
1. In your newly created folder, open the `.Rproj` file. This should open RStudio. 
2. Open the `data_reduction.R` script.
3. Edit line 16 of the script with the correct file name for your raw isodat csv file. 
```r
# file path -------------------------------------------------------------------
file_name = './GCC20221012.csv' # put your file name here
```
4. Run the script line by line. Make sure to read the warnings and inspect each plot before moving on. **Do not just run the script all at once without paying attention**. 
5. The script should output two files. 
   * `run_ID_corrected.csv` has the corrected isotope compositions in it. 
   * `run_ID_summary.csv` has summary statistics for all samples in it. 