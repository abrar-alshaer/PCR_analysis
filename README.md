Analysis of quantitative real-time PCR data
Code written by Abrar Al-Shaer and Eric Helfrich to analyze data generated by Kaylee Helfrich and others from the Susan Smith lab.

1. Introduction- The following manual provides a pipeline to efficiently analyze qPCR data from the raw Ct values to the final fold changes.

2. Expected Results- The expected results from this analysis are the delta Ct values (also called Cq values) for each sample normalized to the reference (needed for statistical analysis and individual data point plotting), standard deviations of the delta Ct values, overall fold changes compared back to a reference sample, SEMs, error ranges, and all other data needed to plot qPCR results. 

3. Methodology & Tools
    a. Required tools
        1. BASH
        2. Python (version 3), pandas, numpy, scipy
    b. Required program code files
        1. PCR_analysis.py
        2. PCR_functions.py
    c. For final analysis of data, see Kaylee Helfrich's GitHub page (kayleehelfrich) for R code to run statistics and graph the final results.
    
4. File Preprocessing Steps- these steps are not automated, because there are some subjective decisions
Note: This program details the instructions for data that is set up by the Bio-Rad CFX Maestro system. Data generated by other qPCR machines may have other formats, and this should be considered when inputting data.
1. Format cells in your Excel file to show ALL decimal places.
2. Save all Excel files as a CSV (using Save As to save a new file without overwriting the original). 
3. Filter 0's from the data (raw Ct's) by highlighting the column with the 0's, then click F5, then select Special, then select Blanks. Afterwards, all the 0's in the column will be highlighted, right click on one of them and select delete then delete entire row.
4. If you wish to remove an entire triplicate set of values be sure to replace all numbers in the set with zeroes, instead of deleting the rows completely. Also make sure to remove the comparison values (ex. of your control gene). 
    a. Example: If the triplicate set of the AIF gene for sample 2103 is bad (with numbers that aren’t close together or are weird), then replace all numeric values in these 3 rows along with the 3 associated reference gene rows with zeroes.
5. Filter any values with high standard deviations you do not want to include in the analysis. When the row with the outlier Cq value is removed make sure to:
    a. Recalculate the new mean among the triplicates
    b. Recalculate the new standard deviation among the triplicates
6. Both steps a) and b) can be done in R to save time and effort. 
    a. For recalculating mean: mean(c(X1, X2)) - X1 = your first number, X2 = your second number, and you may add more numbers in the parentheses if needed.
    b. For recalculating standard deviation: sd(c(X1, X2)) - X1 = your first number, X2 = your second number, and you may add more numbers in the parentheses if needed.
7. Next, isolate each reference and target pair into a separate files. This should only be done if there are multiple genes run on a single plate. For example, if you have GAPDH and AIF and IL6 all in one file, you should create a AIF+GAPDH file and a IL6+GAPDH file. 
    a. To do this, sort the entire file by the Target column and select the reference gene rows along with the first target’s gene rows and paste them into a separate file. MAKE SURE to include the header of the file in the same order in all files. Repeat for each reference and target pair. 
    b. Next, sort the newly created file (reference+target file) by the Sample column in order to restore the file to an alternating row pattern of target and reference.
8. Make sure to save all changes, and ensure that the Python files are in the same folder as the CSV files to be analyzed. 

5. Run Python programs
   1. Check both Python files and make sure that no changes need to be made.
   2. Open BASH and navigate to the correct folder.
   3. Run "python PCR_analysis.py"
   4. Follow the prompts and instructions, and enter the required information.
