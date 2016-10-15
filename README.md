June 8 Version of ThysimModel
Author: Rukan Shao

Author of this README: Alan Chen

This readme is a new document that I have appended to the ThysimModel folder for future reference.  Due to its informal nature there may be parts that are explained poorly.  For example, all line numbers are relative to Rukan's code from 3/8/2016, NOT the code included in this directory.  I have attempted to remedy this by placing comments containing the word "STEP" at appropriate places.  Alan Chen can be contacted at apcucla@g.ucla.edu should it be necessary to contact someone for clarification.

ThysimModel is merely a research testbed and was not built with adaptation for applications in mind.  For a version of the model meant for such a platform, see Simon Han's Thyrosim web application (biocyb1.cs.ucla.edu/thyrosim) and source (bitbucket.org/DistefanoLab/thyrosim/overview).  

The code in this directory has been intentionally stripped down to reduce confusion.  However, this also means that, relative to the code actually produced by Rukan, it lacks some alternative functions and utilities.  Because this readme is intended only for determining how to run the parameter fit, I will not comment on how to use these other files, nor do I know how to use them.  


TABLE OF CONTENTS
1. Running the model for simulation
2. Running a parameter fit
3. Running a custom parameter fit
4. Overview of important files
5. Special notes



1. Running the model for simulation
  1. Open ThysimModel.m
  2. Ensure that variable "searchMode" is equal to 0 (STEP 0)
    2a. If the fit for data other than the blue data is desired, set "fitIndex" to the corresponding number
  3. Run ThysimModel in MATLAB console (open Matlab first if you haven't already)
  4. Screen should show graph of thyroid hormone concentrations in serum, slow pool, and fast pool

2. Running a parameter fit
  1. Open ThysimModel.m
  2. Ensure that variable "searchMode" is equal to 2 (NOTE: searchMode 1 is currently not considered accurate and will probably not work, but searchMode 3 will) (STEP 0)
    2a. If the fit for data other than the blue data is desired, set "fitIndex" to the corresponding number
  3. Select parameters to search over via the "searchPoints" variable
    3a. 5 will search over all parameters which determine flow size
    3b. 11 will search over all parameters which 5 searches over and will search over parameters 
  4. Run ThysimModel in MATLAB console (open Matlab first if you haven't already)
  5. Matlab should spend a few minutes computing fit.  Note that if convergence fails, this may take up to half an hour
   5a. Current cost of fit is displayed each iteration, along with the current parameters being tested
  6. At program end, Matlab will display graph of thyroid hormone concentrations with newly found parameters

3. Running a custom parameter fit
  1. Open ThysimModel.m
  2. Ensure that variable "searchMode" is equal to 2 (STEP 0)
    2a. If the fit for data other than the blue data is desired, set "fitIndex" to the corresponding number
  3. Add lines under both "if fitIndex ==" statements (located at lines 44 and 96 as of this writing) describing a value for searchPoints other than 5 or 11 (STEP 1)
    3a. Format as (where "%yournumber" is whatever number for searchPoints you want):

if searchPoints == %yournumber
  x = [num1 num2 num3 ...];
end

  4. Add lines under "if x ~= 0" located at line 147 (STEP 2)
    4a. Format as:

elseif searchPoints == %yournumber
  param1 = abs(x(1));
  param2 = abs(x(2));
  ...

  5. At the code block starting at line 182 (STEP 3), add the following lines

elseif searchPoints == %yournumber
  startPoint = [param1 param2 ...];

  5. Open CostFunction.m
  6. Add lines under "if searchPoints ==" located at line 27 (STEP 5)
    6a.  Format as (note that this is the same as the format in 4a, just with "input" instead of "x"):

elseif searchPoints == %yournumber
  param1 = abs(input(1));
  param2 = abs(input(2));
  ...

  7. In ThysimModel.m, set searchPoints == %yournumber
  8. Run thyrosim as described in section 2, instruction 4 onwards

4. Overview of important files
  ThysimModel - main function
  CostFunction - generates cost of current fit based on input; search functions will use this function to generate values to be minimized
  InitializeParameters - sets default values for all parameters and initial conditions
  InitializeDataPoints - creates arrays of recorded data values for red and blue fits (other fits not currently used)
  ODEs.m - Differential equations which describe model
  PlotData - plots data points from specified fit
  PlotExperiment - plots thyroid hormone fit according to discovered or preexisting parameter values, depending on whether or not the fit was performed


5. Special notes:

To check returnCost, check the variable in the workspace outliner.  The program will display the two components of returnCost (cost for T4 fit and cost for T3 fit), but the actual returnCost variable will be set to the sum of those two components

To continue a run further, just add another "x = [num1 num2 ...]" line under the one you added under the STEP 1 tag, this time with the values from x in the workspace outliner







