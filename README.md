# CTV-batch-mode-Windows
CTV / ToxValue.org Random Forest Model - batch mode code for Windows

This code allows one to run the ToxValue.org code in batch mode on a Windows machine (get_ad.exe is a Windows executable).  An example input file with CAS numbers, SMILES, and molecular weights is included.  The R script implements the code.

Updated with a new version "Predict_new_chemicals_batch-Ronly.R" that does not rely on Windows executable. The new code get_ad.R implements the applicability domain calculation in R.

Note: One needs to install the original version of rcdk and rcdklibs that were used in model building in order to make valid predictions. The tgz files for those versions are included in "orig packages tgz" folder.
