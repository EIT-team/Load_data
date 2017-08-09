# Load Data Matlab Source code
This folder is divided into the main functions of interest which are found in this directory, and the subfunctions which are not often called by their own.

Further, the functions are categorised by the corresponding EIT system: `ScouseTom_`, `KHU_` and `SwissTom_` etc.

## ScouseTom

1.  `ScouseTom_getHDR.m`






#### KHU EIT system
1.  `KHU_Load.m` Processes the data output from the KHU EIT system. See example usage in example folder

#### SwissTom System
1.  `SwissTom_EIT_reader` Loads the binary `.eit` file into matlab, provided by SwissTom AG
2.  `SwissTom_Load` Converts the output from the above into something a bit more usable and extracts information about the injection protocol etc.
