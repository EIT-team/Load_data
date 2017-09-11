# Load Data Matlab Source code
This folder is divided into the main functions of interest which are found in this directory, and the subfunctions which are not often called by their own.

Further, the functions are categorised by the corresponding EIT system: `ScouseTom_`, `KHU_` and `SwissTom_` etc.

## ScouseTom

-   `ScouseTom_Load.m` - **The most important function!** This reads the infomation stored in the metadata of the EEG file and reads the trigger channels, and then calls the relevant function to demodulate the data. The result is the demodulated magnitude and phase of each measurement across all frames in the dataset.
      - **Conventional EIT recordings** This function calls `ScouseTom_ProcessBV` which
      - **Contact Impedance Checks** This function calls `ScouseTom_ProcessZ` which, after demodulation, estimates the contact impedance of the electrodes based on the voltage on the injection channels. This produces graphs so the user can identify bad channels.
      Example graph to follow soon.


-   `ScouseTom_ProcessBatch.m` - Processes all files within a directory in sequence.  

##### Demodulation
These functions are used during demodulation of the voltages recorded in the EEG files.

-   `ScouseTom_FindFilterSettings.m` - For a given recording, find the optimal filter for use before demodulation. This is based on the carrier frequency, and length of injection. Uses IIR filters where possible, but for short injections it will use FIR.


##### Digital Trigger channels and meta data

-   `ScouseTom_getHDR.m` - This reads (and corrects) the HDR structure created by the `biosig` library, which contains all the metadata in files saved by the EEG systems.

-   `ScouseTom_TrigView.m`- Plots the changes in values on all trigger channels dataset, in something resembling the `strips` function in Matlab. Useful in debugging:

![Trig view example](https://raw.githubusercontent.com/EIT-team/Load_data/master/resources/example_figures/ex_trigview.png)

-   `ScouseTom_TrigReadChn.m` - Reads the information on the digital channels and identifies the starting trigger codes. Also rejects spurious too short pulses. Different EEG systems store this information differently, so the output is common format used in `ScouseTom_TrigProcess`

-   `ScouseTom_TrigProcess.m` - Processes the information encoded in the digital triggers - Conventional/Contact check, single or Multi-Frequency, stimulations, length of injection protocol etc. This is then stored in the `TT` structure which is used for `ProcessBV`












#### KHU EIT system
1.  `KHU_Load.m` Processes the data output from the KHU EIT system. See example usage in example folder

#### SwissTom System
1.  `SwissTom_EIT_reader` Loads the binary `.eit` file into matlab, provided by SwissTom AG
2.  `SwissTom_Load` Converts the output from the above into something a bit more usable and extracts information about the injection protocol etc.
