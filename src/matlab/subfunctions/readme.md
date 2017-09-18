# Processing ScouseTom Data - Subfunctions
These are functions called by the main functions like `ScouseTom_ProcessBV` or `ScouseTom_Load`. These mostly work only in context of those files, so are often unlikely to be of used outside them.

-   `ScouseTom_data_checkfirstinj.m` - Matches the first injection pair found to a line in the current injection protocol. Calls `ScouseTom_data_EstInjPair.m`. Useful if the file was missing the start.

-   `ScouseTom_data_DemodHilbert.m` - Filters and demodulates the data using the hilbert transfrom.

-   `ScouseTom_data_EstInjPair.m` - Estimates the pair of electrodes used for current injection from the largest two voltages.

-   `ScouseTom_data_estZ.m` - Estimates the contact impedance during conventional EIT recordings based on the voltages on the injecting electrodes. If the electrode is used multiple times in a protocol, then the result is averaged.

-   `ScouseTom_data_findprt.m` - Creates the `keep_idx` for valid measurements, and produces the full measurement protocol i.e. `[1 2 1 32; 1 2 2 32]` etc.

-   `ScouseTom_data_getBV.m` - Takes the mean of the demodulated voltage within each injection, to give the standard "boundary voltage".

-   `ScouseTom_data_GetCarrier.m` - Estimates the carrier frequency by finding the largest spectral component in signal

-   `ScouseTom_data_GetFilterTrim.m` - Finds the appropriate bandpass filter given the carrier frequency and injection time. Uses IIR filters if injections are long enough, before switching to FIR. Checks the gain at the carrier frequency is not adversely effected.  
