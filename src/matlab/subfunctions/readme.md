# Processing ScouseTom Data - Subfunctions
These are functions called by the main functions like `ScouseTom_ProcessBV` or `ScouseTom_Load`. These mostly work only in context of those files, so are often unlikely to be of used outside them.

-   `ScouseTom_data_checkfirstinj.m` - Matches the first injection pair found to a line in the current injection protocol. Useful if the data was trunctated at the start.

-
