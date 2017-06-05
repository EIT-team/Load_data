# ScouseTom example data

As the EEG systems used in the [ScouseTom](https://github.com/EIT-team/ScouseTom) are 24bit with a higher sampling rate, the datasets can get very large (> 100 Gb in some cases). So distribution through git gets a bit difficult. So these are some very basic examples

The output of the ScouseTom software for a single recording is:

1. The voltages recorded by the EEG system. Either `###.bdf` for the Biosemi or `###.eeg ###.vmrk ###.vhdr` for the ActiChamp systems.
2. Accompanying matfile `###_log.mat` which stores the `ExpSetup` structure which contains all information about the configuration of the system, as well as storing the order of the frequencies for Multi-Frequency data (as its randomised each injection) and the Phase of the stimulation (as its random with respect to the phase marker)
3. Text file `###_log.txt` which stores all the communications with the system.
