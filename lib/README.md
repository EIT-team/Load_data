# Biosig library
This is a modified version of the [biosig](http://http://biosig.sourceforge.net/) library. The main modification is that it edited to use less memory when creating the `HDR` structure. This was necessary because we can sometimes have 10,000 x the digital triggers than a normal EEG recordings.

### installation
In Matlab run `biosig_installer.m` then `savepath`. Job done
