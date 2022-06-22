# PeakInt1D

Simple script for integrating 1D NMR spectra using lmfit.

Arguments are as follows;

```python proc.py
	-spectra PATH_TO_PROCESSED_DATA
	-output OUTPUT_FILE
	-vclist BRUKER_VC_OR_VDLIST
	-vc TIME_MULTIPLIER
	-initial PEAK_POSITION_GUESSES
	-figs FIG_FOLDER
```

*PATH_TO_PROCESSED_DATA* - path to directory containing processed dataset (eg 2rr files). NMRGlue requires that there by a 2ii and 2ri file present (even if this was recorded as QF in indirect). You can just copy the 2rr file to 2ii and 2ri to trick it into loading anyway.
*OUTPUT_FILE* - file to output the results. Columns are `Time, Background, Background Error, Peak 1, Peak 1 Error, ...`, with peaks ordered as in PEAK_POSITION_GUESSES.
*BRUKER_VC_OR_VDLIST* - path to a Bruker VC or VD list. Used for the time column in OUTPUT_FILE
*TIME_MULTIPLIER* - if a VC list is used, this can be used to set the time increment (eg if a chain of 50 ms pulses is used)
*PEAK_POSITION_GUESSES* - comma separated list of peak positions.
*FIG_FOLDER* - folder to output figures.


