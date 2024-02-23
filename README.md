# dsnd2sfz
Convert *.dsnd audio sample files of 2box to SFZ
---------------------------------------------
# DESCRIPTION

DSND file contains multiple audio samples of drum hits.
It is used by drum modules of company 2box to play electronic drums.
You can convert with this utility a single *.dsnd file into multiple 
WAVs with respective SFZ mapping to use them with LinuxSampler for example.

---------------------------------------------
# USAGE

dsnd2sfz -f dsnd-file-name [options]

   -a dry run write no files
   -b zone boundary volume difference (12 dB)
   -d output directory
   -h help
   -n note (36)
   -i in-threshold (60=-60 dB)
   -o out-threshold (96=-96 dB)
   -v number of velocity layers (5)
   -w window length in samples (200)

---------------------------------------------

# BUILD

gcc -l sndfile -lm -o dsnd2sfz dsnd2sfz.c

---------------------------------------------

# AUTHOR

Miroslav Kovac (mixxoo@gmail.com)
