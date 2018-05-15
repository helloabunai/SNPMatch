SNPMatch: Split PED files into individual chromosomes while maintaining SNP order.
=========================================================

I still need to write a more detailed description of the program haha.
SNPMatch takes PED/Map/FinalReport files from a GIGAMUGA genotyping array, and splits the results for each sample in the *.PED into individual chromosome files.

Usage
=====
General usage is as follows:

    $ snpmatch [-h/--help] [-v] [-p PED FILE] [-r REPORT FILE] [-m MAP FILE] [-t THREADS] [-o OUTPUT]
    e.g.
    $ snpmatch -v -p ~/path/to/pedfile.ped -r ~/path/to/reportfile.txt -m ~/path/to/mapfile.map -t 12 -o ~/path/to/desired/outputfolder

SNPMatch flags are as follows:

    -h/--help:: Simple help message explaining flags in detail
    -v/--verbose:: Enables verbose mode in the terminal (i.e. shows user feedback in the terminal)
    -p/--ped:: A path to your desired *.PED data file, for processing [filepath].
    -r/--report:: A path to your x_FinalReport.txt [filepath].
    -m/--map:: A path to the desired *.MAP data file, to be processed with the above *.PED [filepath].
    -t/--threads:: Number of CPU cores to use in functions that are multi-process compliant [integer].
    -o/--output:: Desired output folder. SNPMatch will create this folder if it does not already exist [directory].
