# Detailed Installation Guide

Installing `Reconstructor` consists of 2-3 main steps:

1. Install DIAMOND aligner (Mac users)
2. Install `Reconstructor` with pip
3. Run tests and download resource files

These steps are discussed in more detail below.

## Download DIAMOND Aligner

### Downloading DIAMOND on MacOS

You must first have the DIAMOND sequence aligner downloaded (**MUST BE >=
VERSION 2.0.15**). Installation instructions can be found on the [DIAMOND Github
page](https://github.com/bbuchfink/diamond).

If you think you already have DIAMOND installed, you can check your version
using the terminal command:

```shell
diamond --version
```

You can also download DIAMOND through homebrew (if you have homebrew installed)
using the terminal command:

```shell
brew install diamond
```

### DIAMOND on Windows

You do not need to install DIAMOND on your windows machine. A Windows executable
function is already pre-packaged within the Reconstructor software.

## Install Reconstructor python package

This can be done via pip in terminal on Mac, or a CMD.exe Prompt launched from
Anaconda Navigator in Widnows:

```shell
pip install reconstructor
```

> [!NOTE]
> *You must be running >= Python 3.8*

To determine your Python version you can use the following command:  

```shell
python --version
```

## Test suite (MUST RUN BEFORE USING RECONSTRUCTOR)

Run the following test to ensure reconstruction was installed correctly and is
functional:

```shell
python -m reconstructor --test yes
```

This command also downloads database files that are necessary for reconstructor
to work. This series of tests should take about an hour to run, dependent on
computer/processor speed. These are runtimes for Reconstructor on a 2020 MacBook
Pro with a 1.4 GHz Quad-Core Intel Core i5 processor.

MAC USERS MAY BE ASKED FOR TERMINAL TO HAVE ACCESS TO DOWNLOADS, CAMERA,
LOCATION, ETC. Please allow terminal to have access to all locations on your
computer. Reconstructor will NOT gather data from your camera, location, or
other sensitive information. Reconstructor is simply searching for the file
titled glpk_interface.py on your local machine (installed when COBRA module is
installed) and replacing it with a newer, functional version.

Use the command below to test reconstructor to ensure correct installation:

*YOU MUST RUN THE TEST SUITE BEFORE PROCEEDING TO USE RECONSTRUCTOR
