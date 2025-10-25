# Detailed Installation Guide

Installing `Reconstructor` consists of 2 (or 3) main steps:

1. Install `Reconstructor` with pip
2. (OPTIONAL) Install the DIAMOND sequence aligner
3. Run tests and download resource files

These steps are discussed in more detail below.

## Install Reconstructor python package

This can be done via pip in terminal on Mac, or a CMD.exe Prompt launched from
Anaconda Navigator in Widnows:

```shell
pip install reconstructor
```

> [!NOTE]
> *You must be running >= Python 3.9*

To determine your Python version you can use the following command:  

```shell
python --version
```

## Install DIAMOND (this step is optional)

If you want Reconstructor to use a version of DIAMOND other than the one that
would be automatically downloaded, you can install DIAMOND manually. See the
[DIAMOND GitHub page](https://github.com/bbuchfink/diamond) for details on how
to install DIAMOND. After installing DIAMOND, you must ensure that DIAMOND is
discoverable on your PATH. Otherwise, Reconstructor will download its own copy.

## Test suite (MUST RUN BEFORE USING RECONSTRUCTOR)

> [!CAUTION]
> YOU MUST RUN THE TEST SUITE BEFORE PROCEEDING TO USE RECONSTRUCTOR

Run the following test to ensure reconstruction was installed correctly and is
functional:

```shell
python -m reconstructor --test yes
```

This command also downloads a database file and a DIAMOND binary that are
necessary for reconstructor to work. This series of tests should take about an
hour to run, dependent on computer/processor speed. These are runtimes for
Reconstructor on a 2020 MacBook Pro with a 1.4 GHz Quad-Core Intel Core i5
processor. You can speed the test suite up somewhat with the optional `--cpu`
argument, which controls how many threads are used for blasting. For example,
the following command would run the test suite with 4 threads for blasting:

```shell
python -m reconstructor --test yes --cpu 4
```

If you have a copy of DIAMOND on your PATH when you run the test suite,
Reconstructor will default to using this copy of DIAMOND rather than
downloading its own. If you would instead like Reconstructor to download and use
its own copy of DIAMOND, you can use the `--diamond` option:

```shell
python -m reconstructor --test yes --diamond
```

This option also allows you specify a specific version of DIAMOND that you would
like Reconstructor to download and use (provided that version is available for
your specific system):

```shell
python -m reconstructor --test yes --diamond 2.1.13
```

If you do not have a version of DIAMOND in your PATH, but you do not want
Reconstructor to download its own version of DIAMOND, you can use the
`--skip-diamond` command:

```shell
python -m reconstructor --test yes --skip-diamond
```
