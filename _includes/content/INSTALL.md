---
layout: default
title: RNA-seq Graph Builder
nav:
  - name: Introduction
    link: "#introduction"
  - name: Installation
    link: "#download_and_installation"
  - name: Usage
    link: "#usage"
  - name: File Formats
    link: "#file_formats"
  - name: License
    link: "#license"
  - name: Contacts
    link: "#contacts"
---

  RNA-seq Graph Builder
==========


## Download and Installation ##

RNA-seq-Graph-Builder is currently distributed only on source form.
It has been developed on Ubuntu Linux machines (v. 10.04 and
10.10) and has been tested on both 32 and 64 bit.
The program requires the C++ library SEQAN available at <http://www.seqan.de>
or it is possible to install the develop package seqan-dev by typing:

    $ sudo apt-get install seqan-dev

### Download ###

RNA-seq-Graph-Builder is developed on the `AlgoLab/RNA-seq-Graph` Git repository hosted by
GitHub.
The repository can be explored using the GitHub web interface at
<https://github.com/AlgoLab/RNA-seq-Graph>.

It is also possible to clone the entire repository using the following
command:

    $ git clone git://github.com/AlgoLab/RNA-seq-Graph.git

The source code is available directly in [zip](https://github.com/AlgoLab/RNA-seq-Graph/zipball/master)
or [tar.gz](https://github.com/AlgoLab/RNA-seq-Graph/tarball/master)

Or, if you have a GitHub account, you can fork the project from the
[repository web page](https://github.com/AlgoLab/RNA-seq-Graph).


### Compilation ###

The program can be compiled by issuing the command at the command
prompt:

    $ make

There is also the possibility to compile with the `low_mem` option in order to reduce the memory consumption.
In this way the memory occupation (i.e. the heap peak) is reduced by ~35% but the time required is increased by ~10%.
Starting from a cleaned repository, the command is:

    $ make low_mem






