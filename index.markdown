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

RNA-seq Graph Builder is a method to reconstruct the Isoform Graph of a gene from RNA-seq
data, without the genome information, where such a graph is a representation
of the variants of alternative splicing of the gene structure.


by [Stefano Beretta](http://bimib.disco.unimib.it/index.php/Beretta_Stefano)


Current release: **2.0.0** (November 3, 2011)

Direct download: [zip](https://github.com/AlgoLab/RNA-seq-Graph/zipball/v2.0.0) - [tar.gz](https://github.com/AlgoLab/RNA-seq-Graph/tarball/v2.0.0)


------------------------------------------------------------------------



## Introduction ##

This program predicts from NGS data the gene structure induced by the different 
full-length isoforms due to alternative splicing. More precisely, it analyzes 
RNA-seq data that have been sampled from the transcripts of a gene, with the goal
of building a graph representation of the variants of alternative splicing
corresponding to those full-length isoforms. The novelty of this method relies on
the fact that it builds such a graph in absence of the genome. A subsequent step 
is to efficiently map the graph to the genome in order to refine the gene structure 
and to compute intron data that, obviously, cannot be present in the isoforms.



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

The source code is available directly in [zip](https://github.com/AlgoLab/RNA-seq-Graph/zipball/v2.0.0) 
or [tar.gz](https://github.com/AlgoLab/RNA-seq-Graph/tarball/v2.0.0)
    
Or, if you have a GitHub account, you can fork the project from the
[repository web page](https://github.com/AlgoLab/RNA-seq-Graph).


### Compilation ###

The program can be compiled by issuing the command at the command
prompt:

    $ make


## Usage ##

The program takes as input a FASTA file with the RNA-seq data of a gene
and returns the RNA-seq Graph. (The file formats are described below.)
The program is executed by typing the following command:

    $ ./bin/build_RNA_seq_graph [options] --reads <RNA-seq_file>

where the possible options are:

    --ref_level {1-5}

 `1` - Standard Algorithm (Default option)

 `2` - Add tiny blocks

 `3` - Add linking edges

 `4` - Add small blocks

 `5` - Refine overlapping nodes

    -o <graphML_out_file> (Default: std output)

    --advanced (For debug)

    --help (Print this screen)

Alternatively it is possible to view the debug options by typing:

    $ ./bin/build_RNA_seq_graph --advanced

that prints the following options:

    --ref_level {1-5}

 `1` - Standard Algorithm (Default option)  

 `2` - Add tiny blocks  

 `3` - Add linking edges  

 `4` - Add small blocks  

 `5` - Refine overlapping nodes  

    -o <graphML_out_file> (Default: std output)

    --debug {1...8}

 `1` - Print left hash table  

 `2` - Print right hash table  

 `3` - Print unspliced RNA-seq sequences  

 `4` - Print spliced RNA-seq sequences

 `5` - Print half spliced RNA-seq sequences

 `6` - Build chains of unspliced reads with half sequence overlap

 `7` - Build chains of unspliced reads with specific overlap

 `8` - Merge chains built of unspliced reads with half sequence overlap

An example of usage is:

    $ ./bin/build_RNA_seq_graph --reads Raw_file.fa -o Out_file

A summary of the available program options can be printed by invoking
`./bin/build_RNA_seq_graph` without parameters.


## File Formats ##

### Input: RNA-seq Dataset ###

RNA-seq file is in FASTA format (.fa, .fas or .fasta).
In particular, each line of the input file describes a single
read and it is composed by at least 2 rows: the first one is the header of
the read (that usually starts with '>') and the second one that contains the
sequence. For example:

    >B2GA004:3:100:1143:752#0/1 /source=region /gene_id=gene /gene_strand=+
    GATGAAATACTACTTCTACCATGGCCTTTCCTGGCCCCAGCTCTCTTACATTGCTGAGGACGAGAATGGGAAGAT

### Output: RNA-seq Graph ###

The program produces as output a file in _txt_ format that contains
a list of nodes and arcs of the RNA-seq graph. It also gives as output 
the same graph in _GDL_ format (<http://www.absint.com/aisee/manual/windows/node58.html>).
It also print on standard output the graph in GraphML format; this latter can be 
redirected into a file in order to visualize or export it.
By default the files are _RNA-seq-graph.txt_ and _RNA-seq-graph.gdl_ .
For example:

    $ ./bin/build_RNA_seq_graph --reads Raw_file.fa > RNA-seq-graph.graphml

If the option `-o` is specified the 3 files will have the specified name:

    $ ./bin/build_RNA_seq_graph --reads Raw_file.fa -o Out-file

generates _Out-file.txt_, _Out-file.gdl_ and _Out-file.graphml_.

## License ##

RNA-seq-Graph-Builder is released under the terms of the GNU General Public License
(GPL) as published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

RNA-seq-Graph-Builder is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

Please refer to file `COPYING` or to the
[GNU website](http://www.gnu.org/licenses/) for a copy of the GNU
General Public License.


## Contacts ##

Please contact *Stefano Beretta* for additional information.  
E-mail:   <beretta@disco.unimib.it>  
Web page: <http://bimib.disco.unimib.it/index.php/Beretta_Stefano>


