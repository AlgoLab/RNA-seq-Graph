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


Current release: **1.0.0** (May 19, 2011)


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

### Download ###

RNA-seq-Graph-Builder is developed on the `AlgoLab/RNA-seq-Graph` Git repository hosted by
GitHub.
The repository can be explored using the GitHub web interface at
<https://github.com/AlgoLab/RNA-seq-Graph>.

It is also possible to clone the entire repository using the following
command:

    $ git clone git://github.com/AlgoLab/RNA-seq-Graph.git
    
Or, if you have a GitHub account, you can fork the project from the
[repository web page](https://github.com/AlgoLab/RNA-seq-Graph).


### Compilation ###

The program can be compiled by issuing the command at the command
prompt:

    $ make


## Usage ##

The program takes as input a FASTA file with the RNA-seq data of a gene
and returns the RNA-seq Graph. (The file formats are described below.)
To program is executed by typing the following command:

    $ ./bin/build_RNA_seq_graph <fasta_file> <option 1...9>

where the possible options are:

`1` - Print left hash table

`2` - Print right hash table

`3` - Print unspliced RNA-seq sequences

`4` - Print spliced RNA-seq sequences

`5` - Print perfectly spliced RNA-seq sequences

`6` - Build chains of unspliced reads with half sequence overlap

`7` - Build chains of unspliced reads with specific overlap

`8` - Merge chains built of unspliced reads with half sequence overlap

`9` - Build the RNA-seq Graph

An example of usage is:

    $ ./bin/build_RNA_seq_graph Raw_file.fa 9

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

The program produces as output a file _RNA-seq-graph.txt_ that contains
a list of nodes and arcs of the RNA-seq graph. It also gives on the
standard output the same graph in GDL format (<http://www.absint.com/aisee/manual/windows/node58.html>).
This latter can be redirected into a file in order to visualize or export it in another format.
For exmaple:

    $ ./bin/build_RNA_seq_graph Raw_file.fa 9 > RNA-seq-graph.gdl

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


