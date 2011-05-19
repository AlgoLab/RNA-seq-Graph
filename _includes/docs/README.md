
  RNA-seq Graph Builder
==========

A SAT-based program to compute a haplotype configuration on pedigrees
with recombinations, genotyping errors, and missing genotypes.

by [Yuri Pirola](http://bimib.disco.unimib.it/index.php/Pirola_Yuri)


Started: September 27, 2010  
Current release: **1.0.0** (May 13, 2011)


------------------------------------------------------------------------



## Introduction ##

This program is based on a reduction of the _Haplotype Configuration
with Recombinations and Errors_ problem to _Boolean Satisfiability_, which
is then solved by a SAT solver.
A haplotype configuration is finally recovered from the satisfying
assignment.



## Download and Installation ##

reHC-* is currently distributed only on source form.
It has been developed on Ubuntu Linux machines (v. 10.04 and
10.10) and has been tested on both 32 and 64 bit.
The program should work on (or should be easily ported to) on MacOS X
but has not been tested and it is not supported on this operating
system.


### Download ###

reHC-* is developed on the `yp/reHCstar` Git repository hosted by
GitHub.
The repository can be explored using the GitHub web interface at
<https://github.com/yp/reHCstar>.

The latest stable version of reHC-* can be downloaded in either
[.tar.gz](https://github.com/yp/reHCstar/tarball/master) or in
[.zip](https://github.com/yp/reHCstar/zipball/master) format.
Previous stable releases can be downloaded from
<https://github.com/yp/reHCstar/archives/master>.

It is also possible to clone the entire repository using the following
command:

    $ git clone git://github.com/yp/reHCstar.git
    
Or, if you have a GitHub account, you can fork the project from the
[repository web page](https://github.com/yp/reHCstar).


### Dependencies ###

- CMake (>= 2.8)
- GNU make
- Boost FileSystem, System, DateTime, ProgramOptions, IOStreams, and
  other include-only libraries (tested with 1.42)
- Apache Log4cxx (tested with 0.10.0)
- optional: Google C++ Testing Framework (tested with 1.30)


### Compilation ###

The program can be compiled by issuing the command at the command
prompt:

    $ make STATUS=Release

The program can be compiled in three different variants:

1. **with** an integrated SAT solver and **without** the ability to
   invoke an external SAT solver.
   This is the default variant, and should be the most efficient on both
   time and memory.
2. **with** an integrated SAT solver and **with** the ability to invoke
   an external SAT solver.
   This is the most flexible variant, but it is also the variant that
   uses more time and memory than the others.
3. **without** an integrated SAT solver but **with** the ability to
   invoke an external SAT solver.
   This variant should be preferred when one wants to use an external
   SAT solver, since it is more memory- and time-efficient than the
   previous one.

Currently, the SAT solvers that can be directly integrated into reHC-*
are [CryptoMiniSat v2.9.1](http://gitorious.org/cryptominisat) by Mate
Soos and [MiniSat v2.2.0](http://www.minisat.se/MiniSat.html) by Niklas
Een and Niklas Sorensson.
reHC-* can interact (in the 2nd and 3rd variants) with any SAT solver
that follow the standard
[requirements](http://www.satcompetition.org/2004/format-solvers2004.html)
of the last SAT competitions.
The variant can be specified by modifying the file `CMakeOptions.txt` in
the root directory.
In particular, two options have to be used to specify the variant:

-  `INTEGRATE_SAT_SOLVER`, which specifies if the SAT solver should be
   integrated into reHC-*;
-  `DISABLE_EXTERNAL_SAT_SOLVERS`, which specifies if the invocation of
   external SAT solvers is allowed.

The following combinations are allowed:

-  `INTEGRATE_SAT_SOLVER=ON` and `DISABLE_EXTERNAL_SAT_SOLVERS=ON`,
   (**default**), which corresponds to the **first** variant;
-  `INTEGRATE_SAT_SOLVER=ON` and `DISABLE_EXTERNAL_SAT_SOLVERS=OFF`,
   which corresponds to the **second** variant;
-  `INTEGRATE_SAT_SOLVER=OFF` and `DISABLE_EXTERNAL_SAT_SOLVERS=OFF`,
   which corresponds to the **third** variant.

The SAT solver that will be integrated (if `INTEGRATE_SAT_SOLVER` is
`ON`) can be specified by setting `USE_CRYPTOMINISAT` or `USE_MINISAT`
to `ON` in the file `CMakeOptions.txt`.

By default, reHC-* generates SAT instances in a augmented CNF format
where the "extended clauses" can also be XORs of literals.
The SAT instance is generally smaller if XORs are allowed,
represents better the "internal structure" of the Boolean formula, and
_should be used_ if it is supported by the SAT solver.
The SAT solver embedded by default, CryptoMiniSat, supports this
functionality, while MiniSat does not.
If MiniSat is chosen instead of CryptoMiniSat by editing the file
`CMakeOptions.txt`, then XOR-clauses are automatically disabled.
To disable the augmented CNF format, reHC-* _must be rebuilt_ specifying
the preprocessor symbol `AVOID_XOR_CLAUSES` in the `CXXFLAGS`.
For example, if the program is compiled in a bash shell in Linux, it
suffices the following command to enable the "pure" CNF format:

    $ CXXFLAGS="-DAVOID_XOR_CLAUSES" make STATUS=Release



## Usage ##

The program takes as input a genotyped pedigree (with missing genotypes)
and returns (if possible) a complete haplotype configuration with at
most *r* recombinations and *e* errors.
(The file formats are described below.)
Depending on the variant that has been compiled, the program works in
four different modes that have to be specified on the command line as
program parameter:

1. `--create` (short form `-1`), that, given a genotyped pedigree,
   creates the associated SAT instance.
   (Available only on variants *2* and *3*.)
2. `--read` (short form `-2`), that, given a genotyped pedigree,
   reads a satisfying model of the associated SAT instance (if such a
   model exists) and computes the associated haplotype configuration.
   (Available only on variants *2* and *3*.)
3. `--create-read` (short form `-3`), that, given a genotyped pedigree,
   creates the associated SAT instance, invokes the external SAT solver,
   reads a satisfying model of the SAT instance (if such a model
   exists), and computes the associated haplotype configuration.
   This mode essentially combines the previous two modes by
   automatically invoking the external SAT solver.
   (Available only on variants *2* and *3*.)
4. `--solve-internal` (short form `-4`), that, given a genotyped pedigree,
   creates the associated SAT instance, uses the integrated SAT solver
   for solving the instance, and, if the SAT instance is satisfiable,
   computes the associated haplotype configuration.
   (Available only on variant *1*.)

The following options are used to specify the input/output files:

-  `--pedigree` (short form `-p`), that specifies the file containing the
   genotyped pedigree (input file);
-  `--sat` (short form `-s`), that specifies the file containing the
   SAT instance associated with the genotyped pedigree (output file);
-  `--result` (short form `-r`), that specifies the file containing the
   results computed by the external SAT solver for the SAT instance
   associated with the genotyped pedigree (input file);
-  `--haplotypes` (short form `-h`), that specifies the file that will
   contain the haplotype configuration of the genotyped pedigree
   computed by reHC-* (output file).

For the `--create-read` mode, the command-line that has to be used to
invoke the external SAT must be specified by using the `--sat-cmdline`
(short form `-c`) program option.
The strings `%%INPUT%%` and `%%OUTPUT%%` are placeholders for,
respectively, the input and the output files of the SAT solver.


### Options for Recombinations and Errors ###

By default, reHC-* search for a haplotype configuration with zero
recombinations and zero errors.
To enable recombinations in the haplotyping process, the program options
`--global-recomb` and `--global-recomb-rate=XX` _must be specified_,
where `XX` is a number between `0.0` and `1.0` that represents the
maximum number of recombinations *r* as a fraction of the total number
of possible recombination loci.

Similarly, to enable genotyping errors in the computed haplotype
configuration, the program options `--global-error` and
`--global-error-rate=XX` _must be specified_, where `XX` is a number
between `0.0` and `1.0` that represents the maximum number of errors *e*
as a fraction of the number of non-missing genotypes.

Other program options allow a finer control over the distribution of
recombinations and errors.
Please refer to the help of the program (that can be obtained by
specifying the `--help` program option) for their presentation and
explanation.


### Other Options ###

reHC-* can also read and write files compressed by GZip.
The GZip compression allows to save some space and, especially for large
instances and when an external SAT solver is used, it could reduce the
running time, since it greatly reduces to time spent for I/O operations.
It is disabled by default since not all the SAT solvers support it.
Three options regulates the GZip compression:

-  `--compress-input`, which enables the GZip compression of some files
   that are read by reHC-* (currently only the `--pedigree` file);
-  `--compress-output`, which enables the GZip compression of some files
   that are written by reHC-* (currently the `--sat` and `--haplotypes`
   files);
-  `--compress` (short form `-z`), which is equivalent to specify both
   `--compress-input` and `--compress-output`.

Temporary files of the `--create-read` mode are automatically removed by
default.
To keep them (for example, for manual inspection), the program option
`--keep` (short form `-k`) has to be specified.

A summary of the available program options can be printed by invoking
reHC-* with the `--help` (short form `-?`) option.


### Example ###

For example, if the genotyped pedigree is described in file
`genotyped-pedigree.txt`, the following commands perform the complete
haplotype inference process (saving the resulting haplotype
configuration in file `haplotype-configuration.txt`).

Using the integrated SAT solver (variant *1* or *2*):

    $ ./bin/reHCstar -4  \
          -p genotyped-pedigree.txt  \
          -h haplotype-configuration.txt

Using an external SAT solver (variant *2* or *3*) with _manual_
invocation of the SAT solver:

    $ ./bin/reHCstar -1  \
          -p genotyped-pedigree.txt  \
          -s instance.cnf
    # ...execution of the external SAT solver, assuming that
    #    it writes the results in file sat-result.txt
    $ ./bin/reHCstar -2  \
          -p genotyped-pedigree.txt  \
          -r sat-result.txt  \
          -h haplotype-configuration.txt

Using an external SAT solver (variant *2* or *3*) with _automatic_
invocation of the SAT solver:

    $ ./bin/reHCstar -3  \
          -p genotyped-pedigree.txt  \
          -h haplotype-configuration.txt  \
          -c "./external-sat-solver %%INPUT%% %%OUTPUT%%"



## File Formats ##

### Input: Genotyped Pedigrees ###

Genotyped pedigrees are described by a single file with the standard PED
format used in
[`plink`](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped).
In particular, each line of the pedigree file fully describes a single
individual and it is composed by at least _six_ whitespace-separated
fields.
The first (mandatory) six fields are:

-  `Family ID` (numeric only)
-  `Individual ID` (numeric only, greater than `0`)
-  `Paternal ID` (the ID of the father, `0` if unknown/not present)
-  `Maternal ID` (the ID of the mother, `0` if unknown/not present)
-  `Sex` (`1`=male, `2`=female)
-  `Phenotype` (ignored, could be any string not containing a whitespace)

**Remark:**
reHC-* currently works only on single-family pedigrees, thus the
`Family ID`s _must be_ the same for all the individuals.

The remaining fields (field 7 onwards) represent the genotype of the
individual, where each field represents a single allele of a single SNP
**biallelic** locus.
Both the alleles of each locus _must be_ specified (they can be
missing alleles), thus the total number of fields of each row _must be_
even.
Major and minor alleles are encoded by the characters `1` and `2`.
Missing genotypes are encoded by the pair `0 0` (i.e. by two fields
containing the missing allele `0`).
The pairs composed by a valid allele (`1` or `2`) and a missing allele
(`0`) _are not valid_.

Rows starting with the character `#` are considered as comments and
ignored.

**Remark:**
The order of the two alleles on each locus is meaningless (i.e., the
pair `2 1` is considered the same as the pair `1 2`).

A simple single-family pedigree composed by 5 individuals genotyped over
5 loci is as follows.

    0 1 0 0 1 phenotype 1 1 2 2 2 2 2 2 1 1
    0 2 0 0 2 phenotype 2 2 1 1 1 1 1 1 1 1
    0 3 1 2 2 phenotype 1 2 0 0 1 2 1 2 1 1
    0 4 0 0 1 phenotype 1 2 1 2 1 1 1 1 0 0
    0 5 4 3 1 phenotype 1 2 1 2 0 0 1 1 1 2


### Output: Haplotype Configuration ###

The haplotype configuration computed by reHC-* is represented in a
PED-like format.
In particular, the first six fields are equal to the PED format.
The remaining fields represent the computed haplotype pair of the
individual, where _each_ field represents the two alleles on a
single locus (separated by the character `|`).
In this case, the order of the two alleles is important and represents
the _phase_ of each locus.
The first allele in each pair is the paternal allele, while the second
one is the maternal allele.

For the example, a zero-recombinant haplotype configuration for the
previous genotyped pedigree is as follows.

    0 1 0 0 1 phenotype 1|1 2|2 2|2 2|2 1|1
    0 2 0 0 2 phenotype 2|2 1|1 1|1 1|1 1|1
    0 3 1 2 2 phenotype 1|2 2|1 2|1 2|1 1|1
    0 4 0 0 1 phenotype 2|1 1|2 1|1 1|1 2|2
    0 5 4 3 1 phenotype 1|2 2|1 1|1 1|1 2|1

where the two (multi-locus) haplotypes of individual `5` are `12112`
(paternal haplotype) and `21111` (maternal haplotype).



## License ##

reHC-* is released under the terms of the GNU General Public License
(GPL) as published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

reHC-* is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

Please refer to file `COPYING` or to the
[GNU website](http://www.gnu.org/licenses/) for a copy of the GNU
General Public License.



## Acknowledgments ##

The template of reHC-* is based on the
[cpp-project-template](http://code.google.com/p/cpp-project-template/)
by Michael Aaron Safyan.

reHC-* incorporates the following SAT solvers:

- [CryptoMiniSat](http://gitorious.org/cryptominisat) version 2.9.1
  (commit bcfdd17915e3, date 7/Apr/2011) by Mate Soos, which is
  distributed under the GNU General Public License version 3;
- [MiniSat](http://www.minisat.se/MiniSat.html) version 2.2.0 by Niklas
  Een and Niklas Sorensson, which is distributed under the MIT license.

We would like to thank Gianluca Della Vedova for useful discussions.



## Contacts ##

Please contact *Yuri Pirola* for additional information.  
E-mail:   <yuri.pirola@gmail.com>  
Web page: <http://bimib.disco.unimib.it/index.php/Pirola_Yuri>


