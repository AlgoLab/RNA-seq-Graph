/**
* RNA-seq-Graph
* Method for reconstructing the Isoform Graph of a gene from RNA-seq data
*
* Copyright (C) 2011 Stefano Beretta <ste.beretta(-at-)gmail.com>
*
* Distributed under the terms of the GNU Affero General Public License (AGPL)
*
*
* This file is part of RNA-seq-Graph.
*
* RNA-seq-Graph is free software: you can redistribute it and/or modify
* it under the terms of the GNU Affero General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* RNA-seq-Graph is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU Affero General Public License for more details.
*
* You should have received a copy of the GNU Affero General Public License
* along with RNA-seq-Graph.  If not, see <http://www.gnu.org/licenses/>.
**/

#ifndef CONFIGURATION_H
#define CONFIGURATION_H
/*****************/
/* INPUT READING */
/*****************/
//Length of the considered reads
#define READ_LEN 64

/*Number of READ_LEN substrings extracted from the input reads
If nothing is specified all the possibile substrings are extracted
File: read_fasta.cpp
Function: read_fasta()
*/
#define TWO_READS

/**************************/
/* GRAPH BUILDING OPTIONS */
/**************************/
/* Gaps length while linking reads
File: join_chains.cpp
Function: link_fragment_chains()
*/
#define GAP_LENGTH 3

/* Merigin procedure of graph nodes
after linking
File: join_chains.cpp
Function: link_fragment_chains()
*/
//#define MERGING

/**********************/
/* REFINEMENT OPTIONS */
/**********************/
/* Minimum length of tiny blocks
File: join_chains.cpp
Function: link_fragment_chains()
*/
#define MIN_LEN_TINY 6

/*********************/
/* OUTPUT GENERATION */
/*********************/
/*Enable GDL output
File: join_chains.cpp
Function: print_graph()
*/
#define GDL_OUT

/*Enable GDL_sign output
in which isolated nodes with length< MIN_SIGN_LENGTH
are omitted 
File: join_chains.cpp
Function: print_graph()
*/
#define SIGNIFICATIVE
#define MIN_SIGN_LENGTH 200

#endif
