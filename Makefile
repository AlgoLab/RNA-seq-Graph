###
#
# RNA-seq-Graph
# Method for reconstructing the Isoform Graph of a gene from RNA-seq data
#
# Copyright (C) 2011 Stefano Beretta <ste.beretta(-at-)gmail.com>
#
# Distributed under the terms of the GNU Affero General Public License (AGPL)
#
#
# This file is part of RNA-seq-Graph.
#
# RNA-seq-Graph is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RNA-seq-Graph is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with RNA-seq-Graph.  If not, see <http://www.gnu.org/licenses/>.
#
###

#Makefile

SRC_DIR:=src
OBJ_DIR:=obj
BIN_DIR:=bin

CFLAGS+= -g -Wall -O2 -UNDEBUG -march=native -Wno-deprecated
CXXFLAGS+= ${CFLAGS}
LIBS = -l boost_graph

.PHONY: all
all:action read_input

.PHONY: action
action:
	@echo "Compiling..."

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp Makefile
	@echo '* Compiling $<'; \
	mkdir -pv $(dir $@) ; \
	$(CXX) $(CXXFLAGS) -o $@ -c $<

read_input_OBJS= \
	${OBJ_DIR}/Main.o \
	${OBJ_DIR}/graph_refinement.o \
	${OBJ_DIR}/join_chains.o \
	${OBJ_DIR}/build_chains.o \
	${OBJ_DIR}/read_fasta.o \
	${OBJ_DIR}/table_entry.o \
	${OBJ_DIR}/RNA_seq.o \

${BIN_DIR}/build_RNA_seq_graph: ${read_input_OBJS}
	@echo 'Linking $@'; \
	mkdir -p ${BIN_DIR}; \
	${CXX} ${CXXFLAGS} -o $@ $^ ${LIBS}

.PHONY: read_input
read_input: ${BIN_DIR}/build_RNA_seq_graph

.PHONY: clean
clean:
	@echo "Cleaning..."; \
	rm -f ${OBJ_DIR}/* ${BIN_DIR}/*
