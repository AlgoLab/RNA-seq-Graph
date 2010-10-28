#Makefile

SRC_DIR:=src
OBJ_DIR:=obj
BIN_DIR:=bin

CFLAGS+= -g -Wall
CXXFLAGS+= ${CFLAGS}
LIBS =

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
	${OBJ_DIR}/build_chains.o \
	${OBJ_DIR}/read_fasta.o \
	${OBJ_DIR}/table_entry.o \
	${OBJ_DIR}/RNA_seq.o \

${BIN_DIR}/read_input: ${read_input_OBJS}
	@echo 'Linking $@'; \
	mkdir -p ${BIN_DIR}; \
	${CXX} ${CXXFLAGS} -o $@ $^ ${LIBS}

.PHONY: read_input
read_input: ${BIN_DIR}/read_input

.PHONY: clean
clean:
	@echo "Cleaning..."; \
	rm -f ${OBJ_DIR}/* ${BIN_DIR}/*
