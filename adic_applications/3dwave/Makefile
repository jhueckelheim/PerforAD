CC=gcc
# CFLAGS=  -DCHECK_ERROR -I..
CFLAGS=  -g -O

all: driver

head.c: adic_output_generator/head.c
	cp $< $@

ifeq ($(N),)
N := 5
endif

CFLAGS+= -DN=$(N) -DDEBUG=1
#CFLAGS+= -DN=$(N) 
RUNTIMEDIR = ./runtime_dense_reverse/
SRC_DRIVER = driver.c head.c			
OBJ_DRIVER = $(subst .c,.lo,$(SRC_DRIVER))
LIBS = -lm 
HIDE_LIBTOOL = @
LIBTOOL = libtool
SRC_RUNTIME_DENSE_REVERSE = runtime_dense_reverse/ad_grad.c runtime_dense_reverse/ad_rev.c runtime_dense_reverse/ad_tape.c runtime_dense_reverse/admm.c
OBJ_RUNTIME_DENSE_REVERSE = $(SRC_RUNTIME_DENSE_REVERSE:.c=.lo)
DENSE_REVERSE_OUTPUT = head.cn.xb.pp.c
DENSE_REVERSE_OUTPUT_OBJECT = $(subst pp.c,pp.lo,$(DENSE_REVERSE_OUTPUT))
SRC_DENSE_REVERSE = $(DENSE_REVERSE_OUTPUT)
OBJ_DENSE_REVERSE = $(subst .cn.xb.pp.c,.cn.xb.pp.lo ,$(SRC_DENSE_REVERSE))
DENSE_REVERSE_DEFINE = -DADIC_DENSE_REVERSE
RUNTIME_INCLUDE = -I$(RUNTIMEDIR)
ADIC_INCLUDE = -I$(ADIC_DIR)/include
ALL_INCLUDES = $(RUNTIME_INCLUDE) $(ADIC_INCLUDE)


head.cn.xb.pp.c: adic_output_generator/head_r.cn.xb.pp.c
	@echo "Creating ADIC2 reverse output"
	cp $< $@
	#cp adic_output_generator/ad_grad_r_saxpy-n_dense.h ./
	#$(ADIC_INST_DIR)/bin/adic2   -mreverse --staticAnnotations --include_dirs "/home/derivs/software/${WHATAMI}/ADIC/rose_inst/include/edg/${CCROSE_BASENAME}-4.8_HEADERS/" -d1  $<
	sed -i 's/"ad_types\.hpp"/"ad_types\.h"/g' $@
	#sed -i 's/"ad_types\.hpp"/"ad_types\.h"/g' ad_grad_r_saxpy-n_dense.h

	
$(OBJ_DENSE_REVERSE): $(SRC_DENSE_REVERSE)
# Runtime library
adic_gradvec_length.h:
	rm adic_gradvec_length.h
	cat "#define ADIC_GRADVEC_LENGTH ($(N)*$(N))" adic_gradvec_length.h

$(RUNTIMEDIR)/adic_gradvec_length.h: adic_gradvec_length.h
	cp ./adic_gradvec_length.h $(RUNTIMEDIR)/adic_gradvec_length.h	

OBJ_RUNTIME_DENSE_REVERSE_INDIVIDUAL: 
	echo "Compiling OBJ_RUNTIME_DENSE_REVERSE_INDIVIDUAL"
	rm -f $(RUNTIMEDIR)/adic_gradvec_length.h
	echo "#define ADIC_GRADVEC_LENGTH ($(N)*$(N)*$(N))" >$(RUNTIMEDIR)/adic_gradvec_length.h
	$(HIDE_LIBTOOL)$(LIBTOOL) --mode=compile --tag=CC $(HIDE_COMPILE) \
				$(CC)  $(CFLAGS) \
			  $(RUNTIME_REVERSE_INCLUDE)  $(DENSE_REVERSE_DEFINE) -c runtime_dense_reverse/ad_grad.c -o runtime_dense_reverse/ad_grad.lo
	$(HIDE_LIBTOOL)$(LIBTOOL) --mode=compile --tag=CC $(HIDE_COMPILE) \
	$(CC)  $(CFLAGS) \
	$(RUNTIME_REVERSE_INCLUDE)  $(DENSE_REVERSE_DEFINE) -c runtime_dense_reverse/ad_rev.c -o runtime_dense_reverse/ad_rev.lo
	$(HIDE_LIBTOOL)$(LIBTOOL) --mode=compile --tag=CC $(HIDE_COMPILE) \
	$(CC)  $(CFLAGS) \
	$(RUNTIME_REVERSE_INCLUDE)  $(DENSE_REVERSE_DEFINE) -c runtime_dense_reverse/ad_tape.c -o runtime_dense_reverse/ad_tape.lo
	$(HIDE_LIBTOOL)$(LIBTOOL) --mode=compile --tag=CC $(HIDE_COMPILE) \
	$(CC)  $(CFLAGS) \
	$(RUNTIME_REVERSE_INCLUDE)  $(DENSE_REVERSE_DEFINE) -c runtime_dense_reverse/admm.c -o runtime_dense_reverse/admm.lo
	echo "\n DONE Compiling OBJ_RUNTIME_DENSE_REVERSE_INDIVIDUAL"

#Rule for the dense reverse runtime library
libruntimedense.la: OBJ_RUNTIME_DENSE_REVERSE_INDIVIDUAL
	@echo -e "\nCompiling libruntimedense.la...\n" 
	$(HIDE_LIBTOOL)$(LIBTOOL) --mode=link --tag=CC $(HIDE_COMPILE) $(CC) $(OBJ_RUNTIME_DENSE_REVERSE) -o $@ 
	@echo -e "\nDONE Compiling libruntimedense.la...\n"

#Compiles the reverse dense src, driver and derivative code
$(DENSE_REVERSE_OUTPUT_OBJECT): $(SRC_RUNTIME_DENSE_REVERSE) $(SRC)
	$(RM) $(DENSE_REVERSE_OUTPUT_OBJECT)
	@echo -e "\nCompiling DENSE_REVERSE_INDIVIDUAL\n"
	$(HIDE_LIBTOOL)$(LIBTOOL) --mode=compile --tag=CC $(HIDE_COMPILE) $(CC) $(CFLAGS) $(ALL_INCLUDES) $(DENSE_REVERSE_DEFINE) -DDIFFERENTIATION -c $(DENSE_REVERSE_OUTPUT) -o $(DENSE_REVERSE_OUTPUT_OBJECT) 

#Rule for the dense reverse executable. This requires grad_length.
driver: libruntimedense.la  $(DENSE_REVERSE_OUTPUT_OBJECT) $(OBJ_DRIVER) $(OBJ_DENSE_REVERSE)
	@echo -e "\n==== Linking AD code and driver routine ====\n"
	$(HIDE_LIBTOOL)$(LIBTOOL) --mode=link --tag=CC  $(HIDE_COMPILE) $(CC) \
	$(CFLAGS) -DDIFFERENTIATION $(ALL_INCLUDES)  $(DENSE_REVERSE_DEFINE) -o $@ $(OBJ_DENSE_REVERSE) $(OBJ_DRIVER) libruntimedense.la -lm

%.o :%.c
	$(CC)  $(CFLAGS) -c $< -o $(<:.c=.o)	

%.lo:%.c
	$(HIDE_LIBTOOL)$(LIBTOOL) --mode=compile --tag=CC $(HIDE_COMPILE) $(CC) $(CFLAGS)  $(DENSE_REVERSE_DEFINE) -I$(RUNTIMEDIR) -c $<
		    
clean:
	rm -Rf *.o *.lo .libs *.la driver *~  head.cn.* head.r.* $(RUNTIMEDIR)/adic_gradvec_length.h ad_grad_*.h
	rm -Rf $(RUNTIMEDIR)/*.o  $(RUNTIMEDIR)/*.lo $(RUNTIMEDIR)/.libs *.la driver $(RUNTIMEDIR)/*~ $(RUNTIMEDIR)/*.cn.* $(RUNTIMEDIR)/*.r.*
	rm head.c
			
.PHONY : clean driver

