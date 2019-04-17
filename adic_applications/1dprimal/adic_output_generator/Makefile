#CC=g++
CXXFLAGS= -I..

all: head_f.cn.xb.pp.c head_r.cn.xb.pp.c 

head_f.cn.xb.pp.c: head.c
	@echo "Creating ADIC2 forward output"
	$(ADIC_INST_DIR)/bin/adic2  -mforward --staticAnnotations --nary --include_dirs "/usr/include" -d1 $<
	mv head.cn.xb.pp.c $@
	mv ad_grad_saxpy-n_dense.h ad_grad_f_saxpy-n_dense.h
	sed -i 's/"ad_grad_saxpy-n_dense\.h"/"ad_grad_f_saxpy-n_dense\.h"/g' $@

ad_grad_f_saxpy-n_dense.h: head_f.cn.xb.pp.c


head_r.cn.xb.pp.c: head.c
	@echo "Creating ADIC2 reverse output"
	$(ADIC_INST_DIR)/bin/adic2  -mreverse --staticAnnotations --admm --include_dirs "/usr/include" -d1  $<
	mv head.cn.xb.pp.c $@
	#mv ad_grad_saxpy-dense.h ad_grad_r_saxpy-n_dense.h
	sed -i 's/"ad_grad_saxpy-n_dense\.h"/"ad_grad_r_saxpy-n_dense\.h"/g' $@

%.o :%.c
	$(CC)  $(CCFLAGS) -c $< -o $(<:.c=.o)	

%.lo:%.c
	$(HIDE_LIBTOOL)$(LIBTOOL) --mode=compile --tag=CC $(HIDE_COMPILE) $(CC) $(CCFLAGS) $(DENSE_REVERSE_DEFINE) -I$(RUNTIMEDIR) -c $<

%.o :%.c
	$(CXX)  $(CCFLAGS) -c $< -o $(<:.c=.o)	

%.lo:%.c
	$(HIDE_LIBTOOL)$(LIBTOOL) --mode=compile --tag=CXX $(HIDE_COMPILE) $(CXX) $(CXXFLAGS) $(DENSE_REVERSE_DEFINE) -I$(RUNTIMEDIR) -c $<

clean:
	rm -Rf *.o *.lo .libs *.la driver *~  
	rm -Rf temp.*  head.r.* ad_grad*
	rm -Rf ad_register_external_functions_internal.* ad_register_external_functions_export.*

.PHONY : clean 

