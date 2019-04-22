CIVL=java -jar /Users/snarayan/Soft/CIVL-1.19_5139/lib/civl-1.19_5139.jar
VERIFY=${CIVL} verify
COMPARE=${CIVL} compare
REPLAY=${CIVL} replay
MYDIR=.
DENSE_REVERSE_DEFINE = -DADIC_DENSE_REVERSE -DDEBUG


all: pad_civl tap_civl compare2

civl: adic_CIVL 

pad_CIVL: driver.c wave3d_perf_b.c pad_driver.c
	${VERIFY} driver.c wave3d_perf_b.c pad_driver.c

tap_CIVL: driver.c wave3d_b.c tap_driver.c
	${VERIFY} mystack.c wave3d_b.c driver.c tap_driver.c
	#${VERIFY} adStack.c adBuffer.c wave3d_b.c driver.c 

compare2: driver.c wave3d_b.c wave3d_perf_b.c mystack.c
	${COMPARE} -impl -D_gAssumption wave3d_perf_b.c driver.c pad_driver.c -spec  mystack.c wave3d_b.c driver.c tap_driver.c


clean:
	rm -f -r CIVLREP adic_out
	rm -Rf -f *~ hand_out tap_out *.o *.lo *.libs *.la
