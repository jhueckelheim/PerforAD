CIVL=java -jar /Users/snarayan/Soft/CIVL-1.19_5139/lib/civl-1.19_5139.jar
VERIFY=${CIVL} verify
COMPARE=${CIVL} compare
REPLAY=${CIVL} replay
MYDIR=.
DENSE_REVERSE_DEFINE = -DADIC_DENSE_REVERSE -DDEBUG


all: pad_civl tap_civl compare2

civl: adic_CIVL 

pad_CIVL: driver.c pad_driver.c burgers1d_perf_b.c
	${VERIFY} driver.c pad_driver.c burgers1d_perf_b.c 

tap_CIVL: driver.c tap_driver.c burgers1d_b.c
	${VERIFY} mystack.c burgers1d_b.c tap_driver.c driver.c 
	#${VERIFY} adStack.c adBuffer.c burgers1d_b.c tap_driver.c driver.c 

compare2: driver.c pad_driver.c tap_driver.c burgers1d_b.c burgers1d_perf_b.c adBuffer.c
	${COMPARE} -impl -D_gAssumption burgers1d_perf_b.c driver.c pad_driver.c -spec  mystack.c burgers1d_b.c tap_driver.c driver.c


clean:
	rm -f -r CIVLREP adic_out
	rm -Rf -f *~ hand_out tap_out *.o *.lo *.libs *.la
