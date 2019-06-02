CC=icc
CFLAGS=-O3 -fopenmp -xhost
EXTRA_SOURCES=adstack/adBuffer.c adstack/adStack.c
INCLUDES=-I adstack/

drivers: driver1dwave burgers1d driver3dwave

driver1dwave: main.c wave1d.c wave1d_b.c wave1d_perf_b.c
	$(CC) $(CFLAGS) -o $@ $^ $(EXTRA_SOURCES) $(INCLUDES)

driver3dwave: driverwave3d.c wave3d.c wave3d_b.c wave3d_perf_b.c
	$(CC) $(CFLAGS) -o $@ $^ $(EXTRA_SOURCES) $(INCLUDES)
diff: wave1d_b wave3d_b burgers1d_b

wave3d_b: wave3d.c
	cpp -DTAPENADE wave3d.c > wave3d_tapenadecpp.c
	tapenade -o wave3d -b wave3d_tapenadecpp.c -head "wave3d(u,u_1,u_2)/(u,u_1,u_2)"

burgers1d_b: burgers1d.c
	cpp -DTAPENADE burgers1d.c > burgers1d_tapenadecpp.c
	tapenade -o burgers1d -b burgers1d_tapenadecpp.c -head "burgers1d(u,u_1)/(u,u_1)"

wave1d_b: wave1d.c
	cpp -DTAPENADE wave1d.c > wave1d_tapenadecpp.c
	tapenade -o wave1d -b wave1d_tapenadecpp.c -head "wave1d(u,u_1,u_2)/(u,u_1,u_2)"

burgers1d: burgers1d.c driver.c burgers1d_b_nostack.c burgers1d_perf_b.c
	$(CC) $(CFLAGS) -o $@ $^ $(EXTRA_SOURCES) $(INCLUDES)
