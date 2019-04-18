import sympy as sp
import perforad
# Define symbols for all examples
c = sp.Function("c")
u_1 = sp.Function("u_1"); u_1_b = sp.Function("u_1_b")
u_2 = sp.Function("u_2"); u_2_b = sp.Function("u_2_b")
u = sp.Function("u")    ; u_b = sp.Function("u_b")
i,j,k,C,D,n = sp.symbols("i,j,k,C,D,n")


######## 1D Wave Equation Example ########
# Build stencil expression
u_xx = u_1(i-1) - 2*u_1(i) + u_1(i+1)
expr = 2.0*u_1(i) - u_2(i) + c(i)*D*u_xx
# Build LoopNest object for this expression
lp = perforad.makeLoopNest(lhs=u(i), rhs=expr, counters = [i], bounds={i:[1,n-2]})
# Output primal and adjoint files
perforad.printfunction(name="wave1d", loopnestlist=[lp])
perforad.printfunction(name="wave1d_perf_b", loopnestlist=lp.diff({u:u_b, u_1:u_1_b, u_2: u_2_b}))


######## 3D Wave Equation Example ########
# Build stencil expression
u_xx = u_1(i-1,j,k) - 2*u_1(i,j,k) + u_1(i+1,j,k)
u_yy = u_1(i,j-1,k) - 2*u_1(i,j,k) + u_1(i,j+1,k)
u_zz = u_1(i,j,k-1) - 2*u_1(i,j,k) + u_1(i,j,k+1)
expr = 2.0*u_1(i,j,k) - u_2(i,j,k) + c(i,j,k)*D*(u_xx + u_yy + u_zz)
# Build LoopNest object for this expression
lp = perforad.makeLoopNest(lhs=u(i,j,k), rhs=expr, counters = [i,j,k], bounds={i:[1,n-2],j:[1,n-2],k:[1,n-2]})
# Output primal and adjoint files
perforad.printfunction(name="wave3d", loopnestlist=[lp])
perforad.printfunction(name="wave3d_perf_b", loopnestlist=lp.diff({u:u_b, u_1:u_1_b, u_2: u_2_b}))


######## 1D Burgers Equation Example ########
# Build stencil expression
ap = sp.functions.Max(u_1(i),0)
am = sp.functions.Min(u_1(i),0)
uxm = u_1(i)-u_1(i-1)
uxp = u_1(i+1)-u_1(i)
ux = ap*uxm+am*uxp
expr = u_1(i) - C * ux + D * (u_1(i+1) + u_1(i-1) - 2.0*u_1(i))
# Build LoopNest object for this expression
lp = perforad.makeLoopNest(lhs=u(i), rhs=expr, counters = [i], bounds={i:[1,n-2]})
# Output primal and adjoint files
perforad.printfunction(name="burgers1d", loopnestlist=[lp])
perforad.printfunction(name="burgers1d_perf_b", loopnestlist=lp.diff({u:u_b, u_1:u_1_b}))
