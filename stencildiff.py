import sympy as sp
import textwrap
from operator import itemgetter

class Loop:
  def __init__(self,body,counter,start,end):
    self.counter = counter
    self.start   = start
    self.end     = end
    self.body    = body

  def __str__(self):
    try:
      # Try and remove loops with 0 iterations.
      if(self.end-self.start<0):
        return ""
    except TypeError:
      # If start or end contain sympy.Symbols, the relation can not be known
      # and this throws an error. We ignore this and assume that the loop may
      # have more than 0 iterations (it does not matter if this turns out to
      # be false at runtime, just adds dead code).
      pass
    try:
      # Try to join the list of statements in the loop body together.
      # If this fails, there is only one statement (no list).
      outlist = []
      for stmt in self.body:
        outlist.append(str(stmt))
      res = "\n".join(outlist)
    except TypeError:
      res = str(self.body)
    try:
      # Try and simplify loops with only one iteration (print the loop
      # body, and assign the correct value to the counter variable).
      if(self.end-self.start==0):
        return "%s=%s\n%s"%(self.counter,self.start,res)
    except TypeError:
      # If that failed, perhaps the loop bounds are again symbolic. Print the
      # loop, everything will be fine at runtime.
      pass
    return "do %s=%s,%s\n%s\nend do"%(self.counter,self.start,self.end,textwrap.indent(str(res),4*" "))

  def diff(self, invar_b, outvar_b):
    body_b = self.body.diff(invar_b, outvar_b)
    body_b.sort(key=itemgetter(0))
    loops_pre = []
    loops_post = []
    for i in range(len(body_b)-1):
      offs_0,stmt_0 = body_b[i]
      offs_1,stmt_1 = body_b[i+1]
      stmts_pre  = list(map(itemgetter(1),body_b))[slice(0,i+1)]
      stmts_post = list(map(itemgetter(1),body_b))[slice(i+1,len(body_b))]
      loops_pre.append(Loop(body=stmts_pre,counter=self.counter,start=self.start+offs_0,end=self.start+offs_1-1))
      loops_post.append(Loop(body=stmts_post,counter=self.counter,start=self.end+offs_0+1,end=self.end+offs_1))
    loops_pre.append(Loop(body=list(map(itemgetter(1),body_b)),counter=self.counter,start=self.start+body_b[-1][0],end=self.end+body_b[0][0]))
    return loops_pre+loops_post
      
# TODO separate presentation/API from logic.
# StencilExpression should only deal with whatever is necessary for the logic,
# and can be extended by a FortranStencilExpression / SympyStencilExpression that
# adds a layer of sugar.
class StencilExpression:
  def __init__(self,outvar,invar,idx_out,offset_in,func):
    self.outvar     = outvar
    self.invar      = invar
    self.idx_out    = idx_out   # should be a loop counter (e.g. i)
    self.offset_in  = offset_in # should be a list of constant offsets (e.g. [i-1,i,i+1])
    self.func       = func

  def __str__(self):
    # Go through the list of input vars and their corresponding list of offsets
    args = []
    for (var,offsets) in list(zip(self.invar,self.offset_in)):
      # Print the var with each offset as an index
      for ofs in offsets:
        # The offset and the array index can have multiple dimensions
        idxlist = []
        for dim,of in list(zip(self.idx_out,ofs)):
          idxlist.append(str( dim+of ))
        args.append("%s[%s]"%(var,",".join(idxlist)))
      lhsargs = ",".join(map(lambda x: str(x),self.idx_out))
      lhs = "%s[%s]"%(self.outvar,lhsargs)
    return "%s = %s(%s)"%(lhs,self.func.__str__(),", ".join(args))

  def diff(self,invar_b,outvar_b):
    # All invars and offsets given to the StencilExpression are
    # zipped so that each offset has the correct invar, like so:
    # [(invar, [(i, -1)]), (invar, [(i, 0)]), (invar, [(i, 1)])]
    inputs = []
    for (var,offsets) in list(zip(self.invar,self.offset_in)):
      # Print the var with each offset as an index
      for ofs in offsets:
        # The offset and the array index can have multiple dimensions
        idxlist = list(zip(self.idx_out,ofs))
        inputs.append((var,idxlist))
    exprs = []
    # zip the list of function arguments and input variables
    for arg,inp in list(zip(self.func.args,inputs)):
      # Differentiate the function wrt. the current input
      func_d = self.func.diff(arg)
      # inpvar is the name of the current input variable.
      # inpidx is a tuple of counter variable (e.g. i) and offset (e.g. -1)
      inpvar, inpidx = inp
      # The output index of the diff'ed expression will be the same as that of
      # the primal expression (that's the whole point of this transformation)
      outidx = self.idx_out
      # We shift all other indices by the offset in inpidx to make this correct
      shifted_idx_in = []
      for (var,offsets) in list(zip(self.invar,self.offset_in)):
        idxlist = []
        for ofs in offsets:
          idxlist.append(list(map(lambda x: (x[1]+x[2][1]),zip(self.idx_out,ofs,inpidx))))
        shifted_idx_in.append(idxlist)
      shifted_idx_in.append([list(map(lambda x: (x[1]),inpidx))])
      expr_d = StencilExpression(outvar = invar_b, invar = self.invar+[outvar_b], idx_out = outidx, offset_in = shifted_idx_in, func = func_d)
      exprs.append((inpidx,expr_d))
    return exprs
      

i, j, n, l, c, r, t, b = sp.symbols('i, j, n, l, c, r, t, b')
outvar, putvar, invar, jnvar = sp.symbols('outvar, putvar, invar, jnvar')
outvar_b, invar_b = sp.symbols('outvar_b, invar_b')
f = sp.Function('f')(l,c,r)
f2d = sp.Function('f2d')(l,b,c,r,t)

stexpr = StencilExpression(outvar, [invar], [i], [[[-1],[0],[1]]],f)
loop1d = Loop(body=stexpr, counter=i, start=2, end=n-1)
print(loop1d)

stexpr2d = StencilExpression(outvar, [invar,jnvar], [i,j], [[[-1,0],[0,-1],[0,0],[1,0],[0,1]],[[0,0]]],f2d)
loop2dinner = Loop(body=stexpr2d, counter=i, start=2, end=n-1)
loop2douter = Loop(body=loop2dinner, counter=i, start=2, end=n-1)
print(loop2douter)

for l,e in (stexpr.diff(invar_b, outvar_b)):
  print(l)
  print(e)
for l,e in (stexpr2d.diff(invar_b, outvar_b)):
  print(l)
  print(e)
#
#stexpr = StencilExpression(outvar, [invar], [i,j], [[-1,0],[0,-1],[0,0],[1,0],[0,1]],f2d)
#lpinner = Loop(body=stexpr, counter=i, start=2, end=n-1)
#lpouter = Loop(body=lpinner, counter=j, start=2, end=n-1)
#print(lpouter)
#
#for l in (lpouter.diff(invar_b, outvar_b)):
#  print(l)

#stexpr = StencilExpression([outvar,putvar], [invar,jnvar], i, [[-1,0,1],[0]],f)
#lpinner = Loop(body=stexpr, counter=i, start=2, end=n-1)
#lpouter = Loop(body=lpinner, counter=j, start=1, end=n)
#print(lpouter)

#def Jacobi(f):
#  return [f.diff(x) for x in f.free_symbols]
#
#def Hessian(f):
#  return [[f.diff(x).diff(y) for x in f.free_symbols] for y in f.free_symbols]
#
#ul = sp.Symbol('ul')
#uc = sp.Symbol('uc')
#ur = sp.Symbol('ur')
#
#f_const = 0.1*ul + 0.2*uc + 0.3*ur
#f_centralself = 0.1*ul + 0.2*uc*uc + 0.3*ur
#f_allself = 0.1*ul*ul + 0.2*uc*uc + 0.3*ur*ur
#f_pairwise = 0.1*ul*uc + 0.2*uc*uc + 0.3*ur*uc
#f_alltoall = 0.1*ul*ul * 0.2*uc*uc * 0.3*ur*ur
#
#
#for func in [f_const, f_centralself, f_allself, f_pairwise, f_alltoall]:
#  print(func)
#  print(Jacobi(func))
#  print(Hessian(func))
#  print()
#
#do i=2,n-1
#  r_i = c1*u_i + c2*u_{i-1} + c3*u_{i+1}
#end do
#
#read write minidx maxidx
#i,   i,    2,     n-1 // c1
#i-1, i,    2,     n-1 // c2
#i+1, i,    2,     n-1 // c3
#------- transform -------
#i,   i,    2,     n-1           // c1
#i,   i+1,  1,     n-2  [+1, -1] // c2
#i,   i-1,  3,     n    [-1, +1] // c3
#------- intersect -------
#do i=1,1
#  r_{i+1} += c2*u_i
#end do
#do i=2,2
#  r_i     += c1*u_i
#  r_{i+1} += c2*u_i
#end do
#do i=3,n-2
#  r_i     += c1*u_i
#  r_{i+1} += c2*u_i
#  r_{i-1} += c3*u_i
#end do
#do i=n-1,n-1
#  r_i     += c1*u_i
#  r_{i-1} += c3*u_i
#end do
#do i=n,n
#  r_{i-1} += c3*u_i
#end do
#-------  unroll   -------
#r_2 += c2*u_2
#r_2 += c1*u_2
#r_3 += c2*u_2
#do i=3,n-2
#  r_i     += c1*u_i
#  r_{i+1} += c2*u_i
#  r_{i-1} += c3*u_i
#end do
#r_{n-1} += c1*u_{n-1}
#r_{n-2} += c3*u_{n-1}
#r_{n-1} += c3*u_{n-1}
#
#

#do i=1,n
#  do j=1,n
#    r(i,j) = 1*u(i-1,j)+2*u(i,j-1)+3*u(i,j)+4*u(i+1,j)+5*u(i,j+1)
#
#do i=1,n
#  do j=1,n
#    u(i-1,j) += 1*rb(i,j)
#    u(i,j-1) += 2*rb(i,j)
#    u(i,j)   += 3*rb(i,j)
#    u(i+1,j) += 4*rb(i,j)
#    u(i,j+1) += 5*rb(i,j)
#
#do i=1,n
#  do j=1,n
#    u(i,j) += 1*rb(i+1,j) [i=0,n-1; j=1,n]
#    u(i,j) += 2*rb(i,j+1) [i=1,n; j=0,n-1]
#    u(i,j) += 3*rb(i,j)   [i=1,n; j=1,n]
#    u(i,j) += 4*rb(i-1,j) [i=2,n+1; j=1,n]
#    u(i,j) += 5*rb(i,j-1) [i=1,n; j=2,n+1]
#
#Outer loop returns:
#  core loop i=2,n-1
#  prequel   i=0,1
#  sequel    i=n,n+1
#
#Inner loop returns:
#  core loop j=2,n-1
#  prequel   j=0,1
#  sequel    j=n,n+1
#do i=2,n-1
#  do j=2,n-1
#    u(i,j) += 1*rb(i+1,j)
#    u(i,j) += 2*rb(i,j+1)
#    u(i,j) += 3*rb(i,j)
#    u(i,j) += 4*rb(i-1,j)
#    u(i,j) += 5*rb(i,j-1)
#do i=0,1
#  do j=0,n+1
#    if [i=0,n-1; j=1,n]: u(i,j) += 1*rb(i+1,j)
#    if [i=1,n; j=0,n-1]: u(i,j) += 2*rb(i,j+1)
#    if [i=1,n; j=1,n]  : u(i,j) += 3*rb(i,j)  
#    if [i=1,n; j=2,n+1]: u(i,j) += 5*rb(i,j-1)
#do i=2,n-1
#  do j=0,1
#    if [i=0,n-1; j=1,n]: u(i,j) += 1*rb(i+1,j)
#    if [i=1,n; j=0,n-1]: u(i,j) += 2*rb(i,j+1)
#    if [i=1,n; j=1,n]  : u(i,j) += 3*rb(i,j)  
#    if [i=2,n+1; j=1,n]: u(i,j) += 4*rb(i-1,j)
#  do j=n,n+1
#    if [i=0,n-1; j=1,n]: u(i,j) += 1*rb(i+1,j)
#    if [i=1,n; j=1,n]  : u(i,j) += 3*rb(i,j)  
#    if [i=2,n+1; j=1,n]: u(i,j) += 4*rb(i-1,j)
#    if [i=1,n; j=2,n+1]: u(i,j) += 5*rb(i,j-1)
#do i=n,n+1
#  do j=0,n+1
#    if [i=1,n; j=0,n-1]: u(i,j) += 2*rb(i,j+1)
#    if [i=1,n; j=1,n]  : u(i,j) += 3*rb(i,j)  
#    if [i=2,n+1; j=1,n]: u(i,j) += 4*rb(i-1,j)
#    if [i=1,n; j=2,n+1]: u(i,j) += 5*rb(i,j-1)
