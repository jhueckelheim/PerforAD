import sympy as sp
import textwrap
from operator import itemgetter

#verboseprint = print if verbose else lambda *a, **k: None

class LoopNest:
  def __init__(self,body,bounds):
    self.body = body
    self.counters = list(bounds.keys())
    self.bounds = bounds
    for counter in self.counters:
      start = bounds[counter][0]
      end = bounds[counter][1]
      body = _Loop_(body,counter,start,end)
    self.loop = body

  def __str__(self):
    return str(self.loop)

  def diff(self, invar_b, outvar_b):
    body_b = self.body.diff(invar_b, outvar_b)
    # A nest is a tuple that contains
    #  - a list containing (offset, statement) tuples
    #  - a dict with {counter: loop bounds}
    # This method takes a nest, and splits it up into several nests
    # such that each nest will iterate over a subset of the original domain,
    # and contains only the statements that are valid in that subset.
    nestlist = [(body_b,self.bounds)]
    # This loop goes over all dimensions. In each dimension, nestlist is replaced
    # with a new nestlist that has been split in that dimension.
    for counter in self.counters:
      print("counter %s"%counter)
      newnestlist = []
      # This loop goes over all nests. Each nest may get split into several new
      # nests, which are appended to newnestlist.
      for nest,nestbound in nestlist:
        print("  nest %s"%str(nestbound))
        nest.sort(key=lambda x: x[0][counter])
        print("  -with offs %s"%list(map(lambda x: x[0],nest)))
        # Multiple offsets may hav the same offset in the current dimension. We need to find
        # the positions in the array where the offset in this dimension changes.
        offsets = list(map(lambda x: x[0][counter],nest))
        uniqueoffsets = list(set(offsets))
        uniqueoffsets.sort()
        chunklimits = [offsets.index(i) for i in uniqueoffsets]
        for i in range(len(chunklimits)-1):
          # Get the range of offsets that this loop nest will contain
          offs_0 = offsets[chunklimits[i]]
          offs_1 = offsets[chunklimits[i+1]]
          # Get all statements that need to be part of the body of this prequel loop
          stmts_pre  = nest[slice(0,chunklimits[i+1])]
          # Get all statements that need to be part of the body of this sequel loop
          stmts_post = nest[slice(chunklimits[i+1],len(nest))]
          # Compute the new loop bounds after applying offsets
          bounds_pre = nestbound.copy()
          bounds_post = nestbound.copy()
          bounds_pre[counter] = nestbound[counter][0]+offs_1-1,nestbound[counter][0]+offs_0
          bounds_post[counter] = nestbound[counter][1]+offs_1,nestbound[counter][1]+offs_0+1
          print("    pre %s"%bounds_pre)
          print("    post %s"%bounds_post)
          # Append the nest to the new list of nests
          newnestlist.append((stmts_pre,bounds_pre))
          newnestlist.append((stmts_post,bounds_post))
        # Finally, create the core loop and append it to the new list of nests
        stmts_core = nest
        bounds_core = nestbound.copy()
        bounds_core[counter] = nestbound[counter][0]+nest[-1][0][counter],nestbound[counter][1]+nest[0][0][counter]
        print("    core %s"%bounds_core)
        newnestlist.append((stmts_core,bounds_core))
      # Replace the old nest list with the refined one, ready for the next iteration
      print(newnestlist)
      nestlist = newnestlist
    # Finally, take all nests and turn them into actual LoopNest objects.
    loops = []
    for body_b,nestbound in nestlist:
      statements = map(itemgetter(1),body_b)
      print(nestbound)
      loops.append(LoopNest(statements,nestbound))
    return loops

class _Loop_:
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
          idxlist.append(list(map(lambda x: (x[1]-x[2][1]),zip(self.idx_out,ofs,inpidx))))
        shifted_idx_in.append(idxlist)
      shifted_idx_in.append([list(map(lambda x: (-x[1]),inpidx))])
      expr_d = StencilExpression(outvar = invar_b, invar = self.invar+[outvar_b], idx_out = outidx, offset_in = shifted_idx_in, func = func_d)
      exprs.append((dict(inpidx),expr_d))
    return exprs
      

i, j, n, l, c, r, t, b, lb, rt, lt, rb = sp.symbols('i, j, n, l, c, r, t, b, lb, rt, lt, rb')
outvar, putvar, invar, jnvar = sp.symbols('outvar, putvar, invar, jnvar')
outvar_b, invar_b = sp.symbols('outvar_b, invar_b')
f = sp.Function('f')(l,c,r)
f2d = sp.Function('f2d')(l,b,c,r,t)
f2df = sp.Function('f2d')(l,b,lb,c,r,t,rt,lt,rb)

stexpr = StencilExpression(outvar, [invar], [i], [[[-1],[0],[1]]],f)
loop1d = LoopNest(body=stexpr, bounds={i:[2,n-1]})
print(loop1d)
for l in (loop1d.diff(invar_b, outvar_b)):
  print(l)

stexpr2d = StencilExpression(outvar, [invar,jnvar], [i,j], [[[-1,0],[0,-1],[0,0],[1,0],[0,1]],[[0,0]]],f2d)
loop2d = LoopNest(body=stexpr2d, bounds={i:[2,n-1],j:[2,n-1]})
print(loop2d)
for l in (loop2d.diff(invar_b, outvar_b)):
  print(l)

stexpr2df = StencilExpression(outvar, [invar], [i,j], [[[-1,0],[0,-1],[-1,-1],[0,0],[1,0],[0,1],[1,1],[-1,1],[1,-1]]],f2df)
loop2df = LoopNest(body=stexpr2df, bounds={i:[2,n-1],j:[2,n-1]})
print(loop2df)
for l in (loop2df.diff(invar_b, outvar_b)):
  print(l)
#for l,e in (stexpr2d.diff(invar_b, outvar_b)):
#  print(l)
#  print(e)
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
