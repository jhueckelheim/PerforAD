import sympy as sp
import textwrap
from operator import itemgetter

verbose = False
verboseprint = print if verbose else lambda *a, **k: None

class LoopNest:
  def __init__(self,body,bounds,counters,arrays,scalars,ints):
    self.body = body
    self.bounds = bounds
    self.counters = counters
    self.arrays = arrays
    self.scalars = scalars
    self.ints = ints

  def __str__(self):
    body = ""
    try:
      # Try to join the list of statements in the loop body together.
      # If this fails, there is only one statement (no list).
      outlist = []
      for stmt in self.body:
        outlist.append(str(stmt))
      body = "\n".join(outlist)
    except TypeError:
      body = str(self.body)
    rcounters = self.counters.copy()
    rcounters.reverse()
    loops = {}
    outermost = None         # keep track which is the outermost loop so we can place an OpenMP pragma in front
    counters_within = []     # keep track of counters that are set within the OpenMP parallel loop, which need to be privatised.
    for counter in rcounters:
      # this builds the loop nest from the innermost to the outermost loop
      counters_within.append(counter)
      start = self.bounds[counter][0]
      end = self.bounds[counter][1]
      try:
        # Try and remove loop nests where any dimension has 0 iterations.
        if(end-start<0):
          return ""
      except TypeError:
        # If start or end contain sympy.Symbols, the relation can not be known
        # and this throws an error. We ignore this and assume that the loop may
        # have more than 0 iterations (it does not matter if this turns out to
        # be false at runtime, just adds dead code).
        pass
      try:
        # Try and simplify loops with only one iteration (print the loop
        # body, and assign the correct value to the counter variable).
        if(end-start==0):
          loops[counter] = False
          continue
      except TypeError:
        # If that failed, perhaps the loop bounds are again symbolic. Print the
        # loop, everything will be fine at runtime.
        pass
      loops[counter] = True
      # whenever we find another nontrivial loop (more than one iteration),
      # save it as the current outermost loop. By the end of this process,
      # this variable will contain the actual outermost loop, and whatever
      # counters we had found within that loop.
      outermost = (counter,counters_within.copy())
    for counter in rcounters:
      start = self.bounds[counter][0]
      end = self.bounds[counter][1]
      is_a_loop = loops[counter]
      # nontrivial loops get a for(...){...} construct.
      if(is_a_loop):
        # the outermost nontrivial loop also gets an OpenMP pragma
        omp = ""
        if(outermost[0] == counter):
          omp = "#pragma omp parallel for private(%s)\n"%",".join(map(str,outermost[1]))
        body = "%sfor ( %s=%s; %s<=%s; %s++ ) {\n%s\n}"%(omp,counter,start,counter,end,counter,textwrap.indent(str(body),4*" "))
      # trivial loops (those with exactly one iteration) just get a statement.
      else:
        body = "%s=%s;\n%s"%(counter,start,body)
    return body

  def diff(self, diffvars):
    body_b = self.body.diff(diffvars)
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
      verboseprint("counter %s"%counter)
      newnestlist = []
      # This loop goes over all nests. Each nest may get split into several new
      # nests, which are appended to newnestlist.
      for nest,nestbound in nestlist:
        verboseprint("  nest %s"%str(nestbound))
        nest.sort(key=lambda x: x[0][counter])
        verboseprint("  -with offs %s"%list(map(lambda x: x[0],nest)))
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
          verboseprint("    pre %s"%bounds_pre)
          verboseprint("    post %s"%bounds_post)
          # Append the nest to the new list of nests
          newnestlist.append((stmts_pre,bounds_pre))
          newnestlist.append((stmts_post,bounds_post))
        # Finally, create the core loop and append it to the new list of nests
        stmts_core = nest
        bounds_core = nestbound.copy()
        bounds_core[counter] = nestbound[counter][0]+nest[-1][0][counter],nestbound[counter][1]+nest[0][0][counter]
        verboseprint("    core %s"%bounds_core)
        newnestlist.append((stmts_core,bounds_core))
      # Replace the old nest list with the refined one, ready for the next iteration
      nestlist = newnestlist
    # Finally, take all nests and turn them into actual LoopNest objects.
    loops = []
    arrays = []
    for v in self.arrays:
      arrays.append(v)
      if(v in diffvars):
        arrays.append(diffvars[v])
    for body_b,nestbound in nestlist:
      statements = map(itemgetter(1),body_b)
      verboseprint(nestbound)
      loops.append(LoopNest(statements,nestbound,self.counters,arrays,self.scalars,self.ints))
    return loops

class SympyFuncStencil:
  def __init__(self, name, args):
    self.name = name
    self.args = args

  def __str__(self):
    return str(self.name)

  def at(self,inputs):
    argstr = ", ".join(map(lambda arg: "%s=%s"%(arg[0],arg[1]), zip(self.args,inputs)))
    return "%s[%s]"%(self.name,argstr)

  def diff(self,wrt):
    diffname = "%s_d"%self.name
    diffargs = self.args + ["%s_d"%wrt]
    return SympyFuncStencil(diffname,diffargs)

class SympyExprStencil:
  def __init__(self, expr, args):
    self.expr = expr
    self.args = args

  def __str__(self):
    return self.expr

  def at(self,inputs):
    subs = dict(zip(self.args,inputs))
    return self.expr.subs(subs)

  def args(self):
    return self.args

  def diff(self,wrt):
    wrtb = sp.Symbol('wrtb')
    return SympyExprStencil(self.expr.diff(wrt)*wrtb,self.args+[wrtb])

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
          idxlist.append(dim+of)
        args.append(var(*idxlist))
      lhsargs = ",".join(map(lambda x: str(x),self.idx_out))
      lhs = "%s(%s)"%(self.outvar,lhsargs)
    return "%s += %s;"%(lhs,self.func.at(args))

  def diff(self,diffvars):
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
      if(inp[0] in diffvars):
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
        expr_d = StencilExpression(outvar = diffvars[inp[0]], invar = self.invar+[diffvars[self.outvar]], idx_out = outidx, offset_in = shifted_idx_in, func = func_d)
        exprs.append((dict(inpidx),expr_d))
    return exprs
      

i, j, n, l, c, r, t, b, lb, rt, lt, rb, a = sp.symbols('i, j, n, l, c, r, t, b, lb, rt, lt, rb, a')
outv = sp.Function('outv')
inv = sp.Function('inv')
vel = sp.Function('vel')
outv_b = sp.Function('outv_b')
inv_b = sp.Function('inv_b')
vel_b = sp.Function('vel_b')

def printfunction(name, loopnestlist):
  counters = loopnestlist[0].counters
  arrays = loopnestlist[0].arrays
  scalars = loopnestlist[0].scalars
  ints = loopnestlist[0].ints
  funcdefs = """
#ifndef TAPENADE
#include <math.h>
#endif
#define Max(x,y) fmax(x,y)
#define Min(x,y) fmin(x,y)
#define Heaviside(x) ((x>=0)?1.0:0.0)
"""
  arrtransformlist = []
  for varname in arrays:
    arglist = list(map(lambda x: x*"x",range(1,len(counters)+1)))
    arrtransformlist.append("#define %s(%s) %s[%s]"%(varname,",".join(arglist),varname,"][".join(arglist)))
  cpp = "%s\n%s\n"%(funcdefs,"\n".join(arrtransformlist))
  args = list(map(lambda x: "double%s %s"%(len(counters)*"*",x),arrays))
  args = args + list(map(lambda x: "double %s"%(x),scalars))
  args = args + list(map(lambda x: "int %s"%(x),ints))
  declarations = "\n".join(list(map(lambda x: "int %s;"%x, counters)))
  body = [textwrap.indent(declarations,4*" ")]
  for loopnest in loopnestlist:
    body.append(textwrap.indent(str(loopnest),4*" "))
  filename = "generated/%s.c"%name
  print("Writing to %s"%filename)
  file = open(filename,"w")
  file.write("%svoid %s(%s) {\n%s\n}"%(cpp,name, ", ".join(args), "\n".join(body)))
  file.close() 

#f = SympyFuncStencil("foo",[l,c,r])
#stexpr = StencilExpression(outv, [inv], [i], [[[-1],[0],[1]]],f)
#loop1d = LoopNest(body=stexpr, bounds={i:[2,n-1]})
#print(loop1d)
#for lp in (loop1d.diff({inv:inv_b, outv:outv_b})):
#  print(lp)

def makeLoopNest(lhs, rhs, counters, bounds):
  functions = list(rhs.atoms(sp.Function))
  functions.sort(key=lambda x: x.func)
  scalars = [s for s in rhs.atoms(sp.Symbol) if s not in counters]
  ints = []
  for b in bounds:
    try:
      ints = ints + list(bounds[b][0].atoms(sp.Symbol))
    except AttributeError:
      pass
    try:
      ints = ints + list(bounds[b][1].atoms(sp.Symbol))
    except AttributeError:
      pass
  ints = list(set(ints))
  funcNames = []
  offsets = []
  funcID = 0
  subs = []
  for func in functions:
    funcName = func.func
    if not (funcName in funcNames):
      funcNames.append(funcName)
      offsets.append([])
    funcArgs = list(map(lambda x: x[0]-x[1], zip(func.args, counters)))
    offsets[-1].append(funcArgs)
    subs.append([func, "perforad_arg_%d"%funcID])
    funcID = funcID + 1
  # TODO check that offsets are const
  # TODO check that lhs args are in correct order
  exprvars = list(map(itemgetter(1),subs))
  f2d = SympyExprStencil(rhs.subs(subs),exprvars)
  stexpr = StencilExpression(lhs.func, funcNames, counters, offsets, f2d)
  loop = LoopNest(body=stexpr, bounds = bounds, counters = counters, arrays = [lhs.func]+funcNames, scalars = scalars, ints = ints)
  return loop
