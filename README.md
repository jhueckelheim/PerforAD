# PerforAD
Automatic Differentiation for high-performance stencil loops

## Why the name?
Stencil computations are a common motif in computer programs. Their name is probably derived from the printing technique. From Wikipedia:
> Stencilling produces an image or pattern by applying pigment to a surface over an intermediate object with designed gaps in it which create the pattern or image by only allowing the pigment to reach some parts of the surface.

When differentiating a stencil, PerforAD (perforate) creates a new stencil with the same shape as the old one.
Also, PerforAD results in really fast, high-performance Automatic Differentiation (AD).
