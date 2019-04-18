/*!#########################################################
! This file is part of OpenAD released under the LGPL.   #
! The full COPYRIGHT notice can be found in the top      #
! level directory of the OpenAD distribution             #
!#########################################################
*/
 #ifndef __AD_TAPE_HPP_
 #define __AD_TAPE_HPP_
 
#   if defined(__cplusplus)
        extern "C" {
#   endif

 #define max_double_tape_size  10000
 #define max_integer_tape_size 10000
 #define max_logical_tape_size 10000
 #define max_character_tape_size 1000
 #define max_stringlength_tape_size 1000
 #define max_address_tape_size 10000

  void __ADIC_TapeInit();
  void __ADIC_Dump();

/* BN: Globals should not be declared in headers */

  extern double __ADIC_double_tape[max_double_tape_size];
  extern int    __ADIC_integer_tape[max_integer_tape_size];
  extern int    __ADIC_logical_tape[max_logical_tape_size];
  extern char*  __ADIC_character_tape[max_character_tape_size];
  extern int    __ADIC_stringlength_tape[max_stringlength_tape_size];
  extern void*  __ADIC_address_tape[max_address_tape_size];

  extern int __ADIC_double_tape_pointer, __ADIC_integer_tape_pointer, __ADIC_logical_tape_pointer, __ADIC_character_tape_pointer, __ADIC_stringlength_tape_pointer, __ADIC_address_tape_pointer; 


//These two macros are a temporary fix to a rather big problem that
  // has to be dealt with. XB is going to change the names of the
  // inlinable functions and rely on Fortran's type resolution
  //mechanism to call the appropriate routine for the appropriate type
  //of arguments (a.k.a. polymorphism). C does not have this capability.

#define push_s0(var)                                               \
{                                                               \
     __ADIC_double_tape[__ADIC_double_tape_pointer]=var;        \
     __ADIC_double_tape_pointer=__ADIC_double_tape_pointer+1;   \
}

#define pop_s0(var)                                                \
{                                                               \
     __ADIC_double_tape_pointer=__ADIC_double_tape_pointer-1;   \
     var = __ADIC_double_tape[__ADIC_double_tape_pointer];      \
}

#define push_i_s0(var)                                               \
{                                                               \
     __ADIC_integer_tape[__ADIC_integer_tape_pointer]=var;        \
     __ADIC_integer_tape_pointer=__ADIC_integer_tape_pointer+1;   \
}

#define pop_i_s0(var)                                                \
{                                                               \
     __ADIC_integer_tape_pointer=__ADIC_integer_tape_pointer-1;   \
     var = __ADIC_integer_tape[__ADIC_integer_tape_pointer];      \
}

#define push_p_s0(var)                                               \
{                                                               \
     __ADIC_address_tape[__ADIC_address_tape_pointer]=var;        \
     __ADIC_address_tape_pointer=__ADIC_address_tape_pointer+1;   \
}

#define pop_p_s0(var)                                                \
{                                                               \
     __ADIC_address_tape_pointer=__ADIC_address_tape_pointer-1;   \
     var = __ADIC_address_tape[__ADIC_address_tape_pointer];      \
}

#   if defined(__cplusplus)
        }
#   endif

#endif
