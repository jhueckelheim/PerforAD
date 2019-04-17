/*!#########################################################
! This file is part of OpenAD released under the LGPL.   #
! The full COPYRIGHT notice can be found in the top      #
! level directory of the OpenAD distribution             #
!#########################################################
*/

  
#include "ad_tape.h"

/* Globals */
  double __ADIC_double_tape[max_double_tape_size];
  int    __ADIC_integer_tape[max_integer_tape_size];
  int    __ADIC_logical_tape[max_logical_tape_size];
  char*  __ADIC_character_tape[max_character_tape_size];
  int    __ADIC_stringlength_tape[max_stringlength_tape_size];
  void*  __ADIC_address_tape[max_address_tape_size];

  int __ADIC_double_tape_pointer, __ADIC_integer_tape_pointer, __ADIC_logical_tape_pointer, __ADIC_character_tape_pointer, __ADIC_stringlength_tape_pointer, __ADIC_address_tape_pointer;


  void __ADIC_TapeInit()
  {
    __ADIC_double_tape_pointer       = 0;
    __ADIC_integer_tape_pointer      = 0;
    __ADIC_logical_tape_pointer      = 0;
    __ADIC_character_tape_pointer    = 0;
    __ADIC_stringlength_tape_pointer = 0;
   }

  
/*
  void dump()
  {
    int i;
    printf("\n double tape");
    for (i=0; i < double_tape_pointer; i++)
      printf("\n %d", double_tape[i]);

    printf("\n integer tape");
    for(i=1; i < integer_tape_pointer-1; i++)
      printf("\n %d", integer_tape[i]);

    printf("\n logical tape");
    for(i=1; i < logical_tape_pointer-1; i++)
      printf("\n %d", logical_tape[i]);

    printf("\n character tape");
    printf("\n %s", character_tape);
    
    printf("\n stringlength tape");
    for(i=1; i < stringlength_tape_pointer-1; i++)
      printf("\n %d", stringlength_tape[i]);
      
  }
*/


