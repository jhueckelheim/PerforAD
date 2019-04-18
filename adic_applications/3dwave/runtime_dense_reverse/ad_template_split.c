/*#########################################################
* This file is part of OpenAD released under the LGPL.   #
* The full COPYRIGHT notice can be found in the top      #
* level directory of the OpenAD distribution             #
* #########################################################
*/
  
#include "ad_tape.hpp"
#include "ad_rev.hpp"
  
void placeholder_replace(int a)
{
 return;
}

void ad_template()
{
  if (our_rev_mode.plain==TRUE){ 
   //original function
   placeholder_replace(1);
  } else if (our_rev_mode.tape == TRUE) {
     //tapingd
   placeholder_replace(2);
  } else if (our_rev_mode.adjoint ==TRUE) { 
    // adjoint
   placeholder_replace(3);  
  } 
}
  
