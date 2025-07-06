/* This index class allows you to explicitly use a part of a tensor.
   If you want to explicitly list all of the indices, just use
   int's. The usual way to do this is to declare a Number like

   Number<0> N; */

#pragma once

namespace FTensor
{
  template <const int N> class Number
  {
  public:
    Number(){};
    operator int() const { return N; }
  };
}
