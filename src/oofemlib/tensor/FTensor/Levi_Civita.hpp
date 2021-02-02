/// Fully Antisymmetric Levi-Civita Tensor

#pragma once

namespace FTensor
{
  /// Levi_Civita Classes
  template <class T = int> class Levi_Civita
  {
  public:
    /// Rank 2
    constexpr T operator()(const int N1, const int N2) const
    {
      return (N1 == N2) ? T(0) : ((N1 == 0) ? T(1) : T(-1));
    }

    template <char i, char j, int Dim0, int Dim1>
    typename std::enable_if<
      (Dim0 <= 2 && Dim1 <= 2),
      Tensor2_Expr<Levi_Civita<T>, T, Dim0, Dim1, i, j>>::type
    operator()(const Index<i, Dim0> &, const Index<j, Dim1> &)
    {
      return Tensor2_Expr<Levi_Civita<T>, T, Dim0, Dim1, i, j>(*this);
    };

    template <char i, int Dim0>
    constexpr auto operator()(const Index<i, Dim0> &, const int &N1)
    {
      static_assert(
        Dim0 <= 2,
        "Index dimension exceeds the maximum for Rank 2 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N1 >= 2 || N1 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()(Index<"
          + std::basic_string<char>{i} + ", " + std::to_string(Dim0) + ">, "
          + std::to_string(N1) + ")\n");
#endif
      auto TensorExpr
        = [this, N1](const int &N0) { return this->operator()(N0, N1); };
      return Tensor1_Expr<decltype(TensorExpr), T, Dim0, i>(TensorExpr);
    };

    template <char j, int Dim1>
    constexpr auto operator()(const int &N0, const Index<j, Dim1> &)
    {
      static_assert(
        Dim1 <= 2,
        "Index dimension exceeds the maximum for Rank 2 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N0 >= 2 || N0 < 0)
        throw std::out_of_range("Bad index in Levi_Civita<T>::operator()("
                                + std::to_string(N0) + ", Index<"
                                + std::basic_string<char>{j} + ", "
                                + std::to_string(Dim1) + ">)\n");
#endif
      auto TensorExpr
        = [this, N0](const int &N1) { return this->operator()(N0, N1); };
      return Tensor1_Expr<decltype(TensorExpr), T, Dim1, j>{TensorExpr};
    };

    /// Rank 3
    constexpr T operator()(const int N1, const int N2, const int N3) const
    {
      return (N1 == N2 || N1 == N3 || N2 == N3)
               ? T(0)
               : (((N1 + 1) % 3 == N2) ? T(1) : T(-1));
    }

    template <char i, char j, char k, int Dim0, int Dim1, int Dim2>
    typename std::enable_if<
      (Dim0 <= 3 && Dim1 <= 3 && Dim2 <= 3),
      Tensor3_Expr<Levi_Civita<T>, T, Dim0, Dim1, Dim2, i, j, k>>::type
    operator()(const Index<i, Dim0> &, const Index<j, Dim1> &,
               const Index<k, Dim2> &)
    {
      return Tensor3_Expr<Levi_Civita<T>, T, Dim0, Dim1, Dim2, i, j, k>(*this);
    };
    
    template <char i, char j, int Dim0, int Dim1>
    auto
    operator()(const Index<i, Dim0> &, const Index<j, Dim1> &, const int &N2)
    {
      static_assert(
        Dim0 <= 3 && Dim1 <= 3,
        "Index dimension exceeds the maximum for Rank 3 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N2 >= 3 || N2 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()("
          "Index<"
          + std::basic_string<char>{i} + ", " + std::to_string(Dim0)
          + ">, "
            "Index<"
          + std::basic_string<char>{j} + ", " + std::to_string(Dim1) + ">, "
          + std::to_string(N2) + ")\n");
#endif
      auto TensorExpr = [this, N2](const int &N0, const int &N1) {
        return this->operator()(N0, N1, N2);
      };
      return Tensor2_Expr<decltype(TensorExpr), T, Dim0, Dim1, i, j>{
        TensorExpr};
    }

    template <char i, char k, int Dim0, int Dim2>
    auto
    operator()(const Index<i, Dim0> &, const int &N1, const Index<k, Dim2> &)
    {
      static_assert(
        Dim0 <= 3 && Dim2 <= 3,
        "Index dimension exceeds the maximum for Rank 3 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N1 >= 3 || N1 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()("
          "Index<"
          + std::basic_string<char>{i} + ", " + std::to_string(Dim0) + ">, "
          + std::to_string(N1)
          + ", "
            "Index<"
          + std::basic_string<char>{k} + ", " + std::to_string(Dim2) + ">)\n");
#endif
      auto TensorExpr = [this, N1](const int &N0, const int &N2) {
        return this->operator()(N0, N1, N2);
      };
      return Tensor2_Expr<decltype(TensorExpr), T, Dim0, Dim2, i, k>{
        TensorExpr};
    }

    template <char j, char k, int Dim1, int Dim2>
    auto
    operator()(const int &N0, const Index<j, Dim1> &, const Index<k, Dim2> &)
    {
      static_assert(
        Dim1 <= 3 && Dim2 <= 3,
        "Index dimension exceeds the maximum for Rank 3 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N0 >= 3 || N0 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()(" + std::to_string(N0)
          + ", "
            "Index<"
          + std::basic_string<char>{j} + ", " + std::to_string(Dim1)
          + ">, "
            "Index<"
          + std::basic_string<char>{k} + ", " + std::to_string(Dim2) + ">)\n");
#endif
      auto TensorExpr = [this, N0](const int &N1, const int &N2) {
        return this->operator()(N0, N1, N2);
      };
      return Tensor2_Expr<decltype(TensorExpr), T, Dim1, Dim2, j, k>{
        TensorExpr};
    }

    template <char i, int Dim0>
    auto operator()(const Index<i, Dim0> &, const int &N1, const int &N2)
    {
      static_assert(
        Dim0 <= 3,
        "Index dimension exceeds the maximum for Rank 3 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N1 >= 3 || N1 < 0 || N2 >= 3 || N2 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()("
          "Index<"
          + std::basic_string<char>{i} + ", " + std::to_string(Dim0) + ">, "
          + std::to_string(N1) + ", " + std::to_string(N2) + ")\n");
#endif
      auto TensorExpr = [this, N1, N2](const int &N0) {
        return this->operator()(N0, N1, N2);
      };
      return Tensor1_Expr<decltype(TensorExpr), T, Dim0, i>{TensorExpr};
    }

    template <char j, int Dim1>
    auto operator()(const int &N0, const Index<j, Dim1> &, const int &N2)
    {
      static_assert(
        Dim1 <= 3,
        "Index dimension exceeds the maximum for Rank 3 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N0 >= 3 || N0 < 0 || N2 >= 3 || N2 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()(" + std::to_string(N0)
          + ", "
            "Index<"
          + std::basic_string<char>{j} + ", " + std::to_string(Dim1) + ">, "
          + std::to_string(N2) + ")\n");
#endif
      auto TensorExpr = [this, N0, N2](const int &N1) {
        return this->operator()(N0, N1, N2);
      };
      return Tensor1_Expr<decltype(TensorExpr), T, Dim1, j>{TensorExpr};
    }

    template <char k, int Dim2>
    auto operator()(const int &N0, const int &N1, const Index<k, Dim2> &)
    {
      static_assert(
        Dim2 <= 3,
        "Index dimension exceeds the maximum for Rank 3 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N0 >= 3 || N0 < 0 || N1 >= 3 || N1 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()(" + std::to_string(N0)
          + ", " + std::to_string(N1)
          + ", "
            "Index<"
          + std::basic_string<char>{k} + ", " + std::to_string(Dim2) + ">)\n");
#endif
      auto TensorExpr = [this, N0, N1](const int &N2) {
        return this->operator()(N0, N1, N2);
      };
      return Tensor1_Expr<decltype(TensorExpr), T, Dim2, k>{TensorExpr};
    }

    /// Rank 4
    constexpr T
    operator()(const int N1, const int N2, const int N3, const int N4) const
    {
      return (N1 == N2 || N1 == N3 || N1 == N4 || N2 == N3 || N2 == N4
              || N3 == N4)
               ? T(0)
               : ((N1 + N2 == 1 || N1 + N2 == 5)
                    ? (((N2 + N3) % 4 == 3) ? T(1) : T(-1))
                    : ((N1 + N2 == 2 || N1 + N2 == 4)
                         ? (((N2 + N3) % 4 == 1) ? T(1) : T(-1))
                         : (N1 + N2 == 3
                              ? (((N2 + N3) % 4 != 1) ? T(1) : T(-1))
                              : T(0))));
    }

    template <char i, char j, char k, char l, int Dim0, int Dim1, int Dim2,
              int Dim3>
    typename std::enable_if<
      (Dim0 <= 4 && Dim1 <= 4 && Dim2 <= 4 && Dim3 <= 4),
      Tensor4_Expr<Levi_Civita<T>, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l>>::type
    operator()(const Index<i, Dim0> &, const Index<j, Dim1> &,
               const Index<k, Dim2> &, const Index<l, Dim3> &)
    {
      return Tensor4_Expr<Levi_Civita<T>, T, Dim0, Dim1, Dim2, Dim3, i, j, k,
                          l>(*this);
    }

    template <char i, char j, char k, int Dim0, int Dim1, int Dim2>
    constexpr auto operator()(const Index<i, Dim0> &, const Index<j, Dim1> &,
                              const Index<k, Dim2> &, const int &N3)
    {
      static_assert(
        Dim0 <= 4 && Dim1 <= 4 && Dim2 <= 4,
        "Index dimension exceeds the maximum for Rank 4 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N3 >= 4 || N3 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()("
          "Index<"
          + std::basic_string<char>{i} + ", " + std::to_string(Dim0)
          + ">, "
            "Index<"
          + std::basic_string<char>{j} + ", " + std::to_string(Dim1)
          + ">, "
            "Index<"
          + std::basic_string<char>{k} + ", " + std::to_string(Dim2) + ">, "
          + std::to_string(N3) + ")\n");
#endif
      auto TensorExpr
        = [this, N3](const int &N0, const int &N1, const int &N2) {
            return this->operator()(N0, N1, N2, N3);
          };
      return Tensor3_Expr<decltype(TensorExpr), T, Dim0, Dim1, Dim2, i, j, k>{
        TensorExpr};
    }

    template <char i, char j, char l, int Dim0, int Dim1, int Dim3>
    constexpr auto operator()(const Index<i, Dim0> &, const Index<j, Dim1> &,
                              const int &N2, const Index<l, Dim3> &)
    {
      static_assert(
        Dim0 <= 4 && Dim1 <= 4 && Dim3 <= 4,
        "Index dimension exceeds the maximum for Rank 4 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N2 >= 4 || N2 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()("
          "Index<"
          + std::basic_string<char>{i} + ", " + std::to_string(Dim0)
          + ">, "
            "Index<"
          + std::basic_string<char>{j} + ", " + std::to_string(Dim1) + ">, "
          + std::to_string(N2)
          + ", "
            "Index<"
          + std::basic_string<char>{l} + ", " + std::to_string(Dim3) + ">)\n");
#endif
      auto TensorExpr
        = [this, N2](const int &N0, const int &N1, const int &N3) {
            return this->operator()(N0, N1, N2, N3);
          };
      return Tensor3_Expr<decltype(TensorExpr), T, Dim0, Dim1, Dim3, i, j, l>{
        TensorExpr};
    }

    template <char i, char k, char l, int Dim0, int Dim2, int Dim3>
    constexpr auto operator()(const Index<i, Dim0> &, const int &N1,
                              const Index<k, Dim2> &, const Index<l, Dim3> &)
    {
      static_assert(
        Dim0 <= 4 && Dim2 <= 4 && Dim3 <= 4,
        "Index dimension exceeds the maximum for Rank 4 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N1 >= 4 || N1 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()("
          "Index<"
          + std::basic_string<char>{i} + ", " + std::to_string(Dim0) + ">, "
          + std::to_string(N1)
          + ", "
            "Index<"
          + std::basic_string<char>{k} + ", " + std::to_string(Dim2)
          + ">, "
            "Index<"
          + std::basic_string<char>{l} + ", " + std::to_string(Dim3) + ">)\n");
#endif
      auto TensorExpr
        = [this, N1](const int &N0, const int &N2, const int &N3) {
            return this->operator()(N0, N1, N2, N3);
          };
      return Tensor3_Expr<decltype(TensorExpr), T, Dim0, Dim2, Dim3, i, k, l>{
        TensorExpr};
    }

    template <char j, char k, char l, int Dim1, int Dim2, int Dim3>
    constexpr auto operator()(const int &N0, const Index<j, Dim1> &,
                              const Index<k, Dim2> &, const Index<l, Dim3> &)
    {
      static_assert(
        Dim1 <= 4 && Dim2 <= 4 && Dim3 <= 4,
        "Index dimension exceeds the maximum for Rank 4 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N0 >= 4 || N0 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()(" + std::to_string(N0)
          + ", "
            "Index<"
          + std::basic_string<char>{j} + ", " + std::to_string(Dim1)
          + ">, "
            "Index<"
          + std::basic_string<char>{k} + ", " + std::to_string(Dim2)
          + ">, "
            "Index<"
          + std::basic_string<char>{l} + ", " + std::to_string(Dim3) + ">)\n");
#endif
      auto TensorExpr
        = [this, N0](const int &N1, const int &N2, const int &N3) {
            return this->operator()(N0, N1, N2, N3);
          };
      return Tensor3_Expr<decltype(TensorExpr), T, Dim1, Dim2, Dim3, j, k, l>{
        TensorExpr};
    }

    template <char i, char j, int Dim0, int Dim1>
    constexpr auto operator()(const Index<i, Dim0> &, const Index<j, Dim1> &,
                              const int &N2, const int &N3)
    {
      static_assert(
        Dim0 <= 4 && Dim1 <= 4,
        "Index dimension exceeds the maximum for Rank 4 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N2 >= 4 || N2 < 0 || N3 >= 4 || N3 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()("
          "Index<"
          + std::basic_string<char>{i} + ", " + std::to_string(Dim0)
          + ">, "
            "Index<"
          + std::basic_string<char>{j} + ", " + std::to_string(Dim1) + ">, "
          + std::to_string(N2) + ", " + std::to_string(N3) + ")\n");
#endif
      auto TensorExpr = [this, N2, N3](const int &N0, const int &N1) {
        return this->operator()(N0, N1, N2, N3);
      };
      return Tensor2_Expr<decltype(TensorExpr), T, Dim0, Dim1, i, j>{
        TensorExpr};
    }

    template <char i, char k, int Dim0, int Dim2>
    constexpr auto operator()(const Index<i, Dim0> &, const int &N1,
                              const Index<k, Dim2> &, const int &N3)
    {
      static_assert(
        Dim0 <= 4 && Dim2 <= 4,
        "Index dimension exceeds the maximum for Rank 4 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N1 >= 4 || N1 < 0 || N3 >= 4 || N3 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()("
          "Index<"
          + std::basic_string<char>{i} + ", " + std::to_string(Dim0) + ">, "
          + std::to_string(N1)
          + ", "
            "Index<"
          + std::basic_string<char>{k} + ", " + std::to_string(Dim2) + ">, "
          + std::to_string(N3) + ")\n");
#endif
      auto TensorExpr = [this, N1, N3](const int &N0, const int &N2) {
        return this->operator()(N0, N1, N2, N3);
      };
      return Tensor2_Expr<decltype(TensorExpr), T, Dim0, Dim2, i, k>{
        TensorExpr};
    }

    template <char j, char k, int Dim1, int Dim2>
    constexpr auto operator()(const int &N0, const Index<j, Dim1> &,
                              const Index<k, Dim2> &, const int &N3)
    {
      static_assert(
        Dim1 <= 4 && Dim2 <= 4,
        "Index dimension exceeds the maximum for Rank 4 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N0 >= 4 || N0 < 0 || N3 >= 4 || N3 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()(" + std::to_string(N0)
          + ", "
            "Index<"
          + std::basic_string<char>{j} + ", " + std::to_string(Dim1)
          + ">, "
            "Index<"
          + std::basic_string<char>{k} + ", " + std::to_string(Dim2) + ">, "
          + std::to_string(N3) + ")\n");
#endif
      auto TensorExpr = [this, N0, N3](const int &N1, const int &N2) {
        return this->operator()(N0, N1, N2, N3);
      };
      return Tensor2_Expr<decltype(TensorExpr), T, Dim1, Dim2, j, k>{
        TensorExpr};
    }

    template <char i, char l, int Dim0, int Dim3>
    constexpr auto operator()(const Index<i, Dim0> &, const int &N1,
                              const int &N2, const Index<l, Dim3> &)
    {
      static_assert(
        Dim0 <= 4 && Dim3 <= 4,
        "Index dimension exceeds the maximum for Rank 4 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N1 >= 4 || N1 < 0 || N2 >= 4 || N2 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()("
          "Index<"
          + std::basic_string<char>{i} + ", " + std::to_string(Dim0) + ">, "
          + std::to_string(N1) + ", " + std::to_string(N2)
          + ", "
            "Index<"
          + std::basic_string<char>{l} + ", " + std::to_string(Dim3) + ">)\n");
#endif
      auto TensorExpr = [this, N1, N2](const int &N0, const int &N3) {
        return this->operator()(N0, N1, N2, N3);
      };
      return Tensor2_Expr<decltype(TensorExpr), T, Dim0, Dim3, i, l>{
        TensorExpr};
    }

    template <char j, char l, int Dim1, int Dim3>
    constexpr auto operator()(const int &N0, const Index<j, Dim1> &,
                              const int &N2, const Index<l, Dim3> &)
    {
      static_assert(
        Dim1 <= 4 && Dim3 <= 4,
        "Index dimension exceeds the maximum for Rank 4 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N0 >= 4 || N0 < 0 || N2 >= 4 || N2 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()(" + std::to_string(N0)
          + ", "
            "Index<"
          + std::basic_string<char>{j} + ", " + std::to_string(Dim1) + ">, "
          + std::to_string(N2)
          + ", "
            "Index<"
          + std::basic_string<char>{l} + ", " + std::to_string(Dim3) + ">)\n");
#endif
      auto TensorExpr = [this, N0, N2](const int &N1, const int &N3) {
        return this->operator()(N0, N1, N2, N3);
      };
      return Tensor2_Expr<decltype(TensorExpr), T, Dim1, Dim3, j, l>{
        TensorExpr};
    }

    template <char k, char l, int Dim2, int Dim3>
    constexpr auto operator()(const int &N0, const int &N1,
                              const Index<k, Dim2> &, const Index<l, Dim3> &)
    {
      static_assert(
        Dim2 <= 4 && Dim3 <= 4,
        "Index dimension exceeds the maximum for Rank 4 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N0 >= 4 || N0 < 0 || N1 >= 4 || N1 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()(" + std::to_string(N0)
          + ", " + std::to_string(N1)
          + ", "
            "Index<"
          + std::basic_string<char>{k} + ", " + std::to_string(Dim2)
          + ">, "
            "Index<"
          + std::basic_string<char>{l} + ", " + std::to_string(Dim3) + ">)\n");
#endif
      auto TensorExpr = [this, N0, N1](const int &N2, const int &N3) {
        return this->operator()(N0, N1, N2, N3);
      };
      return Tensor2_Expr<decltype(TensorExpr), T, Dim2, Dim3, k, l>{
        TensorExpr};
    }

    template <char i, int Dim0>
    constexpr auto operator()(const Index<i, Dim0> &, const int &N1,
                              const int &N2, const int &N3)
    {
      static_assert(
        Dim0 <= 4,
        "Index dimension exceeds the maximum for Rank 4 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N1 >= 4 || N1 < 0 || N2 >= 4 || N2 < 0 || N3 >= 4 || N3 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()("
          "Index<"
          + std::basic_string<char>{i} + ", " + std::to_string(Dim0) + ">, "
          + std::to_string(N1) + ", " + std::to_string(N2) + ", "
          + std::to_string(N3) + ")\n");
#endif
      auto TensorExpr = [this, N1, N2, N3](const int &N0) {
        return this->operator()(N0, N1, N2, N3);
      };
      return Tensor1_Expr<decltype(TensorExpr), T, Dim0, i>{TensorExpr};
    }

    template <char j, int Dim1>
    constexpr auto operator()(const int &N0, const Index<j, Dim1> &,
                              const int &N2, const int &N3)
    {
      static_assert(
        Dim1 <= 4,
        "Index dimension exceeds the maximum for Rank 4 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N0 >= 4 || N0 < 0 || N2 >= 4 || N2 < 0 || N3 >= 4 || N3 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()(" + std::to_string(N0)
          + ", "
            "Index<"
          + std::basic_string<char>{j} + ", " + std::to_string(Dim1) + ">, "
          + std::to_string(N2) + ", " + std::to_string(N3) + ")\n");
#endif
      auto TensorExpr = [this, N0, N2, N3](const int &N1) {
        return this->operator()(N0, N1, N2, N3);
      };
      return Tensor1_Expr<decltype(TensorExpr), T, Dim1, j>{TensorExpr};
    }

    template <char k, int Dim2>
    constexpr auto operator()(const int &N0, const int &N1,
                              const Index<k, Dim2> &, const int &N3)
    {
      static_assert(
        Dim2 <= 4,
        "Index dimension exceeds the maximum for Rank 4 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N0 >= 4 || N0 < 0 || N1 >= 4 || N1 < 0 || N3 >= 4 || N3 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()(" + std::to_string(N0)
          + ", " + std::to_string(N1)
          + ", "
            "Index<"
          + std::basic_string<char>{k} + ", " + std::to_string(Dim2) + ">, "
          + std::to_string(N3) + ")\n");
#endif
      auto TensorExpr = [this, N0, N1, N3](const int &N2) {
        return this->operator()(N0, N1, N2, N3);
      };
      return Tensor1_Expr<decltype(TensorExpr), T, Dim2, k>{TensorExpr};
    }

    template <char l, int Dim3>
    constexpr auto operator()(const int &N0, const int &N1, const int &N2,
                              const Index<l, Dim3> &)
    {
      static_assert(
        Dim3 <= 4,
        "Index dimension exceeds the maximum for Rank 4 Levi-Civita Tensor");
#ifdef FTENSOR_DEBUG
      if(N0 >= 4 || N0 < 0 || N1 >= 4 || N1 < 0 || N2 >= 4 || N2 < 0)
        throw std::out_of_range(
          "Bad index in Levi_Civita<T>::operator()(" + std::to_string(N0)
          + ", " + std::to_string(N1) + ", " + std::to_string(N2)
          + ", "
            "Index<"
          + std::basic_string<char>{l} + ", " + std::to_string(Dim3) + ">)\n");
#endif
      auto TensorExpr = [this, N0, N1, N2](const int &N3) {
        return this->operator()(N0, N1, N2, N3);
      };
      return Tensor1_Expr<decltype(TensorExpr), T, Dim3, l>{TensorExpr};
    }
  };

  /// levi_civita functions to make for easy adhoc use

  /// Normally, this should go in is own file levi_civita.hpp, but
  /// not all filesystems can handle files that differ only in case :(

  /// Rank 2
  template <class T = int, char i, char j, int Dim0, int Dim1>
  constexpr typename std::enable_if<
    (Dim0 <= 2 && Dim1 <= 2),
    Tensor2_Expr<Levi_Civita<T>, T, Dim0, Dim1, i, j>>::type
  levi_civita(const Index<i, Dim0> &, const Index<j, Dim1> &)
  {
    return Levi_Civita<T>()(Index<i, Dim0>(), Index<j, Dim1>());
  }

  template <class T = int, char i, int Dim0>
  constexpr auto levi_civita(const Index<i, Dim0> &, const int &N1)
  {
    return Levi_Civita<T>()(Index<i, Dim0>(), N1);
  }

  template <class T = int, char j, int Dim1>
  constexpr auto levi_civita(const int &N0, const Index<j, Dim1> &)
  {
    return Levi_Civita<T>()(N0, Index<j, Dim1>());
  }

  /// Rank 3
  template <class T = int, char i, char j, char k, int Dim0, int Dim1, int Dim2>
  constexpr typename std::enable_if<
    (Dim0 <= 3 && Dim1 <= 3 && Dim2 <= 3),
    Tensor3_Expr<Levi_Civita<T>, T, Dim0, Dim1, Dim2, i, j, k>>::type
  levi_civita(const Index<i, Dim0> &, const Index<j, Dim1> &,
              const Index<k, Dim2> &)
  {
    return Levi_Civita<T>()(Index<i, Dim0>(), Index<j, Dim1>(),
                            Index<k, Dim2>());
  }

  template <class T = int, char i, char j, int Dim0, int Dim1>
  constexpr auto
  levi_civita(const Index<i, Dim0> &, const Index<j, Dim1> &, const int &N2)
  {
    return Levi_Civita<T>()(Index<i, Dim0>(), Index<j, Dim1>(), N2);
  }

  template <class T = int, char i, char k, int Dim0, int Dim2>
  constexpr auto
  levi_civita(const Index<i, Dim0> &, const int &N1, const Index<k, Dim2> &)
  {
    return Levi_Civita<T>()(Index<i, Dim0>(), N1, Index<k, Dim2>());
  }

  template <class T = int, char j, char k, int Dim1, int Dim2>
  constexpr auto
  levi_civita(const int &N0, const Index<j, Dim1> &, const Index<k, Dim2> &)
  {
    return Levi_Civita<T>()(N0, Index<j, Dim1>(), Index<k, Dim2>());
  }

  template <class T = int, char i, int Dim0>
  constexpr auto
  levi_civita(const Index<i, Dim0> &, const int &N1, const int &N2)
  {
    return Levi_Civita<T>()(Index<i, Dim0>(), N1, N2);
  }

  template <class T = int, char j, int Dim1>
  constexpr auto
  levi_civita(const int &N0, const Index<j, Dim1> &, const int &N2)
  {
    return Levi_Civita<T>()(N0, Index<j, Dim1>(), N2);
  }

  template <class T = int, char k, int Dim2>
  constexpr auto
  levi_civita(const int &N0, const int &N1, const Index<k, Dim2> &)
  {
    return Levi_Civita<T>()(N0, N1, Index<k, Dim2>());
  }

  /// Rank 4
  template <class T = int, char i, char j, char k, char l, int Dim0, int Dim1,
            int Dim2, int Dim3>
  constexpr typename std::enable_if<
    (Dim0 <= 4 && Dim1 <= 4 && Dim2 <= 4 && Dim3 <= 4),
    Tensor4_Expr<Levi_Civita<T>, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l>>::type
  levi_civita(const Index<i, Dim0> &, const Index<j, Dim1> &,
              const Index<k, Dim2> &, const Index<l, Dim3> &)
  {
    return Levi_Civita<T>()(Index<i, Dim0>(), Index<j, Dim1>(),
                            Index<k, Dim2>(), Index<l, Dim3>());
  }

  template <class T = int, char i, char j, char k, int Dim0, int Dim1, int Dim2>
  constexpr auto levi_civita(const Index<i, Dim0> &, const Index<j, Dim1> &,
                             const Index<k, Dim2> &, const int &N3)
  {
    return Levi_Civita<T>()(Index<i, Dim0>(), Index<j, Dim1>(),
                            Index<k, Dim2>(), N3);
  }

  template <class T = int, char i, char j, char l, int Dim0, int Dim1, int Dim3>
  constexpr auto levi_civita(const Index<i, Dim0> &, const Index<j, Dim1> &,
                             const int &N2, const Index<l, Dim3> &)
  {
    return Levi_Civita<T>()(Index<i, Dim0>(), Index<j, Dim1>(), N2,
                            Index<l, Dim3>());
  }

  template <class T = int, char i, char k, char l, int Dim0, int Dim2, int Dim3>
  constexpr auto levi_civita(const Index<i, Dim0> &, const int &N1,
                             const Index<k, Dim2> &, const Index<l, Dim3> &)
  {
    return Levi_Civita<T>()(Index<i, Dim0>(), N1, Index<k, Dim2>(),
                            Index<l, Dim3>());
  }

  template <class T = int, char j, char k, char l, int Dim1, int Dim2, int Dim3>
  constexpr auto levi_civita(const int &N0, const Index<j, Dim1> &,
                             const Index<k, Dim2> &, const Index<l, Dim3> &)
  {
    return Levi_Civita<T>()(N0, Index<j, Dim1>(), Index<k, Dim2>(),
                            Index<l, Dim3>());
  }

  template <class T = int, char i, char j, int Dim0, int Dim1>
  constexpr auto levi_civita(const Index<i, Dim0> &, const Index<j, Dim1> &,
                             const int &N2, const int &N3)
  {
    return Levi_Civita<T>()(Index<i, Dim0>(), Index<j, Dim1>(), N2, N3);
  }

  template <class T = int, char i, char k, int Dim0, int Dim2>
  constexpr auto levi_civita(const Index<i, Dim0> &, const int &N1,
                             const Index<k, Dim2> &, const int &N3)
  {
    return Levi_Civita<T>()(Index<i, Dim0>(), N1, Index<k, Dim2>(), N3);
  }

  template <class T = int, char j, char k, int Dim1, int Dim2>
  constexpr auto levi_civita(const int &N0, const Index<j, Dim1> &,
                             const Index<k, Dim2> &, const int &N3)
  {
    return Levi_Civita<T>()(N0, Index<j, Dim1>(), Index<k, Dim2>(), N3);
  }

  template <class T = int, char i, char l, int Dim0, int Dim3>
  constexpr auto levi_civita(const Index<i, Dim0> &, const int &N1,
                             const int &N2, const Index<l, Dim3> &)
  {
    return Levi_Civita<T>()(Index<i, Dim0>(), N1, N2, Index<l, Dim3>());
  }

  template <class T = int, char j, char l, int Dim1, int Dim3>
  constexpr auto levi_civita(const int &N0, const Index<j, Dim1> &,
                             const int &N2, const Index<l, Dim3> &)
  {
    return Levi_Civita<T>()(N0, Index<j, Dim1>(), N2, Index<l, Dim3>());
  }

  template <class T = int, char k, char l, int Dim2, int Dim3>
  constexpr auto levi_civita(const int &N0, const int &N1,
                             const Index<k, Dim2> &, const Index<l, Dim3> &)
  {
    return Levi_Civita<T>()(N0, N1, Index<k, Dim2>(), Index<l, Dim3>());
  }

  template <class T = int, char i, int Dim0>
  constexpr auto levi_civita(const Index<i, Dim0> &, const int &N1,
                             const int &N2, const int &N3)
  {
    return Levi_Civita<T>()(Index<i, Dim0>(), N1, N2, N3);
  }

  template <class T = int, char j, int Dim1>
  constexpr auto levi_civita(const int &N0, const Index<j, Dim1> &,
                             const int &N2, const int &N3)
  {
    return Levi_Civita<T>()(N0, Index<j, Dim1>(), N2, N3);
  }

  template <class T = int, char k, int Dim2>
  constexpr auto levi_civita(const int &N0, const int &N1,
                             const Index<k, Dim2> &, const int &N3)
  {
    return Levi_Civita<T>()(N0, N1, Index<k, Dim2>(), N3);
  }

  template <class T = int, char l, int Dim3>
  constexpr auto levi_civita(const int &N0, const int &N1, const int &N2,
                             const Index<l, Dim3> &)
  {
    return Levi_Civita<T>()(N0, N1, N2, Index<l, Dim3>());
  }
}
