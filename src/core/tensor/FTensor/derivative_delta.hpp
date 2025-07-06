template <int Dim> class derivative_delta
{
public:
  Tensor1<int, Dim> d_ijk_plus, d_ijk_minus;
  Tensor1<double, Dim> d_xyz;

  derivative_delta(double dx, double dy, double dz, int di, int dj, int dk)
      : d_ijk_plus(di, dj, dk), d_ijk_minus(di, dj, dk), d_xyz(dx, dy, dz)
  {}

  derivative_delta(double dx, double dy, double dz, int di_plus, int dj_plus,
                   int dk_plus, int di_minus, int dj_minus, int dk_minus)
      : d_ijk_plus(di_plus, dj_plus, dk_plus),
        d_ijk_minus(di_minus, dj_minus, dk_minus), d_xyz(dx, dy, dz)
  {}
};
