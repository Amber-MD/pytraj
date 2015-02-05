#ifndef INC_VEC3_H
#define INC_VEC3_H
/// Designed to hold array of size 3 (like XYZ coord etc).
class Vec3 {
  public:
    Vec3() { }
    Vec3(const Vec3& rhs) {
      V_[0] = rhs.V_[0];
      V_[1] = rhs.V_[1];
      V_[2] = rhs.V_[2];
    }
    Vec3(double vx, double vy, double vz) {
      V_[0] = vx;
      V_[1] = vy;
      V_[2] = vz;
    }
    Vec3(double vxyz) { V_[0] = vxyz; V_[1] = vxyz; V_[2] = vxyz; }
    Vec3(const double* XYZ) {
      V_[0] = XYZ[0];
      V_[1] = XYZ[1];
      V_[2] = XYZ[2];
    }
    Vec3(const float* XYZ) {
      V_[0] = (double)XYZ[0];
      V_[1] = (double)XYZ[1];
      V_[2] = (double)XYZ[2];
    }
    Vec3(const int* XYZ) {
      V_[0] = (double)XYZ[0];
      V_[1] = (double)XYZ[1];
      V_[2] = (double)XYZ[2];
    }
    /// Assignment
    Vec3& operator=(const Vec3& rhs) {
      if (this == &rhs) return *this;
      V_[0] = rhs.V_[0];
      V_[1] = rhs.V_[1];
      V_[2] = rhs.V_[2];
      return *this;
    }
    /// Assign const double[3] values to X Y and Z
    void Assign(const double* XYZ) { V_[0] = XYZ[0]; V_[1] = XYZ[1]; V_[2] = XYZ[2]; }
    // Vector OP scalar
    void operator/=(double xIn) {
      V_[0] /= xIn;
      V_[1] /= xIn;
      V_[2] /= xIn;
    }
    Vec3 operator/(double xIn) const {
      return Vec3( V_[0] / xIn, V_[1] / xIn, V_[2] / xIn);
    }
    void operator*=(double xIn) {
      V_[0] *= xIn;
      V_[1] *= xIn;
      V_[2] *= xIn;
    }
    Vec3 operator*(double xIn) const {
      return Vec3( V_[0] * xIn, V_[1] * xIn, V_[2] * xIn);
    }
    void operator+=(double xIn) {
      V_[0] += xIn;
      V_[1] += xIn;
      V_[2] += xIn;
    }
    Vec3 operator+(double xIn) const {
      return Vec3( V_[0] + xIn, V_[1] + xIn, V_[2] + xIn);
    }
    // Vector OP vector
    void operator-=(const Vec3& rhs) {
      V_[0] -= rhs.V_[0];
      V_[1] -= rhs.V_[1];
      V_[2] -= rhs.V_[2];
    }
    Vec3 operator-(const Vec3& rhs) const {
      return Vec3(V_[0]-rhs.V_[0], V_[1]-rhs.V_[1], V_[2]-rhs.V_[2]);
    }
    void operator+=(const Vec3& rhs) {
      V_[0] += rhs.V_[0];
      V_[1] += rhs.V_[1];
      V_[2] += rhs.V_[2];
    }
    Vec3 operator+(const Vec3& rhs) const {
      return Vec3(V_[0]+rhs.V_[0], V_[1]+rhs.V_[1], V_[2]+rhs.V_[2]);
    }
    double operator*(const Vec3& rhs) const { // Dot product
      return ( (V_[0]*rhs.V_[0]) + (V_[1]*rhs.V_[1]) + (V_[2]*rhs.V_[2]) );
    }
    Vec3 operator/(const Vec3& rhs) const {
      return (Vec3(V_[0]/rhs.V_[0], V_[1]/rhs.V_[1], V_[2]/rhs.V_[2]));
    }
    Vec3 Cross(Vec3 const& rhs) const { // Cross product
      return Vec3( (V_[1]*rhs.V_[2]) - (V_[2]*rhs.V_[1]),   // UYVZ+UZVY
                   (V_[2]*rhs.V_[0]) - (V_[0]*rhs.V_[2]),   // UZVX+UXVZ
                   (V_[0]*rhs.V_[1]) - (V_[1]*rhs.V_[0]) ); // UXVY+UYVX
    }
    // TODO: Make const ref only?
    double  operator[](int idx) const { return V_[idx]; }
    double& operator[](int idx)       { return V_[idx]; }
    double Magnitude2() const {
      double x = V_[0] * V_[0];
      double y = V_[1] * V_[1];
      double z = V_[2] * V_[2];
      return (x + y + z);
    }
    void Zero() {
      V_[0] = 0.0;
      V_[1] = 0.0;
      V_[2] = 0.0;
    }
    bool IsZero() const {
      return (V_[0]==0.0 && V_[1]==0.0 && V_[2]==0.0);
    }
    void Neg() {
      V_[0] = -V_[0];
      V_[1] = -V_[1];
      V_[2] = -V_[2];
    }
    void SetVec(double vx, double vy, double vz) {
      V_[0] = vx;
      V_[1] = vy;
      V_[2] = vz;
    }
    double Normalize();
    double Length() const;
    void Print(const char*) const;
    double Angle(Vec3 const&) const;
    double SignedAngle(Vec3 const&, Vec3 const&) const;
    // TODO: Eliminate this routine
    const double* Dptr() const { return V_; }
    double* Dptr() { return V_; }
  private:
    double V_[3];
};
#endif
