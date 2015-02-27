/*
 * Vector class in 3D space
 */
class Vector {
private:
  float r[3];
public:
  Vector();
  Vector(float x, float y, float z);
private:
};

/*
 * 3x3 Matrix class
 */
class Matrix {
private:
  float* r[3];
public:
  Matrix();
  ~Matrix();
private:
};