/*
  Class for three component rank 2 tensor
  including mathematical operations and IO
 */
#include <cmath>   
#include <iostream> 
#include <fstream> 


class ThreeMat
{
 private:
  double coord_[3][3]; // Private data members of 3x3 matrix
  
 public:
  // Default constructor
  ThreeMat() 
    {
      for (int i=0;i<3;++i)
	for (int j=0;j<3;++j)
	  coord_[i][j]=0.0;
    }
 
  // Dyadic constructor
  ThreeMat(ThreeVec a, ThreeVec b)
    {
      for (int i=0;i<3;++i)
	for (int j=0;j<3;++j)
	  coord_[i][j] = a.get(i)*b.get(j);
    }
  // Access function for ith,jth element
  double get(int i,int j)
  {
    return coord_[i][j];
  }
  
  // Modifier method for ith,jth element -> INSERT
  void set(int i,int j,double value)
  {
    coord_[i][j] = value;
  }
   // Alternative modifier method for ith,jth element -> ADD
  void inc(int i,int j,double value)
  {
    coord_[i][j] += value;
  } 
  /*
    Overload the operators +,-,* and to represent matrix operations
  */
  // Addition
  ThreeMat operator+(ThreeMat mat)
  {
    ThreeMat ans;
     for (int i=0;i<3;++i)
	for (int j=0;j<3;++j)
	  ans.set(i,j,coord_[i][j]+mat.get(i,j));
     return ans;
  }
  
  // Subtraction 
  ThreeMat operator-(ThreeMat mat)
  {
    ThreeMat ans;
    for (int i=0;i<3;++i)
      for (int j=0;j<3;++j)
	ans.set(i,j,coord_[i][j]-mat.get(i,j));
    return ans;
  }

  // Transpose matrix
  void transpose()
  {
    double temp_value;
    for (int i=0;i<3;++i)
      for (int j=i;j<3;++j) // note, starts from i, not 0
	{
	  temp_value = coord_[i][j];
	  coord_[i][j]=coord_[j][i];
	  coord_[j][i]=temp_value;
	}
  }
  // Shift the diagonal of the matrix by I*scalar, I-> identity matrix
  void shift(double val)
  {
    for (int i=0;i<3;++i)
      coord_[i][i]+=val;
  }
  // Matrix multiplication with matrix
  ThreeMat operator*(ThreeMat mat)
  {
    ThreeMat ans; // all elements initially 0.0
    for (int i=0;i<3;++i)
      for (int j=0;j<3;++j)
	for (int k=0;k<3;++k) // adding index
	  ans.inc(i,j,coord_[i][k]*mat.get(k,j));
    return ans;
  }
  
  // Matrix multiplication with vector
  ThreeVec operator*(ThreeVec vec)
  {
    ThreeVec ans; // all elements initially 0.0
    for (int i=0;i<3;++i)
      for (int k=0;k<3;++k) // adding index
	ans.inc(i,coord_[i][k]*vec.get(k));
    return ans;
  }
  // Matrix multiplication with scalar
  ThreeMat operator*(double val)
  {
    ThreeMat ans; // all elements initially 0.0
    for (int i=0;i<3;++i)
      for (int j=0;j<3;++j)
	ans.set(i,j,coord_[i][j]*val);
    return ans;
  }
  // Print to screen
  void print()
  {
    for (int i=0;i<3;++i){
      for (int j=0;j<3;++j)
	std::cout << coord_[i][j] << '\t';
      std::cout << '\n';
      }
  }
  
  // Print to file
  void print(std::ofstream &fout)
  {
    for (int i=0;i<3;++i){
      for (int j=0;j<3;++j)
	fout << coord_[i][j] << '\t';
      fout << '\n';
      }
  }
}; 
