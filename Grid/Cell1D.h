/*!\file Cell1D.h
  \brief Header file defining 1D cell types. */

#ifndef _CELL1D_INCLUDED
#define _CELL1D_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

/* Include math macro header file. */

#ifndef _MATH_MACROS_INCLUDED
#include "../Math.h"
#endif // _MATH_MACROS_INCLUDED

/* Define the classes. */

/*!
 * \class: Cell1D_Uniform                                
 * \brief Class definition of the 1D Uniform Cell.
 * 
 * \verbatim
 * Member functions                                     
 *      x       -- Return cell center.                  
 *      dx      -- Return cell length.                  
 *      setsize -- Set cell length.                     
 *      setloc  -- Set cell center location (position).  
 *                                                      
 * Member operators                                     
 *      C -- a cell                                     
 *                                                      
 * C = C;                                               
 * C == C;                                              
 * C != C;                                              
 * cout << C; (output function)                         
 * cin  >> C; (input function)                          
 * \endverbatim                                                     
 */
class Cell1D_Uniform{
private:
public:
  //! @name Cell parameters
  //@{
  double          x;   //!< Cell center.
  static double  dx;   //!< Cell length.
  //@}

  //! @name Constructors
  //@{
  //! Default constructor
  Cell1D_Uniform(void) {
    x = ONE;
  }
  //! Copy constructor
  Cell1D_Uniform(const Cell1D_Uniform &Cell) {
    x = Cell.x;
  }
  //! Value constructor
  Cell1D_Uniform(const double &location) {
    x = location;
  }
  //@}


  //! @name Set cell size
  //@{
  void setsize(void);
  void setsize(const double &length);
  //@}

  //! @name Set cell center location.
  //@{
  void setloc(void);
  void setloc(const Cell1D_Uniform &Cell);
  void setloc(const double &location);
  //@}

  //! Assignment operator.
  Cell1D_Uniform& operator =(const Cell1D_Uniform &Cell1);

  //! @name Relational operators. 
  //@{
  friend int operator ==(const Cell1D_Uniform &Cell1,
			 const Cell1D_Uniform &Cell2);
  friend int operator !=(const Cell1D_Uniform &Cell1,
			 const Cell1D_Uniform &Cell2);
  //@}

  //! @name Input-output operators.
  //@{
  friend ostream &operator << (ostream &out_file,
			       const Cell1D_Uniform &Cell);
  friend istream &operator >> (istream &in_file,
			       Cell1D_Uniform &Cell);
  //@}
};

/******************************************************//**
 * Cell1D_Uniform::setsize -- Set cell size.            *
 ********************************************************/
inline void Cell1D_Uniform::setsize(void) {
    dx = ONE;
}

inline void Cell1D_Uniform::setsize(const double &length) {
    dx = length;
}

/******************************************************//**
 * Cell1D_Uniform::setloc -- Set cell center location.  *
 ********************************************************/
inline void Cell1D_Uniform::setloc(void) {
    x = ONE;
}

inline void Cell1D_Uniform::setloc(const Cell1D_Uniform &Cell) {
    x = Cell.x;
}

inline void Cell1D_Uniform::setloc(const double &location) {
    x = location;
}

/******************************************************//**
 * Cell1D_Uniform -- Relational operators.              *
 ********************************************************/
inline Cell1D_Uniform& Cell1D_Uniform::operator =(const Cell1D_Uniform &Cell1) {
    x = Cell1.x;
    dx = Cell1.dx;
    return (*this);
}

/******************************************************//**
 * Cell1D_Uniform -- Relational operators.              *
 ********************************************************/
inline int operator ==(const Cell1D_Uniform &Cell1,
		       const Cell1D_Uniform &Cell2) {
    return (Cell1.x == Cell2.x);
}

inline int operator !=(const Cell1D_Uniform &Cell1,
		       const Cell1D_Uniform &Cell2) {
    return (Cell1.x != Cell2.x);
}

/******************************************************//**
 * Cell1D_Uniform -- Input-output operators.            *
 ********************************************************/
inline ostream &operator << (ostream &out_file,
                             const Cell1D_Uniform &Cell) {
    out_file.setf(ios::scientific);
    out_file << " " << Cell.x << " " << Cell.dx;
    out_file.unsetf(ios::scientific);
    return (out_file);
}

inline istream &operator >> (istream &in_file,
                             Cell1D_Uniform &Cell) {
    in_file.setf(ios::skipws);
    in_file >> Cell.x;
    in_file.unsetf(ios::skipws);
    return (in_file);
}

/******************************************************//**
 * Useful 1D cell constants.                            *
 ********************************************************/
const Cell1D_Uniform Cell1D_Uniform_ONE(ONE);

#endif /* _CELL1D_INCLUDED  */
