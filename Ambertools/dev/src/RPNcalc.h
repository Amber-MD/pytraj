#ifndef INC_RPNCALC_H
#define INC_RPNCALC_H
#include "DataSetList.h"
/// Reverse Polish notation calculator
class RPNcalc {
  public:
    typedef std::vector<double> Darray;
    RPNcalc();
    void SetDebug(int d) { debug_ = d; }
    /// Set RPNcalc with equation from expression.
    int ProcessExpression(std::string const&);
    /// Evaluate equation; variables = DataSets
    int Evaluate(DataSetList&) const;
    enum AssignType { NO_ASSIGN, YES_ASSIGN, ERR_ASSIGN };
    /// Report on whether equation contains an assignment or not.
    AssignType AssignStatus() const;
    /// Evaluate equation; variables are input parameters.
    int Evaluate(Darray const&, double, double&) const;
    /// \return Name of variable being assigned to.
    std::string const& FirstTokenName() const;
    /// \return Number of A<n> parameter variables, check validity.
    int Nparams() const;
  private:
    class Token;
    class ValType;
    enum TokenType { NONE = 0, NUMBER, VARIABLE,
                     // Left-associative operators 
                     OP_MINUS, OP_PLUS, OP_DIV, OP_MULT, OP_POW,
                     // Right-associative operators
                     OP_NEG, OP_ASSIGN,
                     // Functions
                     FN_SQRT, FN_EXP, FN_LN, FN_ABS,
                     // Trig functions
                     FN_SIN, FN_COS, FN_TAN,
                     // Functions that take a data set
                     FN_SUM, FN_AVG, FN_STDEV, FN_MIN, FN_MAX, 
                     // Parentheses (for infix conversion only)
                     LPAR, RPAR };
    enum Associativity { NO_A = 0, LEFT, RIGHT };
    enum OpClass { NO_C = 0, VALUE, OP, FN };
    /// Hold information about ops
    struct OpType {
      int priority_;
      int nOperands_;
      Associativity assoc_;
      OpClass opclass_;
      int resultIsScalar_;
      const char* description_;
    };

    static inline double DoOperation(double, double, TokenType);

    typedef std::vector<Token> Tarray;
    Tarray tokens_;
    int debug_;
};
/// Hold values/operators for RPN calculator.
class RPNcalc::Token {
  public:
    Token() : type_(NONE), value_(0.0) {}
    /// CONSTRUCTOR - Numerical values
    Token(double val) : type_(NUMBER), value_(val) {}
    /// CONSTRUCTOR - Variables
    Token(std::string const& s) : type_(VARIABLE), value_(0.0), name_(s) {}
    /// CONSTRUCTOR - Operators
    Token(TokenType t) : type_(t), value_(0.0) {}
    /// Set token type. Intended for use with operator.
    void SetType(TokenType t) { type_ = t; }
    /// \return Token type.
    TokenType Type() const { return type_; }
    /// \return Token value.
    double Value() const { return value_; }
    /// \return Token variable name.
    std::string const& Name() const { return name_; }
    const char* name() const { return name_.c_str(); }
    /// \return true if token is a variable or number
    inline bool IsValue() const { return (OpArray_[type_].opclass_ == VALUE); }
    /// \return true if token is an operator
    inline bool IsOperator() const { return (OpArray_[type_].opclass_ == OP); }
    /// \return true if token is a function.
    inline bool IsFunction() const { return (OpArray_[type_].opclass_ == FN); }
    /// \return true if OP is left associative.
    inline bool IsLeftAssociative() const { return (OpArray_[type_].assoc_ == LEFT); }
    /// \return true if OP( dataset ) returns a scalar.
    inline bool ResultIsScalar() const { return (bool)OpArray_[type_].resultIsScalar_; }
    /// \return string indicating token type.
    const char* Description() const { return OpArray_[type_].description_; }
    /// \return operator priority
    int Priority() const { return OpArray_[type_].priority_; }
    /// \return expected number of operands
    int numOperands() const { return OpArray_[type_].nOperands_; }
  private:
    static const OpType OpArray_[];
 
    TokenType type_;
    double value_; ///< Numerical value.
    std::string name_; ///< Variable name
};
/// Hold results from ongoing calculation.
class RPNcalc::ValType {
  public:
    ValType() : ds_(0), val_(0.0), isDataSet_(false) {}
    ValType(double dval) : ds_(0), val_(dval), isDataSet_(false) {}
    ValType(DataSet* dsIn) : ds_(dsIn), val_(0.0), isDataSet_(true) {}
    ValType(ValType const& rhs) : ds_(rhs.ds_), val_(rhs.val_),
                                  isDataSet_(rhs.isDataSet_) {}
    ValType& operator=(ValType const& rhs) {
      if (this != &rhs) {
        ds_ = rhs.ds_;
        val_ = rhs.val_;
        isDataSet_ = rhs.isDataSet_;
      }
      return *this;
    }
    bool IsDataSet() const { return isDataSet_; }
    double Value() const { return val_; }
    DataSet* DS() const { return ds_; }
    void Reset() { ds_=0; val_=0.0; isDataSet_=false; }
    void SetValue(double d) { ds_ = 0; val_ = d; isDataSet_ = false; }
  private:
    DataSet* ds_;
    double val_;
    bool isDataSet_;
};
#endif
