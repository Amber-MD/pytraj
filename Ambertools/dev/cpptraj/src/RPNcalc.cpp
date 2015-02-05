#include <cmath>
#include <sstream>
#include <locale>
#include <stack>
#include "RPNcalc.h"
#include "DataSet_1D.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI

// CONSTRUCTOR
RPNcalc::RPNcalc() {}

static inline bool isOpChar(char cIn) {
  return ( cIn == '(' || cIn == ')' || cIn == '+' || cIn == '-' ||
           cIn == '/' || cIn == '*' || cIn == '^' || cIn == '=');
}

std::string const& RPNcalc::FirstTokenName() const {
  return tokens_.front().Name();
} 

/** Convert infix expression to RPN in tokens_ array. This uses a
  * shunting-yard algorithm which has been slightly modified to
  * recognize unary right-associative operators.
  */
int RPNcalc::ProcessExpression(std::string const& expression) {
  std::locale loc;
  if (expression.empty()) return 1;
  if (debug_ > 0) mprintf("Parsing expression: '%s'\n", expression.c_str());
  tokens_.clear();
  std::stack<Token> op_stack;
  std::string::const_iterator ptr = expression.begin();
  bool lastTokenWasOperator = true;
  while ( ptr != expression.end() )
  {
    //mprintf("DEBUG: Start of loop, char is '%c'\n", *ptr);
    // Skip whitespace
    if (isspace(*ptr, loc)) { ++ptr; continue; }

    // NUMBER ------------------------------------
    if (*ptr == '.' || isdigit(*ptr, loc))
    { // Start of a number
      std::string number;
      std::string exponent;
      int hasExponent = 0;
      bool decimal_point = (*ptr == '.');
      if (decimal_point)
        number.push_back( *(ptr++) );
      while ( ptr != expression.end() && isdigit(*ptr, loc) )
      {
        if (hasExponent != 0)
          exponent.push_back( *(ptr++) );
        else
          number.push_back( *(ptr++) );
        // Check the next character
        if (*ptr == '.')
        {
          if (hasExponent != 0)
          {
            mprinterr("Error: Decimal point not allowed in exponent.\n");
            return 1;
          }
          else if (decimal_point)
          {
            mprinterr("Error: Two decimal points encountered in number: %s\n", number.c_str());
            return 1;
          }
          decimal_point = true;
          number.push_back( *(ptr++) );
        }
        else if (*ptr == 'E' || *ptr == 'e')
        {
          ++ptr;
          if (*ptr == '-')
          {
            hasExponent = -1;
            ++ptr;
          }
          else
            hasExponent = 1;
        } 
      }
      if (debug_ > 0) mprintf("Number detected: %s\n", number.c_str());
      std::istringstream iss(number);
      double val;
      if (!(iss >> val)) {
        mprinterr("Error: Invalid number: %s\n", number.c_str());
        return 1;
      }
      if (hasExponent != 0)
      {
        if (exponent.empty()) {
          mprinterr("Error: Exponent is empty.\n");
          return 1;
        }
        if (debug_ > 0) mprintf("Exponent detected: %s\n", exponent.c_str());
        double eval;
        std::istringstream iss2(exponent);
        if (!(iss2 >> eval)) {
          mprinterr("Error: Invalid exponent: %s\n", exponent.c_str());
          return 1;
        }
        if (hasExponent < 0)
          eval = -eval;
        val *= pow(10, eval);
      }
      tokens_.push_back( Token( val ) );
      lastTokenWasOperator = false;
    }
    // ALPHA (FUNCTION/VAR) ----------------------
    else if ( isalpha(*ptr, loc) )
    { // Look for Function name
      size_t pos = ptr - expression.begin();
      if (expression.compare(pos, 4, "sqrt")==0)
      {
        op_stack.push( Token(FN_SQRT) );
        ptr += 4;
        lastTokenWasOperator = true;
      }
      else if (expression.compare(pos, 3, "exp")==0)
      {
        op_stack.push( Token(FN_EXP) );
        ptr += 3;
        lastTokenWasOperator = true;
      }
      else if (expression.compare(pos, 2, "ln")==0)
      {
        op_stack.push( Token(FN_LN) );
        ptr += 2;
        lastTokenWasOperator = true;
      }
      else if (expression.compare(pos, 3, "abs")==0)
      {
        op_stack.push( Token(FN_ABS) );
        ptr += 3;
        lastTokenWasOperator = true;
      }
      else if (expression.compare(pos, 3, "sin")==0)
      {
        op_stack.push( Token(FN_SIN) );
        ptr += 3;
        lastTokenWasOperator = true;
      }
      else if (expression.compare(pos, 3, "cos")==0)
      {
        op_stack.push( Token(FN_COS) );
        ptr += 3;
        lastTokenWasOperator = true;
      }
      else if (expression.compare(pos, 3, "tan")==0)
      {
        op_stack.push( Token(FN_TAN) );
        ptr += 3;
        lastTokenWasOperator = true;
      }
      else if (expression.compare(pos, 3, "sum")==0)
      {
        op_stack.push( Token(FN_SUM) );
        ptr += 3;
        lastTokenWasOperator = true;
      }
      else if (expression.compare(pos, 3, "avg")==0)
      {
        op_stack.push( Token(FN_AVG) );
        ptr += 3;
        lastTokenWasOperator = true;
      }
      else if (expression.compare(pos, 5, "stdev")==0)
      {
        op_stack.push( Token(FN_STDEV) );
        ptr += 5;
        lastTokenWasOperator = true;
      }
      else if (expression.compare(pos, 3, "min")==0)
      {
        op_stack.push( Token(FN_MIN) );
        ptr += 3;
        lastTokenWasOperator = true;
      } 
      else if (expression.compare(pos, 3, "max")==0)
      {
        op_stack.push( Token(FN_MAX) );
        ptr += 3;
        lastTokenWasOperator = true;
      } 
      // -----------------------------------------
      else if (expression.compare(pos, 2, "PI")==0)
      {
        tokens_.push_back( Token( Constants::PI ) );
        ptr += 2;
        lastTokenWasOperator = false;
      }
      // -----------------------------------------
      else
      { // Assume variable name.
        std::string varname;
        bool has_colon = false; // For index
        enum BracketState { NONE, OPEN, CLOSED };
        BracketState bracket = NONE;
        while (ptr != expression.end() && !isOpChar(*ptr) && !isspace(*ptr,loc))
        {
          //mprintf("DEBUG: Var '%c'\n", *ptr);
          varname.push_back( *(ptr++) );
          // Check for brackets (Aspect)
          if (*ptr == '[' || *ptr == ']')
          {
            if (bracket == CLOSED)
            {
              mprinterr("Error: Multiple bracket sets encountered in set name: %s\n",
                        varname.c_str());
              return 1;
            }
            else if (*ptr == ']')
            {
              if (bracket == NONE)
              {
                mprinterr("Error: ']' encountered before '[' in set name: %s\n", varname.c_str());
                return 1;
              }
              bracket = CLOSED;
            }
            else if (*ptr == '[')
            {
              if (bracket == OPEN)
              {
                mprinterr("Error: Multiple '[' encountered in set name: %s\n", varname.c_str());
                return 1;
              }
              bracket = OPEN;
            }
            varname.push_back( *(ptr++) );
          }
          // Colon (index) can be after bracket or after name
          if (*ptr == ':')
          {
            if (has_colon)
            {
              mprinterr("Error: Two colons encountered in set name: %s\n", varname.c_str());
              return 1;
            }
            has_colon = true;
            varname.push_back( *(ptr++) );
          }
          //mprintf("DEBUG: Post Var '%c'\n", *ptr);
        }
        if (debug_>0) mprintf("Variable name detected: %s\n", varname.c_str());
        tokens_.push_back( Token(varname) );
        lastTokenWasOperator = false;
      }
    }
    // -----------------------------------------
    else if (*ptr == '(')
    {
      op_stack.push( Token(LPAR) );
      ++ptr;
    }
    else if (*ptr == ')')
    {
      if (op_stack.empty()) {
        mprinterr("Error: Mismatched parentheses.\n");
        return 1;
      }
      while (!op_stack.empty() && op_stack.top().Type() != LPAR)
      {
        tokens_.push_back( op_stack.top() );
        op_stack.pop();
      }
      if (op_stack.empty() || op_stack.top().Type() != LPAR) {
        mprinterr("Error: Mismatched parentheses.\n");
        return 1;
      }
      // Pop left parentheses off the stack.
      op_stack.pop();
      // Put function on top of stack in output
      if (!op_stack.empty() && op_stack.top().IsFunction()) {
        tokens_.push_back( op_stack.top() );
        op_stack.pop();
      }
      ++ptr;
    }
    // -----------------------------------------
    else
    { // Operator
      Token O1;
      if (*ptr == '+')
        O1.SetType(OP_PLUS);
      else if (*ptr == '-') {
        if (lastTokenWasOperator)
          O1.SetType(OP_NEG);
        else
          O1.SetType(OP_MINUS);
      } else if (*ptr == '*')
        O1.SetType(OP_MULT);
      else if (*ptr == '/')
        O1.SetType(OP_DIV);
      else if (*ptr == '^')
        O1.SetType(OP_POW);
      else if (*ptr == '=')
        O1.SetType(OP_ASSIGN);
      else {
        mprinterr("Error: Unrecognized character in expression: %c\n", *ptr);
        return 1;
      }
      ++ptr;
      while (!op_stack.empty() && op_stack.top().IsOperator())
      {
        if ((O1.IsLeftAssociative() && O1.Priority() <= op_stack.top().Priority()) ||
            (O1.Priority() < op_stack.top().Priority()))
        {
          tokens_.push_back( op_stack.top() );
          op_stack.pop();
        } else
          break;
      }
      if (debug_>0) mprintf("Operator detected: %s\n", O1.Description());
      lastTokenWasOperator = true;
      op_stack.push( O1 );
    }
  } // END loop over expression
  while (!op_stack.empty()) {
    if (op_stack.top().IsOperator()) {
      tokens_.push_back( op_stack.top() );
      op_stack.pop();
    } else {
      mprinterr("Error: Missing or mismatched parentheses.\n");
      return 1;
    }
  }
  if (debug_ > 0) {
    mprintf("Final RPN expression:\n");
    for (unsigned int i = 0; i != tokens_.size(); i++) {
      mprintf("\t%u: %s", i, tokens_[i].Description());
      if (tokens_[i].IsValue()) {
        if (tokens_[i].Type() == NUMBER)
          mprintf(" %lf", tokens_[i].Value());
        else
          mprintf(" %s", tokens_[i].name());
      }
      mprintf("\n");
    }
  }
  if (tokens_.empty()) {
    mprinterr("Error: No valid tokens detected.\n");
    return 1;
  }
  return 0;
}

double RPNcalc::DoOperation(double d1, double d2, TokenType op_type) {
  switch (op_type) {
    case OP_MINUS: return d2 - d1;
    case OP_PLUS: return d2 + d1;
    case OP_DIV: return d2 / d1;
    case OP_MULT: return d2 * d1;
    case OP_POW: return pow(d2, d1);
    case OP_NEG: return -d1;
    // ---------------------
    case FN_SQRT: return sqrt(d1);
    case FN_EXP: return exp(d1);
    case FN_LN: return log(d1);
    case FN_ABS: return fabs(d1);
    case FN_SIN: return sin(d1);
    case FN_COS: return cos(d1);
    case FN_TAN: return tan(d1);
    case FN_STDEV: return 0.0;
    case FN_SUM:
    case FN_AVG:
    case FN_MIN:
    case FN_MAX:
      return d1;
    default:
      mprinterr("Error: Invalid token type.\n");
  }
  return 0.0;
}

// RPNcalc::Evaluate()
int RPNcalc::Evaluate(DataSetList& DSL) const {
  if (tokens_.empty()) {
    mprinterr("Error: Expression was not set.\n");
    return 1;
  }
  std::stack<ValType> Stack;
  ValType Dval[2]; // NOTE: Must be able to hold max # operands.
  DataSetList LocalList;
  // Are we going to be assigning this?
  bool assigningResult = false;
  AssignType assignStatus = AssignStatus();
  if (assignStatus == ERR_ASSIGN)
    return 1;
  else if (assignStatus == YES_ASSIGN)
    assigningResult = true;
  DataSet* output = 0;
  // Process RPN tokens. 
  for (Tarray::const_iterator T = tokens_.begin(); T != tokens_.end(); ++T)
  {
    if ( T->Type() == NUMBER )
      Stack.push( ValType(T->Value()) );
    else if ( T->Type() == VARIABLE ) {
      DataSet* ds = 0;
      if (assigningResult && T == tokens_.begin())
        ds = output; // NOTE: Will be '0' at this point.
      else {
        ds= DSL.GetDataSet( T->Name() );
        if (ds == 0) {
          mprinterr("Error: Data set with name '%s' not found.\n", T->name());
          return 1;
        }
      }
      Stack.push( ValType( ds ) );
    } else {
      Dval[0].Reset();
      Dval[1].Reset();
      // Operand or function. Get operand(s)
      unsigned int nOps = (unsigned int)T->numOperands();
      if (Stack.size() < nOps) {
        mprinterr("Error: Not enough operands for '%s'.\n", T->Description());
        return 1;
      }
      for (unsigned int i = 0; i != nOps; i++) {
        Dval[i] = Stack.top();
        Stack.pop();
        // Replace 1D datasets of size 1 with the actual value.
        if (Dval[i].IsDataSet() && Dval[i].DS()!=0 && // Probably being assigned to if '0'
            Dval[i].DS()->Ndim()==1 && Dval[i].DS()->Size()==1)
          Dval[i].SetValue(((DataSet_1D*)Dval[i].DS())->Dval(0));
      }
      if (T->Type() == OP_ASSIGN) {
        // Assignment. This should be the last operation.
        if (!assigningResult) {
          mprinterr("Error: Assignment must be the final operation.\n");
          return 1;
        }
        if (!Dval[1].IsDataSet()) {
          mprinterr("Error: Attempting to assign to something that is not a data set.\n");
          return 1;
        }
        if (Dval[1].DS() != output) { // NOTE: Should be '0'
          mprinterr("Internal Error: Assigning to wrong data set!\n");
          return 1;
        }
        output = DSL.AddSet(DataSet::DOUBLE, tokens_.front().Name(), "CALC");
        if (output == 0) return 1;
        if (Dval[0].IsDataSet()) {
          if (debug_>0)
            mprintf("DEBUG: Assigning '%s' to '%s'\n", Dval[0].DS()->Legend().c_str(),
                    output->Legend().c_str());
          // Should be 1D by definition, allocated below in LocalList
          DataSet_1D const& D1 = static_cast<DataSet_1D const&>( *Dval[0].DS() );
          for (unsigned int n = 0; n != D1.Size(); n++) {
            double dval = D1.Dval(n); // TODO: Direct copy
            output->Add(n, &dval);
          }
        } else {
          if (debug_>0)
            mprintf("DEBUG: Assigning %f to '%s'\n", Dval[0].Value(),
                    output->Legend().c_str());
          double dval = Dval[0].Value();
          output->Add(0, &dval); 
        }
        Stack.push(ValType(output));
      } else if (!Dval[0].IsDataSet() && !Dval[1].IsDataSet()) {
        // Neither operand is a data set
        if (debug_>0)
          mprintf("DEBUG: '%f' [%s] '%f'\n", Dval[1].Value(), T->Description(), Dval[0].Value());
        Stack.push(ValType(DoOperation(Dval[0].Value(), Dval[1].Value(), T->Type())));
      } else if (T->numOperands() == 1 && T->ResultIsScalar()) {
        // One operand that is a data set that will be converted to a scalar
        DataSet* ds1 = Dval[0].DS();
        if (debug_ > 0)
          mprintf("DEBUG: [%s] '%s'\n", T->Description(), ds1->Legend().c_str());
        if (ds1->Ndim() != 1) {
          mprinterr("Error: Data set math currently restricted to 1D data sets.\n");
          return 1;
        }
        DataSet_1D const& D1 = static_cast<DataSet_1D const&>( *ds1 );
        if (T->Type() == FN_SUM) {
          double sum = 0.0;
          for (unsigned int n = 0; n != D1.Size(); n++)
            sum += D1.Dval(n);
          Stack.push(ValType(sum));
        } else if (T->Type() == FN_STDEV) {
          double stdev;
          D1.Avg(stdev);
          Stack.push(ValType(stdev));
        } else if (T->Type() == FN_AVG)
          Stack.push(ValType(D1.Avg()));
        else if (T->Type() == FN_MIN)
          Stack.push(ValType(D1.Min()));
        else if (T->Type() == FN_MAX)
          Stack.push(ValType(D1.Max()));
        else {
          mprinterr("Internal Error: OP %s is undefined for data set.\n", T->Description());
          return 1;
        }
      } else {
        // One or both operands is a DataSet. Result is DataSet.
        // Set up temporary data set to hold result.
        DataSet* tempDS = LocalList.AddSetIdx(DataSet::DOUBLE, "TEMP", T-tokens_.begin());
        if (tempDS == 0) return 1;
        // Handle 2 or 1 operand
        if (T->numOperands() == 2) {
          if (Dval[0].IsDataSet() && Dval[1].IsDataSet()) {
            // Both are DataSets. Must have same size.
            DataSet* ds1 = Dval[0].DS();
            DataSet* ds2 = Dval[1].DS();
            if (debug_>0)
              mprintf("DEBUG: '%s' [%s] '%s' => '%s'\n", ds2->Legend().c_str(), T->Description(),
                      ds1->Legend().c_str(), tempDS->Legend().c_str());
            if (ds1->Size() != ds2->Size()) {
              mprinterr("Error: Sets '%s' and '%s' do not have same size, required for %s\n",
                        ds1->Legend().c_str(), ds2->Legend().c_str(), T->name());
              return 1;
            }
            if (ds1->Ndim() != 1 && ds2->Ndim() != 1) {
              mprinterr("Error: Data set math currently restricted to 1D data sets.\n");
              return 1;
            }
            DataSet_1D const& D1 = static_cast<DataSet_1D const&>( *ds1 );
            DataSet_1D const& D2 = static_cast<DataSet_1D const&>( *ds2 );
            for (unsigned int n = 0; n != D1.Size(); n++) {
              double dval = DoOperation(D1.Dval(n), D2.Dval(n), T->Type());
              tempDS->Add(n, &dval);
            }
          } else {
            // DataSet OP Value or Value OP DataSet
            if (Dval[0].IsDataSet()) {
              // DataSet OP Value
              DataSet* ds1 = Dval[0].DS();
              if (debug_ > 0)
                mprintf("DEBUG: %f [%s] '%s' => '%s'\n", Dval[1].Value(), T->Description(),
                        ds1->Legend().c_str(), tempDS->Legend().c_str());
              if (ds1->Ndim() != 1) {
                mprinterr("Error: Data set math currently restricted to 1D data sets.\n");
                return 1;
              }
              DataSet_1D const& D1 = static_cast<DataSet_1D const&>( *ds1 );
              double d2 = Dval[1].Value();
              for (unsigned int n = 0; n != D1.Size(); n++) {
                double dval = DoOperation(D1.Dval(n), d2, T->Type());
                tempDS->Add(n, &dval);
              }
            } else {
              // Value OP DataSet
              DataSet* ds2 = Dval[1].DS();
              if (debug_ > 0)
                mprintf("DEBUG: '%s' [%s] '%f' => '%s'\n", ds2->Legend().c_str(), T->Description(),
                        Dval[0].Value(), tempDS->Legend().c_str());
              if (ds2->Ndim() != 1) {
                mprinterr("Error: Data set math currently restricted to 1D data sets.\n");
                return 1;
              }
              DataSet_1D const& D2 = static_cast<DataSet_1D const&>( *ds2 );
              double d1 = Dval[0].Value();
              for (unsigned int n = 0; n != D2.Size(); n++) {
                double dval = DoOperation(d1, D2.Dval(n), T->Type());
                tempDS->Add(n, &dval);
              }
            }
          }
        } else {
          // Only 1 operand and it is a DataSet
          DataSet* ds1 = Dval[0].DS();
          if (debug_ > 0)
            mprintf("DEBUG: [%s] '%s' => '%s'\n", T->Description(),
                    ds1->Legend().c_str(), tempDS->Legend().c_str());
          if (ds1->Ndim() != 1) {
            mprinterr("Error: Data set math currently restricted to 1D data sets.\n");
            return 1;
          }
          DataSet_1D const& D1 = static_cast<DataSet_1D const&>( *ds1 );
          for (unsigned int n = 0; n != D1.Size(); n++) {
            double dval = DoOperation(D1.Dval(n), 0.0, T->Type());
            tempDS->Add(n, &dval);
          }
        }
        Stack.push(ValType(tempDS));
      }
    }
  }
  if (Stack.size() != 1) {
    mprinterr("Error: Unbalanced expression.\n");
    return 1;
  }
  if (output == 0)
    mprintf("Result: %f\n", Stack.top().Value());
  else
    mprintf("Result stored in '%s'\n", output->Legend().c_str());
  return 0;
}

// RPNcalc::AssignStatus()
RPNcalc::AssignType RPNcalc::AssignStatus() const {
  AssignType assigningResult = NO_ASSIGN;
  if (tokens_.front().IsValue() && tokens_.back().Type() == OP_ASSIGN) {
    if (tokens_.size() < 3) {
      mprinterr("Error: Cannot assign nothing.\n");
      return ERR_ASSIGN;
    }
    if (tokens_.front().Type() != VARIABLE) {
      mprinterr("Error: Must assign to a data set on left hand side.\n");
      return ERR_ASSIGN;
    }
    assigningResult = YES_ASSIGN;
  }
  return assigningResult;
}

// RPNcalc::Nparams()
int RPNcalc::Nparams() const {
  int nparams=0, min_param=-1, max_param=-1;
  bool hasXvar = false;
  for (Tarray::const_iterator T = tokens_.begin(); T != tokens_.end(); ++T)
    if (T->Type() == VARIABLE) {
      if (T->Name()[0] == 'A')
      {
        std::istringstream iss( T->Name().substr(1) );      
        int pnum;
        if (!(iss >> pnum)) {
          mprinterr("Error: Invalid parameter number: %s\n", T->Name().substr(1).c_str());
          return 1;
        }
        if (min_param ==-1 || pnum < min_param) min_param = pnum;
        if (max_param ==-1 || pnum > max_param) max_param = pnum;
        nparams++;
      }
      else if (T->Name() == "X")
        hasXvar = true;
    }
  if (!hasXvar) {
    mprinterr("Error: No X variable in equation.\n");
    return -1;
  }
  if (nparams > 0 && min_param != 0) {
    mprinterr("Error: Minimum paramter is not A0.\n");
    return -1;
  }
  if (nparams > 0 && max_param != nparams-1) {
    mprinterr("Error: %i parameters detected but max parameter is not A%i\n", nparams, max_param);
    return -1;
  }
  return nparams;
}

// RPNcalc::Evaluate()
/** This version of evaluate requires assignment. Should be checked with
  * AssignStatus() prior to call.
  */
int RPNcalc::Evaluate(Darray const& Params, double X, double& Result) const {
  if (tokens_.empty()) {
    mprinterr("Error: Expression was not set.\n");
    return 1;
  }
  std::stack<double> Stack;
  double Dval[2] = {0.0, 0.0}; // NOTE: Must be able to hold max # operands.
    
  for (Tarray::const_iterator T = tokens_.begin(); T != tokens_.end(); ++T)
  {
    if ( T->Type() == NUMBER )
      Stack.push( T->Value() );
    else if ( T->Type() == VARIABLE ) {
      double param = 0.0;
      if (T != tokens_.begin()) { // First var will be output variable, set to 0.0.
        if (T->Name()[0] == 'A') {
          // Find parameter An, where n is parameter position.
          std::istringstream iss( T->Name().substr(1) );
          int nparam;
          if (!(iss >> nparam)) {
            mprinterr("Error: Invalid parameter number: %s\n", T->Name().substr(1).c_str());
            return 1;
          }
          // NOTE: NO CHECK FOR OUT OF BOUNDS.
          param = Params[nparam]; 
        } else if (T->Name()[0] == 'X') {
          param = X;
        } else {
          mprinterr("Error: Invalid variable '%s'. Expect 'X' or 'A<n>'\n", T->Name().c_str());
          return 1;
        }
      }
      Stack.push( param );
    } else {
      // Operand or function. Get operand(s)
      unsigned int nOps = (unsigned int)T->numOperands();
      if (Stack.size() < nOps) {
        mprinterr("Error: Not enough operands for '%s'.\n", T->Description());
        return 1;
      }
      for (unsigned int i = 0; i != nOps; i++) {
        Dval[i] = Stack.top();
        Stack.pop();
      }
      if (T->Type() == OP_ASSIGN) {
        // Assignment. This should be the last operation.
        Stack.push( Dval[0]  );
      } else if (T->numOperands() == 1 && T->ResultIsScalar()) {
        // One operand that is a data set that will be converted to a scalar.
        // Not allowed with this version of Evaluate.
        mprinterr("Error: '%s': Data Set functions not allowed in these equations.\n",
                  T->Description());
        return 1;
      } else {
        // Perform operation. 
        if (debug_ > 0)
          mprintf("DEBUG: '%f' [%s] '%f'\n", Dval[1], T->Description(), Dval[0]);
        Stack.push( DoOperation(Dval[0], Dval[1], T->Type()) );
      }
    }
  }
  if (Stack.size() != 1) {
    mprinterr("Error: Unbalanced expression.\n");
    return 1;
  }
  Result = Stack.top();
  if (debug_ > 0)
    mprintf("Result: Y(%g)= %g\n", X, Result);

  return 0;
}

// -----------------------------------------------------------------------------
/// Priority, #operands, associativity, class, resultIsScalar, description.
const RPNcalc::OpType RPNcalc::Token::OpArray_[] = {
  { 0, 0, NO_A,  NO_C,  0, "None"        }, // NONE
  { 0, 0, NO_A,  VALUE, 0, "Number"      }, // NUMBER
  { 0, 0, NO_A,  VALUE, 0, "Variable"    }, // VARIABLE
  { 1, 2, LEFT,  OP,    0, "Minus"       }, // OP_MINUS
  { 1, 2, LEFT,  OP,    0, "Plus"        }, // OP_PLUS
  { 2, 2, LEFT,  OP,    0, "Divide"      }, // OP_DIV
  { 2, 2, LEFT,  OP,    0, "Multiply"    }, // OP_MULT
  { 3, 2, LEFT,  OP,    0, "Power"       }, // OP_POW
  { 4, 1, RIGHT, OP,    0, "Unary minus" }, // OP_NEG
  { 0, 2, RIGHT, OP,    0, "Assignment"  }, // OP_ASSIGN
  { 0, 1, NO_A,  FN,    0, "Square root" }, // FN_SQRT
  { 0, 1, NO_A,  FN,    0, "Exponential" }, // FN_EXP
  { 0, 1, NO_A,  FN,    0, "Natural log" }, // FN_LN
  { 0, 1, NO_A,  FN,    0, "Abs. value"  }, // FN_ABS
  { 0, 1, NO_A,  FN,    0, "Sine"        }, // FN_SIN
  { 0, 1, NO_A,  FN,    0, "Cosine"      }, // FN_COS
  { 0, 1, NO_A,  FN,    0, "Tangent"     }, // FN_TAN
  { 0, 1, NO_A,  FN,    1, "Summation"   }, // FN_SUM
  { 0, 1, NO_A,  FN,    1, "Average"     }, // FN_AVG
  { 0, 1, NO_A,  FN,    1, "Standard Dev"}, // FN_STDEV
  { 0, 1, NO_A,  FN,    1, "Minimum"     }, // FN_MIN
  { 0, 1, NO_A,  FN,    1, "Maximum"     }, // FN_MAX
  { 0, 0, NO_A,  NO_C,  0, "Left Par"    }, // LPAR
  { 0, 0, NO_A,  NO_C,  0, "Right Par"   }, // RPAR
};
