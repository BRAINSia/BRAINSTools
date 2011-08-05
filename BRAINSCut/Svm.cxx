#include "Svm.h"
typedef float       Qfloat;
typedef signed char schar;
#ifndef USE_OPENCV // Following function conflict with OpenCV
#ifndef min
template <class T>
inline T min(T x, T y)
{
  return ( x < y ) ? x : y;
}

#endif
#ifndef max
template <class T>
inline T max(T x, T y)
{
  return ( x > y ) ? x : y;
}

#endif
template <class T>
inline void swap(T & x, T & y)
{
  T t = x; x = y; y = t;
}

#endif
template <class S, class T>
inline void clone(T * & dst, S const *src, const int n)
{
  dst = new T[n];
  memcpy( (void *)dst, (void const *)src, sizeof( T ) * n );
}

#define INF HUGE_VAL
#define TAU 1e-12
#define Malloc(type, n) (type *)malloc( ( n ) * sizeof( type ) )
class Cache
{
public:
  Cache(int l, int size);
  ~Cache();
  int get_data(const int index, Qfloat * *data, int len);

  void swap_index(int i, int j);  // future_option

private:
  int l;
  int size;
  struct head_t
    {
    head_t *prev, *next;  // a cicular list
    Qfloat *data;
    int len;    // data[0,len) is cached in this entry
    };
  head_t *head;
  head_t  lru_head;
  void lru_delete(head_t *h);

  void lru_insert(head_t *h);
};
Cache::Cache(int l_, int size_) : l(l_), size(size_)
{
  head = (head_t *)calloc( l, sizeof( head_t ) );  // initialized to 0
  size /= sizeof( Qfloat );
  size -= l * sizeof( head_t ) / sizeof( Qfloat );
  size = max( (int)size, 2 * l);  // cache must be large enough for two columns
  lru_head.next = lru_head.prev = &lru_head;
}

Cache::~Cache()
{
  for( head_t *h = lru_head.next; h != &lru_head; h = h->next )
    {
    free(h->data);
    }
  free(head);
}

void Cache::lru_delete(head_t *h)
{
  h->prev->next = h->next;
  h->next->prev = h->prev;
}

void Cache::lru_insert(head_t *h)
{
  h->next = &lru_head;
  h->prev = lru_head.prev;
  h->prev->next = h;
  h->next->prev = h;
}

int Cache::get_data(const int index, Qfloat * *data, int len)
{
  head_t *h = &head[index];

  if( h->len )
    {
    lru_delete(h);
    }
  int more = len - h->len;
  if( more > 0 )
    {
    while( size < more )
      {
      head_t *old = lru_head.next;
      lru_delete(old);
      free(old->data);
      size += old->len;
      old->data = 0;
      old->len = 0;
      }

    h->data = (Qfloat *)realloc(h->data, sizeof( Qfloat ) * len);
    size -= more;
    swap(h->len, len);
    }
  lru_insert(h);
  *data = h->data;
  return len;
}

void Cache::swap_index(int i, int j)
{
  if( i == j )
    {
    return;
    }
  if( head[i].len )
    {
    lru_delete(&head[i]);
    }
  if( head[j].len )
    {
    lru_delete(&head[j]);
    }
  swap(head[i].data, head[j].data);
  swap(head[i].len, head[j].len);
  if( head[i].len )
    {
    lru_insert(&head[i]);
    }
  if( head[j].len )
    {
    lru_insert(&head[j]);
    }
  if( i > j )
    {
    swap(i, j);
    }
  for( head_t *h = lru_head.next; h != &lru_head; h = h->next )
    {
    if( h->len > i )
      {
      if( h->len > j )
        {
        swap(h->data[i], h->data[j]);
        }
      else
        {
        // give up
        lru_delete(h);
        free(h->data);
        size += h->len;
        h->data = 0;
        h->len = 0;
        }
      }
    }
}

class QMatrix
{
public:
  virtual Qfloat * get_Q(int column, int len) const = 0;

  virtual Qfloat * get_QD() const = 0;

  virtual void swap_index(int i, int j) const = 0;

  virtual ~QMatrix()
  {
  }
};
class Kernel : public QMatrix
{
public:
  Kernel(int l, svm_node *const *x, const svm_parameter & param);
  virtual ~Kernel();
  static double k_function(const svm_node *x, const svm_node *y, const svm_parameter & param);

  virtual Qfloat * get_Q(int column, int len) const = 0;

  virtual Qfloat * get_QD() const = 0;

  virtual void swap_index(int i, int j) const
  {
    swap(x[i], x[j]);
    if( x_square )
      {
      swap(x_square[i], x_square[j]);
      }
  }

protected:
  double ( Kernel::*kernel_function )( int i, int j ) const;
private:
  const svm_node * *x;
  double *          x_square;
  const int         kernel_type;
  const double      degree;
  const double      gamma;
  const double      coef0;
  static double dot(const svm_node *px, const svm_node *py);

  double kernel_linear(int i, int j) const
  {
    return dot(x[i], x[j]);
  }

  double kernel_poly(int i, int j) const
  {
    return pow(gamma * dot(x[i], x[j]) + coef0, degree);
  }

  double kernel_rbf(int i, int j) const
  {
    return exp( -gamma * ( x_square[i] + x_square[j] - 2 * dot(x[i], x[j]) ) );
  }

  double kernel_sigmoid(int i, int j) const
  {
    return tanh(gamma * dot(x[i], x[j]) + coef0);
  }
};
Kernel::Kernel(int l, svm_node *const *x_, const svm_parameter & param) : kernel_type(param.kernel_type),
  degree(param.degree),
  gamma(param.gamma), coef0(param.coef0)
{
  kernel_function = &Kernel::kernel_rbf;
  clone(x, x_, l);
  x_square = new double[l];
  for( int i = 0; i < l; i++ )
    {
    x_square[i] = dot(x[i], x[i]);
    }
}

Kernel::~Kernel()
{
  delete[] x;
  delete[] x_square;
}

double Kernel::dot(const svm_node *px, const svm_node *py)
{
  double sum = 0;

  while( px->index != -1 && py->index != -1 )
    {
    if( px->index == py->index )
      {
      sum += px->value * py->value;
      ++px;
      ++py;
      }
    else
      {
      if( px->index > py->index )
        {
        ++py;
        }
      else
        {
        ++px;
        }
      }
    }

  return sum;
}

double Kernel::k_function(const svm_node *x, const svm_node *y,
                          const svm_parameter & param)
{
  double sum = 0;

  while( x->index != -1 && y->index != -1 )
    {
    if( x->index == y->index )
      {
      double d = x->value - y->value;
      sum += d * d;
      ++x;
      ++y;
      }
    else
      {
      if( x->index > y->index )
        {
        sum += y->value * y->value;
        ++y;
        }
      else
        {
        sum += x->value * x->value;
        ++x;
        }
      }
    }

  while( x->index != -1 )
    {
    sum += x->value * x->value;
    ++x;
    }

  while( y->index != -1 )
    {
    sum += y->value * y->value;
    ++y;
    }

  return exp(-param.gamma * sum);
}

class Solver
{
public:
  Solver() :
    active_size(0),
    y(0),
    G(0),            // gradient of objective function
    alpha_status(0), // LOWER_BOUND, UPPER_BOUND, FREE
    alpha(0),
    Q_Solver(0),
    QD_Solver(0),
    eps_Solver(0),
    Cp_Solver(0),
    Cn_Solver(0),
    b(0),
    active_set(0),
    G_bar(0),   // gradient, if we treat free variables as 0
    l_Solver(0),
    unshrinked(false)
  {
  }

  virtual ~Solver()
  {
  }

  struct SolutionInfo
    {
    double obj;
    double rho;
    double upper_bound_p;
    double upper_bound_n;
    //    double r;  // for Solver_NU
    };

  void Solve(int l, const QMatrix & Q, const double *b_, const schar *y_, double *alpha_, double Cp, double Cn,
             double eps, SolutionInfo *si, int shrinking);

protected:
  int     active_size;
  schar * y;
  double *G;    // gradient of objective function
  enum
    {
    LOWER_BOUND, UPPER_BOUND, FREE
    };
  char *         alpha_status; // LOWER_BOUND, UPPER_BOUND, FREE
  double *       alpha;
  const QMatrix *Q_Solver;
  const Qfloat * QD_Solver;
  double         eps_Solver;
  double         Cp_Solver, Cn_Solver;
  double *       b;
  int *          active_set;
  double *       G_bar; // gradient, if we treat free variables as 0
  int            l_Solver;
  bool           unshrinked;

  double get_C(int i)
  {
    return ( y[i] > 0 ) ? this->Cp_Solver : this->Cn_Solver;
  }

  void update_alpha_status(int i)
  {
    if( alpha[i] >= get_C(i) )
      {
      alpha_status[i] = UPPER_BOUND;
      }
    else if( alpha[i] <= 0 )
      {
      alpha_status[i] = LOWER_BOUND;
      }
    else
      {
      alpha_status[i] = FREE;
      }
  }

  bool is_upper_bound(int i)
  {
    return alpha_status[i] == UPPER_BOUND;
  }

  bool is_lower_bound(int i)
  {
    return alpha_status[i] == LOWER_BOUND;
  }

  bool is_free(int i)
  {
    return alpha_status[i] == FREE;
  }

  void swap_index(int i, int j);

  void reconstruct_gradient();

  virtual int select_working_set(int & i, int & j);

  virtual int max_violating_pair(int & i, int & j);

  virtual double calculate_rho();

  virtual void do_shrinking();
};
void Solver::swap_index(int i, int j)
{
  this->Q_Solver->swap_index(i, j);
  swap(y[i], y[j]);
  swap(G[i], G[j]);
  swap(alpha_status[i], alpha_status[j]);
  swap(alpha[i], alpha[j]);
  swap(b[i], b[j]);
  swap(active_set[i], active_set[j]);
  swap(G_bar[i], G_bar[j]);
}

void Solver::reconstruct_gradient()
{
  if( active_size == this->l_Solver )
    {
    return;
    }
  int i;
  for( i = active_size; i < this->l_Solver; i++ )
    {
    G[i] = G_bar[i] + b[i];
    }
  for( i = 0; i < active_size; i++ )
    {
    if( is_free(i) )
      {
      const Qfloat *Q_i = this->Q_Solver->get_Q(i, this->l_Solver);
      double        alpha_i = alpha[i];
      for( int j = active_size; j < this->l_Solver; j++ )
        {
        G[j] += alpha_i * Q_i[j];
        }
      }
    }
}

void Solver::Solve(int local_l,
                   const QMatrix & local_Q,
                   const double *b_,
                   const schar *y_,
                   double *alpha_,
                   double local_Cp,
                   double local_Cn,
                   double local_eps,
                   SolutionInfo *si,
                   int shrinking)
{
  this->l_Solver = local_l;
  this->Q_Solver = &local_Q;
  this->QD_Solver = local_Q.get_QD();
  clone(b, b_, local_l);
  clone(y, y_, local_l);
  clone(alpha, alpha_, local_l);
  this->Cp_Solver = local_Cp;
  this->Cn_Solver = local_Cn;
  this->eps_Solver = local_eps;
  unshrinked = false;
  alpha_status = new char[local_l];
  for( int i = 0; i < local_l; i++ )
    {
    update_alpha_status(i);
    }
  active_set = new int[local_l];
  for( int i = 0; i < local_l; i++ )
    {
    active_set[i] = i;
    }
  active_size = local_l;
  G = new double[local_l];
  G_bar = new double[local_l];
  for( int i = 0; i < local_l; i++ )
    {
    G[i] = b[i];
    G_bar[i] = 0;
    }
  for( int i = 0; i < local_l; i++ )
    {
    if( !is_lower_bound(i) )
      {
      const Qfloat *Q_i = local_Q.get_Q(i, local_l);
      double        alpha_i = alpha[i];
      int           j;
      for( j = 0; j < local_l; j++ )
        {
        G[j] += alpha_i * Q_i[j];
        }
      if( is_upper_bound(i) )
        {
        for( j = 0; j < local_l; j++ )
          {
          G_bar[j] += get_C(i) * Q_i[j];
          }
        }
      }
    }
  int iter = 0;
  int counter = min(local_l, 1000) + 1;
  while( 1 )
    {
    if( --counter == 0 )
      {
      counter = min(local_l, 1000);
      if( shrinking )
        {
        do_shrinking();
        }
      }

    int i, j;
    if( select_working_set(i, j) != 0 )
      {
      reconstruct_gradient();
      active_size = local_l;
      std::cout << "." << std::endl;
      if( select_working_set(i, j) != 0 )
        {
        break;
        }
      else
        {
        counter = 1;
        }
      }
    ++iter;
    const Qfloat *Q_i = local_Q.get_Q(i, active_size);
    const Qfloat *Q_j = local_Q.get_Q(j, active_size);
    double        C_i = get_C(i);
    double        C_j = get_C(j);
    double        old_alpha_i = alpha[i];
    double        old_alpha_j = alpha[j];
    if( y[i] != y[j] )
      {
      double quad_coef = Q_i[i] + Q_j[j] + 2 * Q_i[j];
      if( quad_coef <= 0 )
        {
        quad_coef = TAU;
        }
      double delta = ( -G[i] - G[j] ) / quad_coef;
      double diff = alpha[i] - alpha[j];
      alpha[i] += delta;
      alpha[j] += delta;
      if( diff > 0 )
        {
        if( alpha[j] < 0 )
          {
          alpha[j] = 0;
          alpha[i] = diff;
          }
        }
      else
        {
        if( alpha[i] < 0 )
          {
          alpha[i] = 0;
          alpha[j] = -diff;
          }
        }
      if( diff > C_i - C_j )
        {
        if( alpha[i] > C_i )
          {
          alpha[i] = C_i;
          alpha[j] = C_i - diff;
          }
        }
      else
        {
        if( alpha[j] > C_j )
          {
          alpha[j] = C_j;
          alpha[i] = C_j + diff;
          }
        }
      }
    else
      {
      double quad_coef = Q_i[i] + Q_j[j] - 2 * Q_i[j];
      if( quad_coef <= 0 )
        {
        quad_coef = TAU;
        }
      double delta = ( G[i] - G[j] ) / quad_coef;
      double sum = alpha[i] + alpha[j];
      alpha[i] -= delta;
      alpha[j] += delta;
      if( sum > C_i )
        {
        if( alpha[i] > C_i )
          {
          alpha[i] = C_i;
          alpha[j] = sum - C_i;
          }
        }
      else
        {
        if( alpha[j] < 0 )
          {
          alpha[j] = 0;
          alpha[i] = sum;
          }
        }
      if( sum > C_j )
        {
        if( alpha[j] > C_j )
          {
          alpha[j] = C_j;
          alpha[i] = sum - C_j;
          }
        }
      else
        {
        if( alpha[i] < 0 )
          {
          alpha[i] = 0;
          alpha[j] = sum;
          }
        }
      }
    double delta_alpha_i = alpha[i] - old_alpha_i;
    double delta_alpha_j = alpha[j] - old_alpha_j;
    for( int k = 0; k < active_size; k++ )
      {
      G[k] += Q_i[k] * delta_alpha_i + Q_j[k] * delta_alpha_j;
      }
      {
      bool ui = is_upper_bound(i);
      bool uj = is_upper_bound(j);
      update_alpha_status(i);
      update_alpha_status(j);
      int k;
      if( ui != is_upper_bound(i) )
        {
        Q_i = local_Q.get_Q(i, local_l);
        if( ui )
          {
          for( k = 0; k < local_l; k++ )
            {
            G_bar[k] -= C_i * Q_i[k];
            }
          }
        else
          {
          for( k = 0; k < local_l; k++ )
            {
            G_bar[k] += C_i * Q_i[k];
            }
          }
        }
      if( uj != is_upper_bound(j) )
        {
        Q_j = local_Q.get_Q(j, local_l);
        if( uj )
          {
          for( k = 0; k < local_l; k++ )
            {
            G_bar[k] -= C_j * Q_j[k];
            }
          }
        else
          {
          for( k = 0; k < local_l; k++ )
            {
            G_bar[k] += C_j * Q_j[k];
            }
          }
        }
      }
    }

  si->rho = calculate_rho();
  double v = 0;
  for( int i = 0; i < local_l; i++ )
    {
    v += alpha[i] * ( G[i] + b[i] );
    }
  si->obj = v / 2;
  for( int i = 0; i < local_l; i++ )
    {
    alpha_[active_set[i]] = alpha[i];
    }
  si->upper_bound_p = this->Cp_Solver;
  si->upper_bound_n = this->Cn_Solver;
  delete[] b;
  delete[] y;
  delete[] alpha;
  delete[] alpha_status;
  delete[] active_set;
  delete[] G;
  delete[] G_bar;
}

int Solver::select_working_set(int & out_i, int & out_j)
{
  double Gmax = -INF;
  double Gmax2 = -INF;
  int    Gmax_idx = -1;
  int    Gmin_idx = -1;
  double obj_diff_min = INF;

  for( int t = 0; t < active_size; t++ )
    {
    if( y[t] == +1 )
      {
      if( !is_upper_bound(t) )
        {
        if( -G[t] >= Gmax )
          {
          Gmax = -G[t];
          Gmax_idx = t;
          }
        }
      }
    else
      {
      if( !is_lower_bound(t) )
        {
        if( G[t] >= Gmax )
          {
          Gmax = G[t];
          Gmax_idx = t;
          }
        }
      }
    }
  int           i = Gmax_idx;
  const Qfloat *Q_i = NULL;
  if( i != -1 )
    {
    Q_i = this->Q_Solver->get_Q(i, active_size);
    }
  for( int j = 0; j < active_size; j++ )
    {
    if( y[j] == +1 )
      {
      if( !is_lower_bound(j) )
        {
        double grad_diff = Gmax + G[j];
        if( G[j] >= Gmax2 )
          {
          Gmax2 = G[j];
          }
        if( grad_diff > 0 )
          {
          double obj_diff;
          double quad_coef = Q_i[i] + this->QD_Solver[j] - 2 * y[i] * Q_i[j];
          if( quad_coef > 0 )
            {
            obj_diff = -( grad_diff * grad_diff ) / quad_coef;
            }
          else
            {
            obj_diff = -( grad_diff * grad_diff ) / TAU;
            }
          if( obj_diff <= obj_diff_min )
            {
            Gmin_idx = j;
            obj_diff_min = obj_diff;
            }
          }
        }
      }
    else
      {
      if( !is_upper_bound(j) )
        {
        double grad_diff = Gmax - G[j];
        if( -G[j] >= Gmax2 )
          {
          Gmax2 = -G[j];
          }
        if( grad_diff > 0 )
          {
          double obj_diff;
          double quad_coef = Q_i[i] + this->QD_Solver[j] + 2 * y[i] * Q_i[j];
          if( quad_coef > 0 )
            {
            obj_diff = -( grad_diff * grad_diff ) / quad_coef;
            }
          else
            {
            obj_diff = -( grad_diff * grad_diff ) / TAU;
            }
          if( obj_diff <= obj_diff_min )
            {
            Gmin_idx = j;
            obj_diff_min = obj_diff;
            }
          }
        }
      }
    }
  if( Gmax + Gmax2 < this->eps_Solver )
    {
    return 1;
    }
  out_i = Gmax_idx;
  out_j = Gmin_idx;
  return 0;
}

int Solver::max_violating_pair(int & out_i, int & out_j)
{
  double Gmax1 = -INF;    // max { -y_i * grad(f)_i | i in I_up(\alpha) }
  int    Gmax1_idx = -1;
  double Gmax2 = -INF;    // max { y_i * grad(f)_i | i in I_low(\alpha) }
  int    Gmax2_idx = -1;

  for( int i = 0; i < active_size; i++ )
    {
    if( y[i] == +1 )   // y = +1
      {
      if( !is_upper_bound(i) )   // d = +1
        {
        if( -G[i] >= Gmax1 )
          {
          Gmax1 = -G[i];
          Gmax1_idx = i;
          }
        }
      if( !is_lower_bound(i) )   // d = -1
        {
        if( G[i] >= Gmax2 )
          {
          Gmax2 = G[i];
          Gmax2_idx = i;
          }
        }
      }
    else    // y = -1
      {
      if( !is_upper_bound(i) )   // d = +1
        {
        if( -G[i] >= Gmax2 )
          {
          Gmax2 = -G[i];
          Gmax2_idx = i;
          }
        }
      if( !is_lower_bound(i) )   // d = -1
        {
        if( G[i] >= Gmax1 )
          {
          Gmax1 = G[i];
          Gmax1_idx = i;
          }
        }
      }
    }
  if( Gmax1 + Gmax2 < this->eps_Solver )
    {
    return 1;
    }
  out_i = Gmax1_idx;
  out_j = Gmax2_idx;
  return 0;
}

void Solver::do_shrinking()
{
  int i, j, k;

  if( max_violating_pair(i, j) != 0 )
    {
    return;
    }
  double Gm1 = -y[j] * G[j];
  double Gm2 = y[i] * G[i];
  for( k = 0; k < active_size; k++ )
    {
    if( is_lower_bound(k) )
      {
      if( y[k] == +1 )
        {
        if( -G[k] >= Gm1 )
          {
          continue;
          }
        }
      else
        {
        if( -G[k] >= Gm2 )
          {
          continue;
          }
        }
      }
    else if( is_upper_bound(k) )
      {
      if( y[k] == +1 )
        {
        if( G[k] >= Gm2 )
          {
          continue;
          }
        }
      else
        {
        if( G[k] >= Gm1 )
          {
          continue;
          }
        }
      }
    else
      {
      continue;
      }
    --active_size;
    swap_index(k, active_size);
    --k;  // look at the newcomer
    }
  if( unshrinked || -( Gm1 + Gm2 ) > this->eps_Solver * 10 )
    {
    return;
    }
  unshrinked = true;
  reconstruct_gradient();
  for( k = this->l_Solver - 1; k >= active_size; k-- )
    {
    if( is_lower_bound(k) )
      {
      if( y[k] == +1 )
        {
        if( -G[k] < Gm1 )
          {
          continue;
          }
        }
      else
        {
        if( -G[k] < Gm2 )
          {
          continue;
          }
        }
      }
    else if( is_upper_bound(k) )
      {
      if( y[k] == +1 )
        {
        if( G[k] < Gm2 )
          {
          continue;
          }
        }
      else
        {
        if( G[k] < Gm1 )
          {
          continue;
          }
        }
      }
    else
      {
      continue;
      }
    swap_index(k, active_size);
    active_size++;
    ++k;  // look at the newcomer
    }
}

double Solver::calculate_rho()
{
  double r;
  int    nr_free = 0;
  double ub = INF, lb = -INF, sum_free = 0;

  for( int i = 0; i < active_size; i++ )
    {
    double yG = y[i] * G[i];
    if( is_lower_bound(i) )
      {
      if( y[i] > 0 )
        {
        ub = min(ub, yG);
        }
      else
        {
        lb = max(lb, yG);
        }
      }
    else if( is_upper_bound(i) )
      {
      if( y[i] < 0 )
        {
        ub = min(ub, yG);
        }
      else
        {
        lb = max(lb, yG);
        }
      }
    else
      {
      ++nr_free;
      sum_free += yG;
      }
    }
  if( nr_free > 0 )
    {
    r = sum_free / nr_free;
    }
  else
    {
    r = ( ub + lb ) / 2;
    }
  return r;
}

class SVC_Q : public Kernel
{
public:
  SVC_Q(const svm_problem & prob, const svm_parameter & param, const schar *y_) : Kernel(prob.l_svm, prob.x_svm, param),
    y(0)
  {
    clone(y, y_, prob.l_svm);
    cache = new Cache( prob.l_svm, (int)( param.cache_size * ( 1 << 20 ) ) );
    QD = new Qfloat[prob.l_svm];
    for( int i = 0; i < prob.l_svm; i++ )
      {
      QD[i] = (Qfloat)( this->*kernel_function )(i, i);
      }
  }

  Qfloat * get_Q(int i, int len) const
  {
    Qfloat *data;
    int     start;

    if( ( start = cache->get_data(i, &data, len) ) < len )
      {
      for( int j = start; j < len; j++ )
        {
        data[j] = (Qfloat)( y[i] * y[j] * ( this->*kernel_function )(i, j) );
        }
      }
    return data;
  }

  Qfloat * get_QD() const
  {
    return QD;
  }

  void swap_index(int i, int j) const
  {
    cache->swap_index(i, j);
    Kernel::swap_index(i, j);
    swap(y[i], y[j]);
    swap(QD[i], QD[j]);
  }

  ~SVC_Q()
  {
    delete[] y;
    delete cache;
    delete[] QD;
  }

private:
  schar * y;
  Cache * cache;
  Qfloat *QD;
};
static void solve_c_svc(const svm_problem *prob, const svm_parameter *param,
                        double *alpha, Solver::SolutionInfo *si, double Cp, double Cn)
{
  int     l = prob->l_svm;
  double *minus_ones = new double[l];
  schar * y = new schar[l];
  int     i;

  for( i = 0; i < l; i++ )
    {
    alpha[i] = 0;
    minus_ones[i] = -1;
    if( prob->y_svm[i] > 0 )
      {
      y[i] = +1;
      }
    else
      {
      y[i] = -1;
      }
    }
  Solver s;
  s.Solve(l, SVC_Q(*prob,
                   *param,
                   y), minus_ones, y, alpha, Cp, Cn, param->eps, si, param->shrinking);
/*
  double sum_alpha = 0;
  for( i = 0; i < l; i++ )
    {
    sum_alpha += alpha[i];
    }
*/
  for( i = 0; i < l; i++ )
    {
    alpha[i] *= y[i];
    }
  delete[] minus_ones;
  delete[] y;
}

struct decision_function
  {
  double *alpha;
  double rho;
  };
decision_function svm_train_one(const svm_problem *prob,
                                const svm_parameter *param,
                                double Cp,
                                double Cn)
{
  double *             alpha = Malloc(double, prob->l_svm);
  Solver::SolutionInfo si;

  solve_c_svc(prob, param, alpha, &si, Cp, Cn);
  int nSV = 0;
  int nBSV = 0;
  for( int i = 0; i < prob->l_svm; i++ )
    {
    if( fabs(alpha[i]) > 0 )
      {
      ++nSV;
      if( prob->y_svm[i] > 0 )
        {
        if( fabs(alpha[i]) >= si.upper_bound_p )
          {
          ++nBSV;
          }
        }
      else
        {
        if( fabs(alpha[i]) >= si.upper_bound_n )
          {
          ++nBSV;
          }
        }
      }
    }
  decision_function f;
  f.alpha = alpha;
  f.rho = si.rho;
  return f;
}

struct svm_model
  {
  svm_parameter param;    // parameter
  int l;                  // total #SV
  svm_node * *SV;         // SVs (SV[l])
  double * *sv_coef;      // coefficients for SVs in decision functions
                          // (sv_coef[n-1][l])
  double *rho;            // constants in decision functions (rho[n*(n-1)/2])
  double *probA;          // pariwise probability information
  double *probB;
  int *label;     // label of each class (label[n])
  int *nSV;       // number of SVs for each class (nSV[n]), nSV[0] + nSV[1] +
                  // ... + nSV[n-1] = l
  int free_sv;    // 1 if svm_model is created by svm_load_model, 0 if svm_model
                  // is created by svm_train
  };
void sigmoid_train(int l,
                   const double *dec_values,
                   const double *labels,
                   double & A,
                   double & B)
{
  double prior1 = 0, prior0 = 0;
  int    i;

  for( i = 0; i < l; i++ )
    {
    if( labels[i] > 0 )
      {
      prior1 += 1;
      }
    else
      {
      prior0 += 1;
      }
    }
  int     max_iter = 100;   // Maximal number of iterations
  double  min_step = 1e-10; // Minimal step taken in line search
  double  sigma = 1e-3;     // For numerically strict PD of Hessian
  double  eps = 1e-5;
  double  hiTarget = ( prior1 + 1.0 ) / ( prior1 + 2.0 );
  double  loTarget = 1 / ( prior0 + 2.0 );
  double *t = Malloc(double, l);
  double  fApB, p, q, h11, h22, h21, g1, g2, det, dA, dB, gd, stepsize;
  double  newA, newB, newf, d1, d2;
  int     iter;
  A = 0.0;
  B = log( ( prior0 + 1.0 ) / ( prior1 + 1.0 ) );
  double fval = 0.0;
  for( i = 0; i < l; i++ )
    {
    if( labels[i] > 0 )
      {
      t[i] = hiTarget;
      }
    else
      {
      t[i] = loTarget;
      }
    fApB = dec_values[i] * A + B;
    if( fApB >= 0 )
      {
      fval += t[i] * fApB + log( 1 + exp(-fApB) );
      }
    else
      {
      fval += ( t[i] - 1 ) * fApB + log( 1 + exp(fApB) );
      }
    }
  for( iter = 0; iter < max_iter; iter++ )
    {
    h11 = sigma; // numerically ensures strict PD
    h22 = sigma;
    h21 = 0.0; g1 = 0.0; g2 = 0.0;
    for( i = 0; i < l; i++ )
      {
      fApB = dec_values[i] * A + B;
      if( fApB >= 0 )
        {
        p = exp(-fApB) / ( 1.0 + exp(-fApB) );
        q = 1.0 / ( 1.0 + exp(-fApB) );
        }
      else
        {
        p = 1.0 / ( 1.0 + exp(fApB) );
        q = exp(fApB) / ( 1.0 + exp(fApB) );
        }
      d2 = p * q;
      h11 += dec_values[i] * dec_values[i] * d2;
      h22 += d2;
      h21 += dec_values[i] * d2;
      d1 = t[i] - p;
      g1 += dec_values[i] * d1;
      g2 += d1;
      }
    if( fabs(g1) < eps && fabs(g2) < eps )
      {
      break;
      }
    det = h11 * h22 - h21 * h21;
    dA = -( h22 * g1 - h21 * g2 ) / det;
    dB = -( -h21 * g1 + h11 * g2 ) / det;
    gd = g1 * dA + g2 * dB;
    stepsize = 1;     // Line Search
    while( stepsize >= min_step )
      {
      newA = A + stepsize * dA;
      newB = B + stepsize * dB;
      newf = 0.0;
      for( i = 0; i < l; i++ )
        {
        fApB = dec_values[i] * newA + newB;
        if( fApB >= 0 )
          {
          newf += t[i] * fApB + log( 1 + exp(-fApB) );
          }
        else
          {
          newf += ( t[i] - 1 ) * fApB + log( 1 + exp(fApB) );
          }
        }
      if( newf < fval + 0.0001 * stepsize * gd )
        {
        A = newA;
        B = newB;
        fval = newf;
        break;
        }
      else
        {
        stepsize = stepsize / 2.0;
        }
      }

    if( stepsize < min_step )
      {
      std::cout << "Line search fails in two-class probability estimates\n"
                << std::endl;
      break;
      }
    }
  if( iter >= max_iter )
    {
    std::cout
      << "Reaching maximal iterations in two-class probability estimates\n"
      << std::endl;
    }
  free(t);
}

double sigmoid_predict(double decision_value, double A, double B)
{
  double fApB = decision_value * A + B;

  if( fApB >= 0 )
    {
    return exp(-fApB) / ( 1.0 + exp(-fApB) );
    }
  else
    {
    return 1.0 / ( 1 + exp(fApB) );
    }
}

void multiclass_probability(int k, double * *r, double *p)
{
  int       t, j;
  int       iter = 0, max_iter = max(100, k);
  double * *Q = Malloc(double *, k);
  double *  Qp = Malloc(double, k);
  double    pQp, eps = 0.005 / k;

  for( t = 0; t < k; t++ )
    {
    p[t] = 1.0 / k;  // Valid if k = 1
    Q[t] = Malloc(double, k);
    Q[t][t] = 0;
    for( j = 0; j < t; j++ )
      {
      Q[t][t] += r[j][t] * r[j][t];
      Q[t][j] = Q[j][t];
      }
    for( j = t + 1; j < k; j++ )
      {
      Q[t][t] += r[j][t] * r[j][t];
      Q[t][j] = -r[j][t] * r[t][j];
      }
    }
  for( iter = 0; iter < max_iter; iter++ )
    {
    pQp = 0;
    for( t = 0; t < k; t++ )
      {
      Qp[t] = 0;
      for( j = 0; j < k; j++ )
        {
        Qp[t] += Q[t][j] * p[j];
        }
      pQp += p[t] * Qp[t];
      }
    double max_error = 0;
    for( t = 0; t < k; t++ )
      {
      double error = fabs(Qp[t] - pQp);
      if( error > max_error )
        {
        max_error = error;
        }
      }
    if( max_error < eps )
      {
      break;
      }
    for( t = 0; t < k; t++ )
      {
      double diff = ( -Qp[t] + pQp ) / Q[t][t];
      p[t] += diff;
      pQp =
        ( pQp + diff
          * ( diff * Q[t][t] + 2 * Qp[t] ) ) / ( 1 + diff ) / ( 1 + diff );
      for( j = 0; j < k; j++ )
        {
        Qp[j] = ( Qp[j] + diff * Q[t][j] ) / ( 1 + diff );
        p[j] /= ( 1 + diff );
        }
      }
    }
  if( iter >= max_iter )
    {
    std::cout << "Exceeds max_iter in multiclass_prob\n" << std::endl;
    }
  for( t = 0; t < k; t++ )
    {
    free(Q[t]);
    }
  free(Q);
  free(Qp);
}

void svm_binary_svc_probability(const svm_problem *prob,
                                const svm_parameter *param,
                                double Cp,
                                double Cn,
                                double & probA,
                                double & probB)
{
  int     i;
  int *   perm = Malloc(int, prob->l_svm);
  double *dec_values = Malloc(double, prob->l_svm);

  for( i = 0; i < prob->l_svm; i++ )
    {
    perm[i] = i;
    }
  for( i = 0; i < prob->l_svm; i++ )
    {
    int j = i + rand() % ( prob->l_svm - i );
    swap(perm[i], perm[j]);
    }
  for( i = 0; i < 5; i++ )
    {
    int                begin = i * prob->l_svm / 5;
    int                end = ( i + 1 ) * prob->l_svm / 5;
    int                j, k;
    struct svm_problem subprob;
    subprob.l_svm = prob->l_svm - ( end - begin );
    subprob.x_svm = Malloc(struct svm_node *, subprob.l_svm);
    subprob.y_svm = Malloc(double, subprob.l_svm);
    k = 0;
    for( j = 0; j < begin; j++ )
      {
      subprob.x_svm[k] = prob->x_svm[perm[j]];
      subprob.y_svm[k] = prob->y_svm[perm[j]];
      ++k;
      }
    for( j = end; j < prob->l_svm; j++ )
      {
      subprob.x_svm[k] = prob->x_svm[perm[j]];
      subprob.y_svm[k] = prob->y_svm[perm[j]];
      ++k;
      }
    int p_count = 0, n_count = 0;
    for( j = 0; j < k; j++ )
      {
      if( subprob.y_svm[j] > 0 )
        {
        p_count++;
        }
      else
        {
        n_count++;
        }
      }
    if( p_count == 0 && n_count == 0 )
      {
      for( j = begin; j < end; j++ )
        {
        dec_values[perm[j]] = 0;
        }
      }
    else if( p_count > 0 && n_count == 0 )
      {
      for( j = begin; j < end; j++ )
        {
        dec_values[perm[j]] = 1;
        }
      }
    else if( p_count == 0 && n_count > 0 )
      {
      for( j = begin; j < end; j++ )
        {
        dec_values[perm[j]] = -1;
        }
      }
    else
      {
      svm_parameter subparam = *param;
      subparam.probability = 0;
      subparam.C = 1.0;
      subparam.nr_weight = 2;
      subparam.weight_label = Malloc(int, 2);
      subparam.weight = Malloc(double, 2);
      subparam.weight_label[0] = +1;
      subparam.weight_label[1] = -1;
      subparam.weight[0] = Cp;
      subparam.weight[1] = Cn;
      struct svm_model *submodel = svm_train(&subprob, &subparam);
      for( j = begin; j < end; j++ )
        {
        svm_predict_values( submodel, prob->x_svm[perm[j]],
                            &( dec_values[perm[j]] ) );
        dec_values[perm[j]] *= submodel->label[0];
        }
      svm_destroy_model(submodel);
      svm_destroy_param(&subparam);
      free(subprob.x_svm);
      free(subprob.y_svm);
      }
    }
  sigmoid_train(prob->l_svm, dec_values, prob->y_svm, probA, probB);
  free(dec_values);
  free(perm);
}

void svm_group_classes(const svm_problem *prob,
                       int * *label_ret,
                       int * *start_ret,
                       int * *count_ret,
                       int *perm)
{
  int  l = prob->l_svm;
  int  max_nr_class = 16;
  int  nr_class = 0;
  int *label = Malloc(int, max_nr_class);
  int *count = Malloc(int, max_nr_class);
  int *data_label = Malloc(int, l);
  int  i;

  for( i = 0; i < l; i++ )
    {
    int this_label = (int)prob->y_svm[i];
    int j;
    for( j = 0; j < nr_class; j++ )
      {
      if( this_label == label[j] )
        {
        ++count[j];
        break;
        }
      }
    data_label[i] = j;
    if( j == nr_class )
      {
      if( nr_class == max_nr_class )
        {
        max_nr_class *= 2;
        label = (int *)realloc( label, max_nr_class * sizeof( int ) );
        count = (int *)realloc( count, max_nr_class * sizeof( int ) );
        }
      label[nr_class] = this_label;
      count[nr_class] = 1;
      ++nr_class;
      }
    }
  int *start = Malloc(int, nr_class);
  start[0] = 0;
  for( i = 1; i < nr_class; i++ )
    {
    start[i] = start[i - 1] + count[i - 1];
    }
  for( i = 0; i < l; i++ )
    {
    perm[start[data_label[i]]] = i;
    ++start[data_label[i]];
    }
  start[0] = 0;
  for( i = 1; i < nr_class; i++ )
    {
    start[i] = start[i - 1] + count[i - 1];
    }
  *label_ret = label;
  *start_ret = start;
  *count_ret = count;
  free(data_label);
}

svm_model * svm_train(const svm_problem *prob, const svm_parameter *param)
{
  svm_model *model = Malloc(svm_model, 1);

  model->param = *param;
  model->free_sv = 0;
  int  l = prob->l_svm;
  int *label = NULL;
  int *start = NULL;
  int *count = NULL;
  int *perm = Malloc(int, l);
  svm_group_classes(prob, &label, &start, &count, perm);
  svm_node * *x = Malloc(svm_node *, l);
  int         i;
  for( i = 0; i < l; i++ )
    {
    x[i] = prob->x_svm[perm[i]];
    }
  double *weighted_C = Malloc(double, 2);
  for( i = 0; i < 2; i++ )
    {
    weighted_C[i] = param->C;
    }
  for( i = 0; i < param->nr_weight; i++ )
    {
    int j;
    for( j = 0; j < 2; j++ )
      {
      if( param->weight_label[i] == label[j] )
        {
        break;
        }
      }
    if( j == 2 )
      {
      fprintf(stderr,
              "warning: class label %d specified in weight is not found\n",
              param->weight_label[i]);
      }
    else
      {
      weighted_C[j] *= param->weight[i];
      }
    }
  bool *nonzero = Malloc(bool, l);
  for( i = 0; i < l; i++ )
    {
    nonzero[i] = false;
    }
  decision_function *f = Malloc(decision_function, 1);
  double *           probA = NULL, *probB = NULL;
  if( param->probability )
    {
    probA = Malloc(double, 1);
    probB = Malloc(double, 1);
    }
  svm_problem sub_prob;
  int         si = start[0], sj = start[1];
  int         ci = count[0], cj = count[1];
  sub_prob.l_svm = ci + cj;
  sub_prob.x_svm = Malloc(svm_node *, sub_prob.l_svm);
  sub_prob.y_svm = Malloc(double, sub_prob.l_svm);
  int k;
  for( k = 0; k < ci; k++ )
    {
    sub_prob.x_svm[k] = x[si + k];
    sub_prob.y_svm[k] = +1;
    }
  for( k = 0; k < cj; k++ )
    {
    sub_prob.x_svm[ci + k] = x[sj + k];
    sub_prob.y_svm[ci + k] = -1;
    }
  if( param->probability )
    {
    svm_binary_svc_probability(&sub_prob,
                               param,
                               weighted_C[0],
                               weighted_C[1],
                               probA[0],
                               probB[0]);
    }
  f[0] = svm_train_one(&sub_prob, param, weighted_C[0], weighted_C[1]);
  for( k = 0; k < ci; k++ )
    {
    if( !nonzero[si + k] && fabs(f[0].alpha[k]) > 0 )
      {
      nonzero[si + k] = true;
      }
    }
  for( k = 0; k < cj; k++ )
    {
    if( !nonzero[sj + k] && fabs(f[0].alpha[ci + k]) > 0 )
      {
      nonzero[sj + k] = true;
      }
    }
  free(sub_prob.x_svm);
  free(sub_prob.y_svm);
  model->label = Malloc(int, 2);
  for( i = 0; i < 2; i++ )
    {
    model->label[0] = label[0];
    }
  model->rho = Malloc(double, 1);
  model->rho[0] = f[0].rho;
  if( param->probability )
    {
    model->probA = Malloc(double, 1);
    model->probB = Malloc(double, 1);
    model->probA[0] = probA[0];
    model->probB[0] = probB[0];
    }
  else
    {
    model->probA = NULL;
    model->probB = NULL;
    }
  int  total_sv = 0;
  int *nz_count = Malloc(int, 2);
  model->nSV = Malloc(int, 2);
  for( i = 0; i < 2; i++ )
    {
    int nSV = 0;
    for( int j = 0; j < count[i]; j++ )
      {
      if( nonzero[start[i] + j] )
        {
        ++nSV;
        ++total_sv;
        }
      }
    model->nSV[i] = nSV;
    nz_count[i] = nSV;
    }
  model->l = total_sv;
  model->SV = Malloc(svm_node *, total_sv);
  int p = 0;
  for( i = 0; i < l; i++ )
    {
    if( nonzero[i] )
      {
      model->SV[p++] = x[i];
      }
    }
  int *nz_start = Malloc(int, 2);
  nz_start[0] = 0;
  for( i = 1; i < 2; i++ )
    {
    nz_start[i] = nz_start[i - 1] + nz_count[i - 1];
    }
  model->sv_coef = Malloc(double *, 1);
  model->sv_coef[0] = Malloc(double, total_sv);
  p = 0;
  si = start[0];
  sj = start[1];
  ci = count[0];
  cj = count[1];
  int q = nz_start[0];
  for( k = 0; k < ci; k++ )
    {
    if( nonzero[si + k] )
      {
      model->sv_coef[0][q++] = f[0].alpha[k];
      }
    }
  q = nz_start[1];
  for( k = 0; k < cj; k++ )
    {
    if( nonzero[sj + k] )
      {
      model->sv_coef[0][q++] = f[0].alpha[ci + k];
      }
    }
  free(label);
  free(probA);
  free(probB);
  free(count);
  free(perm);
  free(start);
  free(x);
  free(weighted_C);
  free(nonzero);
  free(f[0].alpha);
  free(f);
  free(nz_count);
  free(nz_start);
  return model;
}

void svm_predict_values(const svm_model *model,
                        const svm_node *x,
                        double *dec_values)
{
  int     i;
  int     l = model->l;
  double *kvalue = Malloc(double, l);

  for( i = 0; i < l; i++ )
    {
    kvalue[i] = Kernel::k_function(x, model->SV[i], model->param);
    }
  int *start = Malloc(int, 2);
  start[0] = 0;
  start[1] = start[0] + model->nSV[0];
  // int p=0;
  int     pos = 0;
  double  sum = 0;
  int     si = start[0];
  int     sj = start[1];
  int     ci = model->nSV[0];
  int     cj = model->nSV[1];
  int     k;
  double *coef1 = model->sv_coef[0];
  double *coef2 = model->sv_coef[0];
  for( k = 0; k < ci; k++ )
    {
    sum += coef1[si + k] * kvalue[si + k];
    }
  for( k = 0; k < cj; k++ )
    {
    sum += coef2[sj + k] * kvalue[sj + k];
    }
  sum -= model->rho[0];
  dec_values[pos++] = sum;
  free(kvalue);
  free(start);
}

double svm_predict_probability(const svm_model *model,
                               const svm_node *x,
                               double *prob_estimates)
{
  int     i;
  double *dec_values = Malloc(double, 1);

  svm_predict_values(model, x, dec_values);
  double    min_prob = 1e-7;
  double * *pairwise_prob = Malloc(double *, 2);
  for( i = 0; i < 2; i++ )
    {
    pairwise_prob[i] = Malloc(double, 2);
    }
  int k = 0;
  for( i = 0; i < 2; i++ )
    {
    for( int j = i + 1; j < 2; j++ )
      {
      pairwise_prob[i][j] =
        min(max(sigmoid_predict(dec_values[k], model->probA[k],
                                model->probB[k]), min_prob), 1 - min_prob);
      pairwise_prob[j][i] = 1 - pairwise_prob[i][j];
      k++;
      }
    }
  multiclass_probability(2, pairwise_prob, prob_estimates);
  int prob_max_idx = 0;
  for( i = 1; i < 2; i++ )
    {
    if( prob_estimates[i] > prob_estimates[prob_max_idx] )
      {
      prob_max_idx = i;
      }
    }
  for( i = 0; i < 2; i++ )
    {
    free(pairwise_prob[i]);
    }
  free(dec_values);
  free(pairwise_prob);
  return model->label[prob_max_idx];
}

int svm_save_model(const char *model_file_name, const svm_model *model)
{
  FILE *fp = fopen(model_file_name, "w");

  if( fp == NULL )
    {
    return -1;
    }
  const svm_parameter & param = model->param;
  fprintf(fp, "gamma %g\n", param.gamma);
  int l = model->l;
  fprintf(fp, "total_sv %d\n", l);
  fprintf(fp, "rho");
  fprintf(fp, " %g", model->rho[0]);
  fprintf(fp, "\n");
  if( model->probA )  // regression has probA only
    {
    fprintf(fp, "probA");
    fprintf(fp, " %g", model->probA[0]);
    fprintf(fp, "\n");
    }
  if( model->probB )
    {
    fprintf(fp, "probB");
    fprintf(fp, " %g", model->probB[0]);
    fprintf(fp, "\n");
    }
  if( model->nSV )
    {
    fprintf(fp, "nr_sv");
    for( int i = 0; i < 2; i++ )
      {
      fprintf(fp, " %d", model->nSV[i]);
      }
    fprintf(fp, "\n");
    }
  fprintf(fp, "SV\n");
  const double *const *  sv_coef = model->sv_coef;
  const svm_node *const *SV = model->SV;
  for( int i = 0; i < l; i++ )
    {
    fprintf(fp, "%.16g ", sv_coef[0][i]);
    const svm_node *p = SV[i];
    while( p->index != -1 )
      {
      fprintf(fp, "%d:%.8g ", p->index, p->value);
      p++;
      }

    fprintf(fp, "\n");
    }
  fclose(fp);
  return 0;
}

svm_model * svm_load_model(const char *model_file_name)
{
  FILE *fp = fopen(model_file_name, "rb");

  if( fp == NULL )
    {
    return NULL;
    }
  svm_model *     model = Malloc(svm_model, 1);
  svm_parameter & param = model->param;
  model->rho = NULL;
  model->probA = NULL;
  model->probB = NULL;
  model->label = NULL;
  model->nSV = NULL;
  model->label = Malloc(int, 2);
  model->label[0] = -1;
  model->label[1] = 1;
  param.svm_type = 0;
  param.kernel_type = 2;
  char cmd[81];
  while( 1 )
    {
    fscanf(fp, "%80s", cmd);
    if( strcmp(cmd, "gamma") == 0 )
      {
      fscanf(fp, "%lf", &param.gamma);
      }
    else if( strcmp(cmd, "total_sv") == 0 )
      {
      fscanf(fp, "%d", &model->l);
      }
    else if( strcmp(cmd, "rho") == 0 )
      {
      model->rho = Malloc(double, 1);
      fscanf(fp, "%lf", &model->rho[0]);
      }
    else if( strcmp(cmd, "probA") == 0 )
      {
      model->probA = Malloc(double, 1);
      fscanf(fp, "%lf", &model->probA[0]);
      }
    else if( strcmp(cmd, "probB") == 0 )
      {
      model->probB = Malloc(double, 1);
      fscanf(fp, "%lf", &model->probB[0]);
      }
    else if( strcmp(cmd, "nr_sv") == 0 )
      {
      model->nSV = Malloc(int, 2);
      for( int i = 0; i < 2; i++ )
        {
        fscanf(fp, "%d", &model->nSV[i]);
        }
      }
    else if( strcmp(cmd, "SV") == 0 )
      {
      while( 1 )
        {
        int c = getc(fp);
        if( c == EOF || c == '\n' )
          {
          break;
          }
        }

      break;
      }
    else
      {
      fprintf(stderr, "unknown text in model file\n");
      free(model->rho);
      free(model->label);
      free(model->nSV);
      free(model);
      return NULL;
      }
    }

  int  elements = 0;
  long pos = ftell(fp);
  bool loop = true;
  while( loop )
    {
    int c = fgetc(fp);

    switch( c )
      {
      case '\n':
      // count the '-1' element
      case ':':
        ++elements;
        break;
      case EOF:
        loop = false;
        break;
      default:
        ;
      }
    }

  fseek(fp, pos, SEEK_SET);
  int m = 1;
  int l = model->l;
  model->sv_coef = Malloc(double *, m);
  int i;
  for( i = 0; i < m; i++ )
    {
    model->sv_coef[i] = Malloc(double, l);
    }
  model->SV = Malloc(svm_node *, l);
  svm_node *x_space = NULL;
  if( l > 0 )
    {
    x_space = Malloc(svm_node, elements);
    }
  int j = 0;
  for( i = 0; i < l; i++ )
    {
    model->SV[i] = &x_space[j];
    for( int k = 0; k < m; k++ )
      {
      // cppcheck-suppress invalidscanf
      fscanf(fp, "%lf", &model->sv_coef[k][i]);
      }
    while( 1 )
      {
      int c;

      do
        {
        c = getc(fp);
        if( c == '\n' )
          {
          goto out2;
          }
        }
      while( isspace(c) );

      ungetc(c, fp);
      fscanf( fp, "%d:%lf", &( x_space[j].index ), &( x_space[j].value ) );
      ++j;
      }

out2:
    x_space[j++].index = -1;
    }
  fclose(fp);
  model->free_sv = 1;  // XXX
  return model;
}

void svm_destroy_model(svm_model *model)
{
  if( model->free_sv && model->l > 0 )
    {
    free( (void *)( model->SV[0] ) );
    }
  free(model->sv_coef[0]);
  free(model->SV);
  free(model->sv_coef);
  free(model->rho);
  free(model->label);
  free(model->probA);
  free(model->probB);
  free(model->nSV);
  free(model);
}

void svm_destroy_param(svm_parameter *param)
{
  free(param->weight_label);
  free(param->weight);
}
