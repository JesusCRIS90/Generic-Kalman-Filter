#include "kalmanfilter.hpp"

#define GET_VARIABLE_NAME(Variable) (#Variable)

void KalmanVerify::Verify(BasicMatrix& mat2verify, arma::uword n)
{
    try {
        if(n <= 0) { throw std::string("Parameter N Invalid"); }
        if(mat2verify.n_rows != mat2verify.n_cols) { throw std::string("Not an SquareMatrix");}
        //if(mat2verify.n_rows != n) { throw std::string("Invalid Dimension for the Matrix");}
    } catch (std::string str) {

        if ( (str.compare("Parameter N Invalid") != 0) || (str.compare("Not an SquareMatrix") != 0) )
        {
            std::cout << "Exception Cacth: " << str <<std::endl;
            throw std::exception();
        }
        else {
            throw std::string(str);
        }
    }
}


void KalmanVerify::Verify(BasicVector& vec2verify, arma::uword n)
{
    try {
        if(n <= 0) { throw std::string("Parameter N Invalid"); }
        if( vec2verify.is_empty() ) { throw std::string("Empty Vector");}
    } catch (std::string str) {
        if ( (str == "Parameter N Invalid") != 0) {throw std::string(str);}
        if ( (str == "Not an SquareMatrix") != 0) {throw std::string(str);}
        std::cout << str << std::endl;
    }
}




/*-------------------------------------------------------------------------------------------------
  -------------------------------------------------------------------------------------------------
  -------------------------------------------------------------------------------------------------*/

KalmanParameters::KalmanParameters(BasicVector init_vector, BasicMatrix q, BasicMatrix r, BasicVector uk)
{
    try {
        std::cout << "Veifiying: " << GET_VARIABLE_NAME(init_vector) << std::endl;
        init_vector.is_empty() ? throw std::string("EmptyVector") : this->N = init_vector.n_rows;
        if(!uk.is_empty()){
            if(uk.n_rows != init_vector.n_rows) { throw std::string("Different dimmension for Zk and Uk");}
            else { this->Uk = uk; set_B().set_D(); }
        }
        else {
            set_B().set_D();
            this->Uk = BasicVector(N, arma::fill::zeros);
        }
        this->Xk_k_corregido = init_vector;

        set_Q(q).set_R(r);
        set_A().set_C().set_Pk();


    } catch (std::string str) {

        std::cout << "Exception Caatch: " << str << std::endl;
        if ( str.compare("EmptyVector") != 0) {throw std::string(str);}
        if ( str.compare("Different dimmension for Zk and Uk") != 0) {throw std::string(str);}
        throw std::exception();
    }
}

KalmanParameters::KalmanParameters(KalmanParameters* ptr)
{
    // Not necessary make a verification because it need to pass it an Correct Object KalmanParameters
    // The normal constructor verify if the object it was correct
    this->A = ptr->A;
    this->B = ptr->B;
    this->C = ptr->C;
    this->D = ptr->D;
    this->Q = ptr->Q;
    this->R = ptr->R;
    this->N = ptr->N;
    this->Uk = ptr->Uk;
    this->KalmanGain = ptr->KalmanGain;
    this->Yk_estimada = ptr->Yk_estimada;
    this->Pk_k_corregida = ptr->Pk_k_corregida;
    this->Xk_k_corregido = ptr->Xk_k_corregido;
    this->Pk_km1_estimada = ptr->Pk_km1_estimada;
    this->Xk_km1_estimado = ptr->Xk_km1_estimado;
}

KalmanParameters::KalmanParameters(KalmanParameters&& ptr)
{
    // Not necessary make a verification because it need to pass it an Correct Object KalmanParameters
    // The normal constructor verify if the object it was correct
    this->A = ptr.A;
    this->B = ptr.B;
    this->C = ptr.C;
    this->D = ptr.D;
    this->Q = ptr.Q;
    this->R = ptr.R;
    this->N = ptr.N;
    this->Uk = ptr.Uk;
    this->KalmanGain = ptr.KalmanGain;
    this->Yk_estimada = ptr.Yk_estimada;
    this->Pk_k_corregida = ptr.Pk_k_corregida;
    this->Xk_k_corregido = ptr.Xk_k_corregido;
    this->Pk_km1_estimada = ptr.Pk_km1_estimada;
    this->Xk_km1_estimado = ptr.Xk_km1_estimado;
}

KalmanParameters::KalmanParameters(KalmanParameters& ptr)
{
    // Not necessary make a verification because it need to pass it an Correct Object KalmanParameters
    // The normal constructor verify if the object it was correct
    this->A = ptr.A;
    this->B = ptr.B;
    this->C = ptr.C;
    this->D = ptr.D;
    this->Q = ptr.Q;
    this->R = ptr.R;
    this->N = ptr.N;
    this->Uk = ptr.Uk;
    this->KalmanGain = ptr.KalmanGain;
    this->Yk_estimada = ptr.Yk_estimada;
    this->Pk_k_corregida = ptr.Pk_k_corregida;
    this->Xk_k_corregido = ptr.Xk_k_corregido;
    this->Pk_km1_estimada = ptr.Pk_km1_estimada;
    this->Xk_km1_estimado = ptr.Xk_km1_estimado;
}

KalmanParameters &KalmanParameters::set_A(BasicMatrix a)
{
    std::cout << "Veifiying: " << GET_VARIABLE_NAME(a) << std::endl;
    KalmanVerify::Verify(a, N);
    a.is_empty() ? this->A = arma::eye(N, N) : this->A = a;
    return *this;
}

KalmanParameters &KalmanParameters::set_B(BasicMatrix b)
{
    std::cout << "Veifiying: " << GET_VARIABLE_NAME(b) << std::endl;
    KalmanVerify::Verify(b, N);
    b.is_empty() ? (  Uk.is_empty() ? B = arma::zeros(N, N) : B = arma::eye(N, N) ) : B = b;
    //b.is_empty() ? ( int(arma::accu(Uk)) == 0 ? B = arma::zeros(N, N) : B = arma::eye(N, N) ) : B = b;
    return *this;
}

KalmanParameters &KalmanParameters::set_C(BasicMatrix c)
{
    std::cout << "Veifiying: " << GET_VARIABLE_NAME(c) << std::endl;
    KalmanVerify::Verify(c, N);
    c.is_empty() ? this->C = arma::eye(N, N) : this->C = c;
    return *this;
}

KalmanParameters &KalmanParameters::set_D(BasicMatrix d)
{
    std::cout << "Veifiying: " << GET_VARIABLE_NAME(d) << std::endl;
    KalmanVerify::Verify(d, N);
    d.is_empty() ? ( Uk.is_empty() ? D = arma::zeros(N, N) : D = arma::eye(N, N) ) : D = d;
    return *this;
}

KalmanParameters &KalmanParameters::set_Q(BasicMatrix q)
{
    std::cout << "Veifiying: " << GET_VARIABLE_NAME(q) << std::endl;
    //KalmanVerify::Verify(q, this->N);
    try {
        if(N <= 0)                  { throw std::string("Parameter N Invalid"); }
        if(q.n_rows != q.n_cols)    { throw std::string("Not an SquareMatrix"); }
        q.is_empty() || (q.n_rows != N) ? throw std::string("Invalid Input Matrix") : this->Q = q;
    } catch (std::string str) {

        if ( (str.compare("Parameter N Invalid") != 0) ||
             (str.compare("Not an SquareMatrix") != 0) ||
             (str.compare("Invalid Input Matrix") != 0) )
        {
            std::cout << "Exception Cacth: " << str <<std::endl;
            throw std::exception();
        }
        else {
            throw std::string(str);
        }

    }
    return *this;
}

KalmanParameters &KalmanParameters::set_R(BasicMatrix r)
{
    std::cout << "Veifiying: " << GET_VARIABLE_NAME(r) << std::endl;
    //KalmanVerify::Verify(r, this->N);
    try {
        if(N <= 0)                       { throw std::string("Parameter N Invalid"); }
        if(r.n_rows != r.n_cols)         { throw std::string("Not an SquareMatrix");}
        r.is_empty() || (r.n_rows != N) ? throw std::string("Invalid Input Matrix") : this->R = r;
    } catch (std::string str) {

        if ( (str.compare("Parameter N Invalid") != 0) ||
             (str.compare("Not an SquareMatrix") != 0) ||
             (str.compare("Invalid Input Matrix") != 0) )
        {
            std::cout << "Exception Cacth: " << str <<std::endl;
            throw std::exception();
        }
        else {
            throw std::string(str);
        }
    }
    return *this;
}

KalmanParameters &KalmanParameters::set_Pk(BasicMatrix p)
{
    std::cout << "Veifiying: " << GET_VARIABLE_NAME(p) << std::endl;
    KalmanVerify::Verify(p, N);
    p.is_empty() ? this->Pk_k_corregida = arma::eye(N, N) : this->Pk_k_corregida = p;
    return *this;
}

void KalmanParameters::copyParameters(KalmanParameters* data2copy)
{
    this->N = data2copy->N;
    set_A(data2copy->A).set_B(data2copy->B).set_C(data2copy->C).set_D(data2copy->D);
    set_Q(data2copy->Q).set_R(data2copy->R).set_Pk(data2copy->Pk_k_corregida);
    this->Uk = data2copy->Uk;
}

void KalmanParameters::copyParameters(KalmanParameters& data2copy)
{
    this->N = data2copy.N;
    set_A(data2copy.A).set_B(data2copy.B).set_C(data2copy.C).set_D(data2copy.D);
    set_Q(data2copy.Q).set_R(data2copy.R).set_Pk(data2copy.Pk_k_corregida);
    this->Xk_k_corregido = data2copy.Xk_k_corregido;
    this->Uk = data2copy.Uk;
}




/* ---------------------------------------------------------------------------------------------
   ---------------------------------------------------------------------------------------------
   ---------------------------------------------------------------------------------------------
   --------------------------------------------------------------------------------------------- */


//KalmanFilter::KalmanFilter()
//{

//}

KalmanFilter::KalmanFilter(BasicVector zk_input, BasicMatrix Q, BasicMatrix R, BasicVector uk_input)
{
    this->Init(zk_input, Q, R, uk_input);
}

KalmanFilter::KalmanFilter(KalmanParameters newparam)
{
    params = new KalmanParameters(newparam);
}

KalmanFilter::~KalmanFilter()
{
    delete params;
}

void KalmanFilter::Init(BasicVector zk_input, BasicMatrix Q, BasicMatrix R, BasicVector uk_input)
{
    params = new KalmanParameters(zk_input, Q, R, uk_input);
}

BasicVector KalmanFilter::update(BasicVector input)
{
    current_input = input;
    PredictionPhase();
    CorrectionPhase();
    return params->Xk_k_corregido;
}

void KalmanFilter::SetConfiguration(KalmanParameters new_params)
{
    params->copyParameters(new_params);
}

KalmanParameters KalmanFilter::Configuration()
{
    // REVISAR ESTA PARTE
    KalmanParameters datos2send(params);
    return datos2send;
}

void KalmanFilter::CorrectionPhase()
{
    params->Yk_estimada = current_input - (params->C*params->Xk_km1_estimado + params->D*params->Uk);

    BasicMatrix inv = arma::inv(params->C*params->Pk_km1_estimada*params->C.t() + params->R);

    params->KalmanGain = params->Pk_km1_estimada*params->C.t()*inv;

    params->Xk_k_corregido = params->Xk_km1_estimado + params->KalmanGain*params->Yk_estimada;

    params->Pk_k_corregida = params->Pk_km1_estimada - params->KalmanGain*params->C*params->Pk_km1_estimada;

}

void KalmanFilter::PredictionPhase()
{
    params->Xk_km1_estimado = (params->A)*(params->Xk_k_corregido) + (params->B)*(params->Uk);
    params->Pk_km1_estimada = params->A*params->Pk_k_corregida*params->A.t() + params->Q;
}



std::ostream &operator<<(std::ostream& salida, const KalmanParameters& param)
{
    std::cout << "Orden del Filtro = " << param.N << std::endl;
    std::cout << "A = " << param.A << std::endl;
    std::cout << "B = " << param.B << std::endl;
    std::cout << "C = " << param.C << std::endl;
    std::cout << "D = " << param.D << std::endl;
    std::cout << "Q = " << param.Q << std::endl;
    std::cout << "R = " << param.R << std::endl;

    std::cout << "----------------------------------------------------"  << std::endl;

    std::cout << "X corregida = " << param.Xk_k_corregido << std::endl;
    std::cout << "X estimado = " << param.Xk_km1_estimado << std::endl;
    std::cout << "Pk corregida = " << param.Pk_k_corregida << std::endl;
    std::cout << "Pk estimada = " << param.Pk_km1_estimada << std::endl;
    std::cout << "Yk estimada = " << param.Yk_estimada << std::endl;
    std::cout << "Ganancia de Kalman = " << param.KalmanGain << std::endl;
    std::cout << "Uk = " << param.Uk << std::endl;

    return salida;
}

std::ostream& operator<<(std::ostream& salida, const KalmanFilter& filter)
{

    std::cout << "Orden del Filtro = " << filter.params->N << std::endl;
    std::cout << "A = " << filter.params->A << std::endl;
    std::cout << "B = " << filter.params->B << std::endl;
    std::cout << "C = " << filter.params->C << std::endl;
    std::cout << "D = " << filter.params->D << std::endl;
    std::cout << "Q = " << filter.params->Q << std::endl;
    std::cout << "R = " << filter.params->R << std::endl;

    std::cout << "----------------------------------------------------"  << std::endl;

    std::cout << "X corregida = " << filter.params->Xk_k_corregido << std::endl;
    std::cout << "X estimado = " << filter.params->Xk_km1_estimado << std::endl;
    std::cout << "Pk corregida = " << filter.params->Pk_k_corregida << std::endl;
    std::cout << "Pk estimada = " << filter.params->Pk_km1_estimada << std::endl;
    std::cout << "Yk estimada = " << filter.params->Yk_estimada << std::endl;
    std::cout << "Ganancia de Kalman = " << filter.params->KalmanGain << std::endl;
    std::cout << "Uk = " << filter.params->Uk << std::endl;

    return salida;
}


