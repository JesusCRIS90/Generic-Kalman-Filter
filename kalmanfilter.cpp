#include "kalmanfilter.hpp"

#define GET_VARIABLE_NAME(Variable) (#Variable)

// Verify that the number of cols and rows of the matrix it is correct
void KalmanVerify::Verify(BasicMatrix& mat2verify, arma::uword rows, arma::uword cols)
{
    //std::cout << "No Rows:" << mat2verify.n_rows << "\t No Cols:" << mat2verify.n_cols << std::endl;
    try {
        if( (rows == 0) || (cols == 0) ) { throw std::string("Invalid Dimension"); }
        if( (mat2verify.n_rows != rows) || (mat2verify.n_cols != cols) ) { throw std::string("Size Invalid");}
    } catch (std::string str) {
        if ( (str.compare("Invalid Dimension") != 0) || (str.compare("Size Invalid") != 0) )
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
        if( vec2verify.n_rows != n ) { throw std::string("Incorrect Size Vector");}
    } catch (std::string str) {
        std::cout << str << std::endl;
        if ( str.compare("Parameter N Invalid") != 0) {throw std::exception();}
        if ( str.compare("Not an SquareMatrix") != 0) {throw std::exception();}
        if ( str.compare("Incorrect Size Vector") != 0) {throw std::exception();}
    }
}




/*-------------------------------------------------------------------------------------------------
  -------------------------------------------------------------------------------------------------
  -------------------------------------------------------------------------------------------------*/

KalmanParameters::KalmanParameters(BasicVector X_0, BasicVector Z_0, BasicMatrix q, BasicMatrix r, BasicVector U_0, BasicMatrix P_0)
{
    try {
        // 1.- Verification that X_0 is not empty
        std::cout << "Veifiying: " << GET_VARIABLE_NAME(X_0) << std::endl;
        X_0.is_empty() ? throw std::string("EmptyVector") : this->Nx = X_0.n_rows;

        // 2.- Verification that Z_0 is not empty
        std::cout << "Veifiying: " << GET_VARIABLE_NAME(Z_0) << std::endl;
        Z_0.is_empty() ? throw std::string("EmptyVector") : this->Nz = Z_0.n_rows;

        // 3.- Verification that U_0
        std::cout << "Veifiying: " << GET_VARIABLE_NAME(U_0) << std::endl;
        //U_0.is_empty() ? this->Nu = this->Nx : this->Nu = U_0.n_cols;
        U_0.is_empty() ? this->Nu = this->Nx : this->Nu = U_0.n_rows;


        // Filling Datas
        U_0.is_empty() ? this->Uk = BasicVector(Nu, arma::fill::zeros) : this->Uk = U_0;
        this->Xk_k_corregido = X_0;

        set_A().set_C().set_B().set_D();
        set_Q(q).set_Pk(P_0).set_R(r);


    } catch (std::string str) {

        std::cout << "Exception Caatch: " << str << std::endl;
        if ( str.compare("EmptyVector") != 0) {throw std::string(str);}
        throw std::exception();
    }
}

KalmanParameters::KalmanParameters(BasicMatrix A, BasicMatrix B, BasicMatrix C, BasicMatrix Q, BasicMatrix R,
                                   BasicVector X0, BasicVector Z0, BasicMatrix Pk0, BasicMatrix D, BasicVector U0)
{
    try {
        // 1.- Checking A
        A.is_empty() ? throw std::string("A_Empty") : this->Nx = A.n_rows;
        set_A(A);

        // 2.- Checking C
        C.is_empty() ? throw std::string("C_Empty") : this->Nz = C.n_rows;
        set_C(C);

        // 3.- Checking B
        B.is_empty() ? throw std::string("B_Empty") : this->Nu = B.n_cols;
        set_B(B);


        std::cout << "Nx:" << Nx << "\t Nu:" << Nu << "\t Nz:" << this->Nz << std::endl;
        // Checking the rest of data
        Q.is_empty() ? throw std::string("Q_Empty") : set_Q(Q);
        R.is_empty() ? throw std::string("R_Empty") : set_R(R);



        if(X0.is_empty())
            throw std::string("X0_Empty");
        else {
            std::cout << "Veifiying: " << GET_VARIABLE_NAME(X0) << std::endl;
            KalmanVerify::Verify(X0, Nx);
            this->Xk_k_corregido = X0;
        }

        if(Z0.is_empty())
            throw std::string("Z0_Empty");
        else {
            std::cout << "Veifiying: " << GET_VARIABLE_NAME(Z0) << std::endl;
            KalmanVerify::Verify(Z0, Nz);
        }

        if(U0.is_empty())
            this->Uk = BasicVector(this->Nu, arma::fill::zeros);
        else {
            std::cout << "Veifiying: " << GET_VARIABLE_NAME(U0) << std::endl;
            KalmanVerify::Verify(U0, Nu);
            this->Uk = U0;
        }

        set_Pk(Pk0).set_D(D);

    } catch (std::string str) {
        std::cout << "Exception Caatch: " << str << std::endl;
        if ( str.compare("A_Empty") != 0) {throw std::string(str);}
        if ( str.compare("C_Empty") != 0) {throw std::string(str);}
        if ( str.compare("B_Empty") != 0) {throw std::string(str);}
        if ( str.compare("Q_Empty") != 0) {throw std::string(str);}
        if ( str.compare("R_Empty") != 0) {throw std::string(str);}
        if ( str.compare("X0_Empty") != 0) {throw std::string(str);}
        if ( str.compare("Z0_Empty") != 0) {throw std::string(str);}
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
    this->Nx = ptr->Nx;
    this->Nu = ptr->Nu;
    this->Nz = ptr->Nz;
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
    this->Nx = ptr.Nx;
    this->Nu = ptr.Nu;
    this->Nz = ptr.Nz;
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
    this->Nx = ptr.Nx;
    this->Nu = ptr.Nu;
    this->Nz = ptr.Nz;
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
    if(a.is_empty()){
        this->A = arma::eye(Nx, Nx);
    }
    else {
        KalmanVerify::Verify(a, Nx, Nx);
        this->A = a;
    }

    return *this;
}

KalmanParameters &KalmanParameters::set_B(BasicMatrix b)
{
    std::cout << "Veifiying: " << GET_VARIABLE_NAME(b) << std::endl;
    if( b.is_empty() ){
        this->B = BasicMatrix(Nx, Nu, arma::fill::zeros);
    }
    else {
        KalmanVerify::Verify(b, Nx, Nu);
        this->B = b;
    }
    return *this;
}

KalmanParameters &KalmanParameters::set_C(BasicMatrix c)
{
    std::cout << "Veifiying: " << GET_VARIABLE_NAME(c) << std::endl;
    if( c.is_empty() ){
        this->C = BasicMatrix(Nz, Nx, arma::fill::ones);
    }
    else {
        KalmanVerify::Verify(c, Nz, Nx);
        this->C = c;
    }
    return *this;
}

KalmanParameters &KalmanParameters::set_D(BasicMatrix d)
{
    std::cout << "Veifiying: " << GET_VARIABLE_NAME(d) << std::endl;
    if( d.is_empty() ){
        this->D = BasicMatrix(Nz, Nu, arma::fill::zeros);
    }
    else {
        KalmanVerify::Verify(d, Nz, Nu);
        this->D = d;
    }
    return *this;
}

KalmanParameters &KalmanParameters::set_Q(BasicMatrix q)
{
    std::cout << "Veifiying: " << GET_VARIABLE_NAME(q) << std::endl;
    if( q.is_empty() ){
        this->Q = arma::eye(Nx, Nx);
    }
    else {
        KalmanVerify::Verify(q, Nx, Nx);
        this->Q = q;
    }
    return *this;
}

KalmanParameters &KalmanParameters::set_R(BasicMatrix r)
{
    std::cout << "Veifiying: " << GET_VARIABLE_NAME(r) << std::endl;
    if( r.is_empty() ){
        this->R = arma::eye(Nz, Nz);
    }
    else {
        std::cout << "No Rows:" << r.n_rows << "\t No Cols:" << r.n_cols << std::endl;
        KalmanVerify::Verify(r, Nz, Nz);
        this->R = r;
    }
    return *this;
}

KalmanParameters &KalmanParameters::set_Pk(BasicMatrix p)
{
    std::cout << "Veifiying: " << GET_VARIABLE_NAME(p) << std::endl;
    if( p.is_empty() ){
        this->Pk_k_corregida = arma::eye(Nx, Nx);
    }
    else {
        KalmanVerify::Verify(p, Nx, Nx);
        this->Pk_k_corregida = p;
    }
    return *this;
}

void KalmanParameters::copyParameters(KalmanParameters* data2copy)
{
    this->Nx = data2copy->Nx; this->Nu = data2copy->Nu; this->Nz = data2copy->Nz;
    set_A(data2copy->A).set_B(data2copy->B).set_C(data2copy->C).set_D(data2copy->D);
    set_Q(data2copy->Q).set_R(data2copy->R).set_Pk(data2copy->Pk_k_corregida);
    this->Uk = data2copy->Uk;
}

void KalmanParameters::copyParameters(KalmanParameters& data2copy)
{
    this->Nx = data2copy.Nx; this->Nu = data2copy.Nu; this->Nz = data2copy.Nz;
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

KalmanFilter::KalmanFilter(BasicVector X_0, BasicVector Z_0, BasicMatrix q, BasicMatrix r, BasicVector U_0, BasicMatrix P_0)
{
    this->Init(X_0, Z_0, q, r, U_0, P_0);
}

KalmanFilter::KalmanFilter(BasicMatrix A, BasicMatrix B, BasicMatrix C, BasicMatrix Q, BasicMatrix R,
                            BasicVector X0, BasicVector Z0, BasicMatrix Pk0, BasicMatrix D, BasicVector U0)
{
    this->Init(A, B, C, Q, R, X0, Z0, Pk0, D, U0);
}


KalmanFilter::KalmanFilter(KalmanParameters newparam)
{
    params = new KalmanParameters(newparam);
}

KalmanFilter::~KalmanFilter()
{
    delete params;
}

void KalmanFilter::Init(BasicVector X_0, BasicVector Z_0, BasicMatrix q, BasicMatrix r, BasicVector U_0, BasicMatrix P_0)
{
//        std::cout << "-------------------------------------------------" << std::endl;
//    std::cout << "No Rows:" << X_0.n_rows << "\t No Cols:" << X_0.n_cols << std::endl;
//    std::cout << "No Rows:" << Z_0.n_rows << "\t No Cols:" << Z_0.n_cols << std::endl;
//    std::cout << "No Rows:" << q.n_rows << "\t No Cols:" << q.n_cols << std::endl;
//    std::cout << "No Rows:" << r.n_rows << "\t No Cols:" << r.n_cols << std::endl;
//    std::cout << "No Rows:" << U_0.n_rows << "\t No Cols:" << U_0.n_cols << std::endl;
//    std::cout << "No Rows:" << P_0.n_rows << "\t No Cols:" << P_0.n_cols << std::endl;
//        std::cout << "-------------------------------------------------" << std::endl;

    params = new KalmanParameters(X_0, Z_0, q, r, U_0, P_0);
}

void KalmanFilter::Init(BasicMatrix A, BasicMatrix B, BasicMatrix C, BasicMatrix Q, BasicMatrix R,
                        BasicVector X0, BasicVector Z0, BasicMatrix Pk0, BasicMatrix D, BasicVector U0)
{
    params = new KalmanParameters(A, B, C, Q, R, X0, Z0, Pk0, D, U0);
}


void KalmanFilter::update(BasicVector input, BasicVector newUk)
{
        KalmanVerify::Verify(input, this->params->Nz);
        current_input = input;
        if( !newUk.is_empty() ) { KalmanVerify::Verify(newUk, this->params->Nu); this->params->Uk = newUk; }
        PredictionPhase();
        CorrectionPhase();
        //return params->Xk_k_corregido;
}

void KalmanFilter::SetConfiguration(KalmanParameters new_params)
{
    params->copyParameters(new_params);
}

KalmanParameters KalmanFilter::Configuration()
{
    KalmanParameters datos2send(params);
    return datos2send;
}

void KalmanFilter::CorrectionPhase()
{
    params->Yk_estimada = current_input - (params->C*params->Xk_km1_estimado/* + params->D*params->Uk*/);

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
    std::cout << "Nx = " << param.Nx << "\t Nu = " << param.Nu << "\t Nz = " << param.Nz<< std::endl;
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

    std::cout << "Nx = " << filter.params->Nx << "\t Nu = " << filter.params->Nu << "\t Nz = " << filter.params->Nz << std::endl;
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


